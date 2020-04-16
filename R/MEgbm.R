#' @export
MEgbm  <- function(X, ...) UseMethod("MEgbm")


#' @title Mixed Effect GBM 
#' @description 
#' Trains a Mixed Effect gradient boosted machine  
#' for longitudinal continuous, binary and count data. A rule based version or these methods 
#' using the \code{inTree} package is also implemented(see [1])  

#' @param form formula  
#' @param dat  data.frame with predictors 
#' @param groups character name of the column containing the group identifier
#' @param family a GLM family, see \code{glm} and \code{family}
#' @param rand.vars random effect variables 
#' @param para named list of gbm training parameters 
#' @param glmer.Control glmer control 
#' @param tol convergence tolerance 
#' @param max.iter maximum number of iterations  
#' @param include.RE (logical) to include random effect Zb as predictor in gbm?  
#' @param verbose verbose for lme4 
#' @param likelihoodCheck (logical) to use log likelihood of glmer to check for convergence? 
#' @param type of predictions of gbm to pass to lme4 as population estimates (these will be used as offset) 
#' @param \dots Further arguments passed to or from other methods.
#' @return An object of class MEgbm; a list with items 
#' \item{gbmfit}{fitted gbm model}
#' \item{glmer.fit}{fitted mixed effect logistic regression model}
#' \item{logLik}{log likelihood of mixed effect logistic regression} 
#' \item{random.effects}{random effect parameter estimates}
#' \item{glmer.form}{lmer4 formula} 
#' \item{glmer.CI}{estimates of mixed effect logistic regression with 
#'     approximate confidence intervals on the logit scale. More accurate values 
#'     can be obtained by bootstrap}
#' \item{fitted.probs}{fitted probabilites for final model}
#' \item{fitted.class}{fitted class labels for final model}
#' \item{train.perf}{various performance measures for final model on training set}
#' \item{threshold}{classification cut-off}
#' \item{predRules}{fitted rules}
#' \item{Y.star}{fitted transform outcome}
#
#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#
#' @references
#' Che Ngufor,  Holly Van Houten, Brian S. Caffo , Nilay D. Shah, Rozalina G. McCoy 
#' Mixed Effect Machine Learning: a framework for predicting longitudinal change in hemoglobin A1c,
#' in Journal of Biomedical Informatics, 2018 
#' @import lme4 caret partykit inTrees gbm
#' 
#' @export
MEgbm <- function(form, dat,  
                       groups = NULL, 
                       family = "binomial", 
                       rand.vars="1", 
                       para = NULL,   
                       tol= 1e-5, 
                       max.iter =100, 
                       include.RE =FALSE, 
                       verbose = FALSE, 
                       maxdepth=5,
                       glmer.Control= glmerControl(optimizer = "bobyqa"), 
                       nAGQ=0, 
                       likelihoodCheck = TRUE,
                       K=3, 
                       decay = 0.05, ...){
  
  if(is.null(groups)) stop("please provide grouping variable")
  rhs.vars <- rhs.form(form)
  resp.vars <- lhs.form(form)
  Y <- dat$Y <- dat[, resp.vars]

  old.lik <- -Inf    
  if(likelihoodCheck == FALSE){
    n.re <- sum(rand.vars != 1)+1 
    b.old <-rep(0, n.re*nlevels(factor(dat[, groups])))  ## initial random effect parameters 
  }
  
  ## family
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  ### initial random effect component: fit a LME model with no fixed effect, ie 
  ### a priori known mean of 0 
  
  form.glmer <- as.formula(paste0("Y ~ ", paste0(c(paste0("1", collapse = "+"), "+", "(", 
                                                   paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))  

  glmer.fit <- myglmer(form = form.glmer, dat = dat,family = family, control = glmer.Control, nAGQ=0, 
                     verbose = as.numeric(verbose))
  pp = predict(glmer.fit, newdata = dat, type = "response") 
  
  if (family$family == "binomial"){
    pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
    w = pp*(1-pp)
    Y.star <- qlogis(pp) + (Y - pp)/w
  }
  else if(family$family == "gaussian") {
    Y.star = Y
    w = rep(1, length(Y))
  } else {
    print(family)
    stop("'family' is not yet implemented")
  }
  
  
  ### get the random effect component 
  Zt <-  getME(glmer.fit, name = "Zt")
  b  <-  getME(glmer.fit, name = "b")	
  Zb <-  as.numeric(cprod(x = Zt, y = b))       
  ### adjusted target 
  Target <- Y.star - Zb 
  dat[, "Target"] <- Target 
  dat[, "Zb"] <- Zb
  
  if(include.RE)   rhs.vars <- unique(c(rhs.vars, "Zb"))
  gbm.form <- as.formula(paste0("Target ~", paste0(rhs.vars, collapse = "+")))
  
  partvars <- c("TreeConditions")
  ### glmer formula 
  glmer.form  <- as.formula(paste0("Y ~ ", paste0(c(paste0(partvars, collapse="+"), "+",
                  "(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             
  
  for(ii in 1:max.iter){    	  	    
    #### fit boosted regression trees 
    if(para$opt.para) {
      
      fitControl <- trainControl(method = para$method, number =  para$number, allowParallel = FALSE)   
      gbm.cv <- train(form=gbm.form, data=dat, method = "gbm", distribution="gaussian", trControl = fitControl,
                      verbose = FALSE, tuneLength = para$tuneLength, metric = "RMSE")
      
      para$n.trees <- gbm.cv$bestTune$n.trees
      para$interaction.depth  <- gbm.cv$bestTune$interaction.depth
      para$shrinkage  <- gbm.cv$bestTune$shrinkage
      para$n.minobsinnode  <- gbm.cv$bestTune$n.minobsinnode                                   
      
      #gbmfit <- gbm.cv$finalModel
      #class(gbmfit) <- "gbm"	
      
      gbmfit <- gbm(form=gbm.form, data= dat, distribution = "gaussian", n.tree = para$n.trees, weights = w, 
                    interaction.depth = para$interaction.depth, shrinkage= para$shrinkage, 
                    n.minobsinnode = para$n.minobsinnode )
      # 
      #      class(gbmfit) <- "gbm"			     
    } else {
      gbmfit <- gbm(form=gbm.form, data= dat, distribution = "gaussian", n.tree = para$n.trees, weights = w, 
                    interaction.depth = para$interaction.depth, shrinkage=para$shrinkage, 
                    n.minobsinnode = para$n.minobsinnode )
    }
    
    treeList <- GBM2List(gbmfit, dat[, rhs.vars]) 
    ruleExec = myextractRules(treeList,dat[, rhs.vars], ntree=floor(0.5*para$n.trees), maxdepth = maxdepth, random=FALSE)    
    ruleExec <- unique(ruleExec)    
    
    
    ruleMetric <- getRuleMetric(ruleExec,dat[,rhs.vars],Target)
    ruleMetric <- pruneRule(ruleMetric,dat[,rhs.vars], Target, decay)
    ruleMetric <- unique(ruleMetric)    
    learner <- buildLearner(ruleMetric,dat[,rhs.vars], Target)
    pred <- myapplyLearner(learner, dat[,rhs.vars])
    dat[, "TreeConditions"] <- factor(pred$predCond)
    levels(dat[, "TreeConditions"]) <- paste0("Rule.", 1:nlevels(dat[, "TreeConditions"]))
    
    if(length(unique(dat[, "TreeConditions"])) < 2) dat[, "TreeConditions"] <- 1
    
    glmer.fit <- myglmer(form = glmer.form, dat= dat,  family = family, control = glmer.Control, 
                       nAGQ= nAGQ, verbose = as.numeric(verbose))   
    
    #  get predicted probabilities and compute transformed response 
    pp <- predict(glmer.fit, newdata = dat, type = "response")
    
    if (family$family == "binomial"){
      pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
      w = pp*(1-pp)
      Y.star <- qlogis(pp) + (Y - pp)/w
    }
    else if(family$family == "gaussian") {
      Y.star = Y
      w = rep(1, length(Y))
    } else {
      print(family)
      stop("'family' is not yet implemented")
    }
    
    Z <-   getME(glmer.fit, name = "Z")
    b  <-  getME(glmer.fit, name = "b")	
    Zb <-  as.numeric(Z%*%b)
    Target  <- Y.star  - Zb
    dat[, "Target"] <- Target 
    dat[, "Zb"] <- Zb
    
    ### test for convergence             
    if(likelihoodCheck){
      new.lik <- as.numeric(logLik(glmer.fit))
      r <- as.numeric(sqrt(t(new.lik - old.lik)%*%(new.lik-old.lik)))
      old.lik <- new.lik	
    } else {
      r <- as.numeric(sqrt(t(b - b.old)%*%(b-b.old)))
      b.old <- b
    } 
    #if(verbose) 
    #    print(paste("Error: ", r))    
    if( r < tol) break     
  } ## for loop 
  
  if(r > tol) warning("EM algorithm did not converge")
  
  
  fitted = predict(gbmfit, newdata = dat, type = "response", n.trees = para$n.trees) + Zb 
  
  threshold = NULL 
  if(family$family == "binomial"){
  pp <- sigmoid(fitted) 
  pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
  perf <- Performance.measures(pp, Y)
  threshold <- perf$threshold	
  }  
  
  rule.levels <- levels(dat[, "TreeConditions"])
  
  ### get confidence intervals for mixed effect logistic regresion: rough estimates using the SEs
  se <- sqrt(diag(as.matrix(vcov(glmer.fit)))) 
  tab <- cbind(Est = fixef(glmer.fit), LL = fixef(glmer.fit) - 1.96 * se, 
               UL = fixef(glmer.fit) + 1.96 *se)
  

  ruleMetric <- getRuleMetric(ruleExec,dat[,rhs.vars],Y)
  ruleMetric <- pruneRule(ruleMetric,dat[,rhs.vars], Y)
  ruleMetric <- unique(ruleMetric)    
  learner <- buildLearner(ruleMetric,dat[,rhs.vars], Y)  
  readableLearner <- data.frame(presentRules(learner, rhs.vars))
  
  res <- list(gbmfit = gbmfit, 
              glmer.fit = glmer.fit,  
              groups = groups, 
              para=para,
              rand.vars=rand.vars, 
              rhs.vars = rhs.vars, 
              logLik=as.numeric(logLik(glmer.fit)), 
              random.effects =ranef(glmer.fit), 
              glmer.form = glmer.form, 
              glmer.CI =tab, 
              Y.star = fitted,
              threshold = threshold, 
              gbm.form = gbm.form,
              include.RE=include.RE, 
              gbmRules = learner, 
              gbmReadableRules = readableLearner, 
              predRules = pred, 
              rule.levels = rule.levels)
  class(res) <- "MEgbm"         
  return(res)         
}






