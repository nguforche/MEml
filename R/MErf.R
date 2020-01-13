#' @title Mixed Effect random forest  
#' @description 
#' Trains a Mixed Effect random forest 
#' for longitudinal continuous, binary and count data. A rule based version or these methods 
#' using the \code{inTree} package is also implemented(see [1])  

#' @name MErf

#' @param form formula  
#' @param dat  data.frame with predictors 
#' @param groups character name of the column containing the group identifier
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
#' \item{rf.fit}{fitted random forest model}
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

#' @import lme4 caret partykit inTrees 
NULL 
#
#
#' @rdname MErf  
#' @export
MErfRules  <- function(X, ...) UseMethod("MErfRules")
#
#' @rdname MErf 
#' @export
#
MErfRules <- function(form, dat,  groups = NULL, rand.vars="1", para = NULL,   
                      tol= 1e-5, max.iter =100, include.RE =FALSE, verbose = FALSE, maxdepth=5,
                      glmer.Control=glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
                      nAGQ=0, likelihoodCheck = TRUE,
                      K=3, decay = 0.05, ...){
  
  
  if(is.null(groups)) stop("please provide grouping variable")
  rhs.vars <- rhs.form(form)
  resp.vars <- lhs.form(form)
  dat$Y <- as.numeric(factor(dat[, resp.vars]))-1
  Y <- dat$Y
  target <- ifelse(Y==1, "Yes", "No")
  
  old.lik <- -Inf    
  if(likelihoodCheck == FALSE){
    n.re <- sum(rand.vars != 1)+1  
    b.old <-rep(0, n.re*nlevels(factor(dat[, groups])))  ## initial random effect parameters 
  }
  
  ### initial random effect component: fit a LME model with no fixed effects, ie 
  ### a priori known mean of 0 
  
  form.glmer <- as.formula(paste0("Y ~ ", paste0(c(paste0("1", collapse = "+"), "+", "(", 
                                                   paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             
  glmer.fit <- glmer(form.glmer, data= dat,family = binomial, control = glmer.Control, nAGQ=0, 
                     verbose = as.numeric(verbose))
  
  pp = predict(glmer.fit, newdata = dat, type = "response")  
  pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
  w = pp*(1-pp)                      
  Y.star <- qlogis(pp) + (Y-pp)/w
  
  ### get the random effect component 
  Zt <-  getME(glmer.fit, name = "Zt")
  b  <-  getME(glmer.fit, name = "b")	
  Zb <-  as.numeric(cprod(x = Zt, y = b))       
  ### adjusted target 
  Target <- Y.star - Zb 
  dat[, "Target"] <- Target 
  dat[, "Zb"] <- Zb
  
  if(include.RE)   rhs.vars <- unique(c(rhs.vars, "Zb"))
  rf.form <- as.formula(paste0("Target ~", paste0(rhs.vars, collapse = "+")))
  
  partvars <- c("TreeConditions")
  ### glmer formula 
  glmer.form  <- as.formula(paste0("Y ~ ", paste0(c(paste0(partvars, collapse="+"), "+",
                                                    "(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             
  
  for(ii in 1:max.iter){    	  	    
    #### fit boosted regression trees 
    #     if(para$opt.para) {            
    #       fitControl <- trainControl(method = para$method, number =  para$number)                                                             	
    #       mod <-  train(form=rf.form, data = dat, method = "RRF", tuneLength = para$tuneLength, trControl = fitControl, metric = "RMSE",
    #                                   ntree = para$ntree, importance = TRUE, keep.inbag = TRUE, replace = FALSE, proximity=FALSE)		 
    #       para$mtry  <-  mod$bestTune$mtry
    #       rf.fit <- mod$finalModel 
    #       class(rf.fit) <- "randomForest"
    #             
    #     } else {
    #       fitControl <- trainControl(method = "none", verboseIter = FALSE)	        
    #       mod <- train(form=rf.form, data = dat,  method = "RRF", trControl = fitControl,
    #                                  tuneGrid = data.frame(mtry = para$mtry, coefReg=para$coefReg, coefImp=para$coefImp), 
    #                    ntree = para$ntree, importance = TRUE, 
    # #                                 proximity=FALSE, keep.inbag=TRUE, replace=FALSE)	
    #      rf.fit <- mod$finalModel
    #       class(rf.fit) <- "randomForest"
    #     
    #     }
    
    rf.fit <- RRF(dat[, rhs.form(rf.form)], dat[, lhs.form(rf.form)], ntree=para$ntree, mtry = floor(sqrt(length(rhs.form(rf.form)))), 
                  importance = FALSE, proximity=FALSE, keep.inbag=TRUE, replace=FALSE)
    
    treeList <- RF2List(rf.fit)    
    ruleExec = extractRules(treeList, dat[, rhs.form(rf.form)], ntree=floor(0.5*para$ntree), maxdepth = maxdepth, random=FALSE)    
    ruleExec <- unique(ruleExec)
    
    #    splt <- SplitVector(Target, K = K)
    #    target <- splt$newV
    
    ruleMetric <- getRuleMetric(ruleExec,dat[,rhs.vars],Target)
    ruleMetric <- pruneRule(ruleMetric,dat[,rhs.vars], Target, decay)
    ruleMetric <- unique(ruleMetric)    
    learner <- buildLearner(ruleMetric,dat[,rhs.vars], Target)
    pred <- myapplyLearner(learner, dat[,rhs.vars])
    dat[, "TreeConditions"] <- factor(pred$predCond)
    levels(dat[, "TreeConditions"]) <- paste0("Rule.", 1:nlevels(dat[, "TreeConditions"]))
    
    if(length(unique(dat[, "TreeConditions"])) < 2) dat[, "TreeConditions"] <- 1
    
    glmer.fit <- glmer(glmer.form, data= dat,  family = binomial, control = glmer.Control, 
                       nAGQ=  nAGQ, verbose = as.numeric(verbose))   
    
    #  get predicted probabilities and compute transformed response 
    pp <- predict(glmer.fit, newdata = dat, type = "response")
    pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
    w = pp*(1-pp)
    Y.star <- qlogis(pp) + (Y - pp)/w
    Z <-   getME(glmer.fit, name = "Z")
    b  <-  getME(glmer.fit, name = "b")	
    Zb <-  as.numeric(Z%*%b)
    Target  <- Y.star  - Zb
    #    target <- round(w*(Target+ Zb -  qlogis(pp) ) + pp) 
    
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
#    if(verbose) 
#    print(paste("Error: ", r))    
    if( r < tol) break     
  } ## for loop 
  
  if(r > tol) warning("EM algorithm did not converge")
  
  
  #  fitted = predict(rf.fit, newdata = dat, type = "raw") + Zb 
  zz <- predict(rf.fit, newdata = dat, type = "response")
  fitted = zz  + Zb 
  pp <- sigmoid(fitted)
  
  pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
  
  perf <- Performance.measures(pp, Y)
  
  threshold <- perf$threshold	
  cls <-ifelse(pp >= threshold, "Yes", "No")  
  names(cls) <- NULL 	
  
  rule.levels <- levels(dat[, "TreeConditions"])
  
  ### get confidence intervals for mixed effect logistic regresion: rough estimates using the SEs
  se <- sqrt(diag(as.matrix(vcov(glmer.fit)))) 
  tab <- cbind(Est = fixef(glmer.fit), LL = fixef(glmer.fit) - 1.96 * se, 
               UL = fixef(glmer.fit) + 1.96 *se)
  
  #   splitV <- splt$splitV
  #   newV <- rep(0, nrow(learner))
  #   labs <- learner[, "pred"]
  #   xx <- unique(labs)
  #   ix <- 1:length(xx)
  #   is <- which(xx == "L1")
  #   newV[which(labs == "L1")] <- splt$splitV[1]
  #   ll <- xx[xx != "L1"]
  #   names(splitV) <- ll
  #   
  #   if (length(splt$splitV) >= 2) {
  #     for (jj in 1:(length(ll)-1)) {    
  #       newV[which(labs == ll[jj])] <- 0.5*(splt$splitV[jj] + splt$splitV[jj+1])
  #     }
  #   }
  #   newV[which(labs == ll[length(ll)])] <- splt$splitV[length(ll)] + 0.5 
  
  
  ruleMetric <- getRuleMetric(ruleExec,dat[,rhs.vars],target)
  ruleMetric <- pruneRule(ruleMetric,dat[,rhs.vars], target)
  ruleMetric <- unique(ruleMetric)    
  learner <- buildLearner(ruleMetric,dat[,rhs.vars], target)  
  readableLearner <- data.frame(presentRules(learner, rhs.vars))
  
  res <- list(rf.fit = rf.fit, 
              glmer.fit = glmer.fit,  
              groups = groups, para=para,
              rand.vars=rand.vars, 
              rhs.vars = rhs.vars, 
              logLik=as.numeric(logLik(glmer.fit)), 
              random.effects =ranef(glmer.fit), 
              glmer.form = glmer.form, 
              glmer.CI =tab,
              Y.star = fitted, 
              fitted.probs = pp, 
              rf.form = rf.form,
              fitted.class = cls, 
              train.perf = perf, 
              threshold = threshold, 
              include.RE=include.RE, 
              rfRules = learner, 
              rfReadableRules = readableLearner, 
              predRules = pred, 
              rule.levels = rule.levels)
  class(res) <- "MErfRules"         
  return(res)         
}


