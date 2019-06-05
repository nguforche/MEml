#' @title Mixture of Mixed Effect GBM 
#' @description 
#' Trains a Mixed Effect gradient boosted machine where the random effects are assumed to follow a 
#' mixture of gaussian distribution.    

#' @name MEmixgbm

#' @param X  data.frame with predictors 
#' @param Y  binary response vector 
#' @param groups character name of the column containing the group identifier
#' @param rand.vars random effect variables 
#' @param gbm.dist gbm loss function   
#' @param para named list of gbm training parameters 
#' @param lme.family glmer control 
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
#' \item{boost.form}{gbm formula for fitted model}
#' \item{glmer.form}{lmer4 formula} 
#' \item{glmer.CI}{estimates of mixed effect logistic regression with 
#'     approximate confidence intervals on the logit scale. More accurate values 
#'     can be obtained by bootstrap}
#' \item{fitted.probs}{fitted probabilites for final model}
#' \item{fitted.class}{fitted class labels for final model}
#' \item{train.perf}{various performance measures for final model on training set}
#' \item{threshold}{classification cut-off}
#
#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#
#' @references
#' Che Ngufor,  Holly Van Houten, Brian S. Caffo , Nilay D. Shah, Rozalina G. McCoy 
#' Mixed Effect Machine Learning: a framework for predicting longitudinal change in hemoglobin A1c,
#' in Journal of Biomedical Informatics, 2018 

#' @import lme4 caret inTrees gbm flexmix
NULL 
#
#
#' @rdname MEmixgbm  
#' @export
MEmixgbm  <- function(X, ...) UseMethod("MEmixgbm")
#
#' @rdname MEmixgbm 
#' @export
#
MEmixgbm <- function(form, dat,  
                     groups = NULL, 
                     rand.vars="1", 
                     para = NULL,   
                     tol= 1e-5, 
                     max.iter =100, 
                     include.RE =FALSE, 
                     verbose = FALSE, 
                     maxdepth=5,
                     glmer.Control=glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
                     nAGQ=0, likelihoodCheck = TRUE,
                     K=3, 
                     krange = 2:5,
                     decay = 0.05, ...){
  
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
  
  ### initial random effect component: fit a LME model with no fixed effect, ie 
  ### a priori known mean of 0 
  
  form.glmer <- as.formula(paste0("Y ~ ", paste0(c(paste0("1", collapse = "+"), "+", "(", 
                                                   paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = ""))) 
  
  glmer.fit <- glmer(form.glmer, data= dat,family = binomial, control = glmer.Control, nAGQ=0, 
                     verbose = as.numeric(verbose))
  
  pp = predict(glmer.fit, newdata = dat, type = "response")  
  pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
  w = pp*(1-pp)      	              
  Y.star <- qlogis(pp) + (Y-pp)/w

  
#   form.mix <- as.formula(paste0("Y ~ ",  paste0(c(paste0(rand.vars, collapse = "+"), "|", id), collapse = "")) )             
#   form.ran <- as.formula(paste0("Y ~",  paste0(rand.vars, collapse = "+")))
#   ### train mixture model   
# #  mixfit <- flexmix(form.mix, k = k, model = FLXMRlmer(random = ~ 1, weighted = FALSE), 
# #                    data = dat, control = list(tolerance = 10^-3))
#   mm <- stepFlexmix(form.mix, k = 2:3, model = FLXMRlmer(form.ran, random = ~1, weighted = FALSE), data = dat, 
#                     control = list(tolerance = 10^-3),verbose = TRUE)
#   mixfit <- getModel(mm, "BIC")
#   post <- mixfit@posterior$scaled

  ### get the random effect component 
  Zt <-  getME(glmer.fit, name = "Zt")
  b  <-  getME(glmer.fit, name = "b")	
  Zb <-  as.numeric(cprod(x = Zt, y = b))       
#  Zb.mix <- apply(post*Zb, 1, sum) 

  ### adjusted target     
  Target <- Y.star - Zb 
#  Target.mix <- Y.star - Zb.mix 
    
  dat[, "Target"] <- Target
#  dat[, "Target.mix"] <- Target.mix
  dat[, "Zb"] <- Zb
#  dat[, "Zb.mix"] <- Zb.mix
  
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
      
      gbmfit <- gbm(form=gbm.form, data= dat, distribution = "gaussian", n.tree = para$n.trees, weights = w, 
                    interaction.depth = para$interaction.depth, shrinkage= para$shrinkage, 
                    n.minobsinnode = para$n.minobsinnode )
    } else {
      gbmfit <- gbm(form=gbm.form, data= dat, distribution = "gaussian", n.tree = para$n.trees, weights = w, 
                    interaction.depth = para$interaction.depth, shrinkage=para$shrinkage, 
                    n.minobsinnode = para$n.minobsinnode )

    }
    
    treeList <- GBM2List(gbmfit, dat[, rhs.vars]) 
    ruleExec = extractRules(treeList,dat[, rhs.vars], ntree=floor(0.5*para$n.trees), maxdepth = maxdepth, random=FALSE)    
    ruleExec <- unique(ruleExec)    

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
#   mixture model  
#   mixfit <- flexmix(form.mix, k = k, model = FLXMRlmer(random = ~1, weighted = FALSE), 
# #                   data = dat, control = list(tolerance = 10^-3))
#     mm <- stepFlexmix(form.mix, k = krange, model = FLXMRlmer(form.ran, random = ~1, weighted = FALSE), data = dat, 
#                       control = list(tolerance = 10^-3), verbose = TRUE)
#     mixfit <- getModel(mm, "BIC")
#     post <- mixfit@posterior$scaled
    
#  get predicted probabilities and compute transformed response 
    pp <- predict(glmer.fit, newdata = dat, type = "response")
    pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
    w = pp*(1-pp)
    Y.star <- qlogis(pp) + (Y - pp)/w
    Z <-   getME(glmer.fit, name = "Z")
    b  <-  getME(glmer.fit, name = "b")	
    Zb <-  as.numeric(Z%*%b)
#    Zb.mix <- apply(post*Zb, 1, sum) 
    
    ### adjusted target     
    Target <- Y.star - Zb 
#    Target.mix <- Y.star - Zb.mix 
    dat[, "Target"] <- Target
#    dat[, "Target.mix"] <- Target.mix
    dat[, "Zb"] <- Zb
#    dat[, "Zb.mix"] <- Zb.mix
    
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

  form.mix <- as.formula(paste0("Target ~ ",  paste0(c(paste0(c(rand.vars), collapse = "+"), "|", id), collapse = "")) )             
  form.ran <- as.formula(paste0("Target ~",  paste0(c(rand.vars), collapse = "+")))
  ### train mixture model   
  #   mixfit <- flexmix(form.mix, k = k, model = FLXMRlmer(random = ~1, weighted = FALSE), 
  #                   data = dat, control = list(tolerance = 10^-3))
  mm <- stepFlexmix(form.mix, k = krange, model = FLXMRlmer(form.ran, random = ~1, weighted = FALSE), data = dat, 
                    control = list(tolerance = 10^-3), verbose = FALSE)
  mixfit <- getModel(mm, "BIC")
  post <- mixfit@posterior$scaled
  Zb.mix <- apply(post*Zb, 1, sum) 
  
  
  fitted = predict(gbmfit, newdata = dat, type = "response", n.trees = para$n.trees) + Zb
  fitted.mix = predict(gbmfit, newdata = dat, type = "response", n.trees = para$n.trees) + Zb.mix
  
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

  ruleMetric <- getRuleMetric(ruleExec,dat[,rhs.vars],target)
  ruleMetric <- pruneRule(ruleMetric,dat[,rhs.vars], target)
  ruleMetric <- unique(ruleMetric)    
  learner <- buildLearner(ruleMetric,dat[,rhs.vars], target)  
  readableLearner <- data.frame(presentRules(learner, rhs.vars))
  
  res <- list(gbmfit = gbmfit, 
              glmer.fit = glmer.fit, 
              mixfit = mixfit,
              groups = groups, 
              para=para,
              rand.vars=rand.vars, 
              rhs.vars = rhs.vars, 
              logLik=as.numeric(logLik(glmer.fit)), 
              random.effects = ranef(glmer.fit), 
              glmer.form = glmer.form, 
              glmer.CI =tab, 
              Y.star = fitted,
              Y.star.mix = fitted.mix,
              fitted.probs = pp, 
              gbm.form = gbm.form,
              fitted.class = cls, 
              train.perf = perf, 
              threshold = threshold, 
              include.RE=include.RE, 
              gbmRules = learner, 
              gbmReadableRules = readableLearner, 
              predRules = pred, 
              rule.levels = rule.levels)
  class(res) <- "MEmixgbm"         
  return(res)         
}






#' @rdname MEmixgbm  
#' @export
MEmixgbm2  <- function(X, ...) UseMethod("MEmixgbm2")
#
#' @rdname MEmixgbm 
#' @export
#
MEmixgbm2 <- function(form, dat,  
                     groups = NULL, 
                     rand.vars="1", 
                     para = NULL,   
                     tol= 1e-5, 
                     max.iter =100, 
                     include.RE =FALSE, 
                     verbose = FALSE, 
                     maxdepth=5,
                     glmer.Control=glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
                     nAGQ=0, likelihoodCheck = TRUE,
                     K=3, 
                     krange = 2:5,
                     decay = 0.05, ...){
  
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
  
  ### initial random effect component: fit a LME model with no fixed effect, ie 
  ### a priori known mean of 0 
  
  form.glmer <- as.formula(paste0("Y ~ ", paste0(c(paste0("1", collapse = "+"), "+", "(", 
                                                   paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = ""))) 
  
  glmer.fit <- glmer(form.glmer, data= dat,family = binomial, control = glmer.Control, nAGQ=0, 
                     verbose = as.numeric(verbose))
  
  pp = predict(glmer.fit, newdata = dat, type = "response")  
  pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
  w = pp*(1-pp)      	              
  Y.star <- qlogis(pp) + (Y-pp)/w
  
  
     # form.mix <- as.formula(paste0("Y ~ ",  paste0(c(paste0(rand.vars, collapse = "+"), "|", id), collapse = "")) )             
     # form.ran <- as.formula(paste0("Y ~",  paste0(rand.vars, collapse = "+")))
  #   ### train mixture model   
  # #  mixfit <- flexmix(form.mix, k = k, model = FLXMRlmer(random = ~ 1, weighted = FALSE), 
  # #                    data = dat, control = list(tolerance = 10^-3))
# #     mm <- stepFlexmix(form.mix, k = 2:3, model = FLXMRlmer(form.ran, random = ~1, weighted = FALSE), data = dat, 
#                        control = list(tolerance = 10^-3),verbose = FALSE)
#      mixfit <- getModel(mm, "BIC")
#      post <- mixfit@posterior$scaled
  
  ### get the random effect component 
  Zt <-  getME(glmer.fit, name = "Zt")
  b  <-  getME(glmer.fit, name = "b")	
  Zb <-  as.numeric(cprod(x = Zt, y = b))       
#  Zb <- apply(post*Zb, 1, sum) 
  
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
  
  form.mix <- as.formula(paste0("Target ~ ",  paste0(c(paste0(c(rand.vars), collapse = "+"), "|", id), collapse = "")) )             
  form.ran <- as.formula(paste0("Target ~",  paste0(c(rand.vars), collapse = "+")))
  
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
      
      gbmfit <- gbm(form=gbm.form, data= dat, distribution = "gaussian", n.tree = para$n.trees, weights = w, 
                    interaction.depth = para$interaction.depth, shrinkage= para$shrinkage, 
                    n.minobsinnode = para$n.minobsinnode )
    } else {
      gbmfit <- gbm(form=gbm.form, data= dat, distribution = "gaussian", n.tree = para$n.trees, weights = w, 
                    interaction.depth = para$interaction.depth, shrinkage=para$shrinkage, 
                    n.minobsinnode = para$n.minobsinnode )
      
    }
    
    treeList <- GBM2List(gbmfit, dat[, rhs.vars]) 
    ruleExec = extractRules(treeList,dat[, rhs.vars], ntree=floor(0.5*para$n.trees), maxdepth = maxdepth, random=FALSE)    
    ruleExec <- unique(ruleExec)    
    
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
    #   mixture model  
    # mixfit <- flexmix(form.mix, k = k, model = FLXMRlmer(random = ~1, weighted = FALSE), 
    #                   data = dat, control = list(tolerance = 10^-3))
         mm <- stepFlexmix(form.mix, k = krange, model = FLXMRlmer(form.ran, random = ~1, weighted = FALSE), data = dat, 
                           control = list(tolerance = 10^-3), verbose = FALSE)
         mixfit <- getModel(mm, "BIC")
         post <- mixfit@posterior$scaled
    
    #  get predicted probabilities and compute transformed response 
    pp <- predict(glmer.fit, newdata = dat, type = "response")
    pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
    w = pp*(1-pp)
    Y.star <- qlogis(pp) + (Y - pp)/w
    Z <-   getME(glmer.fit, name = "Z")
    b  <-  getME(glmer.fit, name = "b")	
    Zb <-  as.numeric(Z%*%b)
    Zb <- apply(post*Zb, 1, sum) 
    
    ### adjusted target     
    Target <- Y.star - Zb 
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
  
  ruleMetric <- getRuleMetric(ruleExec,dat[,rhs.vars],target)
  ruleMetric <- pruneRule(ruleMetric,dat[,rhs.vars], target)
  ruleMetric <- unique(ruleMetric)    
  learner <- buildLearner(ruleMetric,dat[,rhs.vars], target)  
  readableLearner <- data.frame(presentRules(learner, rhs.vars))
  
  res <- list(gbmfit = gbmfit, 
              glmer.fit = glmer.fit, 
              mixfit = mixfit,
              groups = groups, 
              para=para,
              rand.vars=rand.vars, 
              rhs.vars = rhs.vars, 
              logLik=as.numeric(logLik(glmer.fit)), 
              random.effects = ranef(glmer.fit), 
              glmer.form = glmer.form, 
              glmer.CI =tab, 
              Y.star = fitted,
              fitted.probs = pp, 
              gbm.form = gbm.form,
              fitted.class = cls, 
              train.perf = perf, 
              threshold = threshold, 
              include.RE=include.RE, 
              gbmRules = learner, 
              gbmReadableRules = readableLearner, 
              predRules = pred, 
              rule.levels = rule.levels)
  class(res) <- "MEmixgbm2"         
  return(res)         
}







