#' @title Mixed Effect GBM 
#' @description 
#' Trains a Mixed Effect gradient boosted machine  
#' for longitudinal continuous, binary and count data. A rule based version or these methods 
#' using the \code{inTree} package is also implemented(see [1])  
  
#' @name MEgbm

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

#' @import lme4 caret partykit inTrees gbm
NULL 
#
#' @rdname MEgbm  
#' @export
MEgbm  <- function(X, ...) UseMethod("MEgbm")
#
#' @rdname MEgbm 
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mod <- MEgbm(form, rand.form, data)) 
#'
#' }
#
MEgbm <- function(X, Y, groups = NULL, rand.vars="1", para = NULL, lme.family = binomial,  
                 tol= 1e-5, max.iter =100, include.RE =TRUE, verbose = FALSE, 
                 likelihoodCheck = TRUE, ...){
                
     if(is.null(groups)) stop("please provide grouping variable")
     Y <- as.vector(Y) 
     dat <- cbind.data.frame(response = Y, X)   
     resp.vars <- "response"      
     dat[, resp.vars] <- as.numeric(factor(dat[, resp.vars]))-1  
     X[, groups] <- NULL
     Target <- Y 
     Y <- factor(Y); levels(Y) <- c("No", "Yes")	 
         
    old.lik <- -Inf    
	if(likelihoodCheck == FALSE){
	n.re <- sum(rand.vars != 1)+1  
    b.old <-rep(0, n.re*nlevels(factor(dat[, groups])))  ## initial random effect parameters 
    }

### initial random effect component: fit a LME model with no fixed effect, ie 
### a priori known mean of 0     	
form.glmer <- as.formula(paste0("response ~ ", paste0(c("1 + ", 
 "(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             

	glmer.fit <- glmer(form.glmer, data= dat,family = lme.family, 
	             control = glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
		          nAGQ=  0, verbose = as.numeric(verbose))
		          
	pp = predict(glmer.fit, newdata = dat, type = "response")  
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
	w = pp*(1-pp)    		              
  Y.star <- qlogis(pp) + (Target-pp)/w         
### get the random effect component 
   
	Zt <-  getME(glmer.fit, name = "Zt")
	b  <-  getME(glmer.fit, name = "b")	
	Zb <-  as.numeric(cprod(x = Zt, y = b))       
if(include.RE)	X[, "Zb"] <- Zb 

dat[, "tree.fit"] <- rep(0, nrow(dat))

form.glmer <- as.formula(paste0("response ~ ", paste0(c("tree.fit + ", 
             "(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             

for(ii in 1:max.iter){    	  	
#### fit boosted regression trees 
   if(para$opt.para)
	gbmfit <-  train(X, Y.star,  method = "gbm", trControl = 
	           trainControl(method = para$method, number =  para$number), 
	           verbose = verbose,  tuneLength = para$tuneLength, distribution = "gaussian") 
   else 
	gbmfit <-  train(X, Y.star,  method = "gbm", trControl = trainControl(method = "none"),  
	           verbose = verbose, tuneGrid = data.frame(n.trees = para$n.trees, 
               interaction.depth=para$interaction.depth, shrinkage=para$shrinkage), 
             distribution = "gaussian") 
#   if(type == "prob")             
#   pp <-  predict(gbmfit, newdata = X, type = "prob")[, 2]
#   else 
#   pp <-  factor(predict(gbmfit, newdata = X, type = "raw"))
    fitted <-  predict(gbmfit, newdata = X, type = "raw")  
    dat[, "tree.fit"] <- fitted

###  Fit mixed effect logistic regression model with nodes as fixed effect predictors
	glmer.fit <- glmer(form.glmer, data= dat,family = lme.family, 
	              control = glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
		          nAGQ=  0, verbose = 0)

### compute the adjusted response 
### first get the random effect component 

  Zt <-  getME(glmer.fit, name = "Zt")
  b  <-  getME(glmer.fit, name = "b")	
  Zb <-  as.numeric(cprod(x = Zt, y = b))
if(include.RE)   X[,"Zb"] <- Zb 

pp <- sigmoid(fitted + Zb) 
pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  

#   pp = predict(glmer.fit, newdata = dat, type = "response")  
#	  pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
    w = pp*(1-pp)
    Y.star <- qlogis(pp) + (Target-pp)/w

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
	print(paste("Error: ", r))    
	if( r < tol) break 
	
	} ## for loop 
	
	if(r > tol) warning("EM algorithm did not converge")

fitted = predict(gbmfit, newdata = X, type = "raw") + Zb   
pp <- sigmoid(fitted) 
pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  

  perf <- Performance.measures(pp, Y)
	threshold <- perf$threshold	
  cls <-ifelse(pp >= threshold, "Yes", "No")

### get confidence intervals for mixed effect logistic regresion: rough estimates using the SEs
	se <- sqrt(diag(as.matrix(vcov(glmer.fit)))) 
  	tab <- cbind(Est = fixef(glmer.fit), LL = fixef(glmer.fit) - 1.96 * se, 
  	UL = fixef(glmer.fit) + 1.96 *se)

res <- list(gbmfit = gbmfit, glmer.fit = glmer.fit,  groups = groups, 
         rand.vars=rand.vars, logLik=as.numeric(logLik(glmer.fit)), 
         random.effects =ranef(glmer.fit), rand.form = form.glmer, 
         glmer.CI =tab, fitted.probs = pp, 
         fitted.class = cls, train.perf = perf, threshold = threshold, 
         include.RE=include.RE, rhs.vars = colnames(X))
class(res) <- "MEgbm"         
return(res)         
}

#
#' @rdname MEgbm  
#' @export
MEgbmRules  <- function(X, ...) UseMethod("MEgbmRules")
#
#' @rdname MEgbm 
#' @export
#
MEgbmRules <- function(form, dat,  groups = NULL, rand.vars="1", para = NULL,   
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
#      mse <- mod1 <- c() 
#      grid <- para$grid 
#      ix <- sample(nrow(dat), floor(0.75*nrow(dat)))
#      X.trn <- dat[ix, ]; X.val <- dat[-ix, ]
#      for(ii in 1:nrow(grid)){
#        mod <- gbm(form=gbm.form, data= X.trn, distribution = "gaussian", n.tree = 
#                     grid[ii, "n.trees"], interaction.depth = para$interaction.depth, 
#                   shrinkage=grid[ii, "shrinkage"], 
#                   n.minobsinnode = grid[ii, "n.minobsinnode"] )
#        pp <-  predict(mod, newdata = X.val,  type = "response", n.trees = grid[ii, "n.trees"]) 
#        mse <- c(mse, as.numeric(Performance.measures(pred=pp, obs=X.val[, lhs.form(gbm.form)])[, "mse"]))                
#        mod1[[ii]] <- mod 
#      }
#      ix = which.min(mse)
#      gbmfit <- mod1[[ix]]  
#      para$n.trees <- grid[ix, "n.trees"]
#      para$shrinkage <- grid[, "shrinkage"]
#      para$n.minobsinnode <- grid[, "n.minobsinnode"]
#      
#        fitControl <- trainControl(method = para$method, number =  para$number, repeats=para$repeats)                          	                                 		
# 	fit <-  train(form=gbm.form, data=dat, method = "gbm", distribution = "gaussian", trControl = fitControl,
#		                       verbose = FALSE, tuneLength = para$tuneLength, metric = "MSE", weights = w)		                       
#	para$n.trees <- fit$bestTune$n.trees
#	para$interaction.depth  <-  fit$bestTune$interaction.depth
#	para$shrinkage  <- fit$bestTune$shrinkage
#	para$n.minobsinnode  <- fit$bestTune$n.minobsinnode                                   
#        gbmfit <- fit$finalModel 
# n.trees <- para$grid$n.trees
# shrinkage <- para$grid$shrinkage
# gbm.cv <- cv.gbm(X=dat[, rhs.form(gbm.form)],y=dat[, lhs.form(gbm.form)],dist="gaussian", 
# n.trees=n.trees, interaction.depth=para$interaction.depth, n.minobsinnode=para$n.minobsinnode, 
# shrinkage=shrinkage, bag.fraction = .5, nfolds=para$number, seed= para$seed, verbose=FALSE)

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
# 	fitControl <- trainControl(method = "none")
# 	fit <- train(form=gbm.form, data = dat, method = "gbm", distribution = "gaussian", verbose = FALSE, trControl = fitControl, 
# 	              tuneGrid = data.frame(n.trees = para$n.trees,interaction.depth=para$interaction.depth, shrinkage=para$shrinkage), 
#               metric = "MSE")      
#        gbmfit <- fit$finalModel
#       class(gbmfit) <- "GBM"
  
    }

    treeList <- GBM2List(gbmfit, dat[, rhs.vars]) 
    ruleExec = extractRules(treeList,dat[, rhs.vars], ntree=floor(0.5*para$n.trees), maxdepth = maxdepth, random=FALSE)    
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
  class(res) <- "MEgbmRules"         
  return(res)         
}








