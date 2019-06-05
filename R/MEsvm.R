#' @title Mixed Effect support vector machine
#' @description 
#' Train a Mixed Effect support vector machine for binary outcome. 

#' @name MEsvm

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
#' \item{svmfit}{fitted svm model}
#' \item{glmer.fit}{fitted mixed effect logistic regression model}
#' \item{logLik}{log likelihood of mixed effect logistic regression} 
#' \item{random.effects}{random effect parameter estimates}
#' \item{svm.form}{svm formula for fitted model}
#' \item{glmer.form}{lmer4 formula} 
#' \item{glmer.CI}{estimates of mixed effect logistic regression with 
#'     approximate confidence intervals on the logit scale. More accurate values 
#'     can be obtained by bootstrap}
#' \item{fitted.probs}{fitted probabilites for final model}
#' \item{fitted.class}{fitted class labels for final model}
#' \item{fitted.decision}{fitted decision  values for final model} 
#' \item{train.perf}{various performance measures for final model on training set}
#' \item{threshold}{classification cut-off}
#
#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#
#' @references
#' Che Ngufor,  Holly Van Houten, Brian S. Caffo , Nilay D. Shah, Rozalina G. McCoy 
#' Mixed Effect Machine Learning: a framework for predicting longitudinal change in hemoglobin A1c,
#' in Journal of Biomedical Informatics, 2018 

#' @import lme4 caret 
NULL 
#
#' @rdname MEsvm  
#' @export
MEsvm  <- function(X, ...) UseMethod("MEsvm")
#
#' @rdname MEsvm 
#' @export
#
MEsvm <- function(form, dat,  groups = NULL, rand.vars="1", para = NULL,   
                       tol= 1e-5, max.iter =100, include.RE =FALSE, 
                       verbose = FALSE, maxdepth=5,
                       glmer.Control=glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", 
                       check.nobs.vs.nlev="ignore"), 
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
  Target <- Y.star - Zb ## f(x)
  dat[, "Target"] <- Target 
  dat[, "Zb"] <- Zb
  
  if(include.RE)   rhs.vars <- unique(c(rhs.vars, "Zb"))
  svm.form <- as.formula(paste0("Target ~", paste0(rhs.vars, collapse = "+")))
  
  partvars <- c("Target")
  ### glmer formula 
  glmer.form  <- as.formula(paste0("Y ~ ", paste0(c(paste0(partvars, collapse="+"), "+",
                                                    "(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             
  
  for(ii in 1:max.iter){    	  	    
    #### fit boosted regression trees 
    if(para$opt.para) {

      fitControl <- trainControl(method = para$method, number =  para$number, allowParallel = FALSE)  
      
      svmfit <- train(form = svm.form, data = dat, method = Methods()$LaplaceSVM, 
                      trControl = fitControl, tuneLength = para$tuneLength, metric = "RMSE")
      
    } else {
       	  fitControl <- trainControl(method = "none")
       	svmfit <- train(form=svm.form, data = dat, method = Methods()$LaplaceSVM, trControl = fitControl, 
     	              tuneGrid = para$svmGrid,  metric = "RMSE")      
 
    }
    
    dat[, "Target"] <-  predict.caret.ksvm(svmfit, newdata = dat, type = "decision")
    
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
  
  
  fitted = predict.caret.ksvm(svmfit, newdata = dat, type = "decision") + Zb   
  pp <- sigmoid(fitted) 
  pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
  
  perf <- Performance.measures(pp, Y)
  threshold <- perf$threshold	
  cls <-ifelse(pp >= threshold, "Yes", "No")  
  names(cls) <- NULL 	
  
  ### get confidence intervals for mixed effect logistic regresion: rough estimates using the SEs
  se <- sqrt(diag(as.matrix(vcov(glmer.fit)))) 
  tab <- cbind(Est = fixef(glmer.fit), LL = fixef(glmer.fit) - 1.96 * se, 
               UL = fixef(glmer.fit) + 1.96 *se)
  
  res <- list(svmfit = svmfit, 
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
              include.RE=include.RE)
  class(res) <- "MEsvm"         
  return(res)         
}

