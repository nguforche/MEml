#' 
#' Calls the MEml models: MEgbm,  MErf, and MEglm. 
#' The training and test data can be split into lagged training and testing as described in [1]   
#
#' @name  MEml 
#' @param lag time lag between predictors and outcome: e.g if lag = 1, then we use predictors in current  
#'         vist to predict outcome in the next visit.  
#' @param method string name of model. See names(\code{MEml_train}()) for list of possible models. 
#' @param dat  data frame with predictors and outcome  
#' @param id character name of the column containing the group identifier
#' @param rhs.vars  caracter vector of predictors
#' @param order.vars  order variables (usually time variable) 
#' @param rand.vars random effect variables 
#' @param reg.vars reg.vars regressors for MOB  
#' @param part.vars partitioning variables for MOB and predictors  
#' @param para named list of gbm training parameters 
#' @param max.iter maximum number of iterations 
#' @param return.model should the train model be return. Otherwise the return values is only the performance metrics 
#' @return The train MEml model and performance matrics (as data frame) if return.model = TRUE 
#
#' @details
#' \enumerate{
#'  \item \code{MEml_lag} Takes the full data set and calls \code{LongiLagSplit} to split data into lagged 
#'  training and testing.  \code{MEml_lag} also trains the MOB and CTree models (see [1]).  
#'  \item  \code{MEml} is the same as \code{MEml_lag}, except that you pass in the training and test set. So you can 
#'  call \code{LongiLagSlit} and pass the derived training and test sets to \code{MEml}.  
#'  \item  \code{MEml2} is the same as \code{MEml}, except that you don't pass in the test set. 
#'  Also, it is currently implemented only for the GLMER, MEgbm, MErf, and GBM models.  
#' } 
#' 
#' @references
#' Che Ngufor,  Holly Van Houten, Brian S. Caffo , Nilay D. Shah, Rozalina G. McCoy 
#' Mixed Effect Machine Learning: a framework for predicting longitudinal change in hemoglobin A1c, 
#' in Journal of Biomedical Informatics, 2018 

#'# 
#' @author Che Ngufor <Ngufor.Che@@mayo.edu>
#
NULL 
#
### simple train and test splits 
#' @rdname MEml    
#' @export

MEml_lag <- function(lag=NULL, classifier, dat, id, rhs.vars, resp.vars, order.vars, rand.vars=NULL,  
                     reg.vars=NULL, part.vars=NULL, para, max.iter = 10, seed = 1, return.model = TRUE){


  model <- NULL
  ## spilt data into longitudinal training and testing 
  dd <- LongiLagSplit(dat, id, rhs.vars,resp.vars,order.vars,lag=lag)  

  train <- dd$train
  test <- dd$test
  resp.vars <- dd$resp.vars
  rhs.vars <-dd$rhs.vars
  
  names(train)[names(train)%in%c(rhs.vars, resp.vars)] <- make.names(c(rhs.vars, resp.vars))
  names(test)[names(test)%in%c(rhs.vars, resp.vars)] <- make.names(c(rhs.vars, resp.vars))
  rhs.vars <- make.names(dd$rhs.vars)
  resp.vars <- make.names(dd$resp.vars)
  rownames(train) <- NULL
  rownames(test) <- NULL
 
######   
trn <- train 
tst <- test 

mod.res <- tryCatch(
{
lapply(classifier, function(xx) {
##cat("Now Running Classifier:", xx,  "\n")
if(return.model)
  Train.Test()[[xx]](trn=trn, tst=tst, para=para, resp.vars=resp.vars, rand.vars=rand.vars, 
              rhs.vars=rhs.vars, reg.vars=reg.vars, part.vars=part.vars, groups=id)
else 
cbind(Classifier = xx, Trn.Tst[[xx]](trn=trn, tst=tst, para=para, resp.vars=resp.vars, 
                                     rand.vars=rand.vars, rhs.vars=rhs.vars, reg.vars=reg.vars, 
                                     part.vars=part.vars, groups=id)$perf) 
})
}, error=function(e){ 
  cat("Error in the Expression: ",  paste(e$call, collapse= ", "), 
      ": original error message = ", e$message, "\n") 
  list()
}) ## tryCatch
collect.garbage()
names(mod.res) <- classifier
return(mod.res)
}




### simple train and test splits 
#' @rdname MEml    
#' @export
MEml <- function(classifier, dat.trn, dat.tst, id, rhs.vars, resp.vars, rand.vars=NULL,  reg.vars=NULL, part.vars=NULL, 
                             para, max.iter = 10, seed = 1, return.model = FALSE, ...){

mod.res <- tryCatch(
{
  lapply(classifier, function(xx, ...) {
##cat("Now Running Classifier:", xx,  "\n")

trn <- dat.trn[, unique(c(resp.vars, rhs.vars, reg.vars, part.vars, id)), drop = FALSE]
tst <- dat.tst[, unique(c(resp.vars, rhs.vars, reg.vars, part.vars, id)), drop = FALSE]

if(return.model)
Train.Test()[[xx]](trn=trn, tst=tst, para=para, resp.vars=resp.vars, rand.vars=rand.vars, rhs.vars=rhs.vars, 
              reg.vars=reg.vars, part.vars=part.vars, groups=id, ...)
else 
cbind(Classifier = xx, Trn.Tst[[xx]](trn= trn, tst= tst, para=para, resp.vars=resp.vars, rand.vars=rand.vars, 
                                     rhs.vars=rhs.vars, reg.vars=reg.vars, part.vars=part.vars, groups=id, ...)$perf) 
})
}, error=function(e){ 
  cat("Error in the Expression: ",  paste(e$call, collapse= ", "), 
      ": original error message = ", e$message, "\n") 
  list()
}) ## tryCatch
collect.garbage()
names(mod.res) <- classifier
return(mod.res)
}


### simple train and test splits 
#' @rdname MEml    
#' @export
MEml2 <- function(method, data, id,  resp.vars, rhs.vars, rand.vars=NULL, para,  ...){
  
  mod.res <- tryCatch(
    {
        ##cat("Now Running Classifier:", xx,  "\n")
        trn <- data[, unique(c(resp.vars, rhs.vars, rand.vars, id)), drop = FALSE]
        Train_MEml()[[method]](trn=trn, para=para, resp.vars=resp.vars, rand.vars=rand.vars, rhs.vars=rhs.vars, groups = id) 

    }, error=function(e){ 
      cat("Error in the Expression: ",  paste(e$call, collapse= ", "), 
          ": original error message = ", e$message, "\n") 
      list()
    }) ## tryCatch
  collect.garbage()
  return(mod.res)
}

















