#' 
#' Calls the MEml models: MEgbm, MEgbmrules, MErfrules, MEglmtree, MECTree, etc. 
#' The training and test data can be split into lagged training and testing as described in [1]   
#
#' @name  MEml 
#' @param lag time lag between predictors and outcome: e.g if lag = 1, then we use predictors in current  
#'         vist to predict outcome in the next visit.  
#' @param classifier character or character vector with names of classification models. 
#'                    See names(\code{TrainAllModels}()). 
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
#'  \item \code{MEml.lag} Takes the full data set and calls \code{LongiLagSplit} to split data into lagged 
#'  training and testing.  \code{MEml.lag} also trains the MOB and CTree models (see [1]).  
#'  \item  \code{MEml} is the same as \code{MEml.lag}, except that you pass in the training and test set. So you can 
#'  call \code{LongiLagSlit} and pass the derived training and test sets to \code{MEml2}.  
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
#' @examples
#' \dontrun{
#' # parameter list 
#' para <- list(
#' method = "cv", # internal cross-validation method for parameter tuning. See caret package 
#' tuneLength=3, # grid size for parameter search 
#' number = 3, # number of internal cross-validation 
#' n.trees=100, # number of trees in gbm 
#' ntree = 50, # number of trees in random forest 
#' mtry = 5, # mtry in random forest 
#' interaction.depth=4,
#' shrinkage=0.01,
#' n.minobsinnode=10,
#' opt.para= TRUE, # perform parameter tuning through internal cross-validation 
#' coefReg = 0.5, 
#' coefImp=1, 
#' include.RE = FALSE,
#' con.tree = FALSE, 
#' max.iter = 10, alpha=0.05, minsize=20,maxdepth=30,  
#' K = 3, decay = 0.05, tol= 1e-5,
#' seed = 1 # random seed 
#'  )
#' data(heart.valve)
# 
#' dat <- heart.valve
#' dat$id <- as.numeric(dat$id)  ## random effect grouping variable
#' resp.vars <- "inc.lvmi"
#' id <- "id"
#' ## fixed effect variables 
#' rhs.vars <- c("sex", "age", "time", "fuyrs", "grad", "log.grad", "bsa", "lvh", "prenyha", 
#'               "redo", "size", "con.cabg", "creat", "dm", "acei", "lv", "emergenc", 
#'               "hc", "sten.reg.mix", "hs")
#' order.vars = "time"
#' rand.vars= "time"  ## random effect variables 
# 
#' ### split data into lagged training and testing 
#' ### predict two time points in advanced 
#' dd <- LongiLagSplit(dat=dat, id=id, rhs.vars=rhs.vars,resp.vars=resp.vars,
#'                    order.vars=order.vars,lag= 2)  
#' train <- dd$train
#' test <- dd$test
#
#' res <- MEml(classifier="GBM", dat.trn = train, dat.tst=test, id=id, 
#'             rhs.vars=dd$rhs.vars, resp.vars=dd$resp.vars,
#'             rand.vars=rand.vars, para=para, max.iter = 10, seed = 1, 
#'             return.model = FALSE)
#' res$GBM
#' }
#

MEml.lag <- function(lag=NULL, classifier, dat, id, rhs.vars, resp.vars, order.vars, rand.vars=NULL,  
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

Trn.Tst <- lapply(Train.Test(), function(x) x)  

mod.res <- tryCatch(
{
lapply(classifier, function(xx) {
##cat("Now Running Classifier:", xx,  "\n")
if(return.model)
Trn.Tst[[xx]](trn=trn, tst=tst, para=para, resp.vars=resp.vars, rand.vars=rand.vars, 
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

Trn.Tst <- lapply(Train.Test(), function(x) x)  

mod.res <- tryCatch(
{
  lapply(classifier, function(xx, ...) {
##cat("Now Running Classifier:", xx,  "\n")

trn <- dat.trn[, unique(c(resp.vars, rhs.vars, reg.vars, part.vars, id)), drop = FALSE]
tst <- dat.tst[, unique(c(resp.vars, rhs.vars, reg.vars, part.vars, id)), drop = FALSE]

if(return.model)
Trn.Tst[[xx]](trn=trn, tst=tst, para=para, resp.vars=resp.vars, rand.vars=rand.vars, rhs.vars=rhs.vars, 
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




















