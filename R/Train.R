#' Train, Validate and Test various models:   
# 
#' Functions to train, select best hyper-paramters, select best classification thresholds 
#'  and evaluate final model on the test set  using a number of machine 
#' learning classification algorithms 
# 
#' @name Train 
#' @param  classifier character list of classification models. See names(TrainAllModels()). 
#' @param  X.trn,X.val,X.tst matrix of predictors 
#' @param  Y.trn,Y.val,Y.tst matrix of binary \{0,1\} response variable
#' @param  varimp (logical) compute variable importance ?
#' @param  opt.para (logical) tune parameters ?
#' @param  para named list of model parameters 
#' @param  return.model (logical) return trained model ?  
#' @param \dots further arguments passed to or from other methods.
#' @return \code{TrainModels} returns a list of functions for training various algorithms. Currently 
#' 7 machine learning algorithms ELR, GLM, GLMnet, RF, GBM, avNNET, and SVM. Others can be easily added. 
#' \code{Train.Validate.Test} returns a list of performance measures for each classifier traoned:
#' \item{model}{trained model}
#' \item{para}{named list of model hyper-paramters (tunned values if opt.para = TRUE)}
#' \item{run.time}{compute time} 
#' \item{varimp}{variable importance if classifier is GLM, GBM or RF}
#' \item{perf}{a list with two data.frames: val.perf and tst.perf containing performance measures 
#' for validation and test sets }

#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#
NULL 
#' @rdname Train 
#' @export
#
TrainModels <- function(...){
  res = list(
GLMER = function(trn, tst, form, para, rand.vars, ...){
  
      result <- vector("list", length = 3)
      names(result) <-c("perf.tst", "para", "model") 
      class(result) <- "GLMER"
      result$para <- para 			

      result$model <- glmer(form, data=trn,family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=0)
      pp <- predict(result$model, newdata=trn,type="response")
      thresh <- Performance.measures(pp, trn[, lhs.form(form)])$threshold
      pp <- predict(result$model, newdata=tst,type="response")      
      result$perf.tst <- Performance.measures(pp, tst[, lhs.form(form)], threshold = thresh)
      result      
}, 

GLMERboost = function(trn, tst, form, para, rand.vars, ...){
  
  result <- vector("list", length = 3)
  names(result) <-c("perf.tst", "para", "model") 
  class(result) <- "GLMERboost"
  result$para <- para   		
  
  result$model <- glmerLogitBoost(form=form, dat = trn, family = binomial,max.iter = 10, verbose = FALSE)
  pp <- predict(result$model, newdata=tst)
  result$perf.tst <- Performance.measures(pp[,2], tst[, lhs.form(form)], threshold = result$model$threshold)
  result      
}, 

GBMrules = function(trn, tst, form, para, rand.vars, ...){
  
  result <- vector("list", length = 4)
  names(result) <-c("perf.val", "perf.tst", "para", "model") 
  class(result) <- "GBMrules"
  result$para <- para
  
  result$model <- MEgbmRules(form = form, dat=trn,  groups = para$groups, rand.vars= rand.vars, para = para,   
                         tol= 1e-5, max.iter = para$max.iter, include.RE = para$include.RE, 
                         verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=glmerControl(optimizer = "bobyqa"),
                         nAGQ = 0, K = para$K, para$decay)    
  pp <- predict(result$model, newdata=tst)
  result$perf.tst <- Performance.measures(pp$pred[, 2], tst[, lhs.form(form)], threshold = result$model$threshold)
  result      
}, 

MEglmTree = function(X.trn, Y.trn, X.tst, Y.tst, para, rand.vars, part.vars, reg.vars, rhs.vars, ...){
  
  result <- vector("list", length = 4)
  names(result) <-c("perf.val", "perf.tst", "para", "model") 
  class(result) <- "MEglmTree"
  result$para <- para
 
  result$model <- MEglmTree(X=X.trn, Y=Y.trn, part.vars=part.vars, reg.vars=reg.vars, rand.vars=rand.vars, 
                            groups=para$groups,include.RE = para$include.RE, max.iter = 
                              para$max.iter, alpha= para$alpha, minsize=para$minsize, maxdepth= para$maxdepth)    
  pp <- predict(result$model, newdata=X.tst)
  result$perf.tst <- Performance.measures(pp[,2], Y.tst, threshold = result$model$threshold)
  result      
}, 

MECTree = function(X.trn, Y.trn, X.tst, Y.tst, para, rand.vars, part.vars, reg.vars, rhs.vars, ...){
  
  result <- vector("list", length = 4)
  names(result) <-c("perf.val", "perf.tst", "para", "model") 
  class(result) <- "MECTree"
  result$para <- para

  result$model <- MECTree(X=X.trn,Y=Y.trn, con.tree = para$con.tree, rhs.vars=rhs.vars,
                           rand.vars= rand.vars, groups = para$groups,  max.iter=para$max.iter)    
  pp <- predict(result$model, newdata=X.tst)
  result$perf.tst <- Performance.measures(pp[,2], Y.tst, threshold = result$model$threshold)
  result      
} 
)
res
}


#' @rdname Train 
#' @export
#'
Train.Test <- function(...){
res = list(
GLMER = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){
  result <- list()
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(c(paste0(rhs.vars,collapse="+"), "+", "(", 
                     paste0(c(rand.vars),collapse = "+"), "|", groups, ")"), collapse="")))
  result$model <- glmer(form, data= trn,family=binomial,control=glmerControl(optimizer="bobyqa"),nAGQ=0)
 
  result$pred <- predict(result$model, newdata=tst,type="response")
  result$perf <- Performance.measures(result$pred, tst[, resp.vars])
  result      
}, 

MEgbm = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){  
  result <- list()
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 
  result$model <- MEgbmRules(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                         tol= 1e-5, max.iter = para$max.iter, include.RE = para$include.RE, 
                         verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=glmerControl(optimizer = "bobyqa"),
                         nAGQ = 0, K = para$K, para$decay) 

  result$pred <- predict(result$model, newdata= tst) 
  result$perf <- Performance.measures(result$pred$pred[, 2], tst[, resp.vars], threshold = result$model$threshold)
  result      
}, 


MErf = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){  
  result <- list()

  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
  result$model <- MErfRules(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                       tol= para$tol, max.iter = para$max.iter, include.RE = para$include.RE, verbose = FALSE, 
                       maxdepth= para$maxdepth,glmer.Control=glmerControl(optimizer = "bobyqa"), 
                       nAGQ = 0, K = para$K, decay = para$decay)
  result$pred <- predict(result$model, newdata=tst)
  result$perf <- Performance.measures(result$pred$pred[, 2], tst[,resp.vars], threshold = result$model$threshold)
  result      
}, 

MEmixgbm = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){  
  result <- list()
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 
  result$model <- MEmixgbm(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                             tol= 1e-5, max.iter = para$max.iter, include.RE = para$include.RE, 
                             verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=glmerControl(optimizer = "bobyqa"),
                             nAGQ = 0, K = para$K, krang = para$krange, para$decay) 
  
  result$pred <- predict(result$model, newdata= tst) 
  result$perf <- Performance.measures(result$pred$pred[, 2], tst[, resp.vars], threshold = result$model$threshold)
  result      
}, 


MEmixgbm2 = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){  
  result <- list()
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 
  result$model <- MEmixgbm2(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                           tol= 1e-5, max.iter = para$max.iter, include.RE = para$include.RE, 
                           verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=glmerControl(optimizer = "bobyqa"),
                           nAGQ = 0, K = para$K, krang = para$krange, para$decay) 
  
  result$pred <- predict(result$model, newdata= tst) 
  result$perf <- Performance.measures(result$pred$pred[, 2], tst[, resp.vars], threshold = result$model$threshold)
  result      
}, 

MEglmTree = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL,groups, ...){
  result <- list()
  X.trn <- trn
  X <- X.trn[, unique(c(rhs.vars, reg.vars, part.vars, groups)), drop = FALSE]
  Y <- X.trn[, resp.vars]
  
  result$model <- MEglmTree(X=X, Y=Y, part.vars=part.vars, reg.vars=reg.vars, rand.vars=rand.vars, 
                            groups= groups,include.RE = para$include.RE, max.iter = 
                              para$max.iter, alpha= para$alpha, minsize=para$minsize, maxdepth= para$maxdepth, 
                            para=para)

  result$pred <- predict(result$model, newdata= tst)
  result$perf <- Performance.measures(result$pred[,2], tst[, resp.vars], threshold = result$model$threshold)
  result      
}, 

MECTree = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){  
  result <- list()
  X.trn <- trn
  X <- X.trn[, unique(c(rhs.vars, reg.vars, part.vars, groups)), drop = FALSE]  
  Y <- X.trn[, resp.vars]

  result$model <- MECTree(X=X,Y=Y, con.tree = para$con.tree, rhs.vars=rhs.vars,
                           rand.vars=rand.vars, groups = groups, max.iter=para$max.iter)    
  result$pred <- predict(result$model, newdata= tst)
  result$perf <- Performance.measures(result$pred[,2], tst[, resp.vars], threshold = result$model$threshold)
  result      
}, 

GLM = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){
  result <- list()
  xx <- trn 
  xx[, resp.vars] <- factor(ifelse(xx[, resp.vars]==1, "Yes", "No"))
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
  result$model <-  train(form,  data = xx, method = "glm", family = "binomial", control=list(maxit=200))   
  result$pred = predict(result$model, newdata = tst, type = "prob")
  result$perf <- Performance.measures(result$pred[, 2], tst[, resp.vars])
  result      
}, 

GBM = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){
  result <- list()
  xx <- trn 
  xx[, resp.vars] <- factor(ifelse(xx[, resp.vars]==1, "Yes", "No"))
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
if(para$opt.para) {
  fitControl <- trainControl(method = para$method, number =  para$number, 
                             classProbs = TRUE, summaryFunction= twoClassSummary) 
  result$model <-  train(form,  data = xx, method = "gbm", trControl =  fitControl, 
                   verbose = FALSE,  metric = "ROC", tuneLength = para$tuneLength) 
} else {  
  result$model <-  train(form,  data = xx, method = "gbm", trControl = trainControl(method = "none"),  
                   verbose = FALSE, tuneGrid = data.frame(n.trees = para$n.trees, n.minobsinnode = para$n.minobsinnode, 
                    interaction.depth=para$interaction.depth, shrinkage=para$shrinkage))   
}
  result$pred = predict(result$model, newdata = tst, type = "prob")
  result$perf <- Performance.measures(result$pred[, 2], tst[, resp.vars])
  result      
}, 

RF = function(trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){
  result <- list()
  xx <- trn 
  xx[, resp.vars] <- factor(ifelse(xx[, resp.vars]==1, "Yes", "No"))
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
if(para$opt.para) {
  fitControl <- trainControl(method = para$method, number =  para$number, 
                             classProbs = TRUE, summaryFunction= twoClassSummary) 
  result$model <-  train(form,  data = xx, method = "rf", trControl =  fitControl, 
                   verbose = FALSE,  metric = "ROC", tuneLength = para$tuneLength) 
} else {  
  result$model <-  train(form,  data = xx, method = "rf", trControl = trainControl(method = "none"),  
                   verbose = FALSE, tuverbose = FALSE, tuneGrid = data.frame(mtry = floor(length(rhs.vars)/3)), 
                   ntree=para$ntree)   
}
  result$pred <- predict(result$model, newdata = tst, type = "prob")
  result$perf <- Performance.measures(result$pred[, 2], tst[, resp.vars])
  result      
} 
 
)
return(res)
}


## train-test using longitudinal bootstrap 
#' @rdname  Train 
#' @export
#'
Train.Boot <- function(classifier, nboots = 2, XY.dat, resp.vars, rhs.vars, part.vars, reg.vars, 
                    para, parallel=TRUE, n.cores = 2, seed=12345678){
set.seed(seed) 
form.rules <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
form.glm <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(c(paste0(rhs.vars,collapse="+"), "+", "(", 
                paste0(c(rand.vars),collapse = "+"), "|", groups, ")"), collapse="")))

if(parallel) {
  cl <- makeCluster(n.cores)  
  pfun <-  get("parLapply")
} else {
  pfun = get("lapply")
}

names(XY.dat)[names(XY.dat)==para$group] <- "old.id"
tmp <- Multilevel.boot(dat=XY.dat, groups="old.id", nboots = nboots)
trn <- cbind.data.frame(tmp$trn, XY.dat[tmp$trn$trn.ix, ])
tst <- cbind.data.frame(tmp$tst, XY.dat[tmp$tst$tst.ix, ])

#if(parallel) clusterExport(cl, varlist=c("trn", "tst"))
#ix.boot <- unique(trn$Replicate)

res <- pfun(X= unique(trn$Replicate), function(kk, ...){    
    dat.trn <- trn[trn$Replicate==kk, ] 
    dat.tst <- tst[tst$Replicate==kk, ] 
    Y.trn <- dat.trn[, resp.vars]
    X.trn <- dat.trn[, unique(c(reg.vars, part.vars, para$group)), drop = FALSE]
    Y.tst <- dat.tst[, resp.vars]
    X.tst <- dat.tst[, unique(c(reg.vars, part.vars, para$group)), drop = FALSE] 
    
Train.Models <- lapply(TrainModels(), function(x) x)     
mod <- tryCatch(
{
lapply(classifier, function(x) {
  if(any(x%in%c("MEglmTree", "MECTree")))
    Train.Models[[x]](X.trn, Y.trn, X.tst, Y.tst, para, rand.vars, part.vars, reg.vars, rhs.vars, ...)
  else 
  if(x=="GBMrules") form <- form.rules else form <- form.glm  
  Train.Models[[x]](trn=dat.trn, tst=dat.tst, form, para, rand.vars, ...)
})
}, error=function(e){ 
  cat("Error in the Expression: ",  paste(e$call, collapse= ", "), 
      ": original error message = ", e$message, "\n") 
  list()
}) ## tryCatch
cat("Done boot :", kk,  "\n")
names(mod) <- classifier
mod
}, cl = cl)  ## pfun 

stopCluster(cl)
res
}

 ##################################################################
############# Tuning by CV concordance for GBM ###################
##################################################################

#' @rdname  Train 
#' @export
#
cv.gbm <- function(X, y, outcome=NULL, n.trees=seq(150,500,by=10), interaction.depth, n.minobsinnode = 5, 
shrinkage=10^(seq(-3,-1,length=10)), bag.fraction = .5, distribution="bernoulli", foldid=NULL, nfolds=10, seed=220, verbose=FALSE){

  n <- nrow(X)
  if(is.null(foldid))
    foldid <- get.foldid(n,10, seed=seed)
  nfolds <- max(foldid)
  ns <- length(shrinkage)
  nnt <- length(n.trees)
  concord.mat <- foreach(s=1:ns, .combine=rbind)%dopar%{
    concord.s <- rep(0,nnt)
    yhat <- matrix(0,n,nnt)
    for(k in 1:nfolds){
      X.k <- X[foldid==k,,drop=F]
      y.k <- y[foldid==k]
      X.mk <- X[foldid!=k,,drop=F]
      y.mk <- y[foldid!=k]
      if(distribution=="poisson"){
        outcome.mk <- outcome[foldid!=k]
        ans.sk <- gbm.fit(X.mk, outcome.mk, offset=y.mk, n.trees=n.trees[nnt], 
        interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, shrinkage=shrinkage[s], 
        bag.fraction=bag.fraction, distribution=distribution, verbose=verbose)
      }
      else{
        ans.sk <- gbm.fit(X.mk, y.mk, n.trees=n.trees[nnt], interaction.depth=interaction.depth, 
        n.minobsinnode=n.minobsinnode, shrinkage=shrinkage[s], bag.fraction=bag.fraction, distribution=distribution, verbose=verbose)
      }
      if(distribution=="multinomial"){
        for(t in 1:nnt){
	  foo <- predict(ans.sk, X.k, n.trees=n.trees[t],type="response", verbose=verbose)
	  fuh <- foo[cbind(1:length(y.k),as.numeric(y.k),1)]
	  fuh[is.na(fuh)] <- 0
          yhat[foldid==k,t] <- fuh
	}
      }
      else{
        for(t in 1:nnt)
          yhat[foldid==k,t] <- predict(ans.sk, X.k, n.trees=n.trees[t])
      }
    }
    for(t in 1:nnt){
      if(distribution=="coxph")
        concord.s[t] <- survConcordance(y~yhat[,t])$concord
      else if(distribution=="poisson")
        concord.s[t] <- survConcordance(Surv(y,outcome)~yhat[,t])$concord
      else if(distribution=="multinomial")
        concord.s[t] <- sum(log(yhat[,t]))
	#        concord.s[t] <- mean(y==yhat[,t])
      else if(distribution=="gaussian")
        concord.s[t] <- 1-sum((y-yhat[,t])^2)/sum((y-mean(y))^2)
      else
        concord.s[t] <- get.ROC(y, yhat[,t])$auc
    }
    concord.s
  }

  inds <- which(concord.mat==max(concord.mat), arr.ind=TRUE)
  opt.shrink <- shrinkage[inds[1]]
  opt.n.trees <- n.trees[inds[2]]
  opt.concord <- concord.mat[inds[1], inds[2]]
  return(list(shrinkage=opt.shrink, n.trees=opt.n.trees, concord=opt.concord))
}

get.foldid <- function(n, nfolds=10, seed=220){

  replace.seed <- T
  if(missing(seed))
    replace.seed <- F

  if(replace.seed){
   ## set seed to specified value
    if(!any(ls(name='.GlobalEnv', all.names=T)=='.Random.seed')){
      set.seed(1)
    }
    save.seed <- .Random.seed
    set.seed(seed)  
  }

  perm <- sample(1:n, n)
  n.cv <- rep(floor(n/nfolds),nfolds)
  rem <- n - n.cv[1]*nfolds
  n.cv[index(1,rem)] <- n.cv[index(1,rem)]+1
  foldid <- rep(0,n)
  ind2 <- 0
  
  for(i in 1:nfolds){
    ind1 <- ind2+1
    ind2 <- ind2+n.cv[i]
    foldid[perm[ind1:ind2]] <- i
  }
  if(replace.seed){
   ## restore random seed to previous value
    .Random.seed <<- save.seed
  }
  return(foldid)
}

index <- function(m,n){
  if(m<=n) return(m:n)
  else return(numeric(0))
}

which.equal <- function(x, y){

  n <- length(x)
  ans <- rep(0,n)
  for(i in 1:n){
    ans[i] <- any(approx.equal(y,x[i]))
  }
  return(as.logical(ans))
}

approx.equal <- function(x, y, tol=1E-9){

  return(abs(x - y) < tol)
}


   
    
    
    
    
    
    
    
    
    
    
    
    
    
