# Utility functions: 

#' @title collect gabbage 
#
#' @description
#' Collects garbage until the memory is clean
#' @export
collect.garbage = function(){
  #if (exists("collect.garbage") ) rm(collect.garbage)
  ## The following function collects garbage until the memory is clean.
  ## Usage: 1. immediately call this function after you call a function or
  ##        2. rm()
  while (gc()[2,4] != gc()[2,4]){}
}
#
# returns t(x).y 
cprod <- function(x,y=NULL){
 if(is.complex(x) | is.complex(y)){
    if(is.null(y)){
     return(crossprod(Conj(x),x))
    } else {
     return(crossprod(Conj(x),y))
    }
 } else {
    return(crossprod(x,y))
 }
}
###  parse and manipulating mixed-model formulas
## deparse(.) returning \bold{one} string
## @param x,collapse character vector and how to combine them 
## @note Protects against the possibility that results from deparse() will be
##       split after 'width.cutoff' (by default 60, maximally 500)
safeDeparse <- function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)

#

# Random Effects formula only
reOnly <- function(f) {
    reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"))
}
#
##  @param x a random effect (i.e., data frame with rows equal to levels, columns equal to terms
##  @param n vector of new levels
levelfun <- function(x,nl.n,allow.new.levels=FALSE) {
    ## 1. find and deal with new levels
    if (!all(nl.n %in% rownames(x))) {
        if (!allow.new.levels) stop("new levels detected in newdata")
        ## create an all-zero data frame corresponding to the new set of levels ...
        newx <- as.data.frame(matrix(0, nrow=length(nl.n), ncol=ncol(x),
                                     dimnames=list(nl.n, names(x))))
        ## then paste in the matching RE values from the original fit/set of levels
        newx[rownames(x),] <- x
        x <- newx
    }
    ## 2. find and deal with missing old levels
    ## ... these should have been dropped when making the Z matrices
    ##     etc. in mkReTrms, so we'd better drop them here to match ...
    if (!all(r.inn <- rownames(x) %in% nl.n)) {
        x <- x[r.inn,,drop=FALSE]
    }
    return(x)
}

######################################################################################
#################  Train-Test split longitudinal data: train data lags test data by lag 
#################  there must be at least 2 visit times for each observation 
#' @title  Train-Test  longitudinal data split
#
#' @description
#' Split longitudinal data into training and test set, where the training set lags the 
#' test set by an amount "lag". e.g predict an outcome at lag=x hospital visits in advanced.  
#' There must be at least 2 repeated measurements for each unit 
#' @export
LongiLagSplit <- function(dat, id, rhs.vars, resp.vars, order.vars = "time", lag=1){ 

  tab <- data.frame(table(dat[, id]))
  idx <- tab$Var1[tab$Freq >= lag+3]
  dat <- unique(na.omit(dat[dat$id%in%idx, ]))
  dat$id <- factor(dat$id, levels = unique(dat$id))
  
  lag.split <- function(x){
    x <- unique(x) 
    x <- x[order(x[, order.vars], decreasing = FALSE), , drop=FALSE]
    ix1 <- 1:(nrow(x)-lag)
    ix2 <- ix1 + lag 
    dd <- cbind.data.frame(x[ix1, c(id, rhs.vars, resp.vars)], response = x[ix2, resp.vars]) 
    dd[,resp.vars]  <- factor(dd[, resp.vars])
    train <- head(dd, nrow(dd)-1)
    test <- tail(dd, 1) ### test using patients last observation 
    return(list(train = train, test = test))
  }
    
  X <- dlply(dat, .variables = id, .fun = lag.split)
  X[sapply(X, is.null)] <- NULL 
  
  train <- do.call(rbind.data.frame, lapply(X, function(x) x$train))
  test <- do.call(rbind.data.frame, lapply(X, function(x) x$test))
  
  list(train = train, test = test, resp.vars = "response", rhs.vars = c(rhs.vars, resp.vars))
}    

##### do not split into training and testing 
LongiLagSplit.boot <- function(dat, id, rhs.vars, resp.vars, order.vars = "time", lag=1){ 

  tab <- data.frame(table(dat[, id]))
  idx <- tab$Var1[tab$Freq >= lag+3]
  dat <- unique(na.omit(dat[dat$id%in%idx, ]))
  dat$id <- factor(dat$id, levels = unique(dat$id))
  
  lag.split <- function(x){
    x <- unique(x) 
    x <- x[order(x[, order.vars], decreasing = FALSE), , drop=FALSE]
    ix1 <- 1:(nrow(x)-lag)
    ix2 <- ix1 + lag 
    dd <- cbind.data.frame(x[ix1, c(id, rhs.vars, resp.vars)], response = x[ix2, resp.vars]) 
    dd[,resp.vars]  <- factor(dd[, resp.vars])
    dd   
  }
    
  X <- ddply(dat, .variables = id, .fun = lag.split)  
  list(dat = X, resp.vars = "response", rhs.vars = c(rhs.vars, resp.vars))
}    


#### Always predict the next or last visit, i.e use all repeated measure seen so far and predict the next possible one 
### i.e for visit =1 (next visit) data = (x0, y1), (x1, y2), .....
### visit = 2 data (x1, y2), (x2, y3), .... so y0, and y1 must be available before this will work 
LongiLagSplit.visit <- function(dat, id, rhs.vars, resp.vars, order.vars = "time", visit=1){ 
  
  tab <- data.frame(table(dat[, id]))
  idx <- tab$Var1[tab$Freq >= (visit+3)]
  dat <- unique(na.omit(dat[dat$id%in%idx, ]))
  dat$id <- factor(dat$id, levels = unique(dat$id))
  
  lag.split <- function(x){
    x <- unique(x)
    x <- x[order(x[, order.vars], decreasing = FALSE), , drop=FALSE]    
    ix1 <- 1:(nrow(x)-1)
    ix2 <- ix1 + 1    
    dd <- cbind.data.frame(x[ix1, c(id, rhs.vars, resp.vars)], response = x[ix2, resp.vars])  
    train <- head(dd, nrow(dd)-1)
    test <- tail(dd, 1) ### test using patients last observation 
    return(list(train = train, test = test))
  }
  
  X <- dlply(dat, .variables = id, .fun = lag.split)
  X[sapply(X, is.null)] <- NULL 
  
  train <- do.call(rbind.data.frame, lapply(X, function(x) x$train))
  test <- do.call(rbind.data.frame, lapply(X, function(x) x$test))
  
  list(train = train, test = test, resp.vars = "response", rhs.vars = c(rhs.vars, resp.vars))
}    


LongiLagSplit.visit.boot <- function(dat, id, rhs.vars, resp.vars, order.vars = "time", visit=1, parallel=TRUE){ 
  
  tab <- data.frame(table(dat[, id]))
  idx <- tab$Var1[tab$Freq >= (visit+3)]
  dat <- unique(na.omit(dat[dat$id%in%idx, ]))
  dat$id <- factor(dat$id, levels = unique(dat$id))
  
  lag.split <- function(x){
    x <- unique(x)
    x <- x[order(x[, order.vars], decreasing = FALSE), , drop=FALSE]    
    ix1 <- 1:(nrow(x)-1)
    ix2 <- ix1 + 1    
    cbind.data.frame(x[ix1, c(id, rhs.vars, resp.vars)], response = x[ix2, resp.vars])  
   } 
   X <- ddply(dat, .variables = id, .fun = lag.split, .parallel = parallel)     
  list(dat = X, resp.vars = "response", rhs.vars = c(rhs.vars, resp.vars))
}    

LongiLagSplit.current <- function(dat, id, rhs.vars, resp.vars, order.vars = "time"){ 
  
  tab <- data.frame(table(dat[, id]))
  idx <- tab$Var1[tab$Freq >= 3]
  dat <- unique(na.omit(dat[dat$id%in%idx, ]))
  dat$id <- factor(dat$id, levels = unique(dat$id))
  dat$past.HbA1c <- 0
    
  lag.split <- function(x){
    x <- x[order(x[, order.vars], decreasing = FALSE), , drop=FALSE] 
    ix1 <- 1:(nrow(x)-1)
    ix2 <- 2:nrow(x)
    x1 <- x[1, ,drop=FALSE]
    r1 <- x[ix1, resp.vars]
    r2 <- x[ix2, resp.vars]
    dd <- cbind.data.frame(visit = 2:nrow(x), x[ix2, c(id, rhs.vars)], past.HbA1c=r1, response=r2)
    x1 <- cbind.data.frame(visit = 1, x1[, c(id, rhs.vars)], past.HbA1c=x1$past.HbA1c, 
                           response= as.numeric(unlist(x1[, resp.vars])))
    rbind.data.frame(x1, dd)
  }  
X <- dlply(dat, .variables = id, .fun = lag.split)
X[sapply(X, is.null)] <- NULL 
X <- do.call(rbind.data.frame, X)
rhs.vars <- c("past.HbA1c", rhs.vars[rhs.vars!=resp.vars])
list(dd = X, resp.vars = "response", rhs.vars = rhs.vars)
}    



#' @title  Performance measures 
#
#' @description
#' Compute accuracy, AUC, sensitivity, specificity, positive predictive values, etc. 
#' @param pred  predicted probabilities 
#' @param obs  true class 
#' @param threshold class probability threshold (default to NULL)
#' @param prevalence class prevalence in the population (default to NULL)
# 
#' @export
Performance.measures <- function(pred, obs, threshold=NULL, prevalence=NULL){
if(length(unique(obs)) == 2) {
obs <- as.numeric(factor(obs))-1 
## get best cut-off 
if(is.null(threshold))
threshold <- opt.thresh(pred, obs)
### get these performance measures
nme = c("PCC", "PCC.sd", "AUC", "AUC.sd", "sensitivity", "sensitivity.sd", 
"specificity", "specificity.sd")
xx = cbind.data.frame(plotID = 1:length(pred), Observed = obs, Predicted = pred)
accuracy <- presence.absence.accuracy(xx, threshold = threshold, st.dev = TRUE)[, nme]
accuracy$G.mean <- sqrt(as.numeric(accuracy$sensitivity)*as.numeric(accuracy$specificity))
accuracy$BER <- 1 - 0.5*(as.numeric(accuracy$sensitivity) + as.numeric(accuracy$specificity)) 
pred.prev <- predicted.prevalence(DATA=xx,threshold=threshold)[, c("Obs.Prevalence", "Predicted")]


nme <- c("Pos Pred Value", "Neg Pred Value", "Balanced Accuracy")
if(is.null(prevalence)) prevalence <- as.numeric(pred.prev$Obs.Prevalence)
obs <- factor(ifelse(obs==1, "Yes", "No"), levels=c("Yes", "No"))
pred <- factor(ifelse(pred >= threshold, "Yes", "No"), levels = c("Yes","No"))
cmx <- confusionMatrix(data=pred, reference=obs,prevalence=prevalence)$byClass[nme]

#fmeas <- F.measure.single(pred= as.character(pred),labels= as.character(obs))
#accuracy$F.measure <- fmeas["F"]
#accuracy$Accuracy <- fmeas["A"] 
#colnames(cmx) <- c("PPV", "NPV", "BACC")
accuracy$PPV <- cmx[1]; accuracy$NPV = cmx[2]; accuracy$BACC <- cmx[3]
accuracy$F.measure  = 2 *(accuracy$PPV * accuracy$sensitivity) / (accuracy$PPV + accuracy$sensitivity); 

accuracy$threshold = threshold 
} else 
accuracy = data.frame(t(regr.eval(obs, pred, stats= c("mae","mse","rmse"))))  
return(accuracy)
}

sigmoid <- function(x) 1.0 / (1 + exp(-x))
# get left hand side of formula
lhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[response] 
}
# get right hand side of formula 
rhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[-response] 
}


#' @title Optimal threshold
#
#' @description
#' Compute the the optimal classification threshold
#'
#' @param prob Predicted probabilities by a classifier
#' @param obs ground truth (correct) 0-1 labels vector
#' @param opt.method  optima classifcation threshold method see package \code{PresenceAbsence}. Default 
#'  is the minRoc distance: i.e the threshold value at the minimum distance between the ROC curve and the 
#' to left hand corner (0,1)
#' @return threhold 
#' @examples
#' data(cars)
#' logreg <- glm(formula = vs ~ hp + wt,
#'               family = binomial(link = "logit"), data = mtcars)
#' prob <- logreg$fitted.values
#' opt.thresh(prob = pob, obs = mtcars$vs) 
#' @export 
opt.thresh <- function(pred, obs){
thresh = 0.5 
if(length(unique(obs)) > 1){
obs <- as.numeric(as.factor(obs))-1 
SIMDATA = cbind.data.frame(plotID = 1:length(obs), Observed = obs, Predicted = pred)
thresh <- optimal.thresholds(SIMDATA, threshold = 101, which.model = 1, opt.methods = 9)
thresh <- ifelse(length(thresh["Predicted"]) >= 1,as.numeric(thresh["Predicted"]), 0.5)
}
return(thresh)
}


myapplyLearner <- function (learner, X) 
{
  leftIx <- 1:nrow(X)
  predY <- rep("", nrow(X))
  predCond <- rep("", nrow(X))
  for (i in 1:nrow(learner)) {
    ixMatch <- eval(parse(text = paste("which(", learner[i, "condition"], ")")))
    ixMatch <- intersect(leftIx, ixMatch)
    if (length(ixMatch) > 0) {
      predY[ixMatch] <- learner[i, "pred"]
      predCond[ixMatch] <- learner[i, "condition"] 
      leftIx <- setdiff(leftIx, ixMatch)
    }
    if (length(leftIx) == 0) {
      break
    }
  }
  return(list(predY=predY,predCond=predCond))
}


SplitVector <- function (v, K = 3) 
{
  splitV <- quantile(v, probs = seq(0, 1, 1/K), na.rm = FALSE, 
                     names = TRUE, type = 3)
  splitV <- splitV[-c(1, length(splitV))]
  numSplit <- length(splitV)
  if (numSplit == 0) 
    return(v)
  newV <- vector("character", length(v))
  newV[which(v <= splitV[1])] = paste("L1", sep = "")
  if (numSplit >= 2) {
    for (jj in 2:numSplit) {
      newV[which(v > splitV[jj - 1] & v <= splitV[jj])] = paste("L", 
                                                                jj, sep = "")
    }
  }
  newV[which(v > splitV[numSplit])] = paste("L", (numSplit + 
                                                    1), sep = "")
  return(list(newV=newV, splitV = splitV))
}  


predictSplitVector <- function (v, splitV) 
{

  numSplit <- length(splitV)
  if (numSplit == 0) 
    return(v)
  newV <- vector("character", length(v))
  newV[which(v <= splitV[1])] = paste("L1", sep = "")
  if (numSplit >= 2) {
    for (jj in 2:numSplit) {
      newV[which(v > splitV[jj - 1] & v <= splitV[jj])] = paste("L", 
                                                                jj, sep = "")
    }
  }
  newV[which(v > splitV[numSplit])] = paste("L", (numSplit + 
                                                    1), sep = "")
  return(list(newV=newV))
}  



mapLabels <- function (pred, splitV) 
{
  newV <- rep(0, length(pred))
  labs <- pred
  xx <- unique(labs)
  ix <- 1:length(xx)
  is <- which(xx == "L1")
  newV[which(labs == "L1")] <- splitV[1]
  ll <- xx[xx != "L1"]
  names(splitV) <- ll
  
  if (length(splitV) >= 2) {
    for (jj in 1:(length(ll)-1)) {    
      newV[which(labs == ll[jj])] <- 0.5*(splitV[jj] + splitV[jj+1])
    }
  }
  newV[which(labs == ll[length(ll)])] <- splitV[length(ll)] + 0.5     
newV
}  


sparseFactor <- function (dummy, thresh = 0.00016) 
{
  nme <- FALSE     
  lunique = length(unique(dummy))
  
  if (lunique == 1) {
    nme <- TRUE
  }
  if (lunique == 2) {
    levs <- sort(unique(dummy))
    m1 = mean(dummy == levs[1], na.rm = TRUE)
    m0 = mean(dummy == levs[2], na.rm = TRUE)
    if (m1 <= thresh | m0 <= thresh) 
      nme <- TRUE
  }

  return(nme)
}


sortImp <- function (object, top) 
{
  
  best <- "max"
  featureRank <- switch(best, max = rank(-apply(object, 1, max, na.rm = TRUE)), 
                        min = rank(apply(object, 1, min, na.rm = TRUE)), 
                        maxabs = rank(-apply(abs(object), 1, max, na.rm = TRUE)))
  
  tiedRanks <- as.numeric(names(table(featureRank)[table(featureRank) > 1]))
  if (length(tiedRanks) > 0) {
    for (i in seq(along = tiedRanks)) {
      tmp <- featureRank[featureRank == tiedRanks[i]]
      featureRank[featureRank == tiedRanks[i]] <- tmp + 
        runif(length(tmp), min = 0.001, max = 0.999)
    }
  }
  featureOrder <- order(featureRank)
  out <- object[featureOrder, , drop = FALSE]
  out <- out[1:top, , drop = FALSE]
  out
}


varImpPlot.RF <- function (x, sort = TRUE, n.var = min(30, nrow(x$importance)), 
          type = NULL, class = NULL, scale = TRUE, main = "", xlab = "", ...) 
{
  if (!inherits(x, "randomForest")) 
    stop("This function only works for objects of class `randomForest'")
  imp <- importance(x, class = class, scale = scale, type = type, 
                    ...)
  if (ncol(imp) > 2) 
    imp <- imp[, -(1:(ncol(imp) - 2))]
  nmeas <- ncol(imp)
  if (nmeas > 1) {
    op <- par(mfrow = c(1, 2), mar = c(4, 5, 4, 1), mgp = c(2, 
                                                            0.8, 0), oma = c(0, 0, 2, 0), no.readonly = TRUE)
    on.exit(par(op))
  }
  for (i in 1:nmeas) {
    ord <- if (sort) 
      rev(order(imp[, i], decreasing = TRUE)[1:n.var])
    else 1:n.var
    xmin <- if (colnames(imp)[i] %in% c("IncNodePurity", "MeanDecreaseGini")) 
      0
    else min(imp[ord, i])
    dotchart(imp[ord, i], xlab = xlab, ylab = "", 
             main = if (nmeas == 1) 
               main
             else NULL, xlim = c(xmin, max(imp[, i])), color = "blue", ...)
  }
  if (nmeas > 1) 
    mtext(outer = TRUE, side = 3, text = main, cex = 1.2)
  invisible(imp)
}


#
Multilevel.boot <- function(dat, groups, nboots = 10) {
  cid <- unique(dat[, groups])
  ncid <- length(cid)    
  recid <- sample(cid, size = ncid * nboots, replace = TRUE)
  tst <- data.frame()
  trn <-  data.frame()
  for(ii in 1:length(recid)){
    idx <- which(dat[, groups] == recid[ii])
    trn.ix = sample(idx, size = length(idx), replace = TRUE)
    tst.ix =  setdiff(idx, trn.ix)
    if(length(tst.ix) <= 1) next 
    trn = rbind(trn, cbind(id = ii, trn.ix))
    tst = rbind(tst, cbind(id = ii, tst.ix))
  }
  #       rid <- lapply(seq_along(recid), function(i) {
  #           idx <- which(dat[, groups] == recid[i])
  ##           trn.ix = sample(idx, size = length(idx), replace = TRUE)
  #           tst.ix =  setdiff(idx, trn.ix)
  #           list(trn.ix = cbind(id = i, trn.ix), 
  ##           tst.ix = cbind(id = i, tst.ix))})        
  #
  #    trn <- as.data.frame(do.call(rbind, lapply(rid, function(x) x$trn.ix)))
  #    tst <- as.data.frame(do.call(rbind, lapply(rid, function(x) x$tst.ix)))
  trn$Replicate <- cut(trn$id, breaks = c(1, ncid * 1:nboots), include.lowest = TRUE,
                       labels = FALSE) 
  tst$Replicate <- cut(tst$id, breaks = c(1, ncid * 1:nboots), include.lowest = TRUE,
                       labels = FALSE)
  return(list(trn = trn, tst=tst))
}
#
multilevelboot <- function(dat, id = "id", inbag, outbag) {
### for each patient, select time indices that appeared in the inbag sample, sort these and keep aside the last 
### index 
res <- dlply(dat, .variables = id, .fun = function(xx){
p.ix <- xx[, "time.id"]
ix <- sort(inbag[inbag%in%p.ix]) 
last.ix <- ix[which(ix == tail(ix, 1))]  ### select the last index and remove from the inbag, this then goes in to test set 
inbag <- inbag[!inbag%in%last.ix] 
## get all outbags that are greater than last.ix - 1 
ix1 <- sort(inbag[inbag%in%p.ix]) 
last.ix1 <- unique(ix1[which(ix1 == tail(ix1, 1))])
out.ix <- outbag[which(outbag >= last.ix1)] 
test.ix <- unique(last.ix, out.ix)
train.ix <- inbag 
trn <- subset(xx, time.id%in%train.ix)
tst <- subset(xx, time.id%in%test.ix) 
list(trn = trn, tst=tst)
})

trn <- do.call(rbind, lapply(res, function(x) x$trn))
tst <- do.call(rbind, lapply(res, function(x) x$tst))
rownames(trn) <- NULL
rownames(tst) <- NULL
list(trn = trn, tst=tst)
} 


multilevelboot.list <- function(dat.list, inbag, outbag) {
### for each patient, select time indices that appeared in the inbag sample, sort these and keep aside the last 
### index 
res <- llply(dat.list, .fun = function(xx){

p.ix <- xx[, "time.id"]
ix <- sort(inbag[inbag%in%p.ix]) 
last.ix <- ix[which(ix == tail(ix, 1))]  ### select the last index and remove from the inbag, this then goes in to test set 
inbag <- inbag[!inbag%in%last.ix] 
## get all outbags that are greater than last.ix - 1 
ix1 <- sort(inbag[inbag%in%p.ix]) 
last.ix1 <- unique(ix1[which(ix1 == tail(ix1, 1))])
out.ix <- outbag[which(outbag >= last.ix1)] 
test.ix <- unique(last.ix, out.ix)
train.ix <- inbag 
trn <- subset(xx, time.id%in%train.ix)
tst <- subset(xx, time.id%in%test.ix)

list(trn = trn, tst=tst)
})

trn <- do.call(rbind, lapply(res, function(x) x$trn))
tst <- do.call(rbind, lapply(res, function(x) x$tst))
rownames(trn) <- NULL
rownames(tst) <- NULL
list(trn = trn, tst=tst)
} 


createBootSamples <- function (y, times = 10, replace = TRUE, prob = NULL) {        
    trainIndex <- matrix(0, ncol = times, nrow = length(y))
    out <- apply(trainIndex, 2, function(data) {
        index <- seq(along = data)
        out <- sort(sample(index, size = length(index), replace = replace,  prob = prob))
#        out <- sample(index, size = length(index), replace = replace,  prob = prob)
        out
    })
     out
}


multiClassSummary <- function (data, class.names=NULL, eps = 1e-15){  
  #Check data
#  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
#    stop("levels of observed and predicted data do not match")
  
  #Calculate custom one-vs-all stats for each class
  prob_stats <- lapply(class.names, function(class){
    
    #Grab one-vs-all data for the class
    pred <- ifelse(data[, "pred"] == class, 1, 0)
    obs  <- ifelse(data[,  "obs"] == class, 1, 0)
    prob <- data[,class]
    
    #Calculate one-vs-all AUC and logLoss and return
#    cap_prob <- pmin(pmax(prob, .000001), .999999)
    cap_prob <- pmin(pmax(prob, eps), 1-eps)
#    prob_stats <- auc(obs, prob)
#    names(prob_stats) <- 'ROC'
prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
names(prob_stats) <- c('ROC', 'logLoss')
    return(prob_stats) 
  })
  prob_stats <- do.call(rbind, prob_stats)
  rownames(prob_stats) <- paste('Class:', class.names)
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"])
  
  #Aggregate and average class-wise stats
  #Todo: add weights
  class_stats <- cbind(CM$byClass, prob_stats)
  class_stats <- colMeans(class_stats, na.rm=TRUE)
  
  #Aggregate overall stats
  overall_stats <- c(CM$overall)
  
  #Combine overall with class-wise stats and remove some stats we don't want 
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', 
                                       'Prevalence', 'Detection Prevalence')]
  
  #Clean names and return
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  return(stats)
  
}


#' @title normalize data
#
#' @description
#' Normalize matrix or data frame to the min-max range of the variables.
#
#' @param x  matrix or data frame 
#' @return nomalized matrix/data frame with min and max of each variable as attributes 
#' @export
# normalize data to 

normalize <- function(x) { 
  x <- as.matrix(x)
  minAttr=apply(x, 2, min, na.rm=TRUE)
  maxAttr=apply(x, 2, max, na.rm=TRUE)
  x <- sweep(x, 2, minAttr, FUN="-") 
  x=sweep(x, 2,  maxAttr-minAttr, "/") 
  attr(x, 'min') = minAttr
  attr(x, 'max') = maxAttr
  return (x)
} 

normalizeMinusOneToOne <- function(x, a = -1, b = 1) {   
  A= min(x, na.rm=TRUE)
  B <-max(x, na.rm = TRUE)
  (((x - A)*(b-a))/(B-A) + a)
} 

#' @title denomalize data 
#
#' @description 
#' tranform normalized data back to original scale  
#
#' @param  normalized  the normalized data frame/matrix 
#' @param  min,max the min and max of each variable in the data. This can be otained by retrieving the attribute values of 
#' the normalized data, i.e. output of \code{normalize}
#' @export
# denormalize back to original scale

# denormalize back to original scale
denormalize <- function (normalized, min, max) {
  if(dim(normalized)[2]!= length(min)) stop("length of min or max must equal number of columns of data ")
  nme <- colnames(normalized)
  if( !all(nme%in%names(min)) ) stop("colnames of data do not match names of min or max")
  sapply(nme, function(x)   normalized[, x] * (max[x]-min[x]) + min[x] )
}




############################################################
######## 2. F-score ########################################
############################################################
# Function to compute the F-measure for a single class
# Input
# pred : vector of the predicted labels
# labels : vector of the true labels
# Note that 0  stands for negative and 1 for positive.
# In general the first level is negative and the second positive
# Output:
# a named vector with six elements:
# P: precision
# R : recall (sensitivity)
# S : specificity
# F : F measure
# A : 0/1 loss accuracy
# npos : number of positives
F.measure.single <- function(pred,labels) {      
   if (length(pred)!=length(labels))
       stop("F.measure: lengths of true and predicted labels do not match.");
   neg.ex <- which(labels <= 0);	
   np <- length(neg.ex);
   pos.ex <- which(labels > 0);
   npos <- length(pos.ex);	 
   TP <- sum(pred[pos.ex] > 0);
   FN <- sum(pred[pos.ex] <= 0);	
   TN <- sum(pred[neg.ex] <= 0);
   FP <- sum(pred[neg.ex] > 0);	           
   acc <- (TP+TN)/length(labels);
   print(TP)
   print(TN)
   
   if ((TP+FP) == 0)
     precision <- 0
   else 
     precision <- TP/(TP+FP);
   if ((TP+FN) == 0)
     recall <- 0
   else
     recall <- TP/(TP+FN);
   if ((TN+FP) == 0)
     specificity <- 0
   else
     specificity <- TN/(TN+FP);
   if ((precision+recall) == 0)
      F <- 0
   else
      F = 2 *(precision*recall) / (precision+recall); 
   res <- c(precision,recall,specificity,F,acc, npos);
   names(res) <- c("P", "R", "S", "F", "A", "Pos.");
   return (res);
}



################################################################################
################################################################################
#################### predict method for mixed effect model based tree 
################ create a new mixed effect model frame 
## Make new random effect terms from specified object and new data,
## possibly omitting some random effect terms
## @param object fitted model object
## @param newdata (optional) data frame containing new data
## @param re.form,allow.new.levels formula specifying random effect terms to include (NULL=all, ~0)
#
## @note Hidden; _only_ used (twice) in this file
mkNewReTrms <- function(object, newdata, re.form = NULL,  allow.new.levels=FALSE)
{
  ## construct (fixed) model frame in order to find out whether there are
  ## missing data/what to do about them
  ## need rfd to inherit appropriate na.action; need grouping
  ## variables as well as any covariates that are included
  ## in RE terms
  if (is.null(newdata)) 
    rfd <- model.frame(object)
  else 
    rfd <- newdata
  
  ## note: mkReTrms automatically *drops* unused levels
  ReTrms <- mkReTrms(findbars(re.form[[2]]), rfd)
  ## update Lambdat (ugh, better way to do this?)
  ReTrms <- within(ReTrms,Lambdat@x <- unname(getME(object,"theta")[Lind]))
  if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
    stop("NAs are not allowed in prediction data",
         " for grouping variables unless allow.new.levels is TRUE")
  ns.re <- names(re <- ranef(object))
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re))
    stop("grouping factors specified in re.form that were not present in original model")
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  ## fill in/delete levels as appropriate
  re_x <- Map(function(r,n) levelfun(r,n,allow.new.levels=allow.new.levels),
              re[names(new_levels)], new_levels)
  ## pick out random effects values that correspond to
  ##  random effects incorporated in re.form ...
  ## NB: Need integer indexing, as nRnms can be duplicated: (age|Subj) + (sex|Subj) :
  re_new <- lapply(seq_along(nRnms), function(i) {
    rname <- nRnms[i]
    if (!all(Rcnms[[i]] %in% names(re[[rname]])))
      stop("random effects specified in re.form that were not present in original model")
    re_x[[rname]][,Rcnms[[i]]]
  })
  
  re_new <- unlist(lapply(re_new, t))  ## must TRANSPOSE RE matrices before unlisting
  
  Zt <- ReTrms$Zt
  list(Zt=Zt, b=re_new, Lambdat = ReTrms$Lambdat)
}



mkNewReTrmsV2 <- function(object, newdata, re.form=NULL, na.action=na.pass,
                          allow.new.levels=FALSE)
{
  ## construct (fixed) model frame in order to find out whether there are
  ## missing data/what to do about them
  ## need rfd to inherit appropriate na.action; need grouping
  ## variables as well as any covariates that are included
  ## in RE terms
  if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object)
  } else {
    mfnew <- model.frame(delete.response(terms(object,fixed.only=TRUE)),
                         newdata, na.action=na.action)
    ## make sure we pass na.action with new data
    ## it would be nice to do something more principled like
    ## rfd <- model.frame(~.,newdata,na.action=na.action)
    ## but this adds complexities (stored terms, formula, etc.)
    ## that mess things up later on ...
    rfd <- na.action(newdata)
    if (is.null(attr(rfd,"na.action")))
      attr(rfd,"na.action") <- na.action
  }
  if (inherits(re.form, "formula")) {
    ## DROP values with NAs in fixed effects
    if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
      newdata <- newdata[-fit.na.action,]
    }
    ## note: mkReTrms automatically *drops* unused levels
    ReTrms <- mkReTrms(findbars(re.form[[2]]), rfd)
    ## update Lambdat (ugh, better way to do this?)
    ReTrms <- within(ReTrms,Lambdat@x <- unname(getME(object,"theta")[Lind]))
    if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
      stop("NAs are not allowed in prediction data",
           " for grouping variables unless allow.new.levels is TRUE")
    ns.re <- names(re <- ranef(object))
    nRnms <- names(Rcnms <- ReTrms$cnms)
    if (!all(nRnms %in% ns.re))
      stop("grouping factors specified in re.form that were not present in original model")
    new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
    ## fill in/delete levels as appropriate
    re_x <- Map(function(r,n) levelfun(r,n,allow.new.levels=allow.new.levels),
                re[names(new_levels)], new_levels)
    ## pick out random effects values that correspond to
    ##  random effects incorporated in re.form ...
    ## NB: Need integer indexing, as nRnms can be duplicated: (age|Subj) + (sex|Subj) :
    re_new <- lapply(seq_along(nRnms), function(i) {
      rname <- nRnms[i]
      if (!all(Rcnms[[i]] %in% names(re[[rname]])))
        stop("random effects specified in re.form that were not present in original model")
      re_x[[rname]][,Rcnms[[i]]]
    })
    
    re_new <- unlist(lapply(re_new, t))  ## must TRANSPOSE RE matrices before unlisting
    ## FIXME? use vapply(re_new, t, FUN_VALUE=????)
  }
  Zt <- ReTrms$Zt
  attr(Zt, "na.action") <- attr(re_new, "na.action") <- attr(mfnew, "na.action")
  list(Zt=Zt, b=re_new, Lambdat = ReTrms$Lambdat)
}



