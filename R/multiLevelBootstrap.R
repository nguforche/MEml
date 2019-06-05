#' @export
multiLevelBootstrap <- function(lag=NULL, visit = 1, classifier, dat, id, rhs.vars, resp.vars, order.vars, rand.vars,  reg.vars, 
                            groups = "id", para, nBoots = 100, samp.prob=TRUE, max.iter = 10, seed = 1,  return.model = TRUE, ...){
#  if(parallel) {
#     pfun <-  get("mclapply")
# set.seed(seed, kind = "L'Ecuyer-CMRG")
#   } else {
#     pfun = get("lapply")
# set.seed(seed) 
#   }
# 
########  
  split <- NULL  
  model <- NULL
  if(!is.null(lag)){
  split <- lag   
  dd <- LongiLagSplitV2.boot(dat, id, rhs.vars,resp.vars, order.vars,lag=lag)  
  } else {
  split <- visit 
  dd <- LongiLagSplit.visit.boot(dat, id=id, rhs.vars=rhs.vars,resp.vars=resp.vars,order.vars=order.vars,visit=visit)  
  }
 
 XY.dat <- dd$dat 
 resp.vars <- dd$resp.vars
 rhs.vars <-dd$rhs.vars
 names(XY.dat)[names(XY.dat)%in%c(rhs.vars, resp.vars)] <- make.names(c(rhs.vars, resp.vars))
 rhs.vars <-   make.names(rhs.vars)
 resp.vars <- make.names(resp.vars)
 part.vars <- rhs.vars
 
#### get bootstrap data 
## assign time indices to each patient repeated observations 
dat <- ddply(XY.dat, .variables = id, .fun = function(x) {x <- x[order(x[, order.vars]), ]; x$time.id <- 1:nrow(x); x})  
tb <- table(dat$time.id)  ## number of occurences of each time-indexed observation 
ix.boot <- as.numeric(rownames(tb)) ## time-indices that will be bootstrapped. we can also do block bootstrap on these. 

dat.list <- dlply(dat, .variables = id, .fun = function(xx) xx)
dat.list[sapply(dat.list, is.null)] <- NULL 

### bootstrap sample. repeated observations are sampled according to the the proportion of time-index occurence in the data 
if(samp.prob) p <- (tb/sum(tb)+ 0.001) 
else p <- NULL ## sampling probabilities 

BOOT <- createBootSamples(y=ix.boot, times=nBoots, replace = TRUE, prob = p)

Trn.Tst <- lapply(Train.Test(), function(x) x)  

##### do parallel 
Boot.res <- lapply(1:nBoots, function(kk){ 
inbag <-  BOOT[, kk]
outbag <- setdiff(ix.boot, inbag)

#d0 <- multilevelboot(dat= dat,id = id, inbag = inbag, outbag = outbag)
d0 <- multilevelboot.list(dat.list=dat.list, inbag=inbag, outbag=outbag)

trn <- d0$trn 
tst <- d0$tst
mod.res <- tryCatch(
{
lapply(classifier, function(xx) {
##cat("Now Running Classifier:", xx,  "\n")
if(return.model)
Trn.Tst[[xx]](trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars, part.vars, groups)
else 
cbind(Classifier = xx, Lag = split, Trn.Tst[[xx]](trn, tst, para, resp.vars, rand.vars, rhs.vars, reg.vars, part.vars, groups)$perf) 
})
}, error=function(e){ 
  cat("Error in the Expression: ",  paste(e$call, collapse= ", "), 
      ": original error message = ", e$message, "\n") 
  list()
}) ## tryCatch
collect.garbage()
#names(mod.res) <- classifier

cat("Done Bootstrap:", kk,  "\n")
return(mod.res)
}) ## pfun 


return(Boot.res)
}



