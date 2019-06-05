# 
#' Boosted generalized mixed-effect regression models (glmer).  
#  
#' @name glmerBoost
#
#' @param form A formula describing fixed and  random effects (as in lme4)  
#' @param dat Matrix/data frame 
#' @param lme.family family for glmer 
#' @param verbose print iteration outputs ? 
#' @param max.iter maximum iteration 
#' @param coeflearn Adaboot.M1 learning method 
#' @param \dots other arguments 
#' @return a list with items 
#' \item{tree.fit}{fitted lme4 models}
#' \item{fitted.probs}{fitted probabilites for final model}
#' \item{fitted.class}{fitted class labels for final model}
#' \item{train.perf}{various performance measures for final model on training set}
#' \item{threshold}{classification cut-off}
#
#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#' @import lme4 partykit 
NULL 
#
#' @rdname glmerBoost  
#' @export
glmerBoost  <- function(form, ...) UseMethod("glmerBoost")
#
#' @rdname glmerBoost
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mod <- glmerBoost(form, data)) 
#
#' }
#
glmerBoost <- function(form,  dat, lme.family = binomial, max.iter =100, verbose = FALSE, 
                      coeflearn = "Breiman", ...){
    tmp.w <- NULL               
    y <- dat[, lhs.form(form)]
    n <- length(y)
    nclasses <- length(unique(y))
    trees <- list()
    thresh <- alpha <-  rep(0, max.iter)
    w <- rep(1/n, n)
    f.hat <-  0
     

            
     for (m in 1:max.iter) {          
       e <- new.env() 
#       e$tmp.w <- w 
       tmp.w <<- w
        fit <- glmer(form, data= dat,family = lme.family, 
	              control = glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), weights = tmp.w,  
		          nAGQ= 0, verbose = as.numeric(verbose))
	    
	pp = predict(fit, newdata = dat, type = "response")  
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
        thresh[m] <- opt.thresh(pp, y)        
        flearn <- ifelse(pp >= thresh[m], 1, 0)
        ind <- as.numeric(y != flearn)
 #       err <- Performance.measures(pp, y, thresh[m])$BER
       err <- as.numeric(w %*% ind)                            
        c <- log((1 - err)/err)
        
        if (coeflearn == "Breiman")  c <- (1/2) * c        
        if (coeflearn == "Zhu")  c <- c + log(nclasses - 1)
        
        update.vector <- w * exp(c * ind)
        w[ind == 1] <- update.vector[ind == 1]
        w <- w/sum(w)
        maxerror <- 0.5
        eac <- 0.001
        if (coeflearn == "Zhu")  maxerror <- 1 - 1/nclasses
        
        if (err >= maxerror) {
            maxerror <- maxerror - eac
            c <- log((1 - maxerror)/maxerror)
            if (coeflearn == "Breiman") {
                c <- (1/2) * c
            }
            if (coeflearn == "Zhu") {
                c <- c + log(nclasses - 1)
            }
        }
        if (err == 0) {
            c <- log((1 - eac)/eac)
            if (coeflearn == "Breiman") c <- (1/2) * c
            
            if (coeflearn == "Zhu") c <- c + log(nclasses - 1)
           
        }
        
        y.hat <- pp - thresh[m]         
        f.hat <- f.hat + c*y.hat
        alpha[m] <- c           
        trees[[m]] <- fit

        
    }
	H <- f.hat/sum(alpha)
	pp <- sigmoid(H)
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))
	perf <- Performance.measures(pp, y)
	threshold <- perf$threshold	
  cls <-ifelse(pp >= threshold, "OUT.CNTL", "IN.CNTL") 
	

## convert H to interval [0, 1]
	A <- min(H);
	B <- max(H)
	a = 0; b = 1 
	pp2 <- (((H - A)*(b-a))/(B-A) + a)
	perf2 <- Performance.measures(pp2, y)
#rm(tmp.w, envir = e)	      
res <- list(tree.fit = trees, fitted.probs = pp, fitted.class = cls, train.perf = perf, 
train.perf2 = perf2, global.threshold = threshold, alpha = alpha, thresholds = thresh)
class(res) <- "glmerBoost"         
return(res)
}

#
#' @rdname glmerBoost  
#' @export
glmerLogitBoost  <- function(form, ...) UseMethod("glmerLogitBoost")
#
#' @rdname glmerBoost
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mod <- glmerBoost(form, rand.form, data)) 
#
#' }
#
glmerLogitBoost <- function(form,  dat, max.iter =100, verbose = FALSE, ...){
    tmp.w <- NULL                               
    y <- dat[, lhs.form(form)]
    n <- length(y)
    trees <- list()
    thresh <-  rep(0, max.iter)
    f.hat <-  0
    pp <- rep(1/2, n)
    w <- pp*(1-pp)
    z <- (y-pp)/w 
    f.hat <- 0 
    dat[, "z"] <- z 
    
    rhs.var <- as.character(rhs.form(form))
    nn <- length(rhs.var)
    resp.vars <- "z" 
    bar <- tail(rhs.var)[6]
    
    form <- as.formula(paste0(paste0(resp.vars, " ~ "), 
    paste0(c(paste0(head(rhs.var, nn-1),collapse = "+") , "+", paste0(paste0("(", bar),")") ), collapse = "") ))  
           
     for (m in 1:max.iter) {   
       tmp.w <<- w   
#       e <- new.env()     
 #      e$tmp.w <- w 
        fit <- lmer(form, data= dat,control = lmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
                    weights = tmp.w, verbose = as.numeric(verbose))
	    
	fm = predict(fit, newdata = dat, type = "response", allow.new.levels=TRUE)  
	f.hat <- f.hat + 0.5*fm 
	pp <- sigmoid(2.0*f.hat)
	
        w <- pp*(1-pp)
        z <- (y-pp)/w 
        dat[, "z"] <- z 
        trees[[m]] <- fit  	
       }
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))
	perf <- Performance.measures(pp, y)
	threshold <- perf$threshold	
  cls <-ifelse(pp >= threshold, "OUT.CNTL", "IN.CNTL")
#rm(tmp.w, envir = e)		
## convert H to interval [0, 1]
#	A <- min(f.hat);
#	B <- max(f.hat)
#	a = 0; b = 1 
#	pp <- (((f.hat - A)*(b-a))/(B-A) + a)
#	perf2 <- Performance.measures(pp, y)
	      
res <- list(tree.fit = trees, fitted.probs = pp, fitted.class = cls, train.perf = perf, threshold = threshold)
class(res) <- "glmerLogitBoost"         
return(res)
}











