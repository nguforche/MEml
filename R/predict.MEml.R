#' @title Predictions for MEml  
#' Make predictions using fitted MEml models: MEgbm, MEgbmrules, MErfrules, MEglmtree, MECTree, etc.
#
#' @name predict.MEml 
#
#' @param object Fitted model from MECtree.
#' @param newdata A new input data frame.
#' @param type of prediction: "prop" for probabilities and "class" for class labels.
#' @param ... Further arguments passed to or from other methods.
#' @return A list with items 
#' \item{prob}{predicted class probabilities}
#' \item{class}{predicted class memberships obtained by thresholding 
#' class probabilities at the prevalence rate of the positive class} 
#
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @import rpart partykit Matrix
#
NULL 
#
#' @name predict.MEml 
#' @export 
predict.MECTree <- function(object, newdata, 
                 type = c("prob", "class")[1], ...){  
### Get training group identifiers 
	ID <-  getME(object$glmer.fit, name = "flist")[[1]]

	if(sum(ID%in%newdata[, object$groups])==0) {
	Zb <- rep(0, nrow(newdata))	       	   
	} else {            
	re.form <- reOnly(formula(object$glmer.fit)) # RE formula only	
	res <- mkNewReTrms(object=object$glmer.fit, newdata = newdata, re.form = re.form, 
		      allow.new.levels=FALSE)
	b <- res$b 
	Zt <- res$Zt 
	Zb <- as.numeric(cprod(x = Zt, y = b))	
	}

    if(object$con.tree)
	treepred = predict(object$tree.fit, newdata = newdata, type = "response") + Zb 
	else 
	treepred = predict(object$tree.fit, newdata = newdata, type = "vector") + Zb 

	fitted  = treepred + Zb 
	pp <- sigmoid(fitted) 
	pred <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))
	if(type == "class"){     	
	pred <-ifelse(pred >= object$threshold, "OUT.CNTL", "IN.CNTL")
	} else if(type == "prob") {
	pred =  cbind(1-pred, pred)
	colnames(pred) = c("0", "1")
	} else stop("type unknown") 
return(pred)
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
#' @name predict.MEml
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


#' @name predict.MEml
#' @export 
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


#' @name predict.MEml       
#' @export 
predict.MEglmTree <- function(object, newdata, 
                 type = c("prob", "class")[1], ...){
	if (!inherits(object, "MEglmTree")) stop("Object must be a \"MEglmTree \"'")	
### Get training group identifiers 
	ID <-  getME(object$glmer.fit, name = "flist")[[1]]

### If new objects are not in the traing set the available options are  
### 1) set random effect to zero and predict fixed effect only
### 2) TODO use average cluster random effect:   
	if(sum(ID%in%newdata[, object$groups])==0) {
	Zb <- rep(0, nrow(newdata))	       	   
	} else {            
	re.form <- reOnly(formula(object$glmer.fit)) # RE formula only	
	res <- mkNewReTrms(object=object$glmer.fit, newdata = newdata, re.form = re.form, 
		      allow.new.levels=FALSE)
	b <- res$b 
	Zt <- res$Zt 
	Zb <- as.numeric(cprod(x = Zt, y = b))	
	}

 	if(object$include.RE)	
	newdata[,"R.Effect"] <- Zb
	treepred <- predict(object$tree.fit, newdata, type = "response")

	fitted  = treepred + Zb 
	pp <- sigmoid(fitted) 
	pred <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))
	if(type == "class"){     	
	pred <-ifelse(pred >= object$threshold, "OUT.CNTL", "IN.CNTL")
	} else if(type == "prob") {
	pred =  cbind(1-pred, pred)
	colnames(pred) = c("0", "1")
	} else stop("type unknown") 
return(pred)
}
 
#' @name predict.MEml     
#' @export 
predict.MEgbm <- function(object, newdata, 
                 type = c("prob", "raw")[1], ...){
	if (!inherits(object, "MEgbm")) stop("Object must be a \"MEgbm \"'")

#	if(!object$include.RE){
#	 newdata[, object$groups] <- NULL 
#	return(predict(object$gbmfit, newdata, type = type))
#	}
### Get training group identifiers 
	ID <-  getME(object$glmer.fit, name = "flist")[[1]]
    
### If new objects are not in the traing set the available options are  
### 1) set random effect to zero and predict fixed effect only
### 2) TODO use average cluster random effect:   
	if(sum(ID%in%newdata[, object$groups])==0) {
    Zb <- 0
	} else {            
	re.form <- reOnly(formula(object$glmer.fit)) # RE formula only	
	res <- mkNewReTrms(object=object$glmer.fit, newdata = newdata, re.form = re.form, 
	              allow.new.levels=FALSE)
  b <- res$b 
  Zt <- res$Zt        
  Zb <- as.numeric(cprod(x = Zt, y = b))    	   
	}	
if(object$include.RE)  newdata[,"Zb"] <- Zb

sigmoid(predict(object$gbmfit, newdata[, object$rhs.vars], type = "raw") + Zb) 
}

#' @name predict.MEml     
#' @export 
predict.MEgbmRules <- function(object, newdata, 
                          type = c("prob", "class")[1], ...){
  ### Get training group identifiers 
#  ID <-  getME(object$glmer.fit, name = "flist")[[1]]
  
  ### If new objects are not in the traing set the available options are  
  ### 1) set random effect to zero and predict fixed effect only
  ### 2) TODO use average cluster random effect:   
#  if(sum(ID%in%newdata[, object$groups])==0) {
#    newdata[,"Zb"] <- rep(0, nrow(newdata))	   
#     stop("No cases found in the data") 	   
#  } else {            
    re.form <- reOnly(formula(object$glmer.fit)) # RE formula only	
    res <- mkNewReTrms(object=object$glmer.fit, newdata = newdata, re.form = re.form, 
                       allow.new.levels=FALSE)
    b <- res$b 
    Zt <- res$Zt        
    Zb <- as.numeric(cprod(x = Zt, y = b))    
#    newdata[, object$groups] <- NULL 	   
#    newdata[,"Zb"] <- Zb
if(object$include.RE)  newdata[,"Zb"] <- Zb

    Target <- predict(object$gbmfit, newdata = newdata, type = "response", n.trees = object$para$n.trees)    
    fitted = Target + Zb   
    pp <- sigmoid(fitted) 
    prob <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))            
    predRules <- myapplyLearner(object$gbmRules, newdata[,object$rhs.vars])
#     newdata[, "TreeConditions"] <- factor(pred$predCond)    
#     levels(newdata[, "TreeConditions"]) <- object$rule.levels     
# #    pred$Values <- mapLabels(pred$predY, object$splitBins)
# #    pred.Rules <- do.call(cbind.data.frame, pred)  
#     pp <- predict(object$glmer.fit, newdata = newdata, type = "response")   
#     pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
#     w = pp*(1-pp)
#     Y <- ifelse(prob >= object$threshold, 1, 0) 
#     Y.star <- qlogis(pp) + (Y - pp)/w    
#     Target  <- Y.star  - Zb
#     target <- round(w*(Target+ Zb -  qlogis(pp) ) + pp) 
#     target <- ifelse(target==1, "Yes", "No")
#     pred.Rules$status <- target 
    
    
if(type == "class"){     	
  pred <-ifelse(prob >= object$threshold, "Yes", "No")
} else if(type == "prob") {
  pred =  cbind(1-prob, prob)
  colnames(pred) = c("0", "1")
} else stop("type unknown") 
		 
return(list(Y.star = fitted, pred=pred, predRules = predRules))
}

#' @name predict.MEml     
#' @export 
predict.MErfRules <- function(object, newdata, 
                               type = c("prob", "class")[1], ...){
  ### Get training group identifiers 
  #  ID <-  getME(object$glmer.fit, name = "flist")[[1]]
  
  ### If new objects are not in the traing set the available options are  
  ### 1) set random effect to zero and predict fixed effect only
  ### 2) TODO use average cluster random effect:   
  #  if(sum(ID%in%newdata[, object$groups])==0) {
  #    newdata[,"Zb"] <- rep(0, nrow(newdata))     
  #     stop("No cases found in the data") 	   
  #  } else { 
### predict rules 
#  pred <- myapplyLearner(object$learner, newdata[,object$rhs.vars])
#   predRules <- myapplyLearner(object$rfRules, newdata[,object$rhs.vars])
#   
#   print(predRules)
#   stop()
  
  re.form <- reOnly(formula(object$glmer.fit)) # RE formula only

  res <- mkNewReTrms(object=object$glmer.fit, newdata = newdata, re.form = re.form, 
                     allow.new.levels=FALSE)
  b <- res$b 
  Zt <- res$Zt        
  Zb <- as.numeric(cprod(x = Zt, y = b))    
  #    newdata[, object$groups] <- NULL 	   
  #    newdata[,"Zb"] <- Zb
  if(object$include.RE)  newdata[,"Zb"] <- Zb
    
  Target <- predict(object$rf.fit, newdata = newdata, type = "response")    
  fitted = Target + Zb   
  pp <- sigmoid(fitted) 
  prob <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))            
  predRules <- myapplyLearner(object$rfRules, newdata[,object$rhs.vars])
  #     newdata[, "TreeConditions"] <- factor(pred$predCond)    
  #     levels(newdata[, "TreeConditions"]) <- object$rule.levels     
  # #    pred$Values <- mapLabels(pred$predY, object$splitBins)
  # #    pred.Rules <- do.call(cbind.data.frame, pred)  
  #     pp <- predict(object$glmer.fit, newdata = newdata, type = "response")   
  #     pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
  #     w = pp*(1-pp)
  #     Y <- ifelse(prob >= object$threshold, 1, 0) 
  #     Y.star <- qlogis(pp) + (Y - pp)/w    
  #     Target  <- Y.star  - Zb
  #     target <- round(w*(Target+ Zb -  qlogis(pp) ) + pp) 
  #     target <- ifelse(target==1, "Yes", "No")
  #     pred.Rules$status <- target 
  
  
  if(type == "class"){     	
    pred <-ifelse(prob >= object$threshold, "Yes", "No")
  } else if(type == "prob") {
    pred =  cbind(1-prob, prob)
    colnames(pred) = c("0", "1")
  } else stop("type unknown") 
  
  return(list(Y.star = fitted, pred=pred, predRules = predRules))
}

 
 
#' @name predict.MEml     
#' @export 
predict.MEmixgbm <- function(object, newdata, 
                               type = c("prob", "class")[1], ...){

  re.form <- reOnly(formula(object$glmer.fit)) # RE formula only	
  res <- mkNewReTrms(object=object$glmer.fit, newdata = newdata, re.form = re.form, 
                     allow.new.levels=FALSE)
  b <- res$b 
  Zt <- res$Zt        
  Zb <- as.numeric(cprod(x = Zt, y = b))
  #    newdata[, object$groups] <- NULL 	   
  #    newdata[,"Zb"] <- Zb
  if(object$include.RE)  newdata[,"Zb"] <- Zb
  Target <- predict(object$gbmfit, newdata = newdata, type = "response", n.trees = object$para$n.trees)    
  fitted = Target + Zb 

  # predRules <- myapplyLearner(object$gbmRules, newdata[,object$rhs.vars])
  # newdata[, "TreeConditions"] <- factor(pred$predCond)    
  # levels(newdata[, "TreeConditions"]) <- object$rule.levels     
  # pp = predict(object$glmer.fit, newdata = newdata, type = "response")  
  # pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
  # w = pp*(1-pp)  
  # Y.star <- qlogis(pp) + (Y-pp)/w
#  dd <- rbind(newdata, newdata)
  # post <- posterior(object=object$mixfit, newdata = newdata, unscaled = FALSE)
  # Zb.mix <- apply(post*Zb, 1, sum) 
  # fitted = Target + Zb 
  # fitted = fitted[newdata$data==0]
  # fitted.mix = Target + Zb.mix
  # fitted.mix = fitted.mix[newdata$data==0]
  # cluster <- clusters(object=object$mixfit, newdata= newdata)
  # cluster = cluster[newdata$data==0]
  
  
  pp <- sigmoid(fitted) 
  prob <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))            
  predRules <- myapplyLearner(object$gbmRules, newdata[, object$rhs.vars])
  
  if(type == "class"){     	
    pred <-ifelse(prob >= object$threshold, "Yes", "No")
  } else if(type == "prob") {
    pred =  cbind(1-prob, prob)
    colnames(pred) = c("0", "1")
  } else stop("type unknown") 
  
  return(list(Y.star = fitted,  pred=pred, predRules = predRules))
}

 
 
 

#' @name predict.MEml     
#' @export 
predict.MEmixgbm2 <- function(object, newdata, 
                             type = c("prob", "class")[1], ...){
  
  re.form <- reOnly(formula(object$glmer.fit)) # RE formula only	
  res <- mkNewReTrms(object=object$glmer.fit, newdata = newdata, re.form = re.form, 
                     allow.new.levels=FALSE)
  b <- res$b 
  Zt <- res$Zt        
  Zb <- as.numeric(cprod(x = Zt, y = b))
  #    newdata[, object$groups] <- NULL 	   
  #    newdata[,"Zb"] <- Zb
  if(object$include.RE)  newdata[,"Zb"] <- Zb
  Target <- predict(object$gbmfit, newdata = newdata, type = "response", n.trees = object$para$n.trees)    
  fitted = Target + Zb 
  
  # predRules <- myapplyLearner(object$gbmRules, newdata[,object$rhs.vars])
  # newdata[, "TreeConditions"] <- factor(pred$predCond)    
  # levels(newdata[, "TreeConditions"]) <- object$rule.levels     
  # pp = predict(object$glmer.fit, newdata = newdata, type = "response")  
  # pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
  # w = pp*(1-pp)  
  # Y.star <- qlogis(pp) + (Y-pp)/w
  #  dd <- rbind(newdata, newdata)
  # post <- posterior(object=object$mixfit, newdata = newdata, unscaled = FALSE)
  # Zb.mix <- apply(post*Zb, 1, sum) 
  # fitted = Target + Zb 
  # fitted = fitted[newdata$data==0]
  # fitted.mix = Target + Zb.mix
  # fitted.mix = fitted.mix[newdata$data==0]
  # cluster <- clusters(object=object$mixfit, newdata= newdata)
  # cluster = cluster[newdata$data==0]
  
  
  pp <- sigmoid(fitted) 
  prob <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))            
  predRules <- myapplyLearner(object$gbmRules, newdata[, object$rhs.vars])
  
  if(type == "class"){     	
    pred <-ifelse(prob >= object$threshold, "Yes", "No")
  } else if(type == "prob") {
    pred =  cbind(1-prob, prob)
    colnames(pred) = c("0", "1")
  } else stop("type unknown") 
  
  return(list(Y.star = fitted,  pred=pred, predRules = predRules))
}


 
 
      

