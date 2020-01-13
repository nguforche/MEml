#' @title predict MECtree  
#' @description   
#' Make predictions for mixed effect conditional inference trees.
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




#' @title predict MEglmtree   
#' @description   
#' Make predictions for mixed effect model based trees.
#
#' @param object Fitted model from MEglmtree.
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
 




#' @title predict MEmixgbm   
#' @description   
#' Make predictions for mixed effect mixture of GBM.
#
#' @param object Fitted model from MEmixgbm.
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

 
 
 

#' @title predict MEmixgbm2   
#' @description   
#' Make predictions for mixed effect mixture of GBM (version 2).
#
#' @param object Fitted model from MEmixgbm2.
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


 
 
      

