#' @title Predict MEgbm   
#' Make predictions for Mixed Effect GBM.
#
#' @param object Fitted model.
#' @param newdata A new input data frame.
#' @param type of prediction: "prop" for probabilities and "class" for class labels.
#' @param allow.new.levels specifry if new levels of factor variables in the test set should be allowed. Default to FALSE.
#' @return A list with items 
#' \item{Y.star}{predicted transformed outcome}
#' \item{pred}{matrix with predicted class probabilites 
#' \item{predRules}{predicted rules}
#
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @import rpart partykit Matrix

#' @export 
predict.MEgbmRules <- function(object, newdata, 
                               type = c("prob", "class")[1], allow.new.levels = FALSE, ...){
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
                     allow.new.levels=allow.new.levels)
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

#' @title Predict MErf   
#' Make predictions for Mixed Effect random forest.
#
#' @param object Fitted model.
#' @param newdata A new input data frame.
#' @param type of prediction: "prop" for probabilities and "class" for class labels.
#' @param allow.new.levels specifry if new levels of factor variables in the test set should be allowed. Default to FALSE.
#' @return A list with items 
#' \item{Y.star}{predicted transformed outcome}
#' \item{pred}{matrix with predicted class probabilites 
#' \item{predRules}{predicted rules}
#
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @import rpart partykit Matrix

#' @export 
predict.MErfRules <- function(object, newdata, 
                              type = c("prob", "class")[1], allow.new.levels=FALSE, ...){
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
                     allow.new.levels=allow.new.levels)
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

#' @title Predict MEglm  
#' Make predictions for Mixed Effect logistic regression. This is just a wrapper for \code{predict.merMod} in lme4 .
#
#' @param object Fitted model.
#' @param newdata A new input data frame.
#' @param type of prediction: "prop" for probabilities and "class" for class labels.
#' @param allow.new.levels specifry if new levels of factor variables in the test set should be allowed. Default to FALSE.
#' @return A list with items 
#' \item{Y.star}{NULL}
#' \item{pred}{matrix with predicted class probabilites 
#' \item{predRules}{NULL}
#
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @import rpart partykit Matrix

#' @export 
predict.MEglm<- function(object, newdata, 
                              type = c("prob", "class")[1], allow.new.levels=FALSE, ...){
#class(object$MEglm) <- "merMod"
pred = predict(object$MEglm, newdata= newdata, type ="response", allow.new.levels = allow.new.levels)
pred = cbind(1-pred, pred)
colnames(pred) = c("0", "1")
list(Y.star = NULL,  pred = pred, predRules = NULL)  
}


  

