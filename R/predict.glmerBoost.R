#' Make predictions using a fitted MECtree model
#
#' @name predict.glmerBoost
#
#' @param object Fitted model from MECtree.
#' @param newdata A new input data frame.
#' @param type of prediction: "prop" for probabilities and "class" for 
#' class labels.
#' @param ... Further arguments passed to or from other methods.
#' @return A list with items 
#' \item{prob}{predicted class probabilities}
#' \item{class}{predicted class memberships obtained by thresholding 
#' class probabilities at the prevalence rate of the positive class} 
#
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#
#' @name predict.glmerBoost
#' @export 
predict.glmerBoost <- function(object, newdata, 
                 type = c("prob", "class")[1], ...){  

     f.hat <- 0 
     for (m in 1:length(object$tree.fit)) {
     
       pp = predict(object$tree.fit[[m]], newdata = newdata, type = "response")  
       pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
       y.hat <- pp - object$thresh[m]        
       f.hat <- f.hat + object$alpha[m]*y.hat
     }   
	H <- f.hat/sum(object$alpha)
	pp <- sigmoid(H)
	pred <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))
	if(type == "class"){     	
	pred <-ifelse(pred >= object$global.threshold, 1, 0)
	} else if(type == "prob") {
	pred =  cbind(1-pred, pred)
	colnames(pred) = c("0", "1")
	} else stop("type unknown") 
return(pred)
}
#	
#' @name predict.glmerBoost
#' @export 
predict.glmerLogitBoost <- function(object, newdata, 
                 type = c("prob", "class")[1], ...){  

     f.hat <- 0 
     for (m in 1:length(object$tree.fit)) {     
       fm = predict(object$tree.fit[[m]], newdata = newdata, type = "response")  
       f.hat <- f.hat + 0.5*fm     
     }   
	pred <- sigmoid(2.0*f.hat)    
	if(type == "class"){     	
	pred <-ifelse(pred >= object$threshold, 1, 0)
	} else if(type == "prob") {
	pred =  cbind(1-pred, pred)
	colnames(pred) = c("0", "1")
	} else stop("type unknown") 
return(pred)
}
	
 
 
 
 
 
 
 
 
      

