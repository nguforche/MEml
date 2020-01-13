#' @export
MEglm  <- function(X, ...) UseMethod("MEglm")


#' @title Mixed Effect GLM 
#' @description 
#' Trains a Mixed Effect generalized mixed effect models. This function just calls glmer from the Lme4 package. 

#
#' @param form formula  
#' @param data  data.frame with predictors 
#' @param Control glmer control
#' @param nAGQ  see \code{glmer}

#' @export
MEglm <- function(form, data,  control, nAGQ=0 ) {
  res <- glmer(form, data= data,family=binomial,control= control,nAGQ=nAGQ)
  res <- list(model = "MEglm", MEglm= res)
  class(res) <- "MEglm"
  res
  } 
  
