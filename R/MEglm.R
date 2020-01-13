#' @title Mixed Effect GLM 
#' @description 
#' Trains a Mixed Effect generalized mixed effect models. This function just calls glmer from the Lme4 package. 

#' @name MEglm 
#
#' @param form formula  
#' @param data  data.frame with predictors 
#' @param Control glmer control
#'@param nAGQ  see \code{glmer}
#
NULL 
#
#' @rdname MEglm   
#' @export
MEglm  <- function(X, ...) UseMethod("MEglm")
#
#' @rdname MEglm 
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mod <- MEglm(form, data= dat.trn, family=binomial) 
#'
#' }

MEglm <- function(form, data,  control, nAGQ=0 ) {
  res <- glmer(form, data= data,family=binomial,control= control,nAGQ=nAGQ)
  res <- list(model = "MEglm", MEglm= res)
  class(res) <- "MEglm"
  res
  } 
  
