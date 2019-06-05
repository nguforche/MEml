#' @import caret kernlab 
#' 
## #' @importFrom stats predict
## #' @importFrom ksvm predict 
## #' @export
predictionFunctionV2 <- function(method, modelFit, newdata, preProc = NULL, param = NULL)
{
  if(!is.null(newdata) && !is.null(preProc)) newdata <- predict(preProc, newdata)
  out <- method$decision(modelFit = modelFit, 
                         newdata = newdata, 
                         submodels = param)
  ## TODO convert to character with classification
  out 
}

# get predicted decisions from ksvm train with caret. Caret only provide class levels and class probablilities, but 
# we need the decision functions from SVM for each data point. 
#' @export
predict.caret.ksvm <- function(object, newdata = NULL, type = "raw", na.action = na.omit, ...) {
  if(all(names(object) != "modelInfo")) {
    object <- update(object, param = NULL)
  }
#  require(kernlab)
  
  if(!(type %in% c("raw", "prob", "decision"))) stop("type must be either \"raw\" or \"prob\" or \"decision\"  ")
  if(type == "prob") {
    if (is.null(object$modelInfo$prob))
      stop("only classification models that produce probabilities are allowed")
  }
  
  if(!is.null(newdata)) {
    if (inherits(object, "train.formula")) {
      newdata <- as.data.frame(newdata)
      rn <- row.names(newdata)
      Terms <- delete.response(object$terms)
      m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, m)
      keep <- match(row.names(m), rn)
      newdata <- model.matrix(Terms, m, contrasts = object$contrasts)
      xint <- match("(Intercept)", colnames(newdata), nomatch = 0)
      if (xint > 0)
        newdata <- newdata[, -xint, drop = FALSE]
    }
  }
  else if(object$control$method != "oob") {
    if(!is.null(object$trainingData)) {
      if(object$method == "pam") {
        newdata <- object$finalModel$xData
      } else {
        newdata <- object$trainingData
        newdata$.outcome <- NULL
        if("train.formula" %in% class(object) &&
           any(unlist(lapply(newdata, is.factor)))) {
          newdata <- model.matrix(~., data = newdata)[,-1]
          newdata <- as.data.frame(newdata)
        }
      }
    } else stop("please specify data via newdata")
  }
  
  if("xNames" %in% names(object$finalModel) &
     is.null(object$preProcess$method$pca) &
     is.null(object$preProcess$method$ica))
    newdata <- newdata[, colnames(newdata) %in% object$finalModel$xNames, drop = FALSE]
  
  switch(type,
         prob = {
           out <- probFunction(method = object$modelInfo,
                               modelFit = object$finalModel,
                               newdata = newdata,
                               preProc = object$preProcess)
           obsLevels <- levels(object)
           out <- out[, obsLevels, drop = FALSE]
         }, 
         raw = {
           out <- predictionFunction(method = object$modelInfo,
                                     modelFit = object$finalModel,
                                     newdata = newdata,
                                     preProc = object$preProcess)
           if (object$modelType == "Regression") {
             out <- trimPredictions(pred = out,
                                    mod_type =object$modelType,
                                    bounds = object$control$predictionBounds,
                                    limits = object$yLimit)
           } else {
             if(!("levels" %in% names(object)))
               object$levels <- levels(object)
             out <- caret:::outcome_conversion(as.character(out), lv = object$levels)
           }
         }, 
         decision = {
           out <- predictionFunctionV2(method = object$modelInfo,
                                       modelFit = object$finalModel,
                                       newdata = newdata,
                                       preProc = object$preProcess)
         }
  )
  out
}


## Define a custom lsvm model in caret train  
## see https://topepo.github.io/caret/using-your-own-model-in-train.html for a full example 

#' @export
Methods <- function(...){
  res = list(
LaplaceSVM = list(
  type = c("Classification", "Regression"),
  
  library = "kernlab",
  
  parameters = data.frame(parameter = c("C", "sigma"), class = rep("numeric", 2), label = c("Cost", "Sigma")), 

  grid = function(x, y, len = NULL, search = "grid") {
        library(kernlab)
        ## This produces low, middle and high values for sigma 
        ## (i.e. a vector with 3 elements). 
        sigmas <- kernlab::sigest(as.matrix(x), na.action = na.omit, scaled = TRUE)  
        ## To use grid search:
        if(search == "grid") {
          out <- expand.grid(sigma = mean(as.vector(sigmas[-2])),
                             C = 2 ^((1:len) - 3))
        } else {
          ## For random search, define ranges for the parameters then
          ## generate random values for them
          rng <- extendrange(log(sigmas), f = .75)
          out <- data.frame(sigma = exp(runif(len, min = rng[1], max = rng[2])),
                            C = 2^runif(len, min = -5, max = 8))
        }
        out
      }, 
  
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) { 
         kernlab::ksvm(
        x = as.matrix(x), y = y,
        kernel = "rbfdot",
        kpar = list(sigma = param$sigma),
        C = param$C,
        prob.model = classProbs,
        ...
      )
  }, 
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) 
                     kernlab::predict(modelFit, newdata), 
  
  decision = function(modelFit, newdata, preProc = NULL, submodels = NULL)
              kernlab::predict(modelFit, newdata, type = "decision"), 
  
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL)
              kernlab::predict(modelFit, newdata, type = "probabilities"),  
  sort = function(x) x[order(x$C),], 
  
  levels = function(x) kernlab::lev(x), 
  loop = NULL)
)
return(res)
}























