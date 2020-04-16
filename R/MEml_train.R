Train_MEml <- function(...){
  res = list(
    MEglm = function(trn, para, resp.vars, rhs.vars, rand.vars, groups, ...){
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(c(paste0(rhs.vars,collapse="+"), "+", "(", 
                                                                  paste0(c(rand.vars),collapse = "+"), "|", groups, ")"), collapse="")))
      MEglm(form = form, data= trn,control= para$glmer.Control,nAGQ=para$nAGQ)
      
     }, 
    
    MEgbm = function(trn, para, resp.vars, rand.vars, rhs.vars, groups, ...){  
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 
      MEgbm(form = form, dat=trn,  groups = groups, family = para$family, rand.vars= rand.vars, para = para,   
                                 tol= para$tol, max.iter = para$max.iter, include.RE = para$include.RE, 
                                 verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=para$glmer.Control,
                                 nAGQ = para$nAGQ, K = para$K, decay = para$decay, likelihoodCheck = para$likelihoodCheck) 
    }, 
    

    MErf = function(trn, para, resp.vars, rand.vars, rhs.vars,  groups, ...){  
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
      MErf(form = form, dat=trn,  groups = groups, family = para$family, rand.vars= rand.vars, para = para,   
                                tol= para$tol, max.iter = para$max.iter, include.RE = para$include.RE, verbose = FALSE, 
                                maxdepth= para$maxdepth, glmer.Control=para$glmer.Control, 
                                nAGQ = para$nAGQ, K = para$K, decay = para$decay, likelihoodCheck = para$likelihoodCheck)
    }, 
    
    MEmixgbm = function(trn, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){  
      result <- list()
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 
      result$model <- MEmixgbm(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                               tol= 1e-5, max.iter = para$max.iter, include.RE = para$include.RE, 
                               verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=glmerControl(optimizer = "bobyqa"),
                               nAGQ = 0, K = para$K, krang = para$krange, para$decay) 
       result      
    }, 
    
    
    MEmixgbm2 = function(trn, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){  
      result <- list()
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 
      result$model <- MEmixgbm2(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                                tol= 1e-5, max.iter = para$max.iter, include.RE = para$include.RE, 
                                verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=glmerControl(optimizer = "bobyqa"),
                                nAGQ = 0, K = para$K, krang = para$krange, para$decay) 
      result      
    }, 
    
    GBM = function(trn, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){
      result <- list()
      xx <- trn 
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
      
    if(all(unique(dat[, resp.vars])%in%c(0, 1))){
        xx[, resp.vars] <- factor(ifelse(xx[, resp.vars]==1, "Yes", "No"))
    if(para$opt.para) {
      fitControl <- trainControl(method = para$method, number =  para$number, 
                                 classProbs = TRUE, summaryFunction= twoClassSummary) 
      result$model <-  train(form,  data = xx, method = "gbm", trControl =  fitControl, 
                             verbose = FALSE,  metric = "ROC", tuneLength = para$tuneLength) 
    } else {  
      result$model <-  train(form,  data = xx, method = "gbm", trControl = trainControl(method = "none",  classProbs = TRUE, summaryFunction= twoClassSummary),  
                             verbose = FALSE, tuneGrid = data.frame(n.trees = para$n.trees, n.minobsinnode = para$n.minobsinnode, 
                             interaction.depth=para$interaction.depth, shrinkage=para$shrinkage))   
    }
      
  } else {
    
    if(para$opt.para) {
      fitControl <- trainControl(method = para$method, number =  para$number) 
      result$model <-  train(form,  data = xx, method = "gbm", trControl =  fitControl, 
                             verbose = FALSE,   tuneLength = para$tuneLength) 
    } else {  
      result$model <-  train(form,  data = xx, method = "gbm", trControl = trainControl(method = "none"),  
                             verbose = FALSE, tuneGrid = data.frame(n.trees = para$n.trees, n.minobsinnode = para$n.minobsinnode, 
                              interaction.depth=para$interaction.depth, shrinkage=para$shrinkage))   
    }
    
  } 
   
    result      
}, 

RF = function(trn, para, resp.vars, rand.vars, rhs.vars, reg.vars=NULL, part.vars=NULL, groups, ...){
  result <- list()
  xx <- trn 
  form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
  
  if(all(unique(dat[, resp.vars])%in%c(0, 1))){
    xx[, resp.vars] <- factor(ifelse(xx[, resp.vars]==1, "Yes", "No"))
    if(para$opt.para) {
      fitControl <- trainControl(method = para$method, number =  para$number, 
                                 classProbs = TRUE, summaryFunction= twoClassSummary) 
      result$model <-  train(form,  data = xx, method = "rf", trControl =  fitControl, 
                             verbose = FALSE,  metric = "ROC", tuneLength = para$tuneLength) 
    } else {  
      result$model <-  train(form,  data = xx, method = "rf", trControl = trainControl(method = "none",  
                             classProbs = TRUE, summaryFunction= twoClassSummary),  
                             verbose = FALSE, tuneGrid = data.frame(mtry = floor(length(rhs.vars)/3)), 
                             ntree=para$ntree)   
    }
    
  } else {
    
    if(para$opt.para) {
      fitControl <- trainControl(method = para$method, number =  para$number) 
      result$model <-  train(form,  data = xx, method = "rf", trControl =  fitControl, 
                             verbose = FALSE,   tuneLength = para$tuneLength) 
    } else {  
      result$model <-  train(form,  data = xx, method = "rf", trControl = trainControl(method = "none"),  
                             verbose = FALSE, tuneGrid = data.frame(mtry = floor(length(rhs.vars)/3)), 
                             ntree=para$ntree)   
    }
    
  } 
  result      
} 

    )
  return(res)
}

