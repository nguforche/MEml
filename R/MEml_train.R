Train_MEml <- function(...){
  res = list(
    MEglm = function(trn, para, resp.vars, rhs.vars, rand.vars, groups, ...){
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(c(paste0(rhs.vars,collapse="+"), "+", "(", 
                                                                  paste0(c(rand.vars),collapse = "+"), "|", groups, ")"), collapse="")))
      MEglm(form = form, data= trn,control= para$glmer.Control,nAGQ=para$nAGQ)
      
     }, 
    
    MEgbm = function(trn, para, resp.vars, rand.vars, rhs.vars, groups, ...){  
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 
      MEgbmRules(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                                 tol= para$tol, max.iter = para$max.iter, include.RE = para$include.RE, 
                                 verbose = FALSE, maxdepth= para$maxdepth, glmer.Control=para$glmer.Control,
                                 nAGQ = para$nAGQ, K = para$K, decay = para$decay, likelihoodCheck = para$likelihoodCheck) 
    }, 
    
    
    MErf = function(trn, para, resp.vars, rand.vars, rhs.vars,  groups, ...){  
      form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+")))
      MErfRules(form = form, dat=trn,  groups = groups, rand.vars= rand.vars, para = para,   
                                tol= para$tol, max.iter = para$max.iter, include.RE = para$include.RE, verbose = FALSE, 
                                maxdepth= para$maxdepth,glmer.Control=para$glmer.Control, 
                                nAGQ = para$nAGQ, K = para$K, decay = para$decay, likelihoodCheck = para$likelihoodCheck)
    }
    )
  return(res)
}

