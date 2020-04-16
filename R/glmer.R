### fit lmer or glmer

myglmer <- function(form, dat, family, control,  verbose, nAGQ){
  
  if (family$family == "binomial"){
    
    glmer(form, data= dat,  family = family, control = control, 
          nAGQ=  nAGQ, verbose = as.numeric(verbose))   
    
  } else if(family$family == "gaussian") {
    
    lmer(form, data = dat, control = control, verbose = as.numeric(verbose))
  
  }  else {
    print(family)
    stop("'family' is not yet implemented")
  }
}
