#' Trains a Mixed Effect Model Based Recursive partitioning and Mixed Effect 
#'  classification trees for longitudinal continuous, binary and count data. 
#'  These functions can be called directly and more efficiently by the functions \code{MEml} and \code{MEml.lag} 
#'#  
#' @name MEglmTree 
#
#' @param X  data.frame with predictors 
#' @param Y  binary response vector 
#' @param part.vars,rhs.vars partitioning variables for MOB and predictors 
#' @param reg.vars regressors for MOB 
#' @param rand.vars random effect variables 
#' @param groups  character name of the column containing the group identifier 
#' @param initialRandomEffects  [0] a vector of initial values for random effects
#' @param initialProbs [0.5]  a vector of initial conditional probabilites of success  
#' @param tol convergence tolerance 
#' @param max.iter maximum number of iteration for the EM algorithm 
#' @param likelihoodCheck logical: should the likelihood of random effect model be used to check for convergence?  
#' @param verbose logical for printing intermediate trees results 
#' @param include.RE include random effects in glmtree model part?
#' @param con.tree do conditional inference trees instead of rpart?

### These options pertain to the \code{\link[rpart]{rpart} part of estimation
#' @param cv [TRUE] - Should cross-validation be used?
#' @param cpmin [0.0MEBoostedTree001] - complexity parameter used in building a tree before cross-validation
#' @param no.SE [0] - number of standard errors used in pruning (0 if unused)
#' @param minsize,maxdepth,minsplit,minbucket,mincriterion,stump see rpart

#' @param alpha stability parameter in glmtree 
### These options pertain to the \code{\link[lme4]{glmer}} part of estimation
#' @param glmer.Control same as \code{\link[lme4]{glmerControl}} - glmer controls, 
#' default to glmerControl(optimizer = "bobyqa")] 
#' @param nAGQ  as in \code{\link[lme4]{glmer}}, default to 10 
#'		with maximum likelihood or restricted maximum likelihood
#' @param \dots Further arguments passed to or from other methods.
#' @return An object of class \code{MECTree}; a list with items 
#' \item{tree.fit}{fitted classification trees model}
#' \item{glmer.fit}{fitted mixed effect logistic regression model}
#' \item{logLik}{log likelihood of mixed effect logistic regression} 
#' \item{random.effects}{random effect parameter estimates}
#' \item{form}{modified formula for fitted classification trees model}
#' \item{rand.form}{modified formula for random effects} 
#' \item{glmer.Control}{glmer controls}
#' \item{tree.control}{rpart controls}
#' \item{glmer.CI}{estimates of mixed effect logistic regression with 
#'     approximate confidence intervals on the logit scale. More accurate values 
#'     can be obtained by bootstrap}
#' \item{fitted.probs}{fitted probabilites for final model}
#' \item{fitted.class}{fitted class labels for final model}
#' \item{train.perf}{various performance measures for final model on training set}
#' \item{threshold}{classification cut-off}
#
#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#' @import lme4 caret rpart partykit 
NULL 
#
#' @rdname MEglmTree 
#' @export
MEglmTree  <- function(X, ...) UseMethod("MEglmTree")
#'#
#' @rdname MEglmTree 
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mod <- MEMOBTree(form, rand.form, data)) 
#'#
#' }
#
MEglmTree <- function(X, Y, part.vars, reg.vars = "1", rand.vars="1",  
                groups = NULL, initialRandomEffects=rep(0, nrow(X)), 
                initialProbs = rep(0.5, nrow(X)),  tol= 1e-5, max.iter =100, 
                include.RE = TRUE, verbose=FALSE, likelihoodCheck = TRUE, 
                glmer.Control=glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
                nAGQ=0, alpha = 0.05, minsize = 20, maxdepth=50, para = NULL,  ...){
     if(is.null(groups)) stop("please provide grouping variable")
     Y <- as.vector(Y) 
     dat <- cbind.data.frame(response = Y, X)   
     resp.vars <- "response"          
     dat[, resp.vars] <- factor(dat[, resp.vars])
     dat[, "Y"] <- as.numeric(dat[, resp.vars])-1 
        
    old.lik <- -Inf    
	if(likelihoodCheck == FALSE){
	n.re <- sum(rand.vars != 1)+1  
    b.old <-rep(0, n.re*nlevels(factor(dat[, groups])))  ## initial random effect parameters 
    }
# compute transformed response from initial values
	pp <- initialProbs
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp)) 
	### weights 
    w = pp*(1-pp)    
    Y.star <- qlogis(pp) + (Y - pp)/w   
#### Ajusted Target for MOB regresssion trees     
    Target <- Y.star - initialRandomEffects
    
### Initial formula for MOB trees  
    tmp.form <- as.formula(paste(paste0(resp.vars, " ~"),     
                paste0(c(paste0(reg.vars, collapse = "+"), "|", 
                paste0(part.vars, collapse = "+")), collapse = ""))) 
                
### Final formula for MOB trees 
#	if(length(reg.vars) <= 1) 
#	  reg.vars <- "R.Effect"
#	else reg.vars <- c(reg.vars, "R.Effect") 
#				
#     mob.form <- as.formula(paste(paste0("Target", " ~"),     
#                paste0(c(paste0(reg.vars, collapse = "+"), "|", 
#                paste0(part.vars, collapse = "+")), collapse = ""))) 

### Final formula for MOB trees 
	if(include.RE)			
		 mob.form <- as.formula(paste(paste0("Target", " ~"),     
		            paste0(c(paste0(c("1", reg.vars, "R.Effect"), collapse = "+"), "|", 
		            paste0(c(part.vars), collapse = "+")), collapse = ""))) 
	else 
		 mob.form <- as.formula(paste(paste0("Target", " ~"),     
		            paste0(c(paste0(reg.vars, collapse = "+"), "|", 
		            paste0(part.vars, collapse = "+")), collapse = ""))) 
              

              
### create formula for mixed effect logistic regression model using original binary response
### and classfication tree nodes ( nodes > 1) as fixed effect predictors i.e the mode   
### nu_ij = I(X_ij in terminal node)mu_v + Z_ij b_i 
	form.glmer.multinode <- as.formula(paste0("Y ~ ", paste0(c("Tree.Node +", 
       "(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             
### when there is only a single node in decision tree    
 	form.glmer.onenode <- as.formula(paste0("Y ~ ", paste0(c("Tree.Node +", "(", 
 	paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))    
 	
### Make a new data frame to include all the new variables
    newdata <- dat
    dist <- binomial  
    w <- rep(1, nrow(dat))
    gbmfit <- NULL
    newdata[,"Target"] <- Target

for(ii in 1:max.iter){    	  	
	# Compute MOB tree : 
	# use logistic regression model in MOB for initialization. 
	# subsequently use linear regression model.  

#	 tree <- glmtree(tmp.form, data = newdata, family = dist,  maxit = 50, weights = w, caseweights = FALSE, ...)  
	tree <- glmtree(tmp.form, data = newdata, family = dist, weights = w, maxit = 50,caseweights = FALSE, 
	              epsilon = 1e-8, alpha = alpha, minsize = minsize, maxdepth = maxdepth, ...)  	 
	dist <- gaussian 
	tmp.form <- mob.form
 
     	
	if(width(tree) <= 1){
    warning("Tree with one node encounnted: Using gbm !")
# 	  mob.form <- as.formula(paste(paste0("Target", " ~"),     
# 	                               paste0(c(paste0(c("1", reg.vars), collapse = "+"), "|", 
# 	                                        paste0(c(part.vars), collapse = "+")), collapse = ""))) 
# 	  
	  form <- as.formula(paste0(paste0("Target", " ~"), paste0(part.vars, collapse = "+")))
	  
	  gbmfit <-  train(form,  data = newdata, method = "gbm", trControl = trainControl(method = "none"),  
	                   verbose = FALSE, tuneGrid = data.frame(n.trees = para$n.trees, interaction.depth=para$interaction.depth, 
                           shrinkage=para$shrinkage, n.minobsinnode=para$n.minobsinnode))   
	  zz = predict(gbmfit, newdata = newdata, type = "raw")	



	  newdata[,"Tree.Node"] <- zz
    
#	next  
	} else  {	 
	    Tree.Node <- predict(tree, newdata = newdata, type = "node")        
### Estimate New Random Effects and Errors using LME
### Get variables that identify the node for each observation
	 newdata[,"Tree.Node"] <- factor(Tree.Node)
	}	

###  Fit mixed effect logistic regression model with nodes as fixed effect predictors
### Check that the fitted tree has at least two nodes.
### Check that the fitted tree has at least two nodes.
        if(width(tree) <= 1) 
        form.glmer  <- form.glmer.onenode
        else form.glmer <- form.glmer.multinode		           
		 glmer.fit <- glmer(form.glmer, data=newdata,family = binomial, control = glmer.Control, 
		 nAGQ=  nAGQ, verbose = as.numeric(verbose))
		 
#  get predicted probabilities and compute transformed response 
	pp <- predict(glmer.fit, newdata = newdata, type = "response")
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
### weights 
    w = pp*(1-pp)
    Y.star <- qlogis(pp) + (Y - pp)/w
     
### vompute the adjusted response 
### first get the random effect component 
	Z <-   getME(glmer.fit, name = "Z")
	b  <-  getME(glmer.fit, name = "b")	
	Zb <-  as.numeric(Z%*%b)
	Target  <- Y.star  - Zb

### test for convergence             
   if(likelihoodCheck){
	new.lik <- as.numeric(logLik(glmer.fit))
	r <- as.numeric(sqrt(t(new.lik - old.lik)%*%(new.lik-old.lik)))
	old.lik <- new.lik	
	} else {
	r <- as.numeric(sqrt(t(b - b.old)%*%(b-b.old)))
	b.old <- b
    } 
	#if(verbose) 
	print(paste("Error: ", r))    
	if( r < tol) break 
### update if no convergence 
### update data with surogate response and mixed efffects 
	newdata[,"Target"] <- Target
	newdata[,"R.Effect"] <- Zb 
	 		
	} ## for loop 

	if(r > tol) warning("EM algorithm did not converge")
		 
 if(width(tree) <= 1){
   warning("Tree with one node encounnted: Using gbm !")
   form <- as.formula(paste0(paste0("Target", " ~"), paste0(part.vars, collapse = "+")))
   gbmfit <-  train(form,  data = newdata, method = "gbm", trControl = trainControl(method = "none"),  
                    verbose = FALSE, tuneGrid = data.frame(n.trees = para$n.trees, 
                    interaction.depth=para$interaction.depth, shrinkage=para$shrinkage, n.minobsinnode=para$n.minobsinnode))   
    zz = predict(gbmfit, newdata = newdata, type = "raw")
   
  } else 
   zz <- predict(tree, newdata = newdata , type = "response")

  fitted = zz + Zb 
	
	pp <- sigmoid(fitted) 
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
    
         perf <- Performance.measures(pp, Y)
	threshold <- perf$threshold	
  cls <-ifelse(pp >= threshold,1,0)

	names(cls) <- NULL 	
### get confidence intervals for mixed effect logistic regresion: rough estimates using the SEs
	se <- sqrt(diag(as.matrix(vcov(glmer.fit)))) 
  	tab <- cbind(Est = fixef(glmer.fit), LL = fixef(glmer.fit) - 1.96 * se, 
  	UL = fixef(glmer.fit) + 1.96 *se)

res <- list(tree.fit = tree, glmer.fit = glmer.fit, part.vars=part.vars, groups = groups, 
          reg.vars = reg.vars, rand.vars=rand.vars, logLik=as.numeric(logLik(glmer.fit)), 
         random.effects =ranef(glmer.fit), form = mob.form, rand.form = form.glmer, 
         glmer.Control=	glmer.Control, glmer.CI =tab, fitted.probs = pp, 
         fitted.class = cls, train.perf = perf, threshold = threshold, include.RE=include.RE, gbmfit = gbmfit)
class(res) <- "MEglmTree"         
res         
}
#
#' @rdname MEglmTree  
#' @export
###### Mixed effect conditional inference trees 
MECTree  <- function(X, ...) UseMethod("MECTree")
#
#' @rdname MEglmTree 
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mod <-MECTree(fix.form, rand.form, data)) 
#'
#' }
#
MECTree <- function(X, Y, con.tree = FALSE, rhs.vars,  rand.vars="1",  
                groups = NULL, tol= 1e-5, max.iter =100, verbose=FALSE, 
                likelihoodCheck = TRUE,glmer.Control = glmerControl(optimizer = "bobyqa",check.nobs.vs.nRE="ignore", check.nobs.vs.nlev="ignore"), 
                nAGQ=0, cv = TRUE, cpmin=0.0001, minsplit = 50, minbucket = 10, 
                no.SE =1, mincriterion = 0.975, maxdepth=30, stump = FALSE,...){					

	if(is.null(groups)) stop("please provide grouping variable")
	Y <- as.vector(Y) 
	dat <- cbind.data.frame(response = Y, X)   
	resp.vars <- "response"          
	dat[, resp.vars] <- factor(dat[, resp.vars])
	dat[, "Y"] <- as.numeric(dat[, resp.vars])-1
	Y <- dat[, "Y"]

	old.lik <- -Inf    
	if(likelihoodCheck == FALSE){
	n.re <- sum(rand.vars != 1)+1  
	b.old <-rep(0, n.re*nlevels(factor(dat[, groups])))  ## initial random effect parameters 
	}
# compute transformed response from initial values
	pp <- rep(0.5, nrow(dat))
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp)) 
	### weights 
    w = pp*(1-pp)    
    Y.star <- qlogis(pp) + (Y - pp)/w   
#### Ajusted Target for MOB regresssion trees 
         
    Target <- Y.star   
     
    Mec.form <- as.formula(paste(paste0("Target", " ~"), paste0(rhs.vars, collapse = "+")))                               
### create formula for mixed effect logistic regression model using original binary response
### and classfication tree nodes ( nodes > 1) as fixed effect predictors i.e the mode   
### nu_ij = I(X_ij in terminal node)mu_v + Z_ij b_i 
	form.glmer.multinode <- as.formula(paste0("Y ~ ", paste0(c("Tree.Node +", 
       "(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = ""))) 
                   
### when there is only a single node in decision tree    
 	form.glmer.onenode <- as.formula(paste0("Y ~ ", paste0(c("(", 
 	paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))    
 	
### Make a new data frame to include all the new variables
    newdata <- dat 
    newdata[,"Target"] <- Target
    w <- rep(1, nrow(newdata))
        
for(ii in 1:max.iter){  
 	if(con.tree){
	tree <- ctree(Mec.form, data = newdata, weights = w, 
	            control = ctree_control(mtry = floor(length(rhs.vars)/3), 
                stump = stump,  minsplit = minsplit, mincriterion = mincriterion), ...)     	
	Tree.Node <- predict(tree, newdata = newdata, type = "node")
	} else {
	
         if (cv) {
             tree1 <- rpart(Mec.form, data = newdata, weights = w, 
                 method = "anova", control = rpart.control(cp=cpmin, minsplit=minsplit,
                 minbucket=minbucket, maxdepth=maxdepth))
             if (nrow(tree1$cptable)==1){
               tree <- tree1}
             else {
               cventry <- which.min(tree1$cptable[, "xerror"])
               if (no.SE == 0){
                 cpcv <- tree1$cptable[cventry, "CP"]
                 tree <- prune(tree1, cp=cpcv)}
               else {
                 xerrorcv <- tree1$cptable[cventry, "xerror"]
                 sexerrorcv <- xerrorcv + tree1$cptable[cventry, "xstd"] * no.SE
                 cpcvse <- tree1$cptable[which.max(tree1$cptable[, "xerror"] <= sexerrorcv), "CP"]
                 tree <- prune(tree1, cp=cpcvse)}
         	}
	  } else 
             tree <- rpart(Mec.form, data = newdata, weights = w, method = "anova", 
                       control = rpart.control(minsplit=minsplit,minbucket=minbucket, maxdepth=maxdepth))
                               
		Tree.Node <- tree$where
	} 
	
	newdata[,"Tree.Node"] <- factor(Tree.Node)

###  Fit mixed effect logistic regression model with nodes as fixed effect predictors
### Check that the fitted tree has at least two nodes.
	if(length(unique(Tree.Node)) <= 1) 
	form.glmer  <- form.glmer.onenode
	else form.glmer <- form.glmer.multinode		           
	 glmer.fit <- glmer(form.glmer, data=newdata,family = binomial, control = glmer.Control, 
	 nAGQ=  nAGQ, verbose = as.numeric(verbose))
		 
	pp <- predict(glmer.fit, newdata = newdata, type = "response")
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
    w = pp*(1-pp)
    Y.star <- qlogis(pp) + (Y - pp)/w
     
### compute the adjusted response 
	Z <-   getME(glmer.fit, name = "Zt")
	b  <-  getME(glmer.fit, name = "b")	
	Zb <-  as.numeric(cprod(x=Z,y=b))
	Target  <- Y.star  - Zb	
	newdata[,"Target"] <- Target

### test for convergence             
   if(likelihoodCheck){
	new.lik <- as.numeric(logLik(glmer.fit))
	r <- as.numeric(sqrt(t(new.lik - old.lik)%*%(new.lik-old.lik)))
	old.lik <- new.lik	
	} else {
	r <- as.numeric(sqrt(t(b - b.old)%*%(b-b.old)))
	b.old <- b
    } 
	print(paste("Error: ", r))
	if( r < tol) break 
		
	} ## for loop 
	
	if(r > tol) warning("EM algorithm did not converge")
    if(con.tree)
	fitted = predict(tree, newdata = newdata, type = "response") + Zb 
	else 
	fitted = predict(tree, newdata = newdata, type = "vector") + Zb 
	
	pp <- sigmoid(fitted) 
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
    
    perf <- Performance.measures(pp, Y)
	threshold <- perf$threshold	
  cls <-ifelse(pp >= threshold, "OUT.CNTL", "IN.CNTL")	

### get confidence intervals for mixed effect logistic regresion: rough estimates using the SEs
	se <- sqrt(diag(as.matrix(vcov(glmer.fit)))) 
  	tab <- cbind(Est = fixef(glmer.fit), LL = fixef(glmer.fit) - 1.96 * se, 
  	UL = fixef(glmer.fit) + 1.96 *se)

res <- list(tree.fit = tree, glmer.fit = glmer.fit, groups = groups, 
          rand.vars=rand.vars, logLik=as.numeric(logLik(glmer.fit)), 
         random.effects =ranef(glmer.fit), form = Mec.form, rand.form = form.glmer, 
         glmer.Control=	glmer.Control, glmer.CI =tab, fitted.probs = pp, 
         fitted.class = cls, train.perf = perf, threshold = threshold, con.tree=con.tree)
class(res) <- "MECTree"         
res         
}



















