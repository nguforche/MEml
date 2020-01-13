library(MEml)
library(Formula)
library(gbm)
library(plyr)
library(caret)
library(lme4)
library(RRF)
library(inTrees)
library(ggplot2)
library(parallel)
library(class)
library(Hmisc)
library(kernlab)
require(flexmix)
require(gplots)
require(bayou)

### predict longidutinal profile of left ventricular mass index (increase or normal) 
### in  patients undergoing aortic valve surgery.
#### LVM is considered increased if >134 g/m(2) in male patients and >110 g/m(2) in female patients 
data(heart.valve)
dat <- heart.valve

levels(dat$sex) <- c("Male", "Female")
levels(dat$lvh) <- c("No", "Yes")
levels(dat$acei) <- c("No", "Yes")
levels(dat$hc) <- c("Absent", "Treated", "Untreated")
levels(dat$prenyha) <- c("I-II", "III-IV")
levels(dat$dm) <- c("No", "Yes")



## potential modifiable variables 
## acei = ace inhibitor(drug). ACE inhibitors treat a variety of conditions, such as high blood pressure, scleroderma and migraines. 
## dm = preoperative diabetes
## creat = preoperative serum creatinine
## hc = preoperative high cholesterol
## prenyha= preoperative New York Heart Association (NYHA): 1 = I/II; 3 = III/IV 
## bsa = preoperative body surface area.
## fuyrs= maximum follow up time, with surgery date as the time origin (years)

#### train mixed effect machine learning 

seed = 123 
para <- list(
  method = "cv",
  tuneLength=3,
  number = 3,
  n.trees=100,
  ntree = 50, 
  interaction.depth=4,
  shrinkage=0.01,
  n.minobsinnode=10,
  opt.para= TRUE, 
  include.RE = FALSE,
  con.tree = FALSE, 
  max.iter = 10, alpha=0.05, 
  minsize=20,
  maxdepth=30,
  
  glmer.Control = glmerControl(optimizer = "bobyqa"), 
  likelihoodCheck = TRUE, 
  nAGQ=0, 
  
  decay = 0.05, 
  K = 3, 
  tol= 1e-5,
  seed = seed
)

dat$id <- as.numeric(dat$id)  ## random effect grouping variable
resp.vars <- "inc.lvmi"
id <- "id"

## fixed effect variables 
rhs.vars <- c("sex", 
              "age", 
              "time", 
              "grad", 
              "log.grad", 
              "bsa", 
              "lvh", 
              "prenyha", 
              "redo", 
              "size",
              "con.cabg", 
              "creat", 
              "dm", 
              "acei", 
              "lv", 
              "emergenc", 
              "hc", 
              "sten.reg.mix", 
              "hs")

rand.vars= "time"  ## random effect variables 

######

ix <- sample(nrow(dat), floor(nrow(dat)*0.75))
dat.trn <- dat[ix, ]
dat.tst <- dat[-ix, ]


document(pkg = Rcode) 

classifier = c("MEglm", "MEgbm", "MErf")

res <- MEml2(classifier= classifier[2], data = dat.trn, id=id,  resp.vars= resp.vars, rhs.vars= rhs.vars,
             rand.vars=rand.vars, para=para)


pred <- predict(res$MEgbm, newdata= dat.tst, type="prob",  allow.new.levels = TRUE)
perf <- Performance.measures(pred$pred[,2], dat.tst[, resp.vars])



res <- MEml2(classifier= classifier[1], data = dat.trn, id=id,  resp.vars= resp.vars, rhs.vars= rhs.vars,
             rand.vars=rand.vars, para=para)
pred <- predict(res$MEglm, newdata= dat.tst, type="prob",  allow.new.levels = TRUE)
perf <- Performance.measures(pred$pred[,2], dat.tst[, resp.vars])




res <- MEml2(classifier= classifier[3], data = dat.trn, id=id,  resp.vars= resp.vars, rhs.vars= rhs.vars,
             rand.vars=rand.vars, para=para)
pred <- predict(res$MErf, newdata= dat.tst, type="prob",  allow.new.levels = TRUE)
perf <- Performance.measures(pred$pred[,2], dat.tst[, resp.vars])

















