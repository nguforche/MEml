library(MEml)
library(Formula)
library(gbm)
library(plyr)
library(caret)
library(lme4)
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
  max.iter = 10, alpha=0.05, minsize=20,maxdepth=30,  
  K = 3, 
  krange = 2:5,
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
order.vars = "time"

######
form <- as.formula(paste0(paste0(resp.vars, " ~"), paste0(rhs.vars, collapse = "+"))) 

model <- MEmixgbm(form = form, dat=dat, groups = id,  rand.vars= rand.vars,  para = para,   
                max.iter =20, include.RE =FALSE, maxdepth=5, k=3, krange = 2:5, decay = 0.05)



### plot 
dat$risk <- model$Y.star
dat$cluster <- model$mixfit@cluster
dat$time <- round(dat$time)

### mean of trajectory cluster at follow up time 
mn <- ddply(dat, .variables = c("time", "cluster"), .fun = function(xx) {
  xx$meanRisk = mean(xx$risk, na.rm = TRUE)
  xx
})

sp <- dlply(mn, .variables = "cluster", .fun = function(xx){
  sm <- smooth.spline(xx$time, xx$meanRisk, df=3, all.knots = TRUE)
  tab <- cbind.data.frame(x = sm$x, y = sm$y)
  list(smooth = sm, tab = tab)
})

### data for VR 
dd <- ddply(mn, .variables = "cluster", .fun = function(xx){
  sp1 <- sp[[unique(xx$cluster)]]$tab
  ddply(xx, .variables = "time", .fun = function(yy){
    ix <- which(sp1$x%in%unique(yy$time))
    yy$smoothRisk = sp1$y[ix]
    yy
  })})





d1 <-dd[, c("id", "cluster", "time", "risk")]
d2 <- dd[, c("id", "cluster", "time",  "smoothRisk")]
names(d2) <- c("id", "cluster", "time",  "risk")
d2$cluster <- d2$cluster + 3

d1$id <- factor(d1$id)
d1$cluster <- factor(d1$cluster)
d2$id <- factor(d2$id)
d2$cluster <- factor(d2$cluster)

d <- rbind(d1, d2)
d$id <- factor(d$id)
d$cluster <- factor(d$cluster)

col = col2hex(c("darkblue", "darkred", "darkgreen", "purple"))
col1 <- c(makeTransparent(col[1], alpha = 20), makeTransparent(col[2], alpha=20), 
          makeTransparent(col[3], alpha = 20))
col2 <- c(col1, col)


pp <- ggplot() + 
  geom_line(data = d1, aes(x = time, y = risk, group = id, color = cluster, size = cluster)) +  
  geom_line(data = d2, aes(x = time, y = risk, group = id, color = cluster, size = cluster )) +
  scale_size_manual(values = c(1, 1, 1, 1.5, 1.5, 1.5)) + 
  geom_jitter() + 
  guides(colour = FALSE, size = FALSE) +scale_x_continuous(name="Time (years from surgery)",breaks=0:11) + 
  scale_color_manual(values = col2 )  + ylab("Predicted Risk") + 
  annotate("text", x = c(9.2, 11.2, 9.2), y = c(-4, 1.3, 4.5), label = c("T1", "T2", "T3"), size = 5) + 
  theme(axis.title.x=element_text(size=18,face="bold"), 
        axis.title.y=element_text(size=18,face="bold"),
        legend.text = element_text(size=18,face="bold"), 
        axis.text.x = element_text(size = 18, face="bold",colour = "gray40"),
        legend.title = element_text(size=18,face="bold"),
        axis.text.y = element_text(size = 18, face="bold",colour = "gray40"))
print(pp)





### ass a highlight individual longitudinal profile of a single patient 
ix <- sapply(d1$id, length)
ip <- subset(d1, id==5)
ip$id <- paste0("patient", ip$id)
ip$cluster = "patient5"
col3 <- c(col2, makeTransparent(col2hex("darkblue"), 150))

pp <- ggplot() + 
  geom_line(data = d1, aes(x = time1, y = traj, group = id, color = cluster, size = cluster)) +  
  geom_line(data = d2, aes(x = time1, y = traj, group = id, color = cluster, size = cluster )) +
  geom_line(data = ip, aes(x = time1, y = traj, group = id, color = cluster, size = cluster)) +
  scale_size_manual(values = c(1, 1, 1, 1.5, 1.5, 1.5, 1.4)) + 
#  geom_jitter() + 
  guides(colour = FALSE, size = FALSE) +
  scale_color_manual(values = col3 ) + 
  xlab("Time (years from AVR)") + ylab("Predicted continous risk function") + 
  annotate("text", x = c(9.0, 10.0, 8.8, 9.5), y = c(-4.2, 0.5, 4.7, 1.6), 
           label = c("Improving course", "Stable/Slow progression", "Rapid progression", 
                     "Trajectory of an example patient"), size = 8) + 
  theme(axis.title.x=element_text(size=22,face="bold"), 
        axis.title.y=element_text(size=22,face="bold"),
        legend.text = element_text(size=22,face="bold"), 
        axis.text.x = element_text(size = 22, face="bold",colour = "gray40"),
        legend.title = element_text(size=22,face="bold"),
        axis.text.y = element_text(size = 22, face="bold",colour = "gray40"))
print(pp)




