# COX LASSO v 1.0

source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

library(glmnet)

set.seed(321)

# PATH ----

wd<-setwd("~/Analisi/Manfrini/PHENO")
indir=paste0(wd, "/out/clustering")
outdir=paste0(wd,"/out/lasso")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana <- read.csv2(paste0(indir,"/dset.ana.cluster.gmm.csv"),
                        stringsAsFactors=FALSE) #, na.strings = ""
}

# Load dataset descriptor file
nm<-read.csv2(paste0(indir,"/nm.dset.ana.filtered.lasso.csv"), 
              stringsAsFactors=T)


# SET VAR TYPE ----

setVarType(paste0(indir, "/nm.dset.ana.filtered.cluster.csv"))

if(length(numericList)>0){
  dset.ana <- applyNumeric(dset.ana, numericList)  
} else {
  print("No numeric vars")
}
if(length(integerList)>0){
  dset.ana <- applyInteger(dset.ana, integerList)
} else {
  print("No integer vars")
}
if(length(factorList)>0){
  dset.ana <- applyFactors(dset.ana, factorList)
} else {
  print("No factor vars")
}
if(length(orderedList)>0){
  dset.ana <- applyOrdered(dset.ana, orderedList)
} else {
  print("No ordered vars")
}
if(length(logicalList)>0){
  dset.ana <- applyLogical(dset.ana, logicalList)
} else {
  print("No logical vars")
}
if(length(dates)>0){
  dset.ana <- applyDate(dset.ana, dates, "%d/%m/%Y")
} else {
  print("No date vars")
}

# VISUALIZE ANALYSIS DATASET
vis_dat_wrapper(dset.ana)

# LASSO FIT  ----
dset.cox.lasso<-dset.ana
dset.cox.lasso$death<-ifelse(dset.cox.lasso$Outcome==5,1,0)

ddist<-datadist(dset.cox.lasso)
options(datadist=ddist)

units(dset.cox.lasso$time2death) <- 'Day'
dset.cox.lasso[which(dset.cox.lasso$time2death==0), "time2death"]=1
S<-Surv(dset.cox.lasso$time2death, dset.cox.lasso$death)

confounders<-nm[which(nm$lasso==1), "id"]

# confounders<-confounders[1:length(confounders)-1]

expositionList<-nm[which(nm$exp==1), "id"]

x<-data.matrix(dset.cox.lasso[, confounders])

cvfit <- cv.glmnet(x, S, family = "cox", standardize = TRUE, alpha = 1) 

plot(cvfit)

#s = cvfit$lambda.1se
s = cvfit$lambda.min

c<-coef(cvfit, s)

nz<-which(c!=0)

mod.lasso<-data.frame(var=rownames(c)[nz], coef=c[nz])

write.csv2(mod.lasso, file=paste0(outdir, "mod.lasso.csv"), row.names = F)

fit <- glmnet(x, S, family = "cox", alpha = 1, lambda = s, standardize = TRUE)
 
plot(coef(fit))

# 
coef(fit)

# Performance ----

assess.glmnet(fit, newx = x, newy = S)
