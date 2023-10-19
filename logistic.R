source("/Users/mmanfrini/Code/rbiostatfunbox/rbiostatfunbox/funbox.R")

library(MASS)

# PATH ----

wd<-"/Users/mmanfrini/Analisi/malagu"
setwd(wd)
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana <- read.csv2(paste0(indir,"/dset.ana.csv"),
                        stringsAsFactors=TRUE, strip.white = TRUE, quote = "") #, na.strings = ""
}

# Load dataset descriptor file
nm<-read.csv2(paste0(indir,"/nm.dset.ana.filtered.csv"),
              stringsAsFactors=T)

dset.ana$Funzionediastolica <- factor(dset.ana$Funzionediastolica, levels = c("Normale", 
                                                                              "Alterata", 
                                                                              "Nonvalutabile"))

sel<-which(dset.ana$Funzionediastolica=="Nonvalutabile")
sel
if(length(sel)==0){
  dset.ana$Funzionediastolica<-droplevels(dset.ana$Funzionediastolica)
}

# SET VAR TYPE ----

setVarType(paste0(indir, "/nm.dset.ana.filtered.csv"))


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

# LOGISTIC UNIVAR ----

regressorList<-nm[which(nm$rg==1), "id"]
outcomevarindex=nm[which(nm$oc==1), "id"]

logi.univar(dset.ana, outcomevarindex, regressorList, outdir, "logi.AF.univar")

# VARIABLE SELECTION ----




## MODEL 1 ----

# ETA
# IPOTIROIDISMO 
# IPERTENSIONE PLMONARE
# SPLENECTOMIA, 
# Volumetelesistolicoatrialesinistroml
# Funzionediastolica (REF = NORMALE)
# Atrialesnspessoregrassoinmm

dset.ana.complete<-dset.ana[complete.cases(dset.ana[, c(1, 16, 27, 28, 38, 110, 122, 162)]), 
                            c(1, 16, 27, 28, 38, 110, 122, 162)]

# MASS

full.model<-glm(FA ~ Eta
      + Ipotiroidismo 
      + Ipertensionepolmonare 
      + Splenectomia
      + Volumetelesistolicoatrialesinistroml
      + Funzionediastolica
      + Atrialesnspessoregrassoinmm
      , 
      family = binomial(), data = dset.ana.complete)

sel.model <- full.model %>% stepAIC(trace = FALSE)

summary(sel.model)

confint(sel.model)

# BASE R

sel.model <- step(full.model)

summary(sel.model)

confint(sel.model)

# TABLE

or=(exp(cbind("Odds ratio" = coef(sel.model), confint.default(sel.model))))

multi.1<-data.frame(
  
  #names(logi.fit$coefficients),
  OR=round(or[1:dim(or)[1],1],3),
  CI95=paste0(round(or[1:dim(or)[1],2],3)," - " ,round(or[1:dim(or)[1],3],3)),
  p=round(summary(sel.model)$coefficients[1:dim(summary(sel.model)$coefficients)[1],4] , 3)
  
)

# OUT

write.csv2(multi.1, file=paste0(outdir, "/", "multi.1.csv"), row.names = T)


## MODEL 2 ----

# ETA
# IPOTIROIDISMO 
# PAPsmmHg
# SPLENECTOMIA, 
# Volumetelesistolicoatrialesinistroml
# Funzionediastolica (REF = NORMALE)
# Atrialesnspessoregrassoinmm

dset.ana.complete<-dset.ana[complete.cases(dset.ana[, c(1, 16, 121, 28, 38, 110, 122, 162)]), 
                            c(1, 16, 121, 28, 38, 110, 122, 162)]

# MASS

full.model<-glm(FA ~ Eta
                + Ipotiroidismo 
                + PAPsmmHg 
                + Splenectomia
                + Volumetelesistolicoatrialesinistroml
                + Funzionediastolica
                + Atrialesnspessoregrassoinmm
                , 
                family = binomial(), data = dset.ana.complete)

sel.model <- full.model %>% stepAIC(trace = FALSE)

summary(sel.model)

confint(sel.model)

# BASE R

sel.model <- step(full.model)

summary(sel.model)

confint(sel.model)

# TABLE

or=(exp(cbind("Odds ratio" = coef(sel.model), confint.default(sel.model))))

multi.2<-data.frame(
  
  #names(logi.fit$coefficients),
  OR=round(or[1:dim(or)[1],1],3),
  CI95=paste0(round(or[1:dim(or)[1],2],3)," - " ,round(or[1:dim(or)[1],3],3)),
  p=round(summary(sel.model)$coefficients[1:dim(summary(sel.model)$coefficients)[1],4] , 3)
  
)

# OUT

write.csv2(multi.2, file=paste0(outdir, "/", "multi.2.csv"), row.names = T)



## MODEL 3 ----

# ETA
# IPOTIROIDISMO 
# IPERTENSIONE PLMONARE
# SPLENECTOMIA, 
# Volumetelesistolicoatrialesinistroindicizzatomlm2
# Funzionediastolica (REF = NORMALE)
# Atrialesnspessoregrassoinmm

dset.ana.complete<-dset.ana[complete.cases(dset.ana[, c(1, 16, 27, 28, 38, 111, 122, 162)]), 
                            c(1, 16, 27, 28, 38, 111, 122, 162)]

# MASS

full.model<-glm(FA ~ Eta
                + Ipotiroidismo 
                + Ipertensionepolmonare 
                + Splenectomia
                + Volumetelesistolicoatrialesinistroindicizzatomlm2
                + Funzionediastolica
                + Atrialesnspessoregrassoinmm
                , 
                family = binomial(), data = dset.ana.complete)

sel.model <- full.model %>% stepAIC(trace = FALSE)

summary(sel.model)

confint(sel.model)

# BASE R

sel.model <- step(full.model)

summary(sel.model)

confint(sel.model)

# TABLE

or=(exp(cbind("Odds ratio" = coef(sel.model), confint.default(sel.model))))

multi.3<-data.frame(
  
  #names(logi.fit$coefficients),
  OR=round(or[1:dim(or)[1],1],3),
  CI95=paste0(round(or[1:dim(or)[1],2],3)," - " ,round(or[1:dim(or)[1],3],3)),
  p=round(summary(sel.model)$coefficients[1:dim(summary(sel.model)$coefficients)[1],4] , 3)
  
)

# OUT

write.csv2(multi.3, file=paste0(outdir, "/", "multi.3.csv"), row.names = T)





## MODEL 4 ----

# ETA
# IPOTIROIDISMO 
# PAPsmmHg
# SPLENECTOMIA, 
# Volumetelesistolicoatrialesinistroindicizzatomlm2
# Funzionediastolica (REF = NORMALE)
# Atrialesnspessoregrassoinmm

# dset.ana.complete<-dset.ana[complete.cases(dset.ana[, c(1, 16, 121, 28, 38, 111, 122, 162)]), 
#                            c(1, 16,121, 28, 38, 111, 122, 162)]

# MASS

# full.model<-glm(FA ~ Eta
#                 + Ipotiroidismo 
#                 
#                 + Splenectomia
#                 + Volumetelesistolicoatrialesinistroindicizzatomlm2
#                 + Funzionediastolica
#                 + Atrialesnspessoregrassoinmm
#                 , 
#                 family = binomial(), data = dset.ana.complete)

dset.ana.complete<-dset.ana[, c(1, 16, 28, 38, 111, 122, 162)]

cc=complete.cases(dset.ana.complete)

dset.ana.complete<-dset.ana.complete[cc, ]
                                                       

full.model<-glm(FA ~ Eta
                + Ipotiroidismo 
                + Splenectomia
                + Volumetelesistolicoatrialesinistroindicizzatomlm2
                + Funzionediastolica
                + Atrialesnspessoregrassoinmm
                , 
                family = binomial(), data = dset.ana.complete)

sel.model <- full.model %>% stepAIC(trace = FALSE)

summary(sel.model)

confint(sel.model)

# BASE R

sel.model <- step(full.model)

summary(sel.model)

confint(sel.model)

# TABLE

or=(exp(cbind("Odds ratio" = coef(sel.model), confint.default(sel.model))))

multi.4<-data.frame(
  
  #names(logi.fit$coefficients),
  OR=round(or[1:dim(or)[1],1],3),
  CI95=paste0(round(or[1:dim(or)[1],2],3)," - " ,round(or[1:dim(or)[1],3],3)),
  p=round(summary(sel.model)$coefficients[1:dim(summary(sel.model)$coefficients)[1],4] , 3)
  
)

# TEST MODELLO 
# pvalue = 1 - pchisq(78.161 - 52.304, 70 - 67)

# OUT

write.csv2(multi.4, file=paste0(outdir, "/", "multi.4.01.csv"), row.names = T)

# ROC ---- 
library(pROC)
library(verification)

#### pROC ###

# RICHIEDE IN INPUT OBS (binary), PRED (PROB [0,1])

# un predittore

uni.model<-glm(FA ~ 
                 Atrialesnspessoregrassoinmm
               , 
               family = binomial(), data = dset.ana.complete)

summary(uni.model)

obs=as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]
resp=predict(uni.model, type = "response")

robj<-pROC::roc(obs ~ predict(uni.model, type = "response"),
                plot = T, smooth = F, auc = T, ci = T)
robj

### stat

stats <- wilcox.test(resp[obs == 1], resp[obs == 0], alternative = "great")
stats

# modello multivariabile

obs=as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]
resp=predict(sel.model, type = "response")

robj<-pROC::roc(obs, predict(sel.model, type = "response"),
    plot = T, smooth = F, auc = T, ci = T, algorithm = 1)

robj

### stat

stats <- wilcox.test(resp[obs == 1], resp[obs == 0], alternative = "great")
stats

#### verification ###

# RICHIEDE IN INPUT OBS (binary), PRED (PROB [0,1])


# un predittore

 # Binary indicator
 # dset.ana.complete$arbin<-ifelse(dset.ana.complete$Atrialesnspessoregrassoinmm>=4.1,1,0)

uni.model<-glm(FA ~ 
                
                Atrialesnspessoregrassoinmm
                , 
                family = binomial(), data = dset.ana.complete)

summary(uni.model)

obs=as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]
resp=predict(uni.model, type = "response")

roc.area(obs, 
         predict(uni.model, type = "response"))

# modello multivariabile

obs=as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]
resp=predict(sel.model, type = "response")

roc.area(obs, predict(sel.model, type = "response"))

# CUT OFF 

library(cutpointr)

dset.ana.complete$fFA<-as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]

cp <- cutpointr(dset.ana.complete, Atrialesnspessoregrassoinmm, fFA, 
                method = maximize_metric, metric = sum_sens_spec)
cp
plot(cp)

# NOTE: ----
# rocarea function


  id <- is.finite(obs) & is.finite(resp)
  obs <- obs[id]
  pred <- resp[id]
  n1 <- sum(obs)
  n <- length(obs)
  A.tilda <- (mean(rank(pred)[obs == 1]) - (n1 + 1)/2)/(n - 
                                                          n1)
  stats <- wilcox.test(pred[obs == 1], pred[obs == 0], alternative = "great")
  return(list(A = A.tilda, n.total = n, n.events = n1, 
              n.noevents = sum(obs == 0), p.value = stats$p.value))

