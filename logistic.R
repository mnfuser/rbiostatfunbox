#     LOGISTIC REGRESSION TEMPLATE
#
#     GNU GPLv3
# 
#     Copyright (C) 2023  Marco Manfrini
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
#     info@manfrinistudio.it

source("/Users/mmanfrini/Code/rbiostatfunbox/rbiostatfunbox/funbox.R")

library(MASS)
library(readr)

# PATH ----
# SET WORKING DIR, IN AND OUT PATH

wd<-"/Users/mmanfrini/Analisi/malagu"
setwd(wd)
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out")

# GLOBAL VARS ----
# SET GLOBAL VARS

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----
# LOAD DATASET AND DATA DESCRIPTOR FILE

dset.ana <- read_csv2(paste0(indir,"/dset.ana.csv"))
spec(dset.ana)
View(dset.ana)

# Load data descriptor file
ddescr<-read_csv2(paste0(indir,"/ddescr.csv"))
spec(ddescr)
View(ddescr)

# LOGISTIC UNIVARIABLE REGRESSION ----

oFileName="logi.AF.univar"

regressorList<-ddescr[which(ddescr$rg==1), "id"]
outcomevarindex=ddescr[which(ddescrm$oc==1), "id"]

logi.univar(dset.ana, outcomevarindex, regressorList, outdir, oFileName)

# VARIABLE SELECTION ----
# SELECTION OF VARIABLE BY NESTED MODELS EVALUATION OF AIC BIC AND LR

## MULTIVARIABLE MODEL ----

# DEPENDENT VARIABLE
#
# FA
#
# REGRESSORS LIST
#
# ETA
# IPOTIROIDISMO 
# SPLENECTOMIA, 
# Volumetelesistolicoatrialesinistroml
# Funzionediastolica (REF = NORMALE)
# Atrialesnspessoregrassoinmm

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

## MASS

sel.model <- full.model %>% stepAIC(trace = TRUE)

summary(sel.model)

confint(sel.model)

## BASE R

sel.model <- step(full.model)

summary(sel.model)

confint(sel.model)

## OUTPUT TABLE

or=(exp(cbind("Odds ratio" = coef(sel.model), confint.default(sel.model))))

multi.4<-data.frame(
  
  OR=round(or[1:dim(or)[1],1],3),
  CI95=paste0(round(or[1:dim(or)[1],2],3)," - " ,round(or[1:dim(or)[1],3],3)),
  p=round(summary(sel.model)$coefficients[1:dim(summary(sel.model)$coefficients)[1],4] , 3)
  
)

## TEST MODELLO 
## pvalue = 1 - pchisq(78.161 - 52.304, 70 - 67)

# OUT

write.csv2(multi.4, file=paste0(outdir, "/", "multi.00.csv"), row.names = T)

# ROC ---- 

library(pROC)
library(verification)

## pROC

## INPUT OBS (binary), PRED (PROB [0,1])

## un predittore

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

## stat

stats <- wilcox.test(resp[obs == 1], resp[obs == 0], alternative = "great")
stats

## modello multivariabile

## MODEL SELECTED IN PREVIOUS STEP = sel.model 

obs=as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]
resp=predict(sel.model, type = "response")

robj<-pROC::roc(obs, predict(sel.model, type = "response"),
    plot = T, smooth = F, auc = T, ci = T, algorithm = 1)

robj

## stat

stats <- wilcox.test(resp[obs == 1], resp[obs == 0], alternative = "great")
stats

## verification 

# INPUT OBS (binary), PRED (PROB [0,1])

## un predittore

uni.model<-glm(FA ~ 
                
                Atrialesnspessoregrassoinmm
                , 
                family = binomial(), data = dset.ana.complete)

summary(uni.model)

obs=as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]
resp=predict(uni.model, type = "response")

roc.area(obs, 
         predict(uni.model, type = "response"))

## modello multivariabile

## MODEL SELECTED IN PREVIOUS STEP = sel.model 

obs=as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]
resp=predict(sel.model, type = "response")

roc.area(obs, predict(sel.model, type = "response"))

# CUT OFF ----

library(cutpointr)

## INPUT OUTCOME BINARY VARIABLE AS NUMERIC

dset.ana.complete$fFA<-as.numeric(levels(dset.ana.complete$FA))[dset.ana.complete$FA]

cp <- cutpointr(dset.ana.complete, Atrialesnspessoregrassoinmm, fFA, 
                method = maximize_metric, metric = sum_sens_spec)
cp
plot(cp)

# NOTE: ----

# rocarea function


# id <- is.finite(obs) & is.finite(resp)
# obs <- obs[id]
# pred <- resp[id]
# n1 <- sum(obs)
# n <- length(obs)
# A.tilda <- (mean(rank(pred)[obs == 1]) - (n1 + 1)/2)/(n - 
#                                                         n1)
# stats <- wilcox.test(pred[obs == 1], pred[obs == 0], alternative = "great")
# return(list(A = A.tilda, n.total = n, n.events = n1, 
#             n.noevents = sum(obs == 0), p.value = stats$p.value))

