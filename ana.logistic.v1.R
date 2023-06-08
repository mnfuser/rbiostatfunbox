source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

# PATH ----

wd<-setwd("/Users/mmanfrini/Analisi/ELISA/AAT")
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana <- read.csv2(paste0(indir,"/dset.ana.csv"),
                        stringsAsFactors=FALSE) #, na.strings = ""
}

# Load dataset descriptor file
nm<-read.csv2(paste0(indir,"/nm.dset.ana.filtered.csv"), 
              stringsAsFactors=T)


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


# LOGISTIC UNIVAR OSP.NOT TRAIN ----
regressorList<-nm[which(nm$rg==1), "id"]
outcomevarindex=nm[which(nm$oc==1), "id"]

logi.univar(dset.ana, outcomevarindex, regressorList, outdir, "logi.infection.univar")

# LOGISTIC MULTI ----


  
  logi.fit<-glm(Infezione ~ AGE + GENDER + CENTRO + AAT.CONC, 
                family = binomial(), data = dset.ana)

  summary(logi.fit)

  or=(exp(cbind("Odds ratio" = coef(logi.fit), confint.default(logi.fit))))
  
  logi.out<-data.frame(
    
    #names(logi.fit$coefficients),
    OR=round(or[1:dim(or)[1],1],3),
    CI95=paste0(round(or[1:dim(or)[1],2],3)," - " ,round(or[1:dim(or)[1],3],3)),
    p=round(summary(logi.fit)$coefficients[1:dim(summary(logi.fit)$coefficients)[1],4] , 3)
    
  )

print(kable(logi.out, row.names = FALSE))

write.csv2(logi.out, file=paste0(outdir, "/", "logi.multi.csv"), row.names = T)



# ROC ---- 
library(pROC)
library(verification)

robj<-roc(dset.ana$Infezione, dset.ana$AAT.CONC,
    plot = T, smooth = T, auc = T, ci = T, algorithm = 1)

roc.area(dset.ana$Infezione, predict(logi.fit))
