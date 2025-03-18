# TRAIN TEST v1

source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

# PATH ----

wd<-setwd("~/Analisi/Manfrini/PHENO")
indir=paste0(wd, "/out/varsel")
outdir=paste0(wd,"/out/treeexp")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana <- read.csv2(paste0(indir,"/dset.ana.cluster.csv"),
                        stringsAsFactors=FALSE) #, na.strings = ""
}

# Load dataset descriptor file
nm<-read.csv2(paste0(indir,"/nm.dset.ana.filtered.cluster.csv"), 
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


# DATA SPLITTING ----

library(caret)

set.seed(3456)

trainIndex <- createDataPartition(dset.ana$cluster, p = .7, 
                                  list = FALSE, 
                                  times = 1)

dset.ana.train<-dset.ana[trainIndex,]
dset.ana.test<-dset.ana[-trainIndex,]

write.csv2(dset.ana.train, file=paste0(outdir, "/dset.ana.train.csv"))
write.csv2(dset.ana.test, file=paste0(outdir, "/dset.ana.test.csv"))
