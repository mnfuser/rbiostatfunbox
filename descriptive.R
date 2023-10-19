# DESCRIPTIVE.TEMPLATE v1

source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

# PATH ----

wd<-"/Users/mmanfrini/Analisi/tactic/revisione"
setwd(wd)
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
nm<-read.csv2(paste0(indir,"/preprocess/nm.dset.ana.filtered.imputed.csv"),
              stringsAsFactors=T)


# SET VAR TYPE ----

setVarType(paste0(indir, "/preprocess/nm.dset.ana.filtered.imputed.csv"))


if(length(numericList)>0){
  dset <- applyNumeric(dset, numericList)
} else {
  print("No numeric vars")
}
if(length(integerList)>0){
  dset <- applyInteger(dset, integerList)
} else {
  print("No integer vars")
}
if(length(factorList)>0){
  dset <- applyFactors(dset, factorList)
} else {
  print("No factor vars")
}
if(length(orderedList)>0){
  dset <- applyOrdered(dset, orderedList)
} else {
  print("No ordered vars")
}
if(length(logicalList)>0){
  dset <- applyLogical(dset, logicalList)
} else {
  print("No logical vars")
}
if(length(dates)>0){
  dset <- applyDate(dset, dates, "%d/%m/%Y")
} else {
  print("No date vars")
}

# VISUALIZE ANALYSIS DATASET
vis_dat_wrapper(dset.ana)

# SW NORM TEST ----
#+ tidy=TRUE, tidy.opts = list(blank = FALSE, width.cutoff = 50), echo=FALSE, message=FALSE, warning=FALSE, comment=NA

getNumbersFromDataset(dset.ana)

sw.table<-NULL
for(i in numericList){
    sw<-shapiro.test(dset.ana[,i])
    sw.table<-rbind(sw.table, cbind(
      colnames(dset.ana)[i],
      round(sw$statistic, 3),
      ifelse(sw$p.value<0.001,"<0.001",round(sw$p.value, 3))
    ))
}
colnames(sw.table)<-c("Var", "Statistics", "p.value")
sw.table

pdf(file="QQ.ana.set.pdf", height = 4, width = 4)
for(i in 1:length(numericList)){
  qqnorm(dset.ana[, numericList[i]], main=colnames(dset.ana)[numericList[i]], cex.main=0.5)
  qqline(dset.ana[, numericList[i]])
}
dev.off()

# DESCRIPTIVE STATS  DATASET ----
#+ tidy=TRUE, tidy.opts = list(blank = FALSE, width.cutoff = 50), echo=FALSE, message=FALSE, warning=FALSE, comment=NA

colrange<-nm[which(nm$rg==1 | nm$oc==1), "id"]

#baseline<-nm[which(nm$rg==1 | nm$oc==1), "id"]

regressorList<-nm[which(nm$rg==1), "id"]

dset.descr<-dset.ana[,colrange]

# CORRELATION ----

# getNumbersFromDataset(dset.descr)
#
# pdf(file="corplot.ana.set.pdf", height = 4, width = 4)
#   for(i in 1:length(numericList)){
#     corrplot(cor(dset.descr[, numericList], na.rm = T), method = "number", type = "upper")
#   }
# dev.off()
#
# pairs(dset.descr[, numericList])
#
rc<-rcorr(as.matrix(dset.ana[, numericList]))$P
write.csv2(rc, file=paste0(outdir, "/corr.matrix.p.csv"), row.names = F)

rc<-rcorr(as.matrix(dset.ana[, numericList]))$r
write.csv2(rc, file=paste0(outdir, "/corr.matrix.r.csv"), row.names = F)

# DESCR STATS OUTCOME = MORTE_totale ----

# rms

ddist<-datadist(dset.descr)
options(datadist=ddist)

tab1<-summary(MORTE_totale ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.MORTE_totale.rms.csv"))

# Table1

fact<-getFactorsNamesFromDataset(dset.descr)

nnorm<-getNumbersNamesFromDataset(dset.descr)

# Comparison
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[86],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir,"/descriptive.MORTE_totaleE.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.MORTE_totale.overall.csv"))

# DESCR STATS OUTCOME = HOSP_nonCV ----

# rms

ddist<-datadist(dset.descr)
options(datadist=ddist)

tab1<-summary(HOSP_nonCV ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.HOSP_nonCV.rms.csv"))

# Table1

fact<-getFactorsNamesFromDataset(dset.descr)

nnorm<-getNumbersNamesFromDataset(dset.descr)

# Comparison
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[87],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir,"/descriptive.HOSP_nonCV.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.HOSP_nonCV.overall.csv"))

# DESCR STATS OUTCOME = HOSP_CV ----

# rms

ddist<-datadist(dset.descr)
options(datadist=ddist)

tab1<-summary(HOSP_CV ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.HOSP_CV.rms.csv"))

# Table1

fact<-getFactorsNamesFromDataset(dset.descr)

nnorm<-getNumbersNamesFromDataset(dset.descr)

# Comparison
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[88],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir,"/descriptive.HOSP_CV.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.HOSP_CV.overall.csv"))

# DESCR STATS OUTCOME = ANY_HOSP ----

# rms

ddist<-datadist(dset.descr)
options(datadist=ddist)

tab1<-summary(ANY_HOSP ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.ANY_HOSP.rms.csv"))

# Table1

fact<-getFactorsNamesFromDataset(dset.descr)

nnorm<-getNumbersNamesFromDataset(dset.descr)

# Comparison
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[89],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir,"/descriptive.ANY_HOSP.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.ANY_HOSP.overall.csv"))

# DESCR STATS OUTCOME = MORTE_HOSP ----

# rms

ddist<-datadist(dset.descr)
options(datadist=ddist)

tab1<-summary(MORTE_HOSP ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.MORTE_HOSP.rms.csv"))

# Table1

fact<-getFactorsNamesFromDataset(dset.descr)

nnorm<-getNumbersNamesFromDataset(dset.descr)

# Comparison
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[90],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir,"/descriptive.MORTE_HOSP.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.MORTE_HOSP.overall.csv"))
