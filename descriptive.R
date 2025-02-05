#     DESCRIPTIVE
#
#     Descriptive statistics
#
#     GNU GPLv3
#
#     Copyright (C) 2023  Marco Manfrini
#
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

rm(list = ls())

# Source ---------------------------------------------

source("/Users/mmanfrini/Code/rbiostatfunbox/rbiostatfunbox/funbox.R")

# PATH ----

wd<-"/Users/mmanfrini/Analisi/tactic/revisione"
setwd(wd)
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out")

# GLOBAL VARS ----

# LOAD DATASET PRE PROCESSED ----

dset.ana <- read_csv2(paste0(pcadir,"/dset.ana.cluster.csv"))
problems()
View(dset.ana)

## READING METADATA 

meta <- read_csv2(paste0(indir,"/meta.cluster.csv"))
problems()
View(meta)

## SET VAR TYPE

readMetaVarType(paste0(indir, "/meta.cluster.csv"))
dset.ana<-applyVarType(dset.ana, 
                       nList = numericList, 
                       cList = characterList, 
                       fList = factorList)
View(dset.ana)

# DESCRIPTIVE STATS ----

## SELECTING VARS
indexs<-filter(meta, descr==1)|>
  pull(id)

dset.descr<-select(dset.ana, all_of(indexs))

# SW NORM TEST ----

getNumbersFromDataset(dset.descr)

sw.table<-NULL
for(i in numericList){
  sw<-shapiro.test(dset.descr[,i])
  sw.table<-rbind(sw.table, cbind(
    colnames(dset.descr)[i],
    round(sw$statistic, 3),
    ifelse(sw$p.value<0.001,"<0.001",round(sw$p.value, 3))
  ))
}
colnames(sw.table)<-c("Var", "Statistics", "p.value")
sw.table

pdf(file="QQ.ana.set.pdf", height = 4, width = 4)
for(i in 1:length(numericList)){
  qqnorm(dset.descr[, numericList[i]], main=colnames(dset.descr)[numericList[i]], cex.main=0.5)
  qqline(dset.descr[, numericList[i]])
}
dev.off()

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

# END ----