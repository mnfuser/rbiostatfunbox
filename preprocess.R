# PRE-PROCESSING v4

# Functions ---------------------------------------------
source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

# PATH ----

wd<-setwd("/Users/mmanfrini/Analisi/Manfrini/PHENO")
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/preprocess")

# GLOBAL VARS ----

RAW=T
DOIMP=T

# LOAD DATASET RAW ----
if(RAW==T){
  dset <- read.csv2(paste0(indir,"/dataset.csv"),
                    stringsAsFactors=T)
  }

# FORMAT CONVERSION

# Removing all blank strings " "

#dset<-as.data.frame(apply(dset, 2, function(x)gsub('\\s+', '',x)))

# CREATE VARTYPE FILE ----
varclass<-as.matrix(sapply(dset, class))
if(!"vartype.csv" %in% list.files(indir)){
  vartype<-cbind(colnames(dset), varclass)
  colnames(vartype)<-c("var", "type")
  write.csv2(vartype, file=paste0(outdir, "/vartype.csv"), row.names = F)
}

# SET VAR TYPE ----
setVarType(paste0(indir, "/vartype.csv"))

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

# QC PLOTS RAW ----

vis_dat_wrapper(dset)

pdf(file=paste0(outdir, "/boxplot.pdf"), height = 4, width = 4)
for(i in 1:length(numericList)){
  if(sum(dset[, numericList[i]], na.rm = T)>0){
    boxplot(dset[, numericList[i]], main=colnames(dset)[numericList[i]], cex.main=0.5)
  }
}
dev.off()

pdf(file=paste0(outdir, "/qq.pdf"), height = 4, width = 4)
for(i in 1:length(numericList)){
  if(sum(dset[, numericList[i]], na.rm = T)>0){
    qqnorm(scale(dset[, numericList[i]]), main=colnames(dset)[numericList[i]], cex.main=0.5)
    qqline(scale(dset[, numericList[i]]), col = "steelblue", lwd = 2)
  }
}
dev.off()

pdf(file=paste0(outdir, "/barplot.pdf"), height = 6, width = 4)
for(i in 1:length(factorList)){
  tab<-table(dset[,factorList[i]])
  par(mar=c(10,1,1,1))
  barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset)[factorList[i]], cex.main=0.5, cex.names = 0.5)
}
dev.off()

pdf(file=paste0(outdir, "/integers.pdf"), height = 6, width = 4)
for(i in 1:length(integerList)){
  tab<-table(dset[,integerList[i]])
  par(mar=c(10,1,1,1))
  barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset)[integerList[i]], cex.main=0.5, cex.names = 0.5)
}
dev.off()

nm<-nmiss(dset)
nm$vars<-as.factor(rownames(nm))
pdf(file=paste0(outdir, "/missing.pdf"), height = 0.3*dim(nm)[1], width = 20)

nm %>%
  arrange(desc(pm)) %>%
  mutate(vars = factor(vars, unique(vars))) %>%
  ggplot() +
  aes(x=vars, y=pm) +
  geom_segment( aes(x=vars, xend=vars, y=0, yend=pm), color="skyblue") +
  geom_point( color="blue", size=2, alpha=0.6) +
  #geom_text(aes(label = round(pm),2), hjust = 1, size = 3) +
  theme_light() +
  coord_flip() +
  scale_y_continuous(
    "Missing ratio",
    sec.axis = sec_axis(~ ., name = "Missing ratio")
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

dev.off()

# OUTLIERS CORRECTION ----

getNumbersFromDataset(dset)
for(i in numericList){
  dset[, i]<-ApplyCap2Outliers(dset[, i])
}

# RECODE / BUILD NEW VARS ----

# CENSORING SETTINGS

dset$time2death<-as.numeric(difftime(dset$death.date, dset$hospitalization.date, "%Y-%m-%d", units = "days"))
dset$time2death<-ifelse(is.na(dset$time2death), max(dset$time2death, na.rm = T), dset$time2death)

dset$time2ICU<-as.numeric(difftime(dset$ICU.date, dset$hospitalization.date, "%Y-%m-%d", units = "days"))
dset$time2ICU<-ifelse(is.na(dset$time2ICU), max(dset$time2ICU, na.rm = T), dset$time2ICU)

dset$time2Intubation<-as.numeric(difftime(dset$intubation.date, dset$hospitalization.date, "%Y-%m-%d", units = "days"))
dset$time2Intubation<-ifelse(is.na(dset$time2Intubation), max(dset$time2Intubation, na.rm = T), dset$time2Intubation)

nm<-nmiss(dset, fileName = paste0(outdir,"/nm.dset.csv"))

# ANALYSIS SET ----

nm<-read.csv2(paste0(outdir,"/nm.dset.ana.csv"),
          stringsAsFactors=T)

nm.filtered<-nm[which(nm$rg==1 | nm$oc==1), ]

selected<-nm.filtered$id

dset.ana<-cbind(dset[, c(selected)]
               )

# UPDATE DSET DESCRIPTOR

nm.filtered$id<-seq(1:dim(nm.filtered)[1])

write.csv2(nm.filtered, file=paste0(outdir,"/nm.dset.ana.filtered.csv"), row.names = F)

# VISUALIZE ANALYSIS DATASET
vis_dat_wrapper(dset.ana)

# varclass<-as.matrix(sapply(dset.ana, class))

# SET VAR TYPE ANALYSIS SET----

# SET VAR tYPE
setVarType(paste0(indir, "/nm.dset.ana.filtered.csv")) #vartype.dset.ana.csv

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


# QC PLOTS ANALYSIS SET----

vis_dat_wrapper(dset.ana)

pdf(file=paste0(outdir, "/boxplot.ana.set.pdf"), height = 4, width = 4)
for(i in 1:length(numericList)){
  if(sum(dset.ana[, numericList[i]], na.rm = T)>0){
    boxplot(dset.ana[, numericList[i]], main=colnames(dset.ana)[numericList[i]], cex.main=0.5)
  }
}
dev.off()

pdf(file=paste0(outdir, "/qq.ana.set.pdf"), height = 4, width = 4)
for(i in 1:length(numericList)){
  if(sum(dset.ana[, numericList[i]], na.rm = T)>0){
    qqnorm(scale(dset.ana[, numericList[i]]), main=colnames(dset.ana)[numericList[i]], cex.main=0.5)
    qqline(scale(dset.ana[, numericList[i]]), col = "steelblue", lwd = 2)
  }
}
dev.off()

pdf(file=paste0(outdir, "/barplot.ana.set.pdf"), height = 6, width = 4)
for(i in 1:length(factorList)){
  tab<-table(dset.ana[,factorList[i]])
  par(mar=c(10,1,1,1))
  barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset.ana)[factorList[i]], cex.main=0.5, cex.names = 0.5)
}
dev.off()

pdf(file=paste0(outdir, "/integers.ana.set.pdf"), height = 6, width = 4)
for(i in 1:length(integerList)){
  tab<-table(dset.ana[,integerList[i]])
  par(mar=c(10,1,1,1))
  barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset.ana)[integerList[i]], cex.main=0.5, cex.names = 0.5)
}
dev.off()

nm<-nmiss(dset.ana)

pdf(file=paste0(outdir, "/missing.ana.set.pdf"), height = 5, width = 0.2*dim(nm)[1])
barplot(nm$pm, names.arg = substr(rownames(nm), start=1, stop=10), las=2, width = 0.2 ,main="Missing proportions", cex.main=0.5, cex.names = 0.5)
dev.off()


# IMPUTATION ----

if(DOIMP==T){

  library(mice)
  library(VIM)

  pdf(file="pattern.missing.ana.set.pdf", height = 10, width = 0.5*dim(dset.ana)[2])

  aggr_plot <- aggr(dset.ana,
                    col=c('black','gray'),
                    numbers=TRUE,
                    sortVars=TRUE,
                    labels=colnames(dset.ana),
                    cex.axis=.7,
                    gap=3,
                    ylab=c("Histogram of missing data","Pattern")
  )

  dev.off()

  dset.imp<-dset.ana[, c(1:51)]

  tempData <- mice(dset.imp,
                   pred=quickpred(dset.imp, mincor=.3),
                   maxit = 20,
                   m = 5,
                   meth='pmm',
                   seed=500
                   )

  plot(tempData)

  dset.ana.imp.1<-cbind(complete(tempData,1), dset.ana[, c(52:54)])
  dset.ana.imp.2<-cbind(complete(tempData,2), dset.ana[, c(52:54)])
  dset.ana.imp.3<-cbind(complete(tempData,3), dset.ana[, c(52:54)])
  dset.ana.imp.4<-cbind(complete(tempData,4), dset.ana[, c(52:54)])
  dset.ana.imp.5<-cbind(complete(tempData,5), dset.ana[, c(52:54)])

}

# QC PLOTS IMPUTED ----
# FIRST IMPUTED DATASET

dset.ana<-dset.ana.imp.1

vis_dat_wrapper(dset.ana)

pdf(file=paste0(outdir, "/boxplot.ana.set.imputed.pdf"), height = 4, width = 4)
for(i in 1:length(numericList)){
  if(sum(dset.ana[, numericList[i]], na.rm = T)>0){
    boxplot(dset.ana[, numericList[i]], main=colnames(dset.ana)[numericList[i]], cex.main=0.5)
  }
}
dev.off()

pdf(file=paste0(outdir, "/qq.ana.set.imputed.pdf"), height = 4, width = 4)
for(i in 1:length(numericList)){
  if(sum(dset.ana[, numericList[i]], na.rm = T)>0){
    qqnorm(scale(dset.ana[, numericList[i]]), main=colnames(dset.ana)[numericList[i]], cex.main=0.5)
    qqline(scale(dset.ana[, numericList[i]]), col = "steelblue", lwd = 2)
  }
}
dev.off()

pdf(file=paste0(outdir, "/barplot.ana.set.imputed.pdf"), height = 6, width = 4)
for(i in 1:length(factorList)){
  tab<-table(dset.ana[,factorList[i]])
  par(mar=c(10,1,1,1))
  barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset.ana)[factorList[i]], cex.main=0.5, cex.names = 0.5)
}
dev.off()

pdf(file=paste0(outdir, "/integers.ana.set.imputed.pdf"), height = 6, width = 4)
for(i in 1:length(integerList)){
  tab<-table(dset.ana[,integerList[i]])
  par(mar=c(10,1,1,1))
  barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset.ana)[integerList[i]], cex.main=0.5, cex.names = 0.5)
}
dev.off()

nm<-nmiss(dset.ana)

pdf(file=paste0(outdir, "/missing.ana.set.imputed.pdf"), height = 5, width = 0.2*dim(nm)[1])
barplot(nm$pm, names.arg = substr(rownames(nm), start=1, stop=10), las=2, width = 0.2 ,main="Missing proportions", cex.main=0.5, cex.names = 0.5)
dev.off()


# WRITE OUT DATASET ----

write.csv2(dset.ana.imp.1, file=paste0(outdir,"/dset.ana.imp.1.csv"), row.names = F)
write.csv2(dset.ana.imp.2, file=paste0(outdir,"/dset.ana.imp.2.csv"), row.names = F)
write.csv2(dset.ana.imp.3, file=paste0(outdir,"/dset.ana.imp.3.csv"), row.names = F)
write.csv2(dset.ana.imp.4, file=paste0(outdir,"/dset.ana.imp.4.csv"), row.names = F)
write.csv2(dset.ana.imp.5, file=paste0(outdir,"/dset.ana.imp.5.csv"), row.names = F)
