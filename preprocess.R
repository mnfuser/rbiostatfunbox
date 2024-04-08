# PRE-PROCESSING v4

rm(list=ls())

# Sources ----
source("D:/codice/rbiostatfunbox/funbox.R")

# Packages ----
library(readxl)

# PATH ----

wd<-"H:/Drive condivisi/Lavori in corso/VitC"
setwd(wd)
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/preprocess")

# GLOBAL VARS ----
# LOAD DATASET RAW ----
dset <- read.csv2(paste0(indir,"/dataset.csv"),
                    stringsAsFactors=T)
View(dset)

# SET VARTYPE ----
## Metadata - vartype

varclass<-as.matrix(sapply(dset, class))
if(!"vartype.csv" %in% list.files(indir)){
  vartype<-cbind(colnames(dset), varclass)
  colnames(vartype)<-c("var", "type")
  write.csv2(vartype, file=paste0(outdir, "/vartype.csv"), row.names = F)
}
View(varclass)

## Apply vartype
readVarType(paste0(indir, "/vartype.csv"))
apllyVarType(integerList)

# OUTLIERS CORRECTION ----

getNumbersFromDataset(dset)
for(i in numericList){
  dset[, i]<-ApplyCap2Outliers(dset[, i])
}

# RECODE / BUILD NEW VARS ----

## New vars

dset$vitc_postop_gtoreq_23<-ifelse(dset$Vitamin.C..Post.op.>22,1,0)
dset$vitc_postop_ltoreq_23<-ifelse(dset$Vitamin.C..Post.op.<24,1,0)
dset$vitc_delta_perc<-(dset$Vitamin.C..Post.op.-dset$Vitamin.C..Pre.op.)/dset$Vitamin.C..Pre.op.*100
dset$vitc_delta_perc_gtoreq_50<-ifelse(dset$vitc_delta_perc<(-49),1,0)

## CENSORING SETTINGS
# QC PLOTS ----

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

# DATASET METADATA ----

nm<-nmiss(dset, fileName = paste0(outdir,"/nm.dset.csv"))

# WRITE OUT DATASET ----

write.csv2(dset, file=paste0(outdir,"/dset.ana.csv"), row.names = F)
write.csv2(nm, file=paste0(outdir,"/nm.dset.ana.filtered.csv"), row.names = F)