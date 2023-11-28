#     QCTOOLBOX
#
#     Tools to perform data quality check
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

# Source ---------------------------------------------
source("/Users/mmanfrini/Code/rbiostatfunbox/rbiostatfunbox/funbox.R")

# PATH ----

setwd("/Users/mmanfrini/Analisi/maritati/PVL2")
wd<-getwd()
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/preprocess")

# LOAD DATASET ----

dset <- read.csv2(paste0(indir,"/dataset.csv"),
                    stringsAsFactors=T)
View(dset)

# Load dataset dictionary file

ddf<-read.csv2(paste0(indir,"/data_dictionary.csv"),
              stringsAsFactors=T)
View(ddf)

# SET VAR TYPE ----

numericList<-which(ddf$VARTYPE=="Numeric")
factorList<-which(ddf$VARTYPE=="Factor")
stringList<-which(ddf$VARTYPE=="String")
integerList<-which(ddf$VARTYPE=="Integer")
orderedList<-which(ddf$VARTYPE=="Ordered")
logicalList<-which(ddf$VARTYPE=="Logical")
dateList<-which(ddf$VARTYPE=="Date")

print(paste("Found ",  length(numericList), " numerical variables"))
print(paste("Found ",  length(factorList), " factor variables"))
print(paste("Found ",  length(stringList), " string variables"))
print(paste("Found ",  length(integerList), " integer variables"))
print(paste("Found ",  length(orderedList), " ordered variables"))
print(paste("Found ",  length(logicalList), " logical variables"))
print(paste("Found ",  length(dateList), " date variables"))

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

# continuous
## boxplot
if(length(numericList)>0){
  pdf(file=paste0(outdir, "/boxplot.pdf"), height = 4, width = 4)
  for(i in 1:length(numericList)){
    if(sum(dset[, numericList[i]], na.rm = T)>0){
      boxplot(dset[, numericList[i]], main=colnames(dset)[numericList[i]], cex.main=0.5)
    }
  }
  dev.off()

  ## qqplot
  pdf(file=paste0(outdir, "/qq.pdf"), height = 4, width = 4)
  for(i in 1:length(numericList)){
    if(sum(dset[, numericList[i]], na.rm = T)>0){
      qqnorm(scale(dset[, numericList[i]]), main=colnames(dset)[numericList[i]], cex.main=0.5)
      qqline(scale(dset[, numericList[i]]), col = "steelblue", lwd = 2)
    }
  }
  dev.off()
}

# factor
## barplot
if(length(factorList)>0){
  pdf(file=paste0(outdir, "/barplot.pdf"), height = 6, width = 4)
  for(i in 1:length(factorList)){
    tab<-table(dset[,factorList[i]])
    par(mar=c(10,1,1,1))
    barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset)[factorList[i]], cex.main=0.5, cex.names = 0.5)
  }
  dev.off()
}

# integer
## barplot
if(length(integerList)>0){
  pdf(file=paste0(outdir, "/integers.pdf"), height = 6, width = 4)
  for(i in 1:length(integerList)){
    tab<-table(dset[,integerList[i]])
    par(mar=c(10,1,1,1))
    barplot(tab, names.arg = substr(names(tab), start=1, stop=10), las=2, main=colnames(dset)[integerList[i]], cex.main=0.5, cex.names = 0.5)
  }
  dev.off()
}

# missing
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

# WRITE OUT DATASETS ----

## Save

write.csv2(dset, file=paste0(outdir,"/dset.qc.csv"), row.names = F)
write.csv2(nm, file=paste0(outdir,"/ddf.qc.csv"), row.names = F)
