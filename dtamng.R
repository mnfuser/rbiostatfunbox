#     DTAMNG
#
#     Data management
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

setwd("/Users/mmanfrini/Analisi/maritati/PVL2")
wd<-getwd()
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/preprocess")

# LOAD DATASET ----

dset <- read.csv2(paste0(indir,"/dataset.csv"),
                  stringsAsFactors=F)
View(dset)

# Load dataset dictionary file

ddf<-read.csv2(paste0(indir,"/ddf.qc.csv"),
               stringsAsFactors=F)
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
if(length(stringList)>0){
  dset <- applyCharacter(dset, stringList)
} else {
  print("No numeric vars")
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

View(dset)


# RECODE ----

# NEW VARS ----

# FILTER ----

dset<-dset[!is.na(dset$FAILURE),]

# WRITE OUT DATASETS ----

## Save

write.csv2(dset, file=paste0(outdir,"/dset.qc.pp.csv"), row.names = F)
