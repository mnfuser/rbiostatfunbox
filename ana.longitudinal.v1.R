#+ tidy=TRUE, tidy.opts = list(blank = FALSE, width.cutoff = 50), echo=FALSE, message=FALSE, warning=FALSE, comment=NA

# Functions ---------------------------------------------
source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

# PATH ----

wd<-setwd("/Users/mmanfrini/Analisi/maritati/PVL")
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out/descriptive")

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

# OUTLIERS FIX ----

getNumbersFromDataset(dset.ana)

for(i in numericList){
  dset.ana[, i]<-ApplyCap2Outliers(dset.ana[, i])
}
# LINES PLOT ----

## WBC pre-op | discharge ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 2)),
  TIME=c(rep("0",12), rep("1",12)), 
  WBC=c(dset.ana$WBC_pre_op, dset.ana$WBC_discharge),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 2)))
)

s<-(summary(aov(WBC ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=WBC, group = ID)) +
  geom_line() + #aes(color=PVL)
  geom_point(aes(shape = PVL)) + #aes(color=PVL)
  scale_x_discrete(labels=c("Pre-op", "Discharge")) +
  ggtitle("WBC") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  annotate("text", 1.5, 16000, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 43, 53)])


## NEUTROFILI pre-op | discharge ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 2)),
  TIME=c(rep("0",12), rep("1",12)), 
  N=c(dset.ana$neutrophil_pre_op, dset.ana$neutrophil_discharge),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 2)))
)

s<-(summary(aov(N ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=N, group = ID)) +
  geom_line() +
  geom_point(aes(shape = PVL)) +
  scale_x_discrete(labels=c("Pre-op", "Discharge")) +
  ggtitle("Neutrophil") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  annotate("text", 1.5, 16000, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 43, 53)])


## CRP pre-op | discharge ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 2)),
  TIME=c(rep("0",12), rep("1",12)), 
  CRP=c(dset.ana$CRP_pre_op, dset.ana$CRP_discharge),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 2)))
)

s<-(summary(aov(CRP ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=CRP, group = ID)) +
  geom_line() +
  geom_point(aes(shape = PVL)) +
  scale_x_discrete(labels=c("Pre-op", "Discharge")) +
  ggtitle("CRP") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  annotate("text", 1.5, 250, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 40, 54)])

## ESV pre-op | discharge ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 2)),
  TIME=c(rep("0",12), rep("1",12)), 
  ESV=c(dset.ana$ESV_.pre_op, dset.ana$ESV_discharge),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 2)))
)

s<-(summary(aov(ESV ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=ESV, group = ID)) +
  geom_line() +
  geom_point(aes(shape = PVL)) +
  scale_x_discrete(labels=c("Pre-op", "Discharge")) +
  ggtitle("ESV") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  annotate("text", 1.5, 100, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 41, 55)])


## HB pre-op | discharge ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 2)),
  TIME=c(rep("0",12), rep("1",12)), 
  HB=c(dset.ana$HB_pre_op, dset.ana$HB_discharge),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 2)))
)

s<-(summary(aov(HB ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=HB, group = ID)) +
  geom_line() +
  geom_point(aes(shape = PVL)) +
  scale_x_discrete(labels=c("Pre-op", "Discharge")) +
  ggtitle("HB") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  annotate("text", 1.5, 15, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 44, 56)])


## PROTEINS pre-op | discharge ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 2)),
  TIME=c(rep("0",12), rep("1",12)), 
  PR=c(dset.ana$proteins_pre_op, dset.ana$proteins_discharge),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 2)))
)

s<-(summary(aov(PR ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=PR, group = ID)) +
  geom_line() +
  geom_point(aes(shape = PVL)) +
  scale_x_discrete(labels=c("Pre-op", "Discharge")) +
  ggtitle("Proteins") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  annotate("text", 1.5, 7, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 45, 58)])


## WBC pre-op | discharge 1-2-3 month ind_thr ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  WBC=c(dset.ana$WBC_pre_op, 
        dset.ana$WBC_discharge, 
        dset.ana$WBC_1month_ind_thr, 
        dset.ana$WBC_2month_ind_thr, 
        dset.ana$WBC_3month_ind_thr),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

s<-(summary(aov(WBC ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=WBC, group = ID)) +
  geom_line() +
  geom_point(aes(shape=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("WBC") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        annotate("text", 4, 15000, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

p<-ggplot(dset.line, aes(x=TIME, y=WBC, fill = PVL)) +
  geom_boxplot() +
  scale_fill_grey() +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("WBC") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        annotate("text", 3, 15000, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 42, 52, 64, 68, 72)])


## NEUTROFILI pre-op | discharge 1-2-3 month ind_thr ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  N=c(dset.ana$neutrophil_pre_op, 
        dset.ana$neutrophil_discharge, 
        dset.ana$neutrophil_1month_ind_thr, 
        dset.ana$neutrophil_2month_ind_thr, 
        dset.ana$neutrophil_3month_ind_thr),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

s<-(summary(aov(N ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=N, group = ID)) +
  geom_line() +
  scale_fill_grey() +
  geom_point(aes(shape=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("Neutrophil") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  annotate("text", 4, 15000, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

p<-ggplot(dset.line, aes(x=TIME, y=N, fill=PVL)) +
  geom_boxplot() +
  scale_fill_grey() +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("Neutrophil") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", 4, 15000, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 43, 53, 65, 69, 73)])


## CRP pre-op | discharge 1-2-3 month ind_thr ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  CRP=c(dset.ana$CRP_pre_op, 
      dset.ana$CRP_discharge, 
      dset.ana$CRP_1month_ind_thr, 
      dset.ana$CRP_2month_ind_thr, 
      dset.ana$CRP_3month_ind_thr),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

s<-(summary(aov(CRP ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=CRP, group = ID)) +
  geom_line() +
  geom_point(aes(shape=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("CRP") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        annotate("text", 4, 200, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

p<-ggplot(dset.line, aes(x=TIME, y=CRP, fill=PVL)) +
  geom_boxplot() +
  scale_fill_grey() +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("CRP") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        annotate("text", 4, 200, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 40, 54, 66, 70, 74)])


## ESV pre-op | discharge 1-2-3 month ind_thr ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  ESV=c(dset.ana$ESV_.pre_op, 
        dset.ana$ESV_discharge, 
        dset.ana$ESV_1month_ind_thr, 
        dset.ana$ESV_2month_ind_thr, 
        dset.ana$ESV_3month_ind_thr),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

s<-(summary(aov(ESV ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=ESV, group = ID)) +
  geom_line() +
  geom_point(aes(shape=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("ESV") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        annotate("text", 4, 100, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

p<-ggplot(dset.line, aes(x=TIME, y=ESV, fill=PVL)) +
  geom_boxplot() +
  scale_fill_grey() +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month ind_thr", 
                            "2 month ind_thr", 
                            "3 month ind_thr")) +
  ggtitle("ESV") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        annotate("text", 4, 100, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 41, 55, 67, 71, 75)])



## WBC pre-op | discharge 1-3-6 month fup ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  WBC=c(dset.ana$WBC_pre_op, 
        dset.ana$WBC_discharge, 
        dset.ana$WBC_1month_fup, 
        dset.ana$WBC_3month_fup, 
        dset.ana$WBC_6month_fup),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

s<-(summary(aov(WBC ~ factor(TIME)*factor(PVL), data = dset.line)))
s

p<-ggplot(dset.line, aes(x=TIME, y=WBC, group = ID)) +
  geom_line() +
  geom_point(aes(shape=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "3 month fup", 
                            "6 month fup")) +
  ggtitle("WBC") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        annotate("text", 4, 100, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

p<-ggplot(dset.line, aes(x=TIME, y=WBC, fill=PVL)) +
  geom_boxplot() +
  scale_fill_grey() +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "3 month fup", 
                            "6 month fup")) +
  ggtitle("WBC") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        annotate("text", 4, 100, label=paste0("p=", round(s[[1]][5][3,1], 2)))
p

knitr::kable(dset.ana[, c(1, 3, 42, 52, 77, 82, 87)])





## NEUTROFILI pre-op | discharge 1-3-6 month fup ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  N=c(dset.ana$neutrophil_pre_op, 
      dset.ana$neutrophil_discharge, 
      dset.ana$neutrophil_1month_fup, 
      dset.ana$neutrophil_3month_fup, 
      dset.ana$neutrophil_6month_fup),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

p<-ggplot(dset.line, aes(x=TIME, y=N, group = ID)) +
  geom_line(aes(color=PVL)) +
  geom_point(aes(color=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "3 month fup", 
                            "6 month fup")) +
  ggtitle("Neutrophil") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

p<-ggplot(dset.line, aes(x=TIME, y=N)) +
  geom_boxplot(aes(color=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "2 month fup", 
                            "3 month fup")) +
  ggtitle("Neutrophil") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

knitr::kable(dset.ana[, c(1, 3, 43, 53, 78, 83, 88)])


## CRP pre-op | discharge 1-3-6 month fup ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  CRP=c(dset.ana$CRP_pre_op, 
        dset.ana$CRP_discharge, 
        dset.ana$CRP_1month_fup, 
        dset.ana$CRP_3month_fup, 
        dset.ana$CRP_6month_fup),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

p<-ggplot(dset.line, aes(x=TIME, y=CRP, group = ID)) +
  geom_line(aes(color=PVL)) +
  geom_point(aes(color=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "3 month fup", 
                            "6 month fup")) +
  ggtitle("CRP") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

p<-ggplot(dset.line, aes(x=TIME, y=CRP)) +
  geom_boxplot(aes(color=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "2 month fup", 
                            "3 month fup")) +
  ggtitle("CRP") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

knitr::kable(dset.ana[, c(1, 3, 40, 54, 79, 84, 89)])


## ESV pre-op | discharge 1-3-6 month fup ----

dset.line<-data.frame(
  ID=c(rep(dset.ana$id, 5)),
  TIME=c(rep("0",12), rep("1",12), rep("2",12), rep("3",12), rep("4",12)), 
  ESV=c(dset.ana$ESV_.pre_op, 
        dset.ana$ESV_discharge, 
        dset.ana$ESV_1month_fup, 
        dset.ana$ESV_3month_fup, 
        dset.ana$ESV_6month_fup),
  PVL=factor(c(rep(as.character(dset.ana$pvl_test), 5)))
)

p<-ggplot(dset.line, aes(x=TIME, y=ESV, group = ID)) +
  geom_line(aes(color=PVL)) +
  geom_point(aes(color=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "3 month fup", 
                            "6 month fup")) +
  ggtitle("ESV") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

p<-ggplot(dset.line, aes(x=TIME, y=ESV)) +
  geom_boxplot(aes(color=PVL)) +
  scale_x_discrete(labels=c("Pre-op", 
                            "Discharge", 
                            "1 month fup", 
                            "3 month fup", 
                            "6 month fup")) +
  ggtitle("ESV") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

knitr::kable(dset.ana[, c(1, 3, 41, 55, 80, 85, 90)])




# LONGITUDINAL DATA ANALYSIS GLS ----

## GLS DATASET

#dset.gls<-dplyr::select(dset.ana, ... )

# dset.gls<-dset.ana

# WIDE 2 LONG

# dset.gls<-reshape(dset.gls, 
#                   varying=8:25, 
#                   direction = "long", 
#                   idvar = "ID", 
#                   timevar = "TIME")

# TREAT FIRST TIME POINT AS COVAR

# baseline<-subset(data.frame(dset.gls), TIME == 1, -TIME )
# baseline<-upData(baseline, rename=c(AA = "AA.0", ADP5 ="ADP5.0", ADP20 = "ADP20.0", TRAP = "TRAP.0", sCD40L = "sCD40L.0", PSEL = "PSEL.0"), print =FALSE )
# 
# followup<-subset(data.frame(dset.gls), TIME > 1, c(ID ,TIME , AA, ADP5, ADP20, TRAP, sCD40L, PSEL))
# 
# dset.gls.complete<-merge(baseline, followup, by="ID" )

# PLOTS
# 
# # AA
# 
# ggplot(dset.gls, aes(x=TIME, y=AA, color=factor(ID))) +
#   geom_line() + facet_grid(~MORTO.vs.VIVO) +
#   theme(legend.position = "none")
# 
# # ADP5
# 
# ggplot(dset.gls, aes(x=TIME, y=ADP5, color=factor(ID))) +
#   geom_line() + facet_grid(~MORTO.vs.VIVO) +
#   theme(legend.position = "none")
# 
# # ADP20
# 
# ggplot(dset.gls, aes(x=TIME, y=ADP20, color=factor(ID))) +
#   geom_line() + facet_grid(~MORTO.vs.VIVO) +
#   theme(legend.position = "none")
# 
# # TRAP
# 
# ggplot(dset.gls, aes(x=TIME, y=TRAP, color=factor(ID))) +
#   geom_line() + facet_grid(~MORTO.vs.VIVO) +
#   theme(legend.position = "none")
# 
# # sCD40L
# 
# ggplot(dset.gls, aes(x=TIME, y=sCD40L, color=factor(ID))) +
#   geom_line() + facet_grid(~MORTO.vs.VIVO) +
#   theme(legend.position = "none")
# 
# # PSEL
# 
# ggplot(dset.gls, aes(x=TIME, y=PSEL, color=factor(ID))) +
#   geom_line() + facet_grid(~MORTO.vs.VIVO) +
#   theme(legend.position = "none")
#  

# MODELS

# ddist<-datadist(dset.gls.complete)
# options(datadist='ddist')
# 
# require(nlme)

# AA

# cp<-list(corCAR1, corExp, corCompSymm, corLin, corGaus, corSpher)
# 
# z<-vector('list', length (cp))
# 
# for(k in 1:length(cp)){
#   
#   z[[k]]<-gls(AA ~ MORTO.vs.VIVO*TIME,
#               data=dset.gls,
#               correlation=cp[[k]](form=~TIME|ID),
#               na.action = na.exclude
#               )
# }
# 
# anova(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]], z[[6]])
# 
# fit<-gls(AA ~ MORTO.vs.VIVO*TIME +
#          AA.0,
#          data=dset.gls.complete,
#          correlation=corLin(form=~TIME|ID),
#          na.action = na.omit)
# 
# print(fit, coeff=T)
# plot(anova(fit))
# summary(fit)
# 
# emm <- emmeans(fit, ~ MORTO.vs.VIVO*TIME)
# ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
# write.csv2(ct, file=paste0(outdir, "/AA.ct.death.csv"))
# 
# v<-Variogram(fit, distance = corCAR1)
# plot(v)

# PLOT TIME EFFECT

# ggplot(Predict(fit, TIME , MORTO.vs.VIVO , conf.int = F ),
#         adj.subtitle =FALSE , legend.position = ' top ' ) +
#   geom_line()

# CONTRASTS

# k1<-rms::contrast(fit, list(TIME=c(1,2,3), MORTO.vs.VIVO = "1" ),
#                   list(TIME=c(1,2,3), MORTO.vs.VIVO = "0" ))
# 
# print(k1)
# 
# k1<-as.data.frame (k1[c("TIME", "Contrast", "Lower", "Upper")])
# 
# ggplot(Predict(fit, TIME=c(1,2,3) ,  MORTO.vs.VIVO, conf.int =FALSE ),
#         adj.subtitle =FALSE , legend.position ="top" ) +
#   geom_point() +
#   geom_line()
# 
# 
# jpeg(file=paste0(outdir, "/AA.DEATH.EFFECTS.jpg"), height=4, width=4, unit="in", res=600)
# dat <- ggpredict(fit, terms = c("TIME", "MORTO.vs.VIVO"))
# p<-plot(dat,connect.lines = TRUE)
# p<-p+ggtitle("AA MARGINAL EFFECTS PLOT")
# p<-p+annotate("text", x = 1, y = 1, label = paste0("p=", round(summary(fit)$coefficients[1])), size = 3)
# print(p)
# dev.off()

# LONGITUDINAL DATA ANALYSIS RANDOM EFFECTS MIXED MODEL ----

library(ggeffects)
library(emmeans)

library(lmerTest)

## WBC : B : MLE ----

rnd<-c(`lbl#1` = "RND=1", 
       `lbl#2` = "RND=2")

rnd_labeller <- function(variable,value){
  return(rnd[value])
}

## LMER DATASET

dset.lmer<-cbind(ID=c(1:dim(dset.ana)[1]),
           dset.ana[, c(4,  #RND
                        77, #MEL0
                        78, #MEL1
                        79, #WBC0
                        86, #WBC1
                        80, #PLT0
                        87, #PLT1
                        81, #N0
                        88, #N1
                        82, #L0
                        89, #L1
                        83, #M0
                        90, #M1
                        84, #E0
                        91, #E1
                        85, #B0
                        92 #B1
                        )])

# WIDE 2 LONG



# STDZN

#dset.lmer.long[, c(4:11)]<-scale(dset.lmer.long[, c(4:11)], 
#                                 center = F,
#                                 scale = apply(dset.lmer.long[, c(4:11)], 2, sd, na.rm = TRUE))

dset.lmer.long[, c(4:11)]<-scale(dset.lmer.long[, c(4:11)], center = T, scale = T)

# CONTRASTS

dset.lmer.long$TIME<-as.factor(dset.lmer.long$TIME)

contrasts(dset.lmer.long$TIME) <- rbind(
  cbind(" 2:1" = c(0, 1))
)

contrasts(dset.lmer.long$randomization) <- rbind(
  cbind(" 2:1" = c(0, 1))
)

### WBC ----

fit<-lmer(WBC ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/WBC.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/WBC.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/WBC.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/WBC.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/WBC.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) + 
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("WBC")
p<-p+annotate("text", x = 1, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "WBC"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "WBC"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "WBC"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "WBC"], 
                paired = TRUE, 
                alternative = "two.sided")

pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

f_labels <- data.frame(randomization = c(1, 2), label = c(pvalue.t1, pvalue.t2))

jpeg(file=paste0(outdir, "/WBC.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = WBC))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("L (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()

### PLT ----

fit<-lmer(PLT ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/PLT.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/PLT.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/PLT.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/PLT.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/PLT.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("PLT")
p<-p+annotate("text", x = 1, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "PLT"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "PLT"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "PLT"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "PLT"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

f_labels <- data.frame(randomization = c(1, 2), label = c(pvalue.t1, pvalue.t2))

jpeg(file=paste0(outdir, "/PLT.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = PLT))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("PLT (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()

### N ----

fit<-lmer(N ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/N.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/N.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/N.overall.ct.csv"))
pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/N.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/N.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("N")
p<-p+annotate("text", x = 1, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "N"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "N"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "N"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "N"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

f_labels <- data.frame(randomization = c("1", "2"), label = c(pvalue.t1, pvalue.t2))


jpeg(file=paste0(outdir, "/N.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = N))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("N (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()

### L ----

fit<-lmer(L ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/L.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/L.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/L.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/L.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/L.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("L")
p<-p+annotate("text", x = 1, y = 0.5, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "L"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "L"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "L"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "L"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

f_labels <- data.frame(randomization = c("1", "2"), label = c(pvalue.t1, pvalue.t2))

jpeg(file=paste0(outdir, "/L.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = L))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("L (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()


### M ----

fit<-lmer(M ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/M.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/M.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/M.overall.ct.csv"))
pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/M.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/M.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("M")
p<-p+annotate("text", x = 1, y = 0.5, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "M"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "M"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "M"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "M"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

f_labels <- data.frame(id = c("RND1", "RND2"), label = c(pvalue.t1, pvalue.t2))

jpeg(file=paste0(outdir, "/M.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = M))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("M (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()


### E ----

fit<-lmer(E ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/E.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/E.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/E.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/E.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/E.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("E")
p<-p+annotate("text", x = 1, y = 0.5, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "E"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "E"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "E"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "E"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

f_labels <- data.frame(id = c("1", "2"), label = c(pvalue.t1, pvalue.t2))

jpeg(file=paste0(outdir, "/E.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = E))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("E (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()


### B ----

fit<-lmer(B ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/B.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/B.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/B.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/B.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/B.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("B")
p<-p+annotate("text", x = 1, y = 0.01, label = pvalue, size = 3)
print(p)
dev.off()


## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "B"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "B"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "B"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "B"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p=0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p=0.001", paste0("p=",round(w2$p.value,3)))

f_labels <- data.frame(id = c("RND1", "RND2"), label = c(pvalue.t1, pvalue.t2))

jpeg(file=paste0(outdir, "/B.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = B))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("B (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()


### MLE ----

fit<-lmer(MEL ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/MLE.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/MLE.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/MLE.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[4,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/MLE.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/MLE.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("MLE")
p<-p+annotate("text", x = 1, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "MEL"], 
               dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "MEL"], 
               paired = TRUE, 
               alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "MEL"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "MEL"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",w1$p.value))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",w2$p.value))

f_labels <- data.frame(id = c("1", "2"), label = c(pvalue.t1, pvalue.t2))

jpeg(file=paste0(outdir, "/MLE.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = MEL))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("MLE (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()

## LAB ----

## LMER DATASET

dset.lmer<-cbind(ID=c(1:dim(dset.ana)[1]),
                 dset.ana[, c(4,  #RND
                              56, #IL6.0
                              57, #IL6.1
                              58, #IL6.2
                              59, #FABP2.0
                              60, #FABP2.1
                              61, #FABP2.2
                              62, #CystatinC.0
                              63, #CystatinC.1
                              64, #CystatinC.2
                              65, #TNFalfa.0
                              66, #TNFalfa.1
                              67, #TNFalfa.2
                              68, #IL8.0
                              69, #IL8.1
                              70, #IL8.2
                              71, #LIPOCALIN2.0
                              72, #LIPOCALIN2.1
                              73, #LIPOCALIN2.2
                              74, #Enolase.0
                              75, #Enolase.1
                              76 #Enolase.2
                 )])

# WIDE 2 LONG

dset.lmer.long<-reshape(dset.lmer, 
                   varying=c(3:23),
                   direction = "long", 
                   #idvar = "id", 
                   timevar = "TIME")

dset.lmer.long$TIME<-as.factor(dset.lmer.long$TIME)

# STDZN

dset.lmer.long[, c(4:10)]<-scale(dset.lmer.long[, c(4:10)])

# CONTRASTS

contrasts(dset.lmer.long$TIME) <- rbind(
  cbind(" 2:1" = c(0, 1, 0), " 3:1" = c(0, 0, 1))
)

contrasts(dset.lmer.long$randomization) <- rbind(
  cbind(" 2:1" = c(0, 1))
)

### IL6 ----

fit<-lmer(IL6 ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/IL6.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/IL6.overall.confint.csv"))

#emm <- emmeans(fit, ~ TIME*randomization)
emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/IL6.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[6,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/IL6.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/IL6.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("IL6")
p<-p+annotate("text", x = 2, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "IL6"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "IL6"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "IL6"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "IL6"], 
                paired = TRUE, 
                alternative = "two.sided")

w3<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "IL6"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==2), "IL6"], 
                paired = TRUE, 
                alternative = "two.sided")

w4<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "IL6"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==2), "IL6"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

pvalue.t3=ifelse(w3$p.value<0.001, "p<0.001", paste0("p=",round(w3$p.value,3)))
pvalue.t4=ifelse(w4$p.value<0.001, "p<0.001", paste0("p=",round(w4$p.value,3)))


f_labels <- data.frame(id = c("1", "2"), label = c(pvalue.t3, pvalue.t4))


jpeg(file=paste0(outdir, "/IL6.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = IL6))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("IL6 (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()


### IL8 ----

fit<-lmer(IL8 ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/IL8.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/IL8.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/IL8.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[6,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/IL8.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/IL8.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("IL8")
p<-p+annotate("text", x = 2, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "IL8"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "IL8"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "IL8"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "IL8"], 
                paired = TRUE, 
                alternative = "two.sided")

w3<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "IL8"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==2), "IL8"], 
                paired = TRUE, 
                alternative = "two.sided")

w4<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "IL8"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==2), "IL8"], 
                paired = TRUE, 
                alternative = "two.sided")


pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

pvalue.t3=ifelse(w3$p.value<0.001, "p<0.001", paste0("p=",round(w3$p.value,3)))
pvalue.t4=ifelse(w4$p.value<0.001, "p<0.001", paste0("p=",round(w4$p.value,3)))


f_labels <- data.frame(id = c("1", "2"), label = c(pvalue.t3, pvalue.t4))


jpeg(file=paste0(outdir, "/IL8.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = IL8))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("IL8 (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()

### TNF ----

fit<-lmer(TNFalfa ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/TNF.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/TNF.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/TNF.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[6,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/TNF.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/TNF.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("TNF")
p<-p+annotate("text", x = 2, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "TNFalfa"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "TNFalfa"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "TNFalfa"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "TNFalfa"], 
                paired = TRUE, 
                alternative = "two.sided")

w3<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "TNFalfa"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==2), "TNFalfa"], 
                paired = TRUE, 
                alternative = "two.sided")

w4<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "TNFalfa"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==2), "TNFalfa"], 
                paired = TRUE, 
                alternative = "two.sided")

pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

pvalue.t3=ifelse(w3$p.value<0.001, "p<0.001", paste0("p=",round(w3$p.value,3)))
pvalue.t4=ifelse(w4$p.value<0.001, "p<0.001", paste0("p=",round(w4$p.value,3)))


f_labels <- data.frame(randomization = c("1", "2"), label = c(pvalue.t3, pvalue.t4))


jpeg(file=paste0(outdir, "/TNF.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = TNFalfa))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("TNF alpha (z-score)") +
  geom_text(x = 2, y = 2.2, aes(label = label), data = f_labels, ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()



### LIPO ----

fit<-lmer(LIPOCALIN2 ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/LIPO.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/LIPO.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/LIPO.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[6,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/LIPO.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/LIPO.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("LIPO")
p<-p+annotate("text", x = 2, y = 1, label = pvalue, size = 3)
print(p)
dev.off()


## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "LIPOCALIN2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "LIPOCALIN2"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "LIPOCALIN2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "LIPOCALIN2"], 
                paired = TRUE, 
                alternative = "two.sided")

w3<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "LIPOCALIN2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==2), "LIPOCALIN2"], 
                paired = TRUE, 
                alternative = "two.sided")

w4<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "LIPOCALIN2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==2), "LIPOCALIN2"], 
                paired = TRUE, 
                alternative = "two.sided")

pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

pvalue.t3=ifelse(w3$p.value<0.001, "p<0.001", paste0("p=",round(w3$p.value,3)))
pvalue.t4=ifelse(w4$p.value<0.001, "p<0.001", paste0("p=",round(w4$p.value,3)))


f_labels <- data.frame(time = c(0, 1, 2, 3), label = c(pvalue.t1, pvalue.t2,
                                                             pvalue.t3, pvalue.t4))


jpeg(file=paste0(outdir, "/LIPO.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = LIPOCALIN2))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("LIPOCAIN2 (z-score)") +
  geom_text(x = 1, y = c(2), aes(label = label), data = f_labels[1:2,], ) +
  geom_text(x = 3, y = c(2), aes(label = label), data = f_labels[3:4,], ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()

### FABP ----

fit<-lmer(FABP2 ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/FABP.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/FABP.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/FABP.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[6,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/FABP.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/FABP.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("FABP")
p<-p+annotate("text", x = 1.8, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "FABP2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "FABP2"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "FABP2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "FABP2"], 
                paired = TRUE, 
                alternative = "two.sided")

w3<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "FABP2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==2), "FABP2"], 
                paired = TRUE, 
                alternative = "two.sided")

w4<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "FABP2"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==2), "FABP2"], 
                paired = TRUE, 
                alternative = "two.sided")

pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

pvalue.t3=ifelse(w3$p.value<0.001, "p<0.001", paste0("p=",round(w3$p.value,3)))
pvalue.t4=ifelse(w4$p.value<0.001, "p<0.001", paste0("p=",round(w4$p.value,3)))


f_labels <- data.frame(time = c(0, 1, 2, 3), label = c(pvalue.t1, pvalue.t2,
                                                       pvalue.t3, pvalue.t4))


jpeg(file=paste0(outdir, "/FABP2.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = FABP2))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("FABP2 (z-score)") +
  geom_text(x = 1, y = c(2), aes(label = label), data = f_labels[1:2,], ) +
  geom_text(x = 3, y = c(2), aes(label = label), data = f_labels[3:4,], ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()

### ENO ----

fit<-lmer(Enolase ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/ENO.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/ENO.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/ENO.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[6,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/ENO.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/ENO.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("ENO")
p<-p+annotate("text", x = 1.8, y = 1, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "Enolase"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "Enolase"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "Enolase"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "Enolase"], 
                paired = TRUE, 
                alternative = "two.sided")

w3<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "Enolase"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==2), "Enolase"], 
                paired = TRUE, 
                alternative = "two.sided")

w4<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "Enolase"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==2), "Enolase"], 
                paired = TRUE, 
                alternative = "two.sided")

pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

pvalue.t3=ifelse(w3$p.value<0.001, "p<0.001", paste0("p=",round(w3$p.value,3)))
pvalue.t4=ifelse(w4$p.value<0.001, "p<0.001", paste0("p=",round(w4$p.value,3)))


f_labels <- data.frame(time = c(0, 1, 2, 3), label = c(pvalue.t1, pvalue.t2,
                                                       pvalue.t3, pvalue.t4))


jpeg(file=paste0(outdir, "/Enolase.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = Enolase))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("Enolase (z-score)") +
  geom_text(x=1, y=2, aes(label = label), data = f_labels[1:4,], ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()


### CIST ----

fit<-lmer(CystatinC ~ randomization*TIME +(1|id), data=dset.lmer.long)

summary(fit)

write.csv2(summary(fit)$coefficients, file=paste0(outdir, "/CIST.lmer.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/CIST.overall.confint.csv"))

emm <- emmeans(fit, ~ TIME*randomization)
ct<-as.data.frame(emmeans::contrast(emm, interaction = "pairwise"))
write.csv2(ct, file=paste0(outdir, "/CIST.overall.ct.csv"))

pvalue=round(summary(fit)$coefficients[6,5],3)
pvalue=ifelse(pvalue<0.001, "p<0.001", paste0("p=",pvalue))
jpeg(file=paste0(outdir, "/CIST.jpg"), height=4, width=6, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "randomization"))
write.csv2(dat, file=paste0(outdir, "/CIST.overall.marginal.csv"))
p<-plot(dat,connect.lines = TRUE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2)
p<-p+ggtitle("CIST")
p<-p+annotate("text", x = 1.8, y = 0.6, label = pvalue, size = 3)
print(p)
dev.off()

## line plot

w1<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==0), "CystatinC"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "CystatinC"], 
                paired = TRUE, 
                alternative = "two.sided")

w2<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==0), "CystatinC"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "CystatinC"], 
                paired = TRUE, 
                alternative = "two.sided")

w3<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==1), "CystatinC"], 
                dset.lmer.long[which(dset.lmer.long$randomization==1 & dset.lmer.long$TIME==2), "CystatinC"], 
                paired = TRUE, 
                alternative = "two.sided")

w4<-wilcox.test(dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==1), "CystatinC"], 
                dset.lmer.long[which(dset.lmer.long$randomization==2 & dset.lmer.long$TIME==2), "CystatinC"], 
                paired = TRUE, 
                alternative = "two.sided")

pvalue.t1=ifelse(w1$p.value<0.001, "p<0.001", paste0("p=",round(w1$p.value,3)))
pvalue.t2=ifelse(w2$p.value<0.001, "p<0.001", paste0("p=",round(w2$p.value,3)))

pvalue.t3=ifelse(w3$p.value<0.001, "p<0.001", paste0("p=",round(w3$p.value,3)))
pvalue.t4=ifelse(w4$p.value<0.001, "p<0.001", paste0("p=",round(w4$p.value,3)))


f_labels <- data.frame(time = c(0, 1, 2, 3), label = c(pvalue.t1, pvalue.t2,
                                                       pvalue.t3, pvalue.t4))


jpeg(file=paste0(outdir, "/CystatinC.lp.jpg"), height=4, width=6, unit="in", res=600)
p <- ggplot(data = dset.lmer.long, aes(x = TIME, y = CystatinC))
p <- p + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.10) +
  facet_wrap(.~randomization) +
  theme_linedraw() +
  ylab("CystatinC (z-score)") +
  geom_text(x=1, y=2, aes(label = label), data = f_labels[1:4,], ) +
  ggtitle("Randomization group - 1 CTRL, 2 TREATMENT")
print(p)
dev.off()



### WILCOXON ----

dset.delta<-dset.ana[, c(4, 56:92)]

# STDZN

# dset.delta[, c(2:38)]<-scale(dset.delta[, c(2:38)])

# DELTA

dset.delta$IL6.d10<-(dset.delta$IL6.1-dset.delta$IL6.0)/dset.delta$IL6.0
dset.delta$IL6.d21<-(dset.delta$IL6.2-dset.delta$IL6.1)/dset.delta$IL6.1
dset.delta$IL6.d20<-(dset.delta$IL6.2-dset.delta$IL6.0)/dset.delta$IL6.0

dset.delta$IL8.d10<-(dset.delta$IL8.1-dset.delta$IL8.0)/dset.delta$IL8.0
dset.delta$IL8.d21<-(dset.delta$IL8.2-dset.delta$IL8.1)/dset.delta$IL8.1
dset.delta$IL8.d20<-(dset.delta$IL8.2-dset.delta$IL8.0)/dset.delta$IL8.0

dset.delta$FABP2.d10<-(dset.delta$FABP2.1-dset.delta$FABP2.0)/dset.delta$FABP2.0
dset.delta$FABP2.d21<-(dset.delta$FABP2.2-dset.delta$FABP2.1)/dset.delta$FABP2.1
dset.delta$FABP2.d20<-(dset.delta$FABP2.2-dset.delta$FABP2.0)/dset.delta$FABP2.0

dset.delta$CystatinC.d10<-(dset.delta$CystatinC.1-dset.delta$CystatinC.0)/dset.delta$CystatinC.0
dset.delta$CystatinC.d21<-(dset.delta$CystatinC.2-dset.delta$CystatinC.1)/dset.delta$CystatinC.1
dset.delta$CystatinC.d20<-(dset.delta$CystatinC.2-dset.delta$CystatinC.0)/dset.delta$CystatinC.0

dset.delta$TNFalfa.d10<-(dset.delta$TNFalfa.1-dset.delta$TNFalfa.0)/dset.delta$TNFalfa.0
dset.delta$TNFalfa.d21<-(dset.delta$TNFalfa.2-dset.delta$TNFalfa.1)/dset.delta$TNFalfa.1
dset.delta$TNFalfa.d20<-(dset.delta$TNFalfa.2-dset.delta$TNFalfa.0)/dset.delta$TNFalfa.0

dset.delta$LIPOCALIN2.d10<-(dset.delta$LIPOCALIN2.1-dset.delta$LIPOCALIN2.0)/dset.delta$LIPOCALIN2.0
dset.delta$LIPOCALIN2.d21<-(dset.delta$LIPOCALIN2.2-dset.delta$LIPOCALIN2.1)/dset.delta$LIPOCALIN2.1
dset.delta$LIPOCALIN2.d20<-(dset.delta$LIPOCALIN2.2-dset.delta$LIPOCALIN2.0)/dset.delta$LIPOCALIN2.0

dset.delta$Enolase.d10<-(dset.delta$Enolase.1-dset.delta$Enolase.0)/dset.delta$Enolase.0
dset.delta$Enolase2.d21<-(dset.delta$Enolase.2-dset.delta$Enolase.1)/dset.delta$Enolase.1
dset.delta$Enolase2.d20<-(dset.delta$Enolase.2-dset.delta$Enolase.0)/dset.delta$Enolase.0

dset.delta$MEL.d10<-(dset.delta$MEL.1-dset.delta$MEL.0)/dset.delta$MEL.0

dset.delta$WBC.d10<-(dset.delta$WBC.1-dset.delta$WBC.0)/dset.delta$WBC.0

dset.delta$PLT.d10<-(dset.delta$PLT.1-dset.delta$PLT.0)/dset.delta$PLT.0

dset.delta$N.d10<-(dset.delta$N.1-dset.delta$N.0)/dset.delta$N.0

dset.delta$L.d10<-(dset.delta$L.1-dset.delta$L.0)/dset.delta$L.0

dset.delta$M.d10<-(dset.delta$M.1-dset.delta$M.0)/dset.delta$M.0

dset.delta$E.d10<-(dset.delta$E.1-dset.delta$E.0)/dset.delta$E.0

dset.delta$B.d10<-(dset.delta$B.1-dset.delta$B.0)/dset.delta$B.0


## DELTA
### WBC

wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "WBC.d10"],
            dset.delta[which(dset.delta$randomization==2), "WBC.d10"],
            paired = F)

print(wt)
boxplot(WBC.d10 ~ randomization, data = dset.delta, main = "WBC t0-t1", outline = F, 
        ylab = "Relative change in WBC",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.5, round(wt$p.value, 3))


### MLE
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "MEL.d10"],
            dset.delta[which(dset.delta$randomization==2), "MEL.d10"],
            paired = F)

print(wt)
boxplot(MEL.d10 ~ randomization, data = dset.delta, main = "MEL", outline = F,
        ylab = "Relative change in MLE",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.9, round(wt$p.value, 3))

### PLT
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "PLT.d10"],
            dset.delta[which(dset.delta$randomization==2), "PLT.d10"],
            paired = F)

print(wt)

boxplot(PLT.d10 ~ randomization, data = dset.delta, main = "PLT", outline = F,
        ylab = "Relative change in PLT",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0, round(wt$p.value, 3))

### B
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "B.d10"],
            dset.delta[which(dset.delta$randomization==2), "B.d10"],
            paired = F)
print(wt)
boxplot(B.d10 ~ randomization, data = dset.delta, main = "B", outline = F,
        ylab = "Relative change in B",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.6, y=-0.1, round(wt$p.value, 3))

### E
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "E.d10"],
            dset.delta[which(dset.delta$randomization==2), "E.d10"],
            paired = F)
print(wt)
boxplot(E.d10 ~ randomization, data = dset.delta, main = "E", outline = F,
        ylab = "Relative change in E",
        xlab = "Randomization group: 1=ctrl, 2=exp",
        xlim(-1, 0.4)
        )
text(x=2, y=-0.1, round(wt$p.value, 3))

### L
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "L.d10"],
            dset.delta[which(dset.delta$randomization==2), "L.d10"],
            paired = F)
print(wt)
boxplot(L.d10 ~ randomization, data = dset.delta, main = "L", outline = F,
        ylab = "Relative change in L",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.5, round(wt$p.value, 3))

### N
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "N.d10"],
            dset.delta[which(dset.delta$randomization==2), "N.d10"],
            paired = F)
print(wt)
boxplot(N.d10 ~ randomization, data = dset.delta, main ="N", outline = F,
        ylab = "Relative change in N",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.9, round(wt$p.value, 3))

### M
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "M.d10"],
            dset.delta[which(dset.delta$randomization==2), "M.d10"],
            paired = F)
print(wt)
boxplot(M.d10 ~ randomization, data = dset.delta, main = "M", outline = F,
        ylab = "Relative change in M",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=-0.1, round(wt$p.value, 3))

### IL6.d10
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL6.d10"],
            dset.delta[which(dset.delta$randomization==2), "IL6.d10"],
            paired = F)
print(wt)
boxplot(IL6.d10 ~ randomization, data = dset.delta, main = "IL6 t0-t1", outline = F,
        ylab = "Relative change in IL6",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=120, round(wt$p.value, 3))

### IL6.d21
wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL6.d21"],
            dset.delta[which(dset.delta$randomization==2), "IL6.d21"],
            paired = F)

boxplot(IL6.d21 ~ randomization, data = dset.delta, main = "IL6 t1-t2", outline = F,
        ylab = "Relative change in IL6",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.2, round(wt$p.value, 3))

### IL6.d20
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL6.d20"],
            dset.delta[which(dset.delta$randomization==2), "IL6.d20"],
            paired = F)
print(wt)
boxplot(IL6.d20 ~ randomization, data = dset.delta, main = "IL6 t0-t2", outline = F,
        ylab = "Relative change in IL6",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=48, round(wt$p.value, 3))

### IL8.d10
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL8.d10"],
            dset.delta[which(dset.delta$randomization==2), "IL8.d10"],
            paired = F)
print(wt)
boxplot(IL8.d10 ~ randomization, data = dset.delta, main = "IL8 t0-t1", outline = F,
        ylab = "Relative change in IL8",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=11.8, round(wt$p.value, 3))

### IL8.d21
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL8.d21"],
            dset.delta[which(dset.delta$randomization==2), "IL8.d21"],
            paired = F)
print(wt)
boxplot(IL8.d21 ~ randomization, data = dset.delta, main = "IL8 t1-t2", outline = F,
        ylab = "Relative change in IL8",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=9, round(wt$p.value, 3))

### IL8.d20
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL8.d20"],
            dset.delta[which(dset.delta$randomization==2), "IL8.d20"],
            paired = F)
print(wt)
boxplot(IL8.d20 ~ randomization, data = dset.delta, main = "IL8 t0-t2", outline = F,
        ylab = "Relative change in IL8",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=15.5, round(wt$p.value, 3))

### LIPOCALIN2.d10
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "LIPOCALIN2.d10"],
            dset.delta[which(dset.delta$randomization==2), "LIPOCALIN2.d10"],
            paired = F)
print(wt)
boxplot(LIPOCALIN2.d10 ~ randomization, data = dset.delta, main = "LIPOCALIN2 t0-t1", outline = F,
        ylab = "Relative change in LIPOCALIN",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=1.5, round(wt$p.value, 3))

### LIPOCALIN2.d21
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "LIPOCALIN2.d21"],
            dset.delta[which(dset.delta$randomization==2), "LIPOCALIN2.d21"],
            paired = F)
print(wt)
boxplot(LIPOCALIN2.d21 ~ randomization, data = dset.delta, main = "LIPOCALIN2 t1-t2", outline = F,
        ylab = "Relative change in LIPOCALIN",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=1.3, round(wt$p.value, 3))

### LIPOCALIN2.d20
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "LIPOCALIN2.d20"],
            dset.delta[which(dset.delta$randomization==2), "LIPOCALIN2.d20"],
            paired = F)
print(wt)
boxplot(LIPOCALIN2.d20 ~ randomization, data = dset.delta, main = "LIPOCALIN2 t0-t2", outline = F,
        ylab = "Relative change in LIPOCALIN",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=1.5, round(wt$p.value, 3))

### Enolase.d10
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "Enolase.d10"],
            dset.delta[which(dset.delta$randomization==2), "Enolase.d10"],
            paired = F)
print(wt)
boxplot(Enolase.d10 ~ randomization, data = dset.delta, main = "Enolase t0-t1", outline = F,
        ylab = "Relative change in Enolase",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=3, round(wt$p.value, 3))

### Enolase.d21
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "Enolase2.d21"],
            dset.delta[which(dset.delta$randomization==2), "Enolase2.d21"],
            paired = F)
print(wt)
boxplot(Enolase2.d21 ~ randomization, data = dset.delta, main = "Enolase t1-t2", outline = F,
        ylab = "Relative change in Enolase",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=1.5, round(wt$p.value, 3))

### Enolase.d20
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "Enolase2.d20"],
            dset.delta[which(dset.delta$randomization==2), "Enolase2.d20"],
            paired = F)
print(wt)
boxplot(Enolase2.d20 ~ randomization, data = dset.delta, main = "Enolase t0-t2", outline = F,
        ylab = "Relative change in Enolase",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=4, round(wt$p.value, 3))

### FABP2.d10
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "FABP2.d10"],
            dset.delta[which(dset.delta$randomization==2), "FABP2.d10"],
            paired = F)
print(wt)
boxplot(FABP2.d10 ~ randomization, data = dset.delta, main = "FABP2.d10 t0-t1", outline = F,
        ylab = "Relative change in FABP2",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=3, round(wt$p.value, 3))

### FABP2.d21
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "FABP2.d21"],
            dset.delta[which(dset.delta$randomization==2), "FABP2.d21"],
            paired = F)
print(wt)
boxplot(FABP2.d21 ~ randomization, data = dset.delta, main = "FABP2 t1-t2", outline = F,
        ylab = "Relative change in FABP2",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.2, round(wt$p.value, 3))

### FABP2.d20
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "FABP2.d20"],
            dset.delta[which(dset.delta$randomization==2), "FABP2.d20"],
            paired = F)
print(wt)
boxplot(FABP2.d20 ~ randomization, data = dset.delta, main = "FABP2 t0-t2", outline = F,
        ylab = "Relative change in FABP2",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.1, round(wt$p.value, 3))

### CystatinC.d10
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "CystatinC.d10"],
            dset.delta[which(dset.delta$randomization==2), "CystatinC.d10"],
            paired = F)
print(wt)
boxplot(CystatinC.d10 ~ randomization, data = dset.delta, main = "CystatinC t0-t1", outline = F,
        ylab = "Relative change in CystatinC",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.3, round(wt$p.value, 3))

### CystatinC.d21
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "CystatinC.d21"],
            dset.delta[which(dset.delta$randomization==2), "CystatinC.d21"],
            paired = F)
print(wt)
boxplot(CystatinC.d21 ~ randomization, data = dset.delta, main = "CystatinC t1-t2", outline = F,
        ylab = "Relative change in CystatinC",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.6, round(wt$p.value, 3))
        

### CystatinC.d21
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "CystatinC.d20"],
            dset.delta[which(dset.delta$randomization==2), "CystatinC.d20"],
            paired = F)
print(wt)
boxplot(CystatinC.d20 ~ randomization, data = dset.delta, main = "CystatinC t0-t2", outline = F,
        ylab = "Relative change in CystatinC",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.6, round(wt$p.value, 3))

### TNFalfa.d10
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "TNFalfa.d10"],
            dset.delta[which(dset.delta$randomization==2), "TNFalfa.d10"],
            paired = F)
print(wt)
boxplot(TNFalfa.d10 ~ randomization, data = dset.delta, main = "TNFalfa t0-t1", outline = F,
        ylab = "Relative change in TNFalfa",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.4, round(wt$p.value, 3))


### TNFalfa.d21
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "TNFalfa.d21"],
            dset.delta[which(dset.delta$randomization==2), "TNFalfa.d21"],
            paired = F)
print(wt)
boxplot(TNFalfa.d21 ~ randomization, data = dset.delta, main = "TNFalfa t1-t2", outline = F,
        ylab = "Relative change in TNFalfa",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.5, round(wt$p.value, 3))
        

### TNFalfa.d20
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "TNFalfa.d20"],
            dset.delta[which(dset.delta$randomization==2), "TNFalfa.d20"],
            paired = F)
print(wt)
boxplot(TNFalfa.d20 ~ randomization, data = dset.delta, main = "TNFalfa t0-t2", outline = F,
        ylab = "Relative change in TNFalfa",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.4, round(wt$p.value, 3))

### TEST ON MEAN DIFFERENCES ----

dset.delta<-dset.ana[, c(4, 56:92)]

# STDZN

dset.delta[, c(2:38)]<-scale(dset.delta[, c(2:38)])

### WBC

wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "WBC.1"],
                dset.delta[which(dset.delta$randomization==2), "WBC.1"],
                paired = F)

print(wt)
boxplot(WBC.1 ~ randomization, data = dset.delta, main = "WBC t0-t1", outline = F, 
        ylab = "WBC",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=8, round(wt$p.value, 3))

### MLE
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "MEL.1"],
                dset.delta[which(dset.delta$randomization==2), "MEL.1"],
                paired = F)

print(wt)
boxplot(MEL.1 ~ randomization, data = dset.delta, main = "MEL", outline = F,
        ylab = "MLE",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=30, round(wt$p.value, 3))

### PLT
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "PLT.1"],
                dset.delta[which(dset.delta$randomization==2), "PLT.1"],
                paired = F)

print(wt)
boxplot(PLT.1 ~ randomization, data = dset.delta, main = "PLT", outline = F,
        ylab = "PLT",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=150, round(wt$p.value, 3))

### B
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "B.1"],
                dset.delta[which(dset.delta$randomization==2), "B.1"],
                paired = F)
print(wt)
boxplot(B.1 ~ randomization, data = dset.delta, main = "B", outline = F,
        ylab = "B",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.6, y=0.8, round(wt$p.value, 3))

### E
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "E.1"],
                dset.delta[which(dset.delta$randomization==2), "E.1"],
                paired = F)
print(wt)
boxplot(E.1 ~ randomization, data = dset.delta, main = "E", outline = F,
        ylab = "E",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.18, round(wt$p.value, 3))

### L
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "L.1"],
                dset.delta[which(dset.delta$randomization==2), "L.1"],
                paired = F)
print(wt)
boxplot(L.1 ~ randomization, data = dset.delta, main = "L", outline = F,
        ylab = "L",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=2.5, round(wt$p.value, 3))

### N
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "N.1"],
                dset.delta[which(dset.delta$randomization==2), "N.1"],
                paired = F)
print(wt)
boxplot(N.1 ~ randomization, data = dset.delta, main ="N", outline = F,
        ylab = "N",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=5.5, round(wt$p.value, 3))

### M
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "M.1"],
                dset.delta[which(dset.delta$randomization==2), "M.1"],
                paired = F)
print(wt)
boxplot(M.1 ~ randomization, data = dset.delta, main = "M", outline = F,
        ylab = "M",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=0.25, round(wt$p.value, 3))

### IL6.1
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL6.1"],
                dset.delta[which(dset.delta$randomization==2), "IL6.1"],
                paired = F)
print(wt)
boxplot(IL6.1 ~ randomization, data = dset.delta, main = "IL6 t1", outline = F,
        ylab = "IL6",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=300, round(wt$p.value, 3))

### IL6.2
wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL6.2"],
            dset.delta[which(dset.delta$randomization==2), "IL6.2"],
            paired = F)

boxplot(IL6.2 ~ randomization, data = dset.delta, main = "IL6 t2", outline = F,
        ylab = "IL6",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=120, round(wt$p.value, 3))

### IL8.1
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL8.1"],
                dset.delta[which(dset.delta$randomization==2), "IL8.1"],
                paired = F)
print(wt)
boxplot(IL8.1 ~ randomization, data = dset.delta, main = "IL8 t1", outline = F,
        ylab = "IL8",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=45, round(wt$p.value, 3))

### IL8.2
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "IL8.2"],
                dset.delta[which(dset.delta$randomization==2), "IL8.2"],
                paired = F)
print(wt)
boxplot(IL8.2 ~ randomization, data = dset.delta, main = "IL8 t2", outline = F,
        ylab = "IL8",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=120, round(wt$p.value, 3))

### LIPOCALIN2.1
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "LIPOCALIN2.1"],
                dset.delta[which(dset.delta$randomization==2), "LIPOCALIN2.1"],
                paired = F)
print(wt)
boxplot(LIPOCALIN2.1 ~ randomization, data = dset.delta, main = "LIPOCALIN2 t1", outline = F,
        ylab = "LIPOCALIN",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=1400, round(wt$p.value, 3))

### LIPOCALIN2.2
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "LIPOCALIN2.2"],
                dset.delta[which(dset.delta$randomization==2), "LIPOCALIN2.2"],
                paired = F)
print(wt)
boxplot(LIPOCALIN2.2 ~ randomization, data = dset.delta, main = "LIPOCALIN2 t2", outline = F,
        ylab = "LIPOCALIN",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=1500, round(wt$p.value, 3))

### Enolase.1
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "Enolase.1"],
                dset.delta[which(dset.delta$randomization==2), "Enolase.1"],
                paired = F)
print(wt)
boxplot(Enolase.1 ~ randomization, data = dset.delta, main = "Enolase t1", outline = F,
        ylab = "Enolase",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=4500, round(wt$p.value, 3))

### Enolase.2
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "Enolase.2"],
                dset.delta[which(dset.delta$randomization==2), "Enolase.2"],
                paired = F)
print(wt)
boxplot(Enolase.2 ~ randomization, data = dset.delta, main = "Enolase t2", outline = F,
        ylab = "Relative change in Enolase",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=5000, round(wt$p.value, 3))

### FABP2.1
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "FABP2.1"],
                dset.delta[which(dset.delta$randomization==2), "FABP2.1"],
                paired = F)
print(wt)
boxplot(FABP2.1 ~ randomization, data = dset.delta, main = "FABP2.d10 t1", outline = F,
        ylab = "Relative change in FABP2",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=6, round(wt$p.value, 3))

### FABP2.2
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "FABP2.2"],
                dset.delta[which(dset.delta$randomization==2), "FABP2.2"],
                paired = F)
print(wt)
boxplot(FABP2.2 ~ randomization, data = dset.delta, main = "FABP2 t2", outline = F,
        ylab = "Relative change in FABP2",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=2, round(wt$p.value, 3))

### CystatinC.1
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "CystatinC.1"],
                dset.delta[which(dset.delta$randomization==2), "CystatinC.1"],
                paired = F)
print(wt)
boxplot(CystatinC.1 ~ randomization, data = dset.delta, main = "CystatinC t1", outline = F,
        ylab = "Relative change in CystatinC",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=40000, round(wt$p.value, 3))

### CystatinC.2
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "CystatinC.2"],
                dset.delta[which(dset.delta$randomization==2), "CystatinC.2"],
                paired = F)
print(wt)
boxplot(CystatinC.2 ~ randomization, data = dset.delta, main = "CystatinC t2", outline = F,
        ylab = "Relative change in CystatinC",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=40000, round(wt$p.value, 3))

### TNFalfa.1
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "TNFalfa.1"],
                dset.delta[which(dset.delta$randomization==2), "TNFalfa.1"],
                paired = F)
print(wt)
boxplot(TNFalfa.1 ~ randomization, data = dset.delta, main = "TNFalfa t1", outline = F,
        ylab = "Relative change in TNFalfa",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=30, round(wt$p.value, 3))


### TNFalfa.2
wt<-wilcox.test(dset.delta[which(dset.delta$randomization==1), "TNFalfa.2"],
                dset.delta[which(dset.delta$randomization==2), "TNFalfa.2"],
                paired = F)
print(wt)
boxplot(TNFalfa.2 ~ randomization, data = dset.delta, main = "TNFalfa t2", outline = F,
        ylab = "Relative change in TNFalfa",
        xlab = "Randomization group: 1=ctrl, 2=exp")
text(x=1.5, y=30, round(wt$p.value, 3))


