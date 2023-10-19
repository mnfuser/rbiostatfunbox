# Libraries -----------------------------------------------

# library(epitools)
library(survival)
library(survminer)
#library(RcmdrPlugin.EZR)
# library(MatchIt)
# library(cobalt)
library(knitr)
 library(rms)
 library(tidyverse)
 library(lubridate)
 library(tableone)
# library(lsmeans)
# library(arm)
# library(grid)
# library(ggpubr)
# library(ggrepel)
library(visdat)
 library(corrplot)
# library(RColorBrewer)
# library(survMisc)
# library(survivalROC)
# library(ROCR)
# library(plotROC)
# library(survIDINRI)
# library(nricens)
# library(mice)
# library(fastDummies)


# setVarType --------------------------------------------
# INPUT: file with var type description
# OUTPUT: lists of types

characterList<-NULL
dates<-NULL
factorList<-NULL
orderedList<-NULL
numericList<-NULL
integerList<-NULL
logicalList<-NULL

setVarType<-function(fileName){
  varType<-read.csv2(fileName)
  characterList<<-which(varType$type=="c" | varType$type=="character")
  dates<<-which(varType$type=="d" | varType$type=="date")
  factorList<<-which(varType$type=="f" | varType$type=="factor")
  orderedList<<-which(varType$type=="o" | varType$type=="ordinal")
  numericList<<-which(varType$type=="n" | varType$type=="numeric")
  integerList<<-which(varType$type=="i" | varType$type=="integer")
  logicalList<<-which(varType$type=="l" | varType$type=="logical")
}

# nmiss ----
# INPUT: dataset, columns range
# OUTPUT: na for each var in colrange
nmiss<-function(dset, colrange, fileName){
  # missing
  nm<-apply(dset[, colrange], 
            FUN = function(x){
                   length(which(is.na(x)))
                 }
            ,2)
  # missing percent
  pm<-apply(dset[, colrange], 
            FUN = function(x){
              length(which(is.na(x)))/dim(dset)[1]
            }
            ,2)
  # unique values or constant
  vm<-apply(dset[, colrange], 
            FUN = function(x){
              l=length(unique(x))
              ifelse(l>1,l,"CONST")
            }
            ,2)
  # data type
  cm<-as.matrix(sapply(dset, class))
  # regressor
  rg<-rep(0, dim(dset)[2])
  # outcome
  oc<-rep(0, dim(dset)[2])
  # position in the dataset
  id<-seq(1:dim(dset)[2])
  if(!missing(fileName)){
    write.csv2(data.frame(id=id,nm=nm,pm=pm,vm=vm,type=cm,rg=rg,oc=oc), file=fileName)  
  }
  return(data.frame(id=id,nm=nm, pm=pm,vm=vm,cm=cm,rg=rg,oc=oc))
}

# applyFactors ----
applyFactors<-function(x, c){
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      x[,c[i]]<-as.factor(x[,c[i]])  
    }  
  }
  return(x)
}

# do.panel.quant ----
do.panel.quant<-function(d,xx,s){
  if(str_length(s>50)) { s<-substr(s, 0, 50) }
  if(check.na.limit(d[,xx], FALSE)){
    quant.var.panel<-as.matrix(cbind(min(d[,xx], na.rm = TRUE), t(quantile(d[,xx],na.rm = TRUE))[2], median(d[,xx], na.rm = TRUE), mean(d[,xx], na.rm = TRUE), sd(d[,xx], na.rm = TRUE),  t(quantile(d[,xx],na.rm = TRUE))[4], max(d[,xx],na.rm = TRUE)))
    colnames(quant.var.panel)<-c("Min", "Q2", "Median", "Mean", "SD", "Q3", "Max")
    print(kable(quant.var.panel, "latex", longtable = T, booktabs = F, caption=paste0("Variable: ", s))%>%
            kable_styling(full_width = T, latex_options = c("hold_position")))
  }
}

# do.plot.quant ----
do.plot.quant<-function(d,xx,s){
  if(check.na.limit(d[,xx], FALSE)){
    print(ggplot(d, 
                 aes(x=1, y=d[,xx]))+
            geom_boxplot()+
            theme(axis.title.x=element_blank(),
                  axis.title.y = element_text(size=8),
                  axis.text.x=element_blank(),
                  axis.ticks.x =element_blank(),
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "grey98"),
                  plot.background = element_rect(fill = "lightblue"))+
            labs(y = s))
    return(1)
  } else return(0)
}

# do.panel.qual ----
do.panel.qual<-function(d,xx,s){
  if(str_length(s>50)) { s<-substr(s, 0, 50) }
  if(check.na.limit(d[,xx], FALSE)){
    qual.var.panel<-as.data.frame(table(d[,xx]))
    percent<-qual.var.panel[,2]/sum(qual.var.panel[,2])
    qual.var.panel<-cbind(qual.var.panel, percent)
    colnames(qual.var.panel)<-c("Levels", "Counts", "%")
    print(kable(qual.var.panel, "latex", longtable = T, booktabs = T, caption=paste0("Variable: ", s))%>%
            kable_styling(full_width = T, latex_options = c("repeat_header", "hold_position"))%>%
            column_spec(1, width = "8cm"))
  }
}

# do.plot.qual ----
do.plot.qual<-function(d,xx,s){
  if(check.na.limit(d[,xx], FALSE)){
    
    if(length(levels(as.factor(xx)))<=6){
      scaledLabels<-c(1:length(levels(as.factor(d[,xx]))))
    } else {
      scaledLabels<-""
    }
    
    print(ggplot(data=subset(d, !is.na(d[,xx])), aes(na.omit(d[,xx])))+
            geom_bar()+
            theme(axis.title.x=element_text(size=8),
                  axis.title.y = element_text(size=8),
                  axis.text.x=element_blank(),
                  axis.ticks.x =element_blank(),
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "grey98"),
                  plot.background = element_rect(fill = "lightblue"))+
            #scale_x_discrete(labels=scaledLabels) +
            labs(y = s, x="Levels"))
    return(1)  
  } else return(0)
}

# check.var.type ----
check.var.type<-function(d,x,s,p){
  r=0
  if(class(d[,x])=="integer" || class(d[,x])=="numeric") {
    if(p==FALSE) {do.panel.quant(d,x,s)}
    else {
      r=r+do.plot.quant(d,x,s)} 
  }
  else {
    if(p==FALSE) {do.panel.qual(d,x,s)}
    else {
      r=r+do.plot.qual(d,x,s)  
    }
  } 
  return(r)
}

# check.na.limit ----
check.na.limit<-function(x, s){
  v = round(length(which(is.na(x)))/length(x), 3)
  if(s==TRUE) print(v+"\n")
  if(v<=na.cutoff) return (TRUE)
  else return (FALSE)
}

# String encoding
encToUTF8<-function(x){
  if(class(x)=="character") {
    x<-stri_enc_toutf8(x, is_unknown_8bit = TRUE, validate = TRUE)
  }
  return(x)
}

# String trimming
trimStrings<-function(x){
  if(class(x)=="character") {
    x<-str_trim(x)
  }
  return(x)
}

# char to factor conversion
convertCharToFactor<-function(x, c){
  if(class(x)=="character") {
    x<-as.factor(x)
  }
  return(x)  
}

# blank to NA conversion
convertBlankToNA<-function(x){
  if(class(x)=="character") {
    x<-gsub("^$", NA, x)
  }
  return(x)
}

# applyCharacter ----
applyCharacter<-function(x, c){
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      x[,c[i]]<-as.character(x[,c[i]])  
    }  
  }
  return(x)
}

# applyFactors ----
applyFactors<-function(x, c){
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      x[,c[i]]<-as.factor(x[,c[i]])  
    }  
  }
  return(x)
}

# applyOrdered ----
applyOrdered<-function(x, c){
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      x[,c[i]]<-as.ordered(x[,c[i]])  
    }  
  }
  return(x)
}

# applyNumeric ----
applyNumeric<-function(x, c){
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      tryCatch(
        x[,c[i]]<-as.numeric(x[,c[i]])
        , warning=function(w){print(paste0("Warning on column: ", i))})
    }  
  }
  return(x)
}

# applyInteger----
applyInteger<-function(x, c){
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      x[,c[i]]<-as.integer(x[,c[i]])  
    }  
  }
  return(x)
}

# applyLogical ----
applyLogical<-function(x, c){
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      x[,c[i]]<-as.logical(x[,c[i]])  
    }  
  }
  return(x)
}

# applyDate ----
applyDate<-function(x, c, f){
  # f=format
  # f=1->"%d/%m/%Y"
  # f=2->"%Y-%m-%d"
  if(missing(c)){
    return(x)
  } else {
    for(i in 1:length(c)){
      ifelse(f==1,
             x[,c[i]]<-dmy(x[,c[i]]),
             x[,c[i]]<-ymd(x[,c[i]]))
    }  
  }
  return(x)
}

# findLongFactor ----
findLongFactor<-function(x){
  return(names(which(sapply(sapply(x[,sapply(x, is.factor)], levels), length)>30)))
}


# findDates ----
is.Date <- function(x) inherits(x, 'Date')

findDates<-function(x){
  return(names(which(sapply(x, is.Date))))
}

# sanitizeColNames ----
sanitizeColNames<-function(d){
  colnames(d)<-gsub("_", ".", colnames(d))
  for(i in 1:length(colnames(d))){
    if(str_length(colnames(d)[i])>50){
      colnames(d)[i]<-paste0(substr(colnames(d)[i], 0, 50),"[...]")
    }
  }
  return(d)
}

# vis_dat_wrapper ----
vis_dat_wrapper<-function(dset){
  vis_dat(dset, sort_type = FALSE, warn_large_data = FALSE)+
    theme(axis.title.x=element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.x=element_blank()
    )+xlab("Variables")+
    ggtitle("Post processing var types")  
}

# getFactors ----
getFactors<-function(fileName){
  varType<-read.csv2(paste0(getwd(),"/", fileName))
  factorList<<-which(varType$type=="f")
}

# getFactorsNamesFromDataset ----
getFactorsNamesFromDataset<-function(DatasetName){
  factorNames<<-DatasetName %>% Filter(f = is.factor) %>% names
}

# getNumbersNamesFromDataset ----
getNumbersNamesFromDataset<-function(DatasetName){
  numericNames<<-DatasetName %>% Filter(f = is.numeric) %>% names
}

# getFactorsFromDataset ----
getFactorsFromDataset<-function(DatasetName){
  factorNames<-DatasetName %>% Filter(f = is.factor) %>% names
  factorList<<-which(colnames(DatasetName) %in% factorNames)
}

# getNumbersFromDataset -----

getNumbersFromDataset<-function(DatasetName){
  numericNames<-DatasetName %>% Filter(f = is.numeric) %>% names
  numericList<<-which(colnames(DatasetName) %in% numericNames)
}

# GPS (glaswgow prognostic score) -----

GPS<-function(x,y){
  gps=0
  if(is.na(x) | is.na(y)) return (NA)
  if(x<=1 & y<3.5) return (1)
  if(x<=1 & y>=3.5) return (0)
  if(x>1 & y<3.5) return (2)
  if(x>1 & y>=3.5) return (1)  
}

# MPI (multidimensional prognostic index)  MNA-SF ----
## Input vector vars SPMSQ, ESS, ADL, IADL, CIRS, MNA, NUMBER OF MEDICATIONS, SOCIAL SUPPORT NETWORK
## Score matrix vars SPMSQ, ESS, ADL, IADL, CIRS, MNA, NUMBER OF MEDICATIONS, SOCIAL SUPPORT NETWORK
## Score matrix cutoff LOW, MILD, HIGH
mpi<-function(in.v){
  score.levels=NULL
  score.tmp=0
  score.m<-matrix(c(3,7,10,
                    20,15,9,
                    6,4,2,
                    8,5,3,
                    0,2,3,
                    14,11,7,
                    3,6,7,
                    0,1,2), 
                    byrow=F, 
                    ncol=8, 
                    dimnames = list(
                      c("LOW","MILD", "HIGH"), 
                      c("SPMSQ", "ESS", "ADL", "IADL", "CIRS", "MNA", "NUMBER OF MEDICATIONS", "SOCIAL SUPPORT NETWORK")))
  # SPMSQ 
  score.tmp<-score.tmp+ifelse(in.v[1]<=score.m[1,1],0,
                    ifelse(in.v[1]<=score.m[2,1],0.5,1))
  # ESS 
  score.tmp<-score.tmp+ifelse(in.v[2]<=score.m[3,2],1,
                    ifelse(in.v[2]<=score.m[2,2],0.5,0))
  # ADL 
  score.tmp<-score.tmp+ifelse(in.v[3]<=score.m[3,3],1,
                    ifelse(in.v[3]<=score.m[2,3],0.5,0))
  # IADL 
  score.tmp<-score.tmp+ifelse(in.v[4]<=score.m[3,4],1,
                    ifelse(in.v[4]<=score.m[2,4],0.5,0))
  # CIRS
  score.tmp<-score.tmp+ifelse(in.v[5]==0,0,
                    ifelse(in.v[5]<=score.m[2,5],0.5,1))
  # MNA 
  score.tmp<-score.tmp+ifelse(in.v[6]<=score.m[3,6],1,
                    ifelse(in.v[6]<=score.m[2,6],0.5,0))
  # NUMBER OF MEDICATION
  score.tmp<-score.tmp+ifelse(in.v[7]<=score.m[1,7],0,
                    ifelse(in.v[7]<=score.m[2,7],0.5,1))
  # SOCIAL SUPPORT NETWORK 
  score.tmp<-score.tmp+ifelse(in.v[8]<=score.m[1,8],0,
                    ifelse(in.v[8]<=score.m[2,8],0.5,1))
  
  score.final=score.tmp/8
  
  score.levels<-ifelse(score.final<=0.33,0,
         ifelse(score.final<=0.66,1,2))
  
  #return (c(score.final, score.levels))  
  return (list(score=score.final, class=score.levels))  
}

# GIC SCORE ----
gic<-function(dset){
  GIC<-rep(NA, dim(dset)[1])
  for(i in 1:dim(dset)[1]){
    v<-c(
      dset[i, "rf_ischemic_or_org_heartfa"],
      dset[i, "rf_primary_arrhythmiasfa"],
      dset[i, "rf_heart_diseases_w_nonifa"],
      dset[i, "rf_hypertension_gicfa"],
      dset[i, "rf_strokefa"],
      dset[i, "rf_periph_vascular_diseasfa"],
      dset[i, "rf_diabete_mellitusfa"],
      dset[i, "rf_anemiafa"],
      dset[i, "rf_gastrointesti_diseasesfa"],
      dset[i, "rf_hepatobiliary_diseasesfa"],
      dset[i, "rf_renal_diseasesfa"],
      dset[i, "rf_respiratory_diseasesfa"],
      dset[i, "rf_musculoskelet_disordersfa"],
      dset[i, "rf_malignanciesfa"]
    )
    print(v)
    g1<-which(v==1)
    g2<-which(v==2)
    g3<-which(v==3)
    g4<-which(v==4)
    ifelse(length(g4)>0, GIC[i]<-4, 
           ifelse(length(g3)>1, GIC[i]<-4, 
           ifelse(length(g3)==1, GIC[i]<-3, 
           ifelse(length(g2)>0, GIC[i]<-2, 
           ifelse(length(g1)>0, GIC[i]<-1, GIC[i]<-0)))))
  }
  return(GIC)
}

# EMT ----
# Calculate EMT score
EMT<-function(sex, # 0 = female; 1 = male
              chairRises, # <15s = 0; >= 15s = 1; 
              spmsq, # (congnitive impairment) < 8 = 1; >8 = 0; 
              hemoglobin, # Male >= 13g/dl = 0 (female 12g/dl); Male < 13g/dl = < (female 12g/dl) 
              serumAlbumine # >= 3.5g/dl = 0; < 3.5g/dl = 1  
              ){
  emtScore = rep(0, length(sex))
  emtScore=ifelse(chairRises == 0, emtScore + 2,
                  ifelse(chairRises >= 15, emtScore + 1, emtScore))
  emtScore=ifelse(spmsq < 8, emtScore + 1, emtScore)
  emtScore=ifelse((sex = 0 & hemoglobin < 12) | 
                    (sex = 1 & hemoglobin < 13), emtScore + 1, emtScore)
  emtScore=ifelse(serumAlbumine < 3.5, emtScore + 1, emtScore)
  return (emtScore)  
}

# get model stats from rms -----

get_model_stats = function(x) {
  cap = capture.output(print(x))
  
  #model stats
  stats = c()
  stats$R2.adj = str_match(cap, "R2 adj\\s+ (\\d\\.\\d+)") %>% na.omit() %>% .[, 2] %>% as.numeric()
  
  #coef stats lines
  coef_lines = cap[which(str_detect(cap, "Coef\\s+S\\.E\\.")):(length(cap) - 1)]
  
  #parse
  coef_lines_table = suppressWarnings(readr::read_table(coef_lines %>% stringr::str_c(collapse = "\n")))
  colnames(coef_lines_table)[1] = "Predictor"
  
  list(
    stats = stats,
    coefs = coef_lines_table
  )
}

# CalculateAucFromDxy -----

CalculateAucFromDxy <- function(validate) {
  ## Test if the object is correct
  stopifnot(class(validate) == "validate")
  
  ## Calculate AUCs from Dxy's
  aucs <- (validate["Dxy", c("index.orig","training","test","optimism","index.corrected")])/2 + 0.5
  
  ## Get n
  n <- validate["Dxy", c("n")]
  
  ## Combine as result
  res <- rbind(validate, AUC = c(aucs, n))
  
  ## Fix optimism
  res["AUC","optimism"] <- res["AUC","optimism"] - 0.5
  
  ## Return results
  res
}

# ApplyCap2Outliers -----

ApplyCap2Outliers<-function(dset.col){
  x <- dset.col
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  caps <- quantile(x, probs=c(.05, .95), na.rm = T)
  H <- 1.5 * IQR(x, na.rm = T)
  x[x < (qnt[1] - H)] <- caps[1]
  x[x > (qnt[2] + H)] <- caps[2]
  return(x)
}


# Cox.univar.rms ----------------------------------------------------------
# PERFORM UNIVARIABLE COX REGRESSION OF A SET OF VARIABLES AGAINST SELECTED OUTCOME VARIABLE 
# dset: dataset
# event.time: time of event vector
# event: event vector
# time.units: 'Day' or other
# var.interval: range for variables to be tested
# outdir: output directory
# f.out: output file name
cox.univar.rms<-function(
  dset,
  event.times,
  event,
  time.units,
  var.interval,
  outdir,
  f.out
) {
  ddist<<-datadist(dset)
  options(datadist='ddist')
  
  units(event.times) <- time.units
  S<-Surv(event.times, event==1)
  
  summary_list=vector("list", length(var.interval))
  l=1
  
  cox.out=data.frame()
  for(i in var.interval){
    fit<-cph(as.formula(paste("S~", colnames(dset)[i])),
             data = dset,
             y=TRUE, x=TRUE)
    print(fit, coef = TRUE)
    if(is.factor(dset[,i])){
      refl<-ifelse(!is.na(min(levels(dset[,i]), na.rm = T))
                   , min(levels(dset[,i]), na.rm = T)
                   , 1)
      s<-eval(parse(text = paste0("summary(fit, ",colnames(dset)[i], "='", refl,"')")))
    } else {
      s<-summary(fit)
    }
    pval<-get_model_stats(fit)$coefs[dim(get_model_stats(fit)$coefs)[2]]
    cox.out<-rbind(cox.out,
                   cbind(rownames(s)[seq(1, dim(s)[1], 2)],
                         round(s[seq(2, dim(s)[1], 2), 4],3),
                         paste0(round(s[seq(2, dim(s)[1], 2),6],3), " - ", round(s[seq(2, dim(s)[1], 2),7],3)),
                         pval))
    summary_list[[l]]<-summary(fit)
    l=l+1
  }
  
  colnames(cox.out)<-c("Variable", "HR", "95ci", "p")
  #cox.out$Pvalue<-as.numeric(levels(cox.out$`Pvalue`))[cox.out$`Pvalue`]
  print(kable(cox.out, row.names = FALSE))
  
  write.csv2(cox.out, file=paste0(outdir, "/", f.out, ".csv"))
  
  cox.out<-cox.out[which(cox.out[,4]<0.05),]
  print(kable(cox.out, row.names = FALSE))
  
  write.csv2(cox.out, file=paste0(outdir, "/", f.out, ".sig.csv"))
             
  return(summary_list)
}


# Cox.univar.rms2 ----------------------------------------------------------
# PERFORM UNIVARIABLE COX REGRESSION OF A SET OF VARIABLES AGAINST SELECTED OUTCOME VARIABLE 
# dset: dataset
# event.time: time of event vector
# event: event vector
# time.units: 'Day' or other
# var.interval: range for variables to be tested
# outdir: output directory
# f.out: output file name
cox.univar.rms2<-function(
  dset,
  event.times,
  event,
  time.units,
  var.interval,
  outdir,
  f.out
) {
  ddist<<-datadist(dset)
  options(datadist='ddist')
  
  units(event.times) <- time.units
  S<-Surv(event.times, event==1)
  
  summary_list=vector("list", length(var.interval))
  l=1
  
  cox.out=data.frame()
  for(i in var.interval){
    fit<-cph(as.formula(paste("S~", colnames(dset)[i])),
             data = dset,
             y=TRUE, x=TRUE)
    print(fit, coef = TRUE)
    s<-summary(fit)
    pval<-get_model_stats(fit)$coefs[dim(get_model_stats(fit)$coefs)[2]]
    cox.out<-rbind(cox.out,
                   cbind(rownames(s)[seq(1, dim(s)[1], 2)],
                         round(s[seq(2, dim(s)[1], 2), 4],3),
                         paste0(round(s[seq(2, dim(s)[1], 2),6],3), " - ", round(s[seq(2, dim(s)[1], 2),7],3)),
                         pval))
    summary_list[[l]]<-summary(fit)
    l=l+1
  }
  
  colnames(cox.out)<-c("Variable", "HR", "95ci", "p")
  #cox.out$Pvalue<-as.numeric(levels(cox.out$`Pvalue`))[cox.out$`Pvalue`]
  print(kable(cox.out, row.names = FALSE))
  
  write.csv2(cox.out, file=paste0(outdir, "/", f.out, ".csv"))
  
  cox.out<-cox.out[which(cox.out[,4]<0.05),]
  print(kable(cox.out, row.names = FALSE))
  
  write.csv2(cox.out, file=paste0(outdir, "/", f.out, ".sig.csv"))
  
  return(summary_list)
}


# Cox.univar.robust.rms ----------------------------------------------------------
# PERFORM UNIVARIABLE COX REGRESSION OF A SET OF VARIABLES AGAINST SELECTED OUTCOME VARIABLE 
# dset: dataset
# event.time: time of event vector
# event: event vector
# time.units: 'Day' or other
# var.interval: range for variables to be tested
# outdir: output directory
# f.out: output file name
cox.univar.robust.rms<-function(
  dset,
  Id,
  event.times,
  event,
  time.units,
  var.interval,
  outdir,
  f.out
) {
  ddist<<-datadist(dset)
  options(datadist='ddist')
  
  units(event.times) <- time.units
  S<-Surv(event.times, event==1)
  
  summary_list=vector("list", length(var.interval))
  l=1
  
  cox.out=data.frame()
  for(i in var.interval){
    fit<-cph(as.formula(paste("S~", colnames(dset)[i], "+cluster(", Id, ")")),
             robust = T,
             data = dset,
             y=TRUE, x=TRUE)
    print(fit, coef = TRUE)
    if(is.factor(dset[,i])){
      s<-eval(parse(text = paste("summary(fit,",colnames(dset)[i], "=0",")")))
    } else {
      s<-summary(fit)
    }
    pval<-get_model_stats(fit)$coefs[dim(get_model_stats(fit)$coefs)[2]]
    cox.out<-rbind(cox.out,
                   cbind(rownames(s)[seq(1, dim(s)[1], 2)],
                         round(s[seq(2, dim(s)[1], 2), 4],3),
                         paste0(round(s[seq(2, dim(s)[1], 2),6],3), " - ", round(s[seq(2, dim(s)[1], 2),7],3)),
                         pval))
    summary_list[[l]]<-summary(fit)
    l=l+1
  }
  
  colnames(cox.out)<-c("Variable", "HR", "95ci", "p")
  #cox.out$Pvalue<-as.numeric(levels(cox.out$`Pvalue`))[cox.out$`Pvalue`]
  print(kable(cox.out, row.names = FALSE))
  
  write.csv2(cox.out, file=paste0(outdir, "/", f.out, ".csv"))
  
  cox.out<-cox.out[which(cox.out[,4]<0.05),]
  print(kable(cox.out, row.names = FALSE))
  
  write.csv2(cox.out, file=paste0(outdir, "/", f.out, ".sig.csv"))
  
  return(summary_list)
}


# logi.bayes.univar -------------------------------------------------------
# PERFORM UNIVARIABLE BAYESAN LOGISTIC REGRESSION OF A SET OF VARIABLES AGAINST SELECTED OUTCOME VARIABLE 
# dset: dataset
# oucome.var: colum index for outcome
# var.interval: range for variables to be tested
# outdir: output directory
# f.out: output file name
logi.bayes.univar<-function(
  dset,
  outcome.var,
  var.interval,
  outdir,
  f.out  
){
  summary_list=vector("list", length(var.interval))
  l=1
  logi.out=data.frame()
  for(i in var.interval){
    
    print(colnames(dset)[i])
    
    logi.fit<-bayesglm(
      as.formula(
        paste(
          outcome.var,"~", colnames(dset)[i]
               )
        ), 
      family = binomial(), 
      #prior.df=Inf,
      #prior.mean = 0.5,
      data = dset)
  
    or=(exp(cbind("Odds ratio" = coef(logi.fit), confint.default(logi.fit))))
    logi.out<-rbind(logi.out, cbind(
      names(logi.fit$coefficients),
      round(or[1:dim(or)[1],1],3),
      paste0 (round(or[1:dim(or)[1],2],3)," - " ,round(or[1:dim(or)[1],3],3)),
      round(summary(logi.fit)$coefficients[1:dim(summary(logi.fit)$coefficients)[1],4] , 3)
    ))
    summary_list[[l]]<-summary(logi.fit)
    l=l+1
  }
  colnames(logi.out)<-c("Variable", "OR", "95ci", "p")
  logi.out$`OR`<-as.numeric(levels(logi.out$`OR`))[logi.out$`OR`]
  logi.out$p<-as.numeric(levels(logi.out$`p`))[logi.out$`p`]
  logi.out$Variable<-gsub("dset[, i]", "", logi.out$Variable, fixed=T)
  logi.out<-logi.out[which(logi.out$Variable!="(Intercept)"),]
  print(kable(logi.out, row.names = FALSE))
  
  write.csv2(logi.out, file=paste0(outdir, "/", f.out, ".csv"))
  
  logi.out<-logi.out[which(logi.out[,4]<0.05),]
  print(kable(logi.out, row.names = FALSE))
  
  write.csv2(logi.out, file=paste0(outdir, "/", f.out, ".sig.csv"))
  
  return(summary_list)
  
}

# logi.univar ----
# PERFORM UNIVARIABLE LOGISTIC REGRESSION OF A SET OF VARIABLES AGAINST SELECTED OUTCOME VARIABLE 
# dset: dataset
# oucome.var: colum index for outcome
# var.interval: range for variables to be tested
# outdir: output directory
# f.out: output file name
logi.univar<-function(
  dset,
  outcome.var,
  var.interval,
  outdir,
  f.out  
){
  summary_list=vector("list", length(var.interval))
  l=1
  logi.out=data.frame()
  for(i in var.interval){
    print(colnames(dset)[i])
    logi.fit<-glm(dset[,outcome.var] ~ dset[,i], family = binomial(), data = dset)
    or=(exp(cbind("Odds ratio" = coef(logi.fit), confint.default(logi.fit))))
    Variable=ifelse(is.factor(dset[,i]), 
           paste0(colnames(dset)[i], ", ", levels(dset[,i])[length(levels(dset[,i]))]),
           colnames(dset)[i])
    logi.out<-rbind(logi.out, 
                    data.frame(
                      Col=i,
                      Variable,
                      OR=round(or[1:dim(or)[1],1],3),
                      LCI=round(or[1:dim(or)[1],2],3),
                      UCI=round(or[1:dim(or)[1],3],3),
                      pvalue=round(summary(logi.fit)$coefficients[1:dim(summary(logi.fit)$coefficients)[1],4] , 3)
                    )
    )
    summary_list[[l]]<-summary(logi.fit)
    l=l+1
  }
  sel<-grep("Intercept", rownames(logi.out))
  logi.out<-logi.out[-sel,]
  print(kable(logi.out, row.names = FALSE))
  
  write.csv2(logi.out, file=paste0(outdir, "/", f.out, ".csv"), row.names = F)
  
  logi.out<-logi.out[which(logi.out[,6]<0.05),]
  print(kable(logi.out, row.names = FALSE))
  
  write.csv2(logi.out, file=paste0(outdir, "/", f.out, ".sig.csv"), row.names = F)
  
  return(logi.out)
  
}



# lm.univar -------------------------------------------------------
# PERFORM UNIVARIABLE LINEAR REGRESSION OF A SET OF VARIABLES AGAINST SELECTED OUTCOME VARIABLE 
# dset: dataset
# oucome.var: colum index for outcome
# var.interval: range for variables to be tested
# outdir: output directory
# f.out: output file name
lm.univar<-function(
  dset,
  outcome.var,
  var.interval,
  outdir,
  f.out  
){
  summary_list=vector("list", length(var.interval))
  l=1
  lm.out=data.frame()
  for(i in var.interval){
    print(colnames(dset)[i])
    lm.fit<-lm(dset[, outcome.var] ~ dset[,i], data = dset)
    Variable=ifelse(is.factor(dset[,i]), 
                    paste0(colnames(dset)[i], ", ", levels(dset[,i])[length(levels(dset[,i]))]),
                    colnames(dset)[i])
    lm.out<-rbind(lm.out, 
                  data.frame(
                    Col=i,
                    Variable,
                    B=round(summary(lm.fit)$coefficients[1:dim(summary(lm.fit)$coefficients)[1],1], 3), #beta
                    CI=round(confint(lm.fit), 3),# confint
                    pvalue=round(summary(lm.fit)$coefficients[1:dim(summary(lm.fit)$coefficients)[1],4], 3) #pval
    ))
    summary_list[[l]]<-summary(lm.fit)
    l=l+1
  }
  sel<-grep("Intercept", rownames(lm.out))
  lm.out<-lm.out[-sel,]
  colnames(lm.out)<-c("Col", "Variable", "B", "LCI", "UCI", "p-value")
  print(kable(lm.out, row.names = FALSE))
  
  write.csv2(lm.out, file=paste0(outdir, "/", f.out, ".csv"), row.names = F)
  
  lm.out<-lm.out[which(lm.out[,6]<0.05),]
  print(kable(lm.out, row.names = FALSE))
  
  write.csv2(lm.out, file=paste0(outdir, "/", f.out, ".sig.csv"), row.names = F)
  
  return(summary_list)
  
}

idi.out<-function (data, outcome, predrisk1, predrisk2, cutoff){
    c1 <- cut(predrisk1, breaks = cutoff, include.lowest = TRUE, 
              right = FALSE)
    c2 <- cut(predrisk2, breaks = cutoff, include.lowest = TRUE, 
              right = FALSE)
    tabReclas <- table(`Initial Model` = c1, `Updated Model` = c2)
    cat(" _________________________________________\n")
    cat(" \n     Reclassification table    \n")
    cat(" _________________________________________\n")
    ta <- table(c1, c2, data[, cOutcome])
    cat("\n Outcome: absent \n  \n")
    TabAbs <- ta[, , 1]
    tab1 <- cbind(TabAbs, ` % reclassified` = round((rowSums(TabAbs) - 
                                                       diag(TabAbs))/rowSums(TabAbs), 2) * 100)
    names(dimnames(tab1)) <- c("Initial Model", "Updated Model")
    print(tab1)
    cat("\n \n Outcome: present \n  \n")
    TabPre <- ta[, , 2]
    tab2 <- cbind(TabPre, ` % reclassified` = round((rowSums(TabPre) - 
                                                       diag(TabPre))/rowSums(TabPre), 2) * 100)
    names(dimnames(tab2)) <- c("Initial Model", "Updated Model")
    print(tab2)
    cat("\n \n Combined Data \n  \n")
    Tab <- tabReclas
    tab <- cbind(Tab, ` % reclassified` = round((rowSums(Tab) - 
                                                   diag(Tab))/rowSums(Tab), 2) * 100)
    names(dimnames(tab)) <- c("Initial Model", "Updated Model")
    print(tab)
    cat(" _________________________________________\n")
    c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
    c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
    x <- improveProb(x1 = as.numeric(c11) * (1/(length(levels(c11)))), 
                     x2 = as.numeric(c22) * (1/(length(levels(c22)))), y = data[, 
                                                                                cOutcome])
    y <- improveProb(x1 = predrisk1, x2 = predrisk2, y = data[, 
                                                              cOutcome])
    cat("\n NRI(Categorical) [95% CI]:", round(x$nri, 4), "[", 
        round(x$nri - 1.96 * x$se.nri, 4), "-", round(x$nri + 
                                                        1.96 * x$se.nri, 4), "]", "; p-value:", round(2 * 
                                                                                                        pnorm(-abs(x$z.nri)), 5), "\n")
    cat(" NRI(Continuous) [95% CI]:", round(y$nri, 4), "[", round(y$nri - 
                                                                    1.96 * y$se.nri, 4), "-", round(y$nri + 1.96 * y$se.nri, 
                                                                                                    4), "]", "; p-value:", round(2 * pnorm(-abs(y$z.nri)), 
                                                                                                                                 5), "\n")
    cat(" IDI [95% CI]:", round(y$idi, 4), "[", round(y$idi - 
                                                        1.96 * y$se.idi, 4), "-", round(y$idi + 1.96 * y$se.idi, 
                                                                                        4), "]", "; p-value:", round(2 * pnorm(-abs(y$z.idi)), 
                                                                                                                     5), "\n")
    
    return(cbind(idi=round(y$idi, 4),
                 lower=round(y$idi - 1.96 * y$se.idi, 4),
                 upper=round(y$idi + 1.96 * y$se.idi, 4),
                 pval=round(2 * pnorm(-abs(y$z.idi)), 3)
           ))
  }



# pairewise chisquare ----
pairwise.chisq.test <- function(x, g, p.adjust.method = p.adjust.methods, ...) {
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  p.adjust.method <- match.arg(p.adjust.method)
  
  compare.levels <- function(i, j) {
    xi <- x[as.integer(g) == i]
    xj <- x[as.integer(g) == j]
    chisq.test(table(xi, xj), ...)$p.value
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  ans <- list(method = "chi-squared test", data.name = DNAME, p.value = PVAL, 
              p.adjust.method = p.adjust.method)
  class(ans) <- "pairwise.htest"
  ans
}
# write.coxph.table ----
write.coxph.table<-function(fit, d, fileName){
  cox.out=data.frame()
  s<-summary(fit)
  cox.out<-data.frame(
                Variables=rownames(s$coefficients),
                HR=round(s$coefficients[,2],d),
                CI95=paste0(round(s$conf.int[,3],d), " - ", round(s$conf.int[,4],d)),
                pvalue=round(s$coefficients[,6],d)
  )
  
  cox.out
  
  write.csv2(cox.out, file=paste0(outdir, "/", fileName, ".csv"), row.names = F)
}
# print.cox ----
print.cox<-function(fit){
  cox.out<-data.frame()
  s<-summary(fit)
  cox.out<-rbind(cox.out,
                 cbind(rownames(s$coefficients),
                       round(s$coefficients[,2],3),
                       paste0(round(exp(confint(fit1)[,1]),3), " - ", round(exp(confint(fit1)[,2]),3)),
                       round(s$coefficients[,5],3))
  )
  colnames(cox.out)<-c("Variable", "HR", "95% c.i.", "p value")
  return(cox.out)
}
