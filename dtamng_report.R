# DATA MANAGEMENT AND TABLE EXTRACTION FROM EXPORTED EUROCAT DATASET 2016-2020 ALL ANOMALIES
# 2023 March 24
# Marco Manfrini Ph.D.
# University of Ferrara - Dept Moedical Sciences
# Center for Clinical and Epidemiological Research
# IMER Registry - Emilia Romagna Registry of Birth Defects

# PACKAGES ----

library(data.table)
library(epiR)
library(ggplot2)
library(tableone)
library(rms)
library(readxl)


source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

# FUNCTIONS ----

findICDCases<-function(codes, data){
  selected<- NULL
  for(c in 1:dim(data)[2]){
    for(code in codes){
      selected<-c(selected, grep(paste0("^",code), data[,c]))
    }
  }
  # apply(data,
  #       MARGIN = 2,
  #       FUN = function(x, codes, selected){
  #         for (code in codes){
  #           selected<-c(selected, grep(code, x))
  #         }
  #       },
  #       codes,
  #       selected
  #     )
  return(selected)
}

smrCalc = function(tmp, tmp_c, index){
  oo=as.numeric(tmp["Totali", index])
  if(oo!=0){
    mm<-ifelse(oo<=5, "mid.p", "vandenbroucke")
    etn=as.numeric(tmp_c["Totali", index])
    ep=as.numeric(tmp_c["Popolazione", index])
    otn=as.numeric(tmp["Popolazione", index])
    ee=ee=etn/ep*otn
    esmr<-epi.smr(o=oo, e=ee ,method = mm, conf.level = 0.95)
    return(paste0(round(esmr$est,2), 
                  " (", 
                  round(esmr$lower,2), 
                  " - ", 
                  round(esmr$upper,2), 
                  ")"))  
  } else {
    return(0)
  }
}

getPPoints<-function(tmp, group, index){
  tt<-as.numeric(tmp["Totali", index])
  pp<-as.numeric(tmp["Popolazione", index])
  PrevTot=round(tt/pp*10000, 2)
  conf.ll<-round((1.96/2 - sqrt(tt + 0.02))^2/pp*10000,2)
  conf.ul<-round((1.96/2 + sqrt(tt + 0.96))^2/pp*10000,2)
  pps<-data.frame(Anno=tmp["Anno", index],
             PrevTot=PrevTot,
             Pll=conf.ll,
             Pul=conf.ul,
             group=group)
  return(pps)
}

makePlot<-function(pPoints, anom){
  
  jpeg(filename = paste0(plotDir, "/", anom, ".jpg"), 
       width = 8, height = 4, units = "in", pointsize = 12,
       res = 300)
  
  p<-ggplot(pPoints, aes(x = Anno, y = PrevTot, group = group)) + 
  geom_line(aes(color = group)) + 
  geom_point(shape=21, size=2, aes(color = group, fill = group)) + 
  geom_ribbon(aes(ymin = Pll, ymax = Pul, fill = group), 
              alpha=0.1) +
  #theme_light() +
  theme(# axis.line.y.right = element_line(linewidth = 1, colour = "grey80"),
  #       axis.line.x.top = element_line(linewidth = 1, colour = "grey80"),
  #       axis.line.y.left = element_line(linewidth = 1, colour = "grey80"),
  #       axis.line.x.bottom = element_line(linewidth = 1, colour = "grey80"),
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.title=element_blank()) +
  xlab("Anno") +
  ylab("Prevalenza totale (x10000)")
  #annotate("text", x=2, y=min(pPoints$Pll), label= anom)
  
  print(p)
  
  dev.off()
}

# GLOBAL VARS DATA STRUCTURES ----

colRegistry=1
colAnomaly=2

# SET WORKING DIR ----

wd="/Users/mmanfrini/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/tabelle_eurocat"

outdir="/Users/mmanfrini/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023"

setwd(wd)

anomalyGroupsDir=paste0(wd,"/anomaly_groups")

otherSyndromeGroupsDir=paste0(wd,"/other_syndrome")

geneticAnomaliesDir=paste0(wd,"/genetic_anomalies")

selectedAnomaliesDir=paste0(wd,"/selected_anomalies")

plotDir=paste0(wd,"/plot")

# DATA IMPORT ----

dataset_2016_2020_combined_exc_gen_con <- read.csv("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/tabelle_eurocat/dataset_2016_2020_combined_exc_gen_con.csv", sep=";")

dataset_2016_2020_combined_inc_gen_con <- read.csv("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/tabelle_eurocat/dataset_2016_2020_combined_inc_gen_con.csv", sep=";")

dataset_2016_2020_exc_gen_con <- read.csv("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/tabelle_eurocat/dataset_2016_2021_exc_gen_con.csv", sep=";")

dataset_2016_2020_inc_gen_con <- read.csv("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/tabelle_eurocat/dataset_2016_2021_inc_gen_con.csv", sep=";")

imer_cedap_2020 <- read_excel("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/database/2020/completo/Imer_cedap_2020.xlsx")

mamme_2021_RER <- read_excel("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/cedap/mamme_2021_RER.xlsx", na = "-")
mamme_2021_RER <- mamme_2021_RER[,
                                 c(
                                   12, # Età della madre
                                   13, # Età del padre
                                   17, # Cittadinanza madre
                                   26, # Stato civile madre                                 
                                   27, # Cittadinanza padre 
                                   29, # Titolo studio madre  
                                   33, # Titolo studio padre 
                                   37, # Concepimenti precedenti                            
                                   38, # Parti precedenti                                   
                                   39, # Aborti spontanei                                   
                                   40, # I.V.G.                                             
                                   41, # Nati vivi                                          
                                   42, # Nati morti                                         
                                   43, # Tagli cesarei  
                                   45, # Abitudine fumo                                     
                                   46, # Altezza madre                                      
                                   47, # Peso pregravidico                                  
                                   48, # Peso al parto                                      
                                   49, # Variazione ponderale                               
                                   50, # BMI                                                
                                   51, # Consanguineità      
                                   56, # Corso preparto                                     
                                   57, # Test combinato                                     
                                   58, # Villocentesi                                       
                                   59, # Amniocentesi                                       
                                   60, # Fetoscopia/Funicolocentesi                         
                                   61, # Ecografia >22 sett   
                                   64, # Procreazione assistita                             
                                   65, # Metodo proc assist 
                                   67, # Età gestazionale  
                                   69, # Metodo calcolo dur grav 
                                   86, # Genere parto                                       
                                   87, # Nati maschi                                        
                                   88 # Nati femmine    
                                 )
                                 ]

Cedap2016 <- read.csv2("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/database/Cedap2016-2020/Cedap2016.csv", 
                        na="-")

Cedap2016 <- Cedap2016[,
                       c(
                         2,
                         16,
                         17,
                         18,
                         19,
                         29,
                         32,
                         34,
                         36,
                         68,
                         69,
                         73,
                         82,
                         85,
                         89,
                         93,
                         94,
                         95,
                         96,
                         97,
                         98,
                         101,
                         106,
                         107,
                         112,
                         113,
                         114,
                         115,
                         116,
                         117,
                         120,
                         121,
                         123,
                         142,
                         143,
                         144
                       )
                       ]

Cedap2017 <- read.csv2("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/database/Cedap2016-2020/Cedap2017.csv", 
                        na="-")

Cedap2017 <- Cedap2017[,
                       c(2,
                         113,
                         114,
                         115,
                         116,
                         126,
                         129,
                         131,
                         133,
                         18,
                         19,
                         23,
                         32,
                         35,
                         39,
                         43,
                         44,
                         45,
                         46,
                         47,
                         48,
                         51,
                         56,
                         57,
                         62,
                         63,
                         64,
                         65,
                         66,
                         67,
                         70,
                         71,
                         73,
                         92,
                         93,
                         94
                       )
]

Cedap2018 <- read.csv2("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/database/Cedap2016-2020/Cedap2018.csv", 
                        na="-")

Cedap2018 <- Cedap2018[,
                       c(1,
                         110,
                         111,
                         112,
                         113,
                         123,
                         125,
                         127,
                         129,
                         14,
                         15,
                         19,
                         28,
                         31,
                         35,
                         39,
                         40,
                         41,
                         42,
                         43,
                         44,
                         47,
                         52,
                         53,
                         58,
                         59,
                         60,
                         61,
                         62,
                         63,
                         66,
                         67,
                         69,
                         87,
                         88,
                         89  
                       )]

Cedap2019 <- read.csv2("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/database/Cedap2016-2020/Cedap2019.csv", 
                        na="-")

Cedap2019 <- Cedap2019[,
                       c(2,
                         17,
                         18,
                         19,
                         20,
                         31,
                         34,
                         36,
                         38,
                         70,
                         71,
                         76,
                         85,
                         88,
                         92,
                         96,
                         97,
                         98,
                         99,
                         100,
                         101,
                         104,
                         109,
                         110,
                         115,
                         116,
                         117,
                         118,
                         119,
                         120,
                         123,
                         124,
                         126,
                         145,
                         146,
                         147
                       )]


Cedap2020 <- read.csv2("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/database/Cedap2016-2020/Cedap2020.csv", 
                        na="-")

Cedap2020 <- Cedap2020[,
                       c(2,
                         117,
                         118,
                         119,
                         120,
                         131,
                         134,
                         136,
                         138,
                         22,
                         23,
                         28,
                         37,
                         40,
                         44,
                         48,
                         49,
                         50,
                         51,
                         52,
                         53,
                         56,
                         61,
                         62,
                         67,
                         68,
                         69,
                         70,
                         71,
                         72,
                         75,
                         76,
                         78,
                         97,
                         98,
                         99
                       )]

# c2016=colnames(Cedap2016)
# c2017=colnames(Cedap2017)
# c2018=colnames(Cedap2018)
# c2019=colnames(Cedap2019)
# c2020=colnames(Cedap2020)
# 
# length(c2016)<-length(c2020)
# length(c2017)<-length(c2020)
# length(c2018)<-length(c2020)
# 
# 
# cn<-data.frame(c2016=c2016,
#               c2017=c2017,
#               c2018=c2018,
#               c2019=c2019,
#               c2020=c2020
#               
#               )
# 
# write.csv2(cn, file=paste0(outdir, "/cn.csv"), row.names = F)

colnames(Cedap2016)<-tolower(colnames(Cedap2016))
colnames(Cedap2017)<-tolower(colnames(Cedap2017))
colnames(Cedap2018)<-tolower(colnames(Cedap2018))
colnames(Cedap2019)<-tolower(colnames(Cedap2019))
colnames(Cedap2020)<-tolower(colnames(Cedap2020))
colnames(Cedap2019)[16]<-"concepimenti_pr"
colnames(Cedap2020)[16]<-"concepimenti_pr"

cedap.2016.2020<-rbind(Cedap2016, Cedap2017, Cedap2018, Cedap2019, Cedap2020)
rm(Cedap2016, Cedap2017, Cedap2018, Cedap2019, Cedap2020)

Imer1978.2020 <- read.csv("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/database/Imer1978-2020.csv", sep=";")
Imer2016.2020 <- Imer1978.2020[which(Imer1978.2020$Anno>2015),]
rm(Imer1978.2020)

# SET COUNTRY ----
## dataset_2016_2020_combined_exc_gen_com----

country=''
for(i in 1:dim(dataset_2016_2020_combined_exc_gen_con)[1]){
  ifelse(dataset_2016_2020_combined_exc_gen_con[i, colRegistry]!='', 
         country<-dataset_2016_2020_combined_exc_gen_con[i, colRegistry],
         dataset_2016_2020_combined_exc_gen_con[i, colRegistry]<-country)
}

## dataset_2016_2020_combined_inc_gen_com----

country=''
for(i in 1:dim(dataset_2016_2020_combined_inc_gen_con)[1]){
  ifelse(dataset_2016_2020_combined_inc_gen_con[i, colRegistry]!='', 
         country<-dataset_2016_2020_combined_inc_gen_con[i, colRegistry],
         dataset_2016_2020_combined_inc_gen_con[i, colRegistry]<-country)
}

## dataset_2016_2020_exc_gen_con ----

country=''
for(i in 1:dim(dataset_2016_2020_exc_gen_con)[1]){
  ifelse(dataset_2016_2020_exc_gen_con[i, colRegistry]!='', 
         country<-dataset_2016_2020_exc_gen_con[i, colRegistry],
         dataset_2016_2020_exc_gen_con[i, colRegistry]<-country)
}

## dataset_2016_2020_inc_gen_con ----

country=''
for(i in 1:dim(dataset_2016_2020_inc_gen_con)[1]){
  ifelse(dataset_2016_2020_inc_gen_con[i, colRegistry]!='', 
         country<-dataset_2016_2020_inc_gen_con[i, colRegistry],
         dataset_2016_2020_inc_gen_con[i, colRegistry]<-country)
}


# SET ANOMALY ----
## dataset_2016_2020_combined_exc_gen_con ----

anomaly=''
for(i in 1:dim(dataset_2016_2020_combined_exc_gen_con)[1]){
  ifelse(dataset_2016_2020_combined_exc_gen_con[i, colAnomaly]!='', 
         anomaly<-dataset_2016_2020_combined_exc_gen_con[i, colAnomaly],
         dataset_2016_2020_combined_exc_gen_con[i, colAnomaly]<-anomaly)
}

## dataset_2016_2020_combined_inc_gen_con ----

anomaly=''
for(i in 1:dim(dataset_2016_2020_combined_inc_gen_con)[1]){
  ifelse(dataset_2016_2020_combined_inc_gen_con[i, colAnomaly]!='', 
         anomaly<-dataset_2016_2020_combined_inc_gen_con[i, colAnomaly],
         dataset_2016_2020_combined_inc_gen_con[i, colAnomaly]<-anomaly)
}

## dataset_2016_2020_exc_gen_con ----

anomaly=''
for(i in 1:dim(dataset_2016_2020_exc_gen_con)[1]){
  ifelse(dataset_2016_2020_exc_gen_con[i, colAnomaly]!='', 
         anomaly<-dataset_2016_2020_exc_gen_con[i, colAnomaly],
         dataset_2016_2020_exc_gen_con[i, colAnomaly]<-anomaly)
}

## dataset_2016_2020_inc_gen_con ----

anomaly=''
for(i in 1:dim(dataset_2016_2020_inc_gen_con)[1]){
  ifelse(dataset_2016_2020_inc_gen_con[i, colAnomaly]!='', 
         anomaly<-dataset_2016_2020_inc_gen_con[i, colAnomaly],
         dataset_2016_2020_inc_gen_con[i, colAnomaly]<-anomaly)
}

# FIX SLASH IN ANOMALY ----

dataset_2016_2020_combined_exc_gen_con[, colAnomaly]<-gsub("/", " - ", dataset_2016_2020_combined_exc_gen_con[, colAnomaly])
dataset_2016_2020_combined_inc_gen_con[, colAnomaly]<-gsub("/", " - ", dataset_2016_2020_combined_inc_gen_con[, colAnomaly])
dataset_2016_2020_exc_gen_con[, colAnomaly]<-gsub("/", " - ", dataset_2016_2020_exc_gen_con[, colAnomaly])
dataset_2016_2020_inc_gen_con[, colAnomaly]<-gsub("/", " - ", dataset_2016_2020_inc_gen_con[, colAnomaly])

dataset_2016_2020_combined_exc_gen_con[, colAnomaly]<-gsub("\xd0", " - ", dataset_2016_2020_combined_exc_gen_con[, colAnomaly])
dataset_2016_2020_combined_inc_gen_con[, colAnomaly]<-gsub("\xd0", " - ", dataset_2016_2020_combined_inc_gen_con[, colAnomaly])
dataset_2016_2020_exc_gen_con[, colAnomaly]<-gsub("\xd0", " - ", dataset_2016_2020_exc_gen_con[, colAnomaly])
dataset_2016_2020_inc_gen_con[, colAnomaly]<-gsub("\xd0", " - ", dataset_2016_2020_inc_gen_con[, colAnomaly])

dataset_2016_2020_combined_exc_gen_con[, colAnomaly]<-gsub("\xd5s", " ", dataset_2016_2020_combined_exc_gen_con[, colAnomaly])
dataset_2016_2020_combined_inc_gen_con[, colAnomaly]<-gsub("\xd5s", " ", dataset_2016_2020_combined_inc_gen_con[, colAnomaly])
dataset_2016_2020_exc_gen_con[, colAnomaly]<-gsub("\xd5s", " ", dataset_2016_2020_exc_gen_con[, colAnomaly])
dataset_2016_2020_inc_gen_con[, colAnomaly]<-gsub("\xd5s", " ", dataset_2016_2020_inc_gen_con[, colAnomaly])


# REGS ----

rgs<-levels(as.factor(dataset_2016_2020_inc_gen_con[,colRegistry]))

# ANOMALIES  ----

l<-levels(as.factor(dataset_2016_2020_inc_gen_con[,colAnomaly]))

# ANALISYS EXCLUDING GENTICAL ANOMALIES
# ANOMALY GROUPS ----

ags<-l[c(3, 70, 47, 43, 32, 83, 74, 48, 1, 28, 51, 61)]

# OTHER ANOMALIES/SYNDROMES ----

oass<-l[c(37, 89, 35, 103, 77, 20, 88, 84, 102, 60, 93, 101, 67)]

# GENETIC ANOMALIES ----

lg<-levels(as.factor(dataset_2016_2020_inc_gen_con[,colAnomaly]))

ga<-l[c(90, 41, 75, 45, 99, 97)]

# SINGLE SELECTED ANOMALIES ----

ssa<-l[-(which(l %in% ags))]
ssa<-ssa[-(which(ssa %in% oass))]
ssa<-ssa[-(which(ssa %in% ga))]
ssa<-ssa[-38] # GENETIC DISORDERS

# SUBTABLES EXTRACTION AND PLOTS ANOMALY GROUPS ----

total.prev<-data.frame()

print("SUBTABLES EXTRACTION AND PLOTS ANOMALY GROUPS")
tabname=''
tabname_c=''
colNames <- c("Registro",
              "Anomalia",
              "Anno",
              "Popolazione",
              "Nati vivi",
              "Nati morti",
              "IVG",
              "Totali",
              "Prevalenza",
              "Prevalenza nati vivi",
              "Prevalenza nati morti",
              "Prevalenza IVG")

colNames_c <- colNames

for(k in ags){
  
  pPoints<-data.frame()
  
  esmr<-vector()
  tmp<-t(
    dataset_2016_2020_exc_gen_con[which(dataset_2016_2020_exc_gen_con[,colRegistry]==rgs[7] & 
                                               dataset_2016_2020_exc_gen_con[,colAnomaly]==k), c(1:12)]
  )
  tmp_c<-t(
    dataset_2016_2020_combined_exc_gen_con[which(dataset_2016_2020_combined_exc_gen_con[,colAnomaly]==k), c(1:12)]
  )
  if(dim(tmp)[2]>0){
    rownames(tmp)<-colNames
    rownames(tmp_c)<-colNames
    
    for(j in 1:dim(tmp)[2]){
      esmr<-c(esmr, smrCalc(tmp, tmp_c, j))
      pPoints<-rbind(pPoints, getPPoints(tmp, "IMER", j))
    }
    
    for(j in 1:dim(tmp_c)[2]){
      pPoints<-rbind(pPoints, getPPoints(tmp_c, "EUROCAT", j))
    }
    
    #makePlot(pPoints, k)
    
    tmp<-rbind(tmp, esmr)
    row.names(tmp)<-c(rownames(tmp)[1:12], "SMR")
    
    tabname=paste0(rgs[7], "_", k)
    tabname_c=paste0("Total registries", "_", k, "_cmb")
    
    tpop.local<-sum(as.numeric(tmp[4,]))
    tcases.local<-sum(as.numeric(tmp[8,]))
    
    tprev.local<-round(tcases.local/tpop.local*10000,2)
    
    tprev.local.conf.ll<-round((1.96/2 - sqrt(tcases.local + 0.02))^2/tpop.local*10000,2)
    tprev.local.conf.ul<-round((1.96/2 + sqrt(tcases.local + 0.96))^2/tpop.local*10000,2)
    
    trow<-data.frame(Anomaly=tabname,
                     Population=tpop.local,
                     Cases=tcases.local,
                     Prevalence=tprev.local,
                     Lower=tprev.local.conf.ll,
                     Upper=tprev.local.conf.ul)
    
    total.prev=rbind(total.prev, trow)
    
    #write.csv2(tmp, paste0(anomalyGroupsDir, "/", tabname, ".csv"), col.names = F)
    #write.csv2(tmp_c, paste0(anomalyGroupsDir, "/", tabname_c, ".csv"), col.names = F)  
  }
  
  write.csv2(total.prev, paste0(anomalyGroupsDir, "/total.prev.csv"), row.names = F)
  
}

# SUBTABLES EXTRACTION AND PLOTS OTHER ANOMALIES/SYNDROMES ----

total.prev<-data.frame()

print("OTHER ANOMALIES/SYNDROMES")
tabname=''
tabname_c=''
colNames <- c("Registro",
              "Anomalia",
              "Anno",
              "Popolazione",
              "Nati vivi",
              "Nati morti",
              "IVG",
              "Totali",
              "Prevalenza",
              "Prevalenza nati vivi",
              "Prevalenza nati morti",
              "Prevalenza IVG")

colNames_c <- colNames

for(k in oass){
  
  pPoints<-data.frame()
  
  esmr<-vector()
  tmp<-t(
    dataset_2016_2020_exc_gen_con[which(dataset_2016_2020_exc_gen_con[,colRegistry]==rgs[7] & 
                                          dataset_2016_2020_exc_gen_con[colAnomaly]==k), c(1:12)]
  )
  tmp_c<-t(
    dataset_2016_2020_combined_exc_gen_con[which(dataset_2016_2020_combined_exc_gen_con[,colAnomaly]==k), c(1:12)]
  )
  if(dim(tmp)[2]>0){
    rownames(tmp)<-colNames
    rownames(tmp_c)<-colNames
    
    for(j in 1:dim(tmp)[2]){
      esmr<-c(esmr, smrCalc(tmp, tmp_c, j))
      pPoints<-rbind(pPoints, getPPoints(tmp, "IMER", j))
    }
    
    for(j in 1:dim(tmp_c)[2]){
      pPoints<-rbind(pPoints, getPPoints(tmp_c, "EUROCAT", j))
    }
    
    #makePlot(pPoints, k)
    
    tmp<-rbind(tmp, esmr)
    row.names(tmp)<-c(rownames(tmp)[1:12], "SMR")
    
    tabname=paste0(rgs[7], "_", k)
    tabname_c=paste0("Total registries", "_", k, "_cmb")
    
    tpop.local<-sum(as.numeric(tmp[4,]))
    tcases.local<-sum(as.numeric(tmp[8,]))
    
    tprev.local<-round(tcases.local/tpop.local*10000,2)
    
    tprev.local.conf.ll<-round((1.96/2 - sqrt(tcases.local + 0.02))^2/tpop.local*10000,2)
    tprev.local.conf.ul<-round((1.96/2 + sqrt(tcases.local + 0.96))^2/tpop.local*10000,2)
    
    trow<-data.frame(Anomaly=tabname,
                     Population=tpop.local,
                     Cases=tcases.local,
                     Prevalence=tprev.local,
                     Lower=tprev.local.conf.ll,
                     Upper=tprev.local.conf.ul)
    
    total.prev=rbind(total.prev, trow)
    
    #write.csv2(tmp, paste0(otherSyndromeGroupsDir, "/", tabname, ".csv"), col.names = F)
    #write.csv2(tmp_c, paste0(otherSyndromeGroupsDir, "/", tabname_c, ".csv"), col.names = F) 
  }
  
  write.csv2(total.prev, paste0(otherSyndromeGroupsDir, "/total.prev.csv"), row.names = F)
  
}

# SUBTABLES EXTRACTION AND PLOTS SINGLE SELECTED ANOMALIES ----

total.prev<-data.frame()

tabname=''
tabname_c=''
colNames <- c("Registro",
              "Anomalia",
              "Anno",
              "Popolazione",
              "Nati vivi",
              "Nati morti",
              "IVG",
              "Totali",
              "Prevalenza",
              "Prevalenza nati vivi",
              "Prevalenza nati morti",
              "Prevalenza IVG")

colNames_c <- colNames

for(k in ssa){
  
  pPoints<-data.frame()
  
  esmr<-vector()
  tmp<-t(
    dataset_2016_2020_exc_gen_con[which(dataset_2016_2020_exc_gen_con[, colRegistry]==rgs[7] & 
                                          dataset_2016_2020_exc_gen_con[, colAnomaly]==k), c(1:12)]
  )
  tmp_c<-t(
    dataset_2016_2020_combined_exc_gen_con[which(dataset_2016_2020_combined_exc_gen_con[, colAnomaly]==k), c(1:12)]
  )
  if(dim(tmp)[2]>0){
    rownames(tmp)<-colNames
    rownames(tmp_c)<-colNames
    
    for(j in 1:dim(tmp)[2]){
      esmr<-c(esmr, smrCalc(tmp, tmp_c, j))
      pPoints<-rbind(pPoints, getPPoints(tmp, "IMER", j))
    }
    
    for(j in 1:dim(tmp_c)[2]){
      pPoints<-rbind(pPoints, getPPoints(tmp_c, "EUROCAT", j))
    }
    
    #makePlot(pPoints, k)
    
    tmp<-rbind(tmp, esmr)
    row.names(tmp)<-c(rownames(tmp)[1:12], "SMR")
    
    tabname=paste0(rgs[7], "_", k)
    tabname_c=paste0("Total registries", "_", k, "_cmb")
    
    tpop.local<-sum(as.numeric(tmp[4,]))
    tcases.local<-sum(as.numeric(tmp[8,]))
    
    tprev.local<-round(tcases.local/tpop.local*10000,2)
    
    tprev.local.conf.ll<-round((1.96/2 - sqrt(tcases.local + 0.02))^2/tpop.local*10000,2)
    tprev.local.conf.ul<-round((1.96/2 + sqrt(tcases.local + 0.96))^2/tpop.local*10000,2)
    
    trow<-data.frame(Anomaly=tabname,
                     Population=tpop.local,
                     Cases=tcases.local,
                     Prevalence=tprev.local,
                     Lower=tprev.local.conf.ll,
                     Upper=tprev.local.conf.ul)
    
    total.prev=rbind(total.prev, trow)
    
    #write.csv2(tmp, paste0(selectedAnomaliesDir, "/", tabname, ".csv"), col.names = F)
    #write.csv2(tmp_c, paste0(selectedAnomaliesDir, "/", tabname_c, ".csv"), col.names = F)  
  }
  
  write.csv2(total.prev, paste0(selectedAnomaliesDir, "/total.prev.csv"), row.names = F)
  
}


# SUBTABLES EXTRACTION AND PLOTS GENETIC DISORDERS ----

total.prev<-data.frame()

tabname=''
tabname_c=''
colNames <- c("Registro",
              "Anomalia",
              "Anno",
              "Popolazione",
              "Nati vivi",
              "Nati morti",
              "IVG",
              "Totali",
              "Prevalenza",
              "Prevalenza nati vivi",
              "Prevalenza nati morti",
              "Prevalenza IVG")

colNames_c <- colNames

for(k in ga){
  
  pPoints<-data.frame()
  
  esmr<-vector()
  tmp<-t(
    dataset_2016_2020_inc_gen_con[which(dataset_2016_2020_inc_gen_con[, colRegistry]==rgs[7] & 
                                          dataset_2016_2020_inc_gen_con[, colAnomaly]==k), c(1:12)]
  )
  tmp_c<-t(
    dataset_2016_2020_combined_inc_gen_con[which(dataset_2016_2020_combined_inc_gen_con[, colAnomaly]==k), c(1:12)]
  )
  if(dim(tmp)[2]>0){
    rownames(tmp)<-colNames
    rownames(tmp_c)<-colNames
    
    for(j in 1:dim(tmp)[2]){
      esmr<-c(esmr, smrCalc(tmp, tmp_c, j))
      pPoints<-rbind(pPoints, getPPoints(tmp, "IMER", j))
    }
    
    for(j in 1:dim(tmp_c)[2]){
      pPoints<-rbind(pPoints, getPPoints(tmp_c, "EUROCAT", j))
    }
    
    # makePlot(pPoints, k)
    
    tmp<-rbind(tmp, esmr)
    row.names(tmp)<-c(rownames(tmp)[1:12], "SMR")
    
    tabname=paste0(rgs[7], "_", k)
    tabname_c=paste0("Total registries", "_", k, "_cmb")
    
    tpop.local<-sum(as.numeric(tmp[4,]))
    tcases.local<-sum(as.numeric(tmp[8,]))
    
    tprev.local<-round(tcases.local/tpop.local*10000,2)
    
    tprev.local.conf.ll<-round((1.96/2 - sqrt(tcases.local + 0.02))^2/tpop.local*10000,2)
    tprev.local.conf.ul<-round((1.96/2 + sqrt(tcases.local + 0.96))^2/tpop.local*10000,2)
    
    trow<-data.frame(Anomaly=tabname,
                     Population=tpop.local,
                     Cases=tcases.local,
                     Prevalence=tprev.local,
                     Lower=tprev.local.conf.ll,
                     Upper=tprev.local.conf.ul)
   
    total.prev=rbind(total.prev, trow)
    
  #  write.csv2(tmp, paste0(geneticAnomaliesDir, "/", tabname, ".csv"))
  #  write.csv2(tmp_c, paste0(geneticAnomaliesDir, "/", tabname_c, ".csv")) 
  }
  
  write.csv2(total.prev, paste0(geneticAnomaliesDir, "/total.prev.csv"), row.names = F)
}

# PRENATAL DIAGNOSIS BY MATERNAL AGE PLOT ----

Prenatal_detection_rates_By_maternal_age_plot <- read_excel("prenatal/Prenatal.detection.rates.By.maternal.age.plot.xlsx")
u35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Under 35`, " "), "[[", 1))))
eoo35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`35 and over`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_By_maternal_age_plot$Registry, 3)
a<-c(rep("1", length(u35)), rep("2", length(eoo35)), rep("3", length(monk)))
p<-c(u35,eoo35,monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_By_maternal_age_plot$Registry,
                      Age=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_maternal_age", ".jpg"), 
     width = 6, height = 4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Age, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=31.19), color="goldenrod1", lty=2) +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Classe di età",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Meno di 35 anni", "35 anni o più", "Non nota"))

dev.off()


# PRENATAL DIAGNOSIS BY OUTCOME PLOT ----

Prenatal_detection_rates_by_outcome_plot <- read_excel("prenatal/Prenatal.detection.rates.by.outcome.2016.2020.plot.xlsx")
nv<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Live births`, " "), "[[", 1))))
nm<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Fetal deaths`, " "), "[[", 1))))
ivg<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$TOPFAs, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_outcome_plot$Registry, 3)
a<-c(rep("1", length(nv)), rep("2", length(nm)), rep("3", length(ivg)))
p<-c(nv,nm,ivg)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_outcome_plot$Registry,
                      Outcome=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_outcome", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Outcome, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=31.19), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Esito",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("Nati vivi", "Nati morti", "IVG"))

dev.off()


# PRENATAL DIAGNOSIS BY GESTATIONAL AGE PLOT ----

Prenatal_detection_rates_by_gest_age_plot <- read_excel("prenatal/Prenatal.detection.rates.by.gest.age.2016.2020.plot.xlsx")
l14<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Less 14 weeks`, " "), "[[", 1))))
b1423<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`14-23 weeks`, " "), "[[", 1))))
m23<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`More 23 weeks`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_gest_age_plot$Registry, 3)
a<-c(rep("1", length(l14)), rep("2", length(b1423)), rep("3", length(m23)), rep("4", length(monk)))
p<-c(l14,b1423,m23, monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_gest_age_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_gest_age", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=31.19), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Età gestazionale",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("< 14 settimane", "14 - 23 settimane", "> 23 settimane", "Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY INDICATION PLOT ----

Prenatal_detection_rates_by_indication_plot <- read_excel("prenatal/Prenatal.detection.rates.by.indication.2016.2020.plot.xlsx")
s<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Screening, " "), "[[", 1))))
u<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Ultrasound, " "), "[[", 1))))
ao<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Age and other`, " "), "[[", 1))))
m<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_indication_plot$Registry, 4)
a<-c(rep("1", length(s)), rep("2", length(u)), rep("3", length(ao)), rep("4", length(m)))
p<-c(s,u,ao,m)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_indication_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_indication", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Indicazione",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Screening", "Ultrasuoni", "Età", "Mancante/Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY MATERNAL AGE PLOT DOWN SYNDROME ----

Prenatal_detection_rates_By_maternal_age_plot <- read_excel("prenatal/Prenatal_detection_rates_By_maternal_age_2016_2020_down.plot.xlsx")
u35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Under 35`, " "), "[[", 1))))
eoo35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`35 and over`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_By_maternal_age_plot$Registry, 3)
a<-c(rep("1", length(u35)), rep("2", length(eoo35)), rep("3", length(monk)))
p<-c(u35,eoo35,monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_By_maternal_age_plot$Registry,
                      Age=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_maternal_age_Down", ".jpg"), 
     width = 6, height = 4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Age, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=81.03), color="goldenrod1", lty=2) +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Classe di età",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Meno di 35 anni", "35 anni o più", "Non nota"))

dev.off()


# PRENATAL DIAGNOSIS BY OUTCOME PLOT DOWN SYNDROME ----

Prenatal_detection_rates_by_outcome_plot <- read_excel("prenatal/Prenatal_detection_rates_By_outcome_2016_2020_Down_syn.plot.xlsx")
nv<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Live births`, " "), "[[", 1))))
nm<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Fetal deaths`, " "), "[[", 1))))
ivg<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$TOPFAs, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_outcome_plot$Registry, 3)
a<-c(rep("1", length(nv)), rep("2", length(nm)), rep("3", length(ivg)))
p<-c(nv,nm,ivg)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_outcome_plot$Registry,
                      Outcome=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_outcome_Down", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Outcome, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=81.03), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Esito",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("Nati vivi", "Nati morti", "IVG"))

dev.off()


# PRENATAL DIAGNOSIS BY GESTATIONAL AGE PLOT DOWN SYNDROME ----

Prenatal_detection_rates_by_gest_age_plot <- read_excel("prenatal/Prenatal_detection_rates_By_gestational_age_at_discovery_2016_2020_down.plot.xlsx")
l14<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Less 14 weeks`, " "), "[[", 1))))
b1423<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`14-23 weeks`, " "), "[[", 1))))
m23<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`More 23 weeks`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_gest_age_plot$Registry, 3)
a<-c(rep("1", length(l14)), rep("2", length(b1423)), rep("3", length(m23)), rep("4", length(monk)))
p<-c(l14,b1423,m23, monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_gest_age_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_gest_age_Down", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=81.03), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Età gestazionale",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("< 14 settimane", "14 - 23 settimane", "> 23 settimane", "Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY INDICATION PLOT DOWN SYNDROME ----

Prenatal_detection_rates_by_indication_plot <- read_excel("prenatal/Prenatal.detection.rates.by.indication.2016.2020.plot.xlsx")
s<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Screening, " "), "[[", 1))))
u<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Ultrasound, " "), "[[", 1))))
ao<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Age and other`, " "), "[[", 1))))
m<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_indication_plot$Registry, 4)
a<-c(rep("1", length(s)), rep("2", length(u)), rep("3", length(ao)), rep("4", length(m)))
p<-c(s,u,ao,m)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_indication_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_indication", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Indicazione",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Screening", "Ultrasuoni", "Età", "Mancante/Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY MATERNAL AGE PLOT EDWARD SYNDROME ----

Prenatal_detection_rates_By_maternal_age_plot <- read_excel("prenatal/Prenatal_detection_rates_By_maternal_age_2016_2020.Edward.plot.xlsx")
u35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Under 35`, " "), "[[", 1))))
eoo35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`35 and over`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_By_maternal_age_plot$Registry, 3)
a<-c(rep("1", length(u35)), rep("2", length(eoo35)), rep("3", length(monk)))
p<-c(u35,eoo35,monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_By_maternal_age_plot$Registry,
                      Age=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_maternal_age_Edward", ".jpg"), 
     width = 6, height = 4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Age, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=100.00), color="goldenrod1", lty=2) +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Classe di età",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Meno di 35 anni", "35 anni o più", "Non nota"))

dev.off()


# PRENATAL DIAGNOSIS BY OUTCOME PLOT EDWARD SYNDROME ----

Prenatal_detection_rates_by_outcome_plot <- read_excel("prenatal/Prenatal_detection_rates_By_outcome_2016_2020_edward.plot.xlsx")
nv<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Live births`, " "), "[[", 1))))
nm<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Fetal deaths`, " "), "[[", 1))))
ivg<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$TOPFAs, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_outcome_plot$Registry, 3)
a<-c(rep("1", length(nv)), rep("2", length(nm)), rep("3", length(ivg)))
p<-c(nv,nm,ivg)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_outcome_plot$Registry,
                      Outcome=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_outcome_Edward", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Outcome, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=100.00), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Esito",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("Nati vivi", "Nati morti", "IVG"))

dev.off()


# PRENATAL DIAGNOSIS BY GESTATIONAL AGE PLOT EDWARD SYNDROME ----

Prenatal_detection_rates_by_gest_age_plot <- read_excel("prenatal/Prenatal_detection_rates_By_gestational_age_at_discovery_2016_2020.Edward.plot.xlsx")
l14<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Less 14 weeks`, " "), "[[", 1))))
b1423<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`14-23 weeks`, " "), "[[", 1))))
m23<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`More 23 weeks`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_gest_age_plot$Registry, 3)
a<-c(rep("1", length(l14)), rep("2", length(b1423)), rep("3", length(m23)), rep("4", length(monk)))
p<-c(l14,b1423,m23, monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_gest_age_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_gest_age_Edward", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=100.00), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Età gestazionale",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("< 14 settimane", "14 - 23 settimane", "> 23 settimane", "Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY INDICATION PLOT EDWARD SYNDROME ----

Prenatal_detection_rates_by_indication_plot <- read_excel("prenatal/Prenatal.detection.rates.by.indication.2016.2020.plot.xlsx")
s<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Screening, " "), "[[", 1))))
u<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Ultrasound, " "), "[[", 1))))
ao<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Age and other`, " "), "[[", 1))))
m<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_indication_plot$Registry, 4)
a<-c(rep("1", length(s)), rep("2", length(u)), rep("3", length(ao)), rep("4", length(m)))
p<-c(s,u,ao,m)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_indication_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_indication", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Indicazione",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Screening", "Ultrasuoni", "Età", "Mancante/Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY MATERNAL AGE PLOT PATAU SYNDROME ----

Prenatal_detection_rates_By_maternal_age_plot <- read_excel("prenatal/Prenatal_detection_rates_By_maternal_age_2016_2020.patau.plot.xlsx")
u35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Under 35`, " "), "[[", 1))))
eoo35<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`35 and over`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_By_maternal_age_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_By_maternal_age_plot$Registry, 3)
a<-c(rep("1", length(u35)), rep("2", length(eoo35)), rep("3", length(monk)))
p<-c(u35,eoo35,monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_By_maternal_age_plot$Registry,
                      Age=a,
                      Perc=p)

sel<-which(is.na(plotTable$Perc))
plotTable<-plotTable[-sel,]

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_maternal_age_Patau", ".jpg"), 
     width = 6, height = 4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Age, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=100.00), color="goldenrod1", lty=2) +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Classe di età",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Meno di 35 anni", "35 anni o più", "Non nota"))

dev.off()


# PRENATAL DIAGNOSIS BY OUTCOME PLOT PATAU SYNDROME ----

Prenatal_detection_rates_by_outcome_plot <- read_excel("prenatal/Prenatal_detection_rates_By_outcome_2016_2020_Patau.plot.xlsx")
nv<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Live births`, " "), "[[", 1))))
nm<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Fetal deaths`, " "), "[[", 1))))
ivg<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$TOPFAs, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_outcome_plot$Registry, 3)
a<-c(rep("1", length(nv)), rep("2", length(nm)), rep("3", length(ivg)))
p<-c(nv,nm,ivg)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_outcome_plot$Registry,
                      Outcome=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_outcome_Patau", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Outcome, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=100.00), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Esito",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("Nati vivi", "Nati morti", "IVG"))

dev.off()


# PRENATAL DIAGNOSIS BY GESTATIONAL AGE PLOT PATAU SYNDROME ----

Prenatal_detection_rates_by_gest_age_plot <- read_excel("prenatal/Prenatal_detection_rates_By_gestational_age_at_discovery_2016_2020.patau.plot.xlsx")
l14<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Less 14 weeks`, " "), "[[", 1))))
b1423<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`14-23 weeks`, " "), "[[", 1))))
m23<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`More 23 weeks`, " "), "[[", 1))))
monk<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_gest_age_plot$`Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_gest_age_plot$Registry, 3)
a<-c(rep("1", length(l14)), rep("2", length(b1423)), rep("3", length(m23)), rep("4", length(monk)))
p<-c(l14,b1423,m23, monk)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_gest_age_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_gest_age_Patau", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(aes(yintercept=100.00), color="goldenrod1", lty=2) +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Età gestazionale",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("< 14 settimane", "14 - 23 settimane", "> 23 settimane", "Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY INDICATION PLOT PATAU SYNDROME ----

Prenatal_detection_rates_by_indication_plot <- read_excel("prenatal/Prenatal.detection.rates.by.indication.2016.2020.plot.xlsx")
s<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Screening, " "), "[[", 1))))
u<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$Ultrasound, " "), "[[", 1))))
ao<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Age and other`, " "), "[[", 1))))
m<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_indication_plot$`Missing/Not known`, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_indication_plot$Registry, 4)
a<-c(rep("1", length(s)), rep("2", length(u)), rep("3", length(ao)), rep("4", length(m)))
p<-c(s,u,ao,m)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_indication_plot$Registry,
                      GestAge=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_indication", ".jpg"), 
     width = 6, height =4, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=GestAge, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  #coord_flip() +
  ylab("N(%)") +
  xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Indicazione",
                    values=c("#56B4E9", "darkgreen", "black", "salmon"),
                    labels=c("Screening", "Ultrasuoni", "Età", "Mancante/Non noto"))

dev.off()


# PRENATAL DIAGNOSIS BY OUTCOME PLOT MALFORMAZIONI SELEZIONATE ----

Prenatal_detection_rates_by_outcome_plot <- read_excel("prenatal/Prenatal_detecetion_rates_by_outcome_ER_total.plot.xlsx")
nv<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Live births`, " "), "[[", 1))))
nm<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$`Fetal deaths`, " "), "[[", 1))))
ivg<-as.numeric(gsub("%","",unlist(lapply(strsplit(Prenatal_detection_rates_by_outcome_plot$TOPFAs, " "), "[[", 1))))

r<-rep(Prenatal_detection_rates_by_outcome_plot$Anomaly, 3)
a<-c(rep("1", length(nv)), rep("2", length(nm)), rep("3", length(ivg)))
p<-c(nv,nm,ivg)

plotTable<-data.frame(Registry=Prenatal_detection_rates_by_outcome_plot$Anomaly,
                      Outcome=a,
                      Perc=p)

jpeg(filename = paste0(plotDir, "/Prenatal_detection_rate_by_outcome_selected_malformations", ".jpg"), 
     width = 6, height = 3, units = "in", pointsize = 12,
     res = 300)

ggplot(plotTable, aes(fill=Outcome, y=Perc, x=Registry)) + 
  geom_bar(position="stack", stat="identity") +
  #geom_hline(aes(yintercept=100.00), color="goldenrod1", lty=2) +
  coord_flip() +
  ylab("N(%)") +
  xlab("") +
  #theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(#colour = c(rep("gray50", 5), "firebrick3", rep("grey50", 24)),
                                   angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="Esito",
                    values=c("#56B4E9", "darkgreen", "salmon", "black"),
                    labels=c("Nati vivi", "Nati morti", "IVG"))

dev.off()


# IMER CEDAP 2016.2020 TABLES ----
# rms

# ddist<-datadist(dset.descr)
# options(datadist=ddist)
# 
# tab1<-summary(Cittadinanza ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
# tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
# write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.BY_GENERE.rms.csv"))

# Table1

dset.descr<-Imer2016.2020[,c(
  "ChiaveIMER",
  "TipoDiNascita",
  "NatoSingoloOGemello",
  "Sesso",
  "StatoCivileMadre",
  "ScolaritaMadre",
  "AlcoolMadre",
  "FumoMadre",
  "AlcoolInGravidanza",
  "FumoInGravidanza",
  "FarmaciMadre",
  "DrogheMadre",
  "RadiazioniMadre",
  "AccertamentiSierologici",
  "Vaccinazione",
  "VirosiMadre",
  "TotaleGravidanzeATermine",
  "TotaleGravidanzePretermine",
  "TotaleNatiMorti",
  "TotaleGravidanzeExtraUterine",
  "TotaleAbortiSpontanei",
  "TotaleAbortiProvocati",
  "TotaleMalformati",
  "TotaleGemelli",
  "DiagnosiPrenatale",
  "Amnio",
  "Villo",
  "AltraDP",
  "Cittadinanza",
  "Gemellarita",
  "SettimaneDiGestazione",
  "TipoParto",
  "PresentazioneParto",
  "LiquidoAmniotico",
  "Placenta",
  "PesoNeonato",
  "EtaMadre"
  #"SettimaneDiagn"
)]

sel<-which(dset.descr$ChiaveIMER=="H")
dset.descr[sel, "ChiaveIMER"]<-"I"
sel<-which(dset.descr$Cittadinanza!=100)
dset.descr[sel, "Cittadinanza"]<-"extraIT"
sel<-which(dset.descr$Cittadinanza==100)
dset.descr[sel, "Cittadinanza"]<-"IT"
dset.descr$Cittadinanza<-as.factor(dset.descr$Cittadinanza)

fact<-colnames(dset.descr)[c(1:16, 17:24, 25:30, 32:35)]

nnorm<-colnames(dset.descr)[c(31, 36:37)]

# Comparison chiave imer
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[1],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, nonnormal=nnorm) #exact = fact
write.csv2(tabMat, file=paste0(outdir,"/descriptive.chiave.imer.2016.2020.csv"))

# Comparison cittadinanza
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[29],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, nonnormal=nnorm) #exact = fact
write.csv2(tabMat, file=paste0(outdir,"/descriptive.cittadinanza.imer.2016.2020.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.imer.overall.2016.2020.csv"))

# IMER CEDAP 2016 CHIAVE IMER  ----
# rms

# ddist<-datadist(dset.descr)
# options(datadist=ddist)
# 
# tab1<-summary(Cittadinanza ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
# tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
# write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.BY_GENERE.rms.csv"))

# Table1

dset.descr<-Imer2016.2020[,c(
  "ChiaveIMER",
  "Anno"
)]

sel<-which(dset.descr$ChiaveIMER=="H")
dset.descr[sel, "ChiaveIMER"]<-"I"

fact<-colnames(dset.descr)[c(1, 2)]

# Comparison chiave imer
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[1],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE) #exact = fact
write.csv2(tabMat, file=paste0(outdir,"/descriptive.chiave.imer.Anno.2016.2020.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.imer.Anno.overall.2016.2020.csv"))

# IMER CEDAP 2016 CHIAVE IMER  DIAGNOSI PRENATALE ----
# rms

# ddist<-datadist(dset.descr)
# options(datadist=ddist)
# 
# tab1<-summary(Cittadinanza ~ ., method = 'rev', overall = T, test = T, data = dset.descr)
# tab1.printed<-print(tab1, prtest = c('P'), npct = c('numerator'), prmsd = TRUE, digits = 2, prType = "latex")
# write.csv2(tab1.printed, file=paste0(outdir, "/descriptive.BY_GENERE.rms.csv"))

# Table1

dset.descr<-Imer2016.2020[,c(
  "ChiaveIMER",
  "DiagnosiPrenatale"
)]

sel<-which(dset.descr$ChiaveIMER=="H")
dset.descr[sel, "ChiaveIMER"]<-"I"

fact<-colnames(dset.descr)[c(1, 2)]

# Comparison chiave imer
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    strata = colnames(dset.descr)[1],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE) #exact = fact
write.csv2(tabMat, file=paste0(outdir,"/descriptive.chiave.imer.diag.pren.Anno.2016.2020.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr),
                    factorVars = fact,
                    data = dset.descr,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.imer.diag.pren.Anno.overall.2016.2020.csv"))

# NO SDO

dset.descr.no.sdo<-Imer2016.2020[which(Imer2016.2020$Numero<2000 | Imer2016.2020$Numero>4999)
  ,c(
  "ChiaveIMER",
  "DiagnosiPrenatale"
)]

sel<-which(dset.descr.no.sdo$ChiaveIMER=="H")
dset.descr.no.sdo[sel, "ChiaveIMER"]<-"I"

fact<-colnames(dset.descr.no.sdo)[c(1, 2)]

# Comparison chiave imer
tab<-CreateTableOne(vars = colnames(dset.descr.no.sdo),
                    factorVars = fact,
                    data = dset.descr.no.sdo,
                    strata = colnames(dset.descr.no.sdo)[1],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE) #exact = fact
write.csv2(tabMat, file=paste0(outdir,"/descriptive.chiave.imer.diag.pren.no.sdo.Anno.2016.2020.csv"))

# Overall
tab<-CreateTableOne(vars = colnames(dset.descr.no.sdo),
                    factorVars = fact,
                    data = dset.descr.no.sdo,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0(outdir, "/descriptive.imer.diag.pren.no.sdo.Anno.overall.2016.2020.csv"))

# CEDAP MAMME 2021 TABLES ----

fact_mamme_2021_RER<-c(3:15, 21:34 )
nnorm_mamme_2021_RER<-c(1, 2, 16:20)

# Overall
tab<-CreateTableOne(vars = colnames(mamme_2021_RER),
                    factorVars = fact_mamme_2021_RER,
                    data = mamme_2021_RER,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm_mamme_2021_RER)
write.csv2(tabMat, file=paste0("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/cedap/descriptive.overall.mamme.cedap.csv"))

# CEDAP 2016 - 2021 TABLES ----

fact<-colnames(cedap.2016.2020)[c(1,6:9, 12:22, 24:32, 34:36)]
nnorm<-colnames(cedap.2016.2020)[c(2:5, 10, 11, 23, 33)]

# Overall
tab<-CreateTableOne(vars = colnames(cedap.2016.2020),
                    factorVars = fact,
                    data = cedap.2016.2020,
                    test = FALSE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, exact = fact, nonnormal=nnorm)
write.csv2(tabMat, file=paste0("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/cedap/descriptive.overall.csv"))

# Comparison
tab<-CreateTableOne(vars = colnames(cedap.2016.2020),
                    factorVars = fact,
                    data = cedap.2016.2020,
                    strata = colnames(cedap.2016.2020)[1],
                    test = TRUE,
                    includeNA = FALSE)
tabMat <- print(tab, quote = FALSE, noSpaces = TRUE, printToggle = TRUE, nonnormal=nnorm) #exact = fact
write.csv2(tabMat, file=paste0("~/Library/CloudStorage/GoogleDrive-mnfmrc@unife.it/Drive condivisi/IMER/documenti/report/2023/cedap/descriptive.byYear.csv"))

# CONTRIBUTO FONTI DATI 2021 ----

library(ggrepel)

# Schede IMER <1000
# CedAP 1000 < 2000
# SDO 2000 < 4000
# SDO7+ 4000 < 5000
# SDO confermate con Scheda IMER dai referenti 5000 < 6000
# Rare >=6000

# IMER = 4
# CEDAP = 3
# SDO = 1
# SDO+7 = 2
# SDO confermate con Scheda IMER dai referenti = 1
# RARE = 6
# IVG = 5

Imer2016.2020$SOURCE=NA

sel<-which(Imer2016.2020$Numero<1000)
Imer2016.2020[sel, "SOURCE"]<-4

sel<-which(Imer2016.2020$Numero>=1000 & Imer2016.2020$Numero<2000)
Imer2016.2020[sel, "SOURCE"]<-3

sel<-which(Imer2016.2020$Numero>=2000 & Imer2016.2020$Numero<4000)
Imer2016.2020[sel, "SOURCE"]<-1

sel<-which(Imer2016.2020$Numero>=4000 & Imer2016.2020$Numero<5000)
Imer2016.2020[sel, "SOURCE"]<-2

sel<-which(Imer2016.2020$Numero>=5000 & Imer2016.2020$Numero<6000)
Imer2016.2020[sel, "SOURCE"]<-1

sel<-which(Imer2016.2020$Numero>=6000)
Imer2016.2020[sel, "SOURCE"]<-6

sel<-which(Imer2016.2020$TipoDiNascita==4)
Imer2016.2020[sel, "SOURCE"]<-5
table(Imer2016.2020$SOURCE)

Imer2016.2020$SOURCE<-factor(Imer2016.2020$SOURCE, levels=c(1,2,3,4,5,6))

relTab<-data.frame(prop.table(table(Imer2016.2020$SOURCE)))
#relTab<-relTab[-3,]
relTab$Var1<-droplevels(relTab$Var1)

relTablabs <- relTab %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

jpeg(filename = paste0(plotDir, "/contributo.fonti.dati", ".jpg"), 
     width = 4, height = 4, units = "in", pointsize = 12,
     res = 300)

p<-ggplot(relTab, aes(x="", y=Freq, fill=Var1)) + 
  geom_bar(stat="identity", width=1, color="white") +
  geom_label_repel(data = relTablabs,
                   aes(y = pos, label = paste0(round(Freq*100,2), "%")),
                   segment.colour="black", size = 4, nudge_x = 1, show.legend = FALSE) +
  coord_polar("y", start=0) +
  scale_fill_manual(labels = c("SDO", "SDO >7d", "CedAP", "IMER", "IVG", "Registro malattie rare"), 
                    values = c("lightskyblue3", "lightblue", "orange", "mediumpurple2", "darkolivegreen3", "indianred1")) +
  ggtitle("Contributi delle fonti di accertamento", subtitle = "2016 - 2020") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        # #panel.grid.major = element_blank(),  
        # panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom"
  )
p

dev.off()



