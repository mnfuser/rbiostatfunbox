# GMM 1.0

source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")


library(ClusterR)
library(dplyr)
library(ggplot2)
library(readr)
library(Rtsne)

set.seed(1234)

# PATH ----

wd<-setwd("~/Analisi/Manfrini/PHENO")
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out/pca")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana <- read.csv2(paste0(indir,"/dset.ana.imp.1.csv"),
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

# GMM ----

sel<-c(1, 26:35, 41, 43, 52:54)
dset.gmm<-(dset.ana[,-sel])

# VISUALIZE ANALYSIS DATASET
vis_dat_wrapper(dset.gmm)

getFactorsFromDataset(dset.gmm)
dset.gmm[,factorList]<-apply(dset.gmm[, factorList], 2, as.numeric)

# VISUALIZE ANALYSIS DATASET
vis_dat_wrapper(dset.gmm)

# scaling numeric vars
getNumbersFromDataset(dset.gmm)
dset.gmm[, numericList]<-scale(dset.gmm[, numericList])

# optimal cluster

opt_gmm = Optimal_Clusters_GMM(dset.gmm, max_clusters = 3, criterion = "BIC", 
                               
                               dist_mode = "maha_dist", seed_mode = "random_subset",
                               
                               km_iter = 10, em_iter = 10, var_floor = 1e-10, 
                               
                               plot_data = T)

gmm = GMM(dset.gmm, 3, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10,
          em_iter = 10, verbose = F)          

# predict centroids, covariance matrix and weights
pr = predict(gmm, newdata = dset.gmm)

dset.ana$cluster<-pr

write.csv2(dset.ana, file=paste0(outdir, "/dset.ana.cluster.gmm.csv"))

# Perf

res = external_validation(as.numeric(dset.ana$Outcome), pr, 
                          
                          method = "adjusted_rand_index", summary_stats = T)

res



# CLUSTER DESCRIPTION ----

fviz_dend(hcpc,
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = T, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8,      # Augment the room for labels
          show_labels = F,
          type = "rectangle",
          horiz=T,
          main="Patients grouped by cluster",
          xlab = "Cluster"
)

fviz_cluster(hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

plot(hcpc, choice = "3D.map")

library(ggplot2)
library(tidytext)

# Basic barplot
cluster.names <- as_labeller(c(`1` = "Cluster 1", `2` = "Cluster 2", `3` = "Cluster 3", `4` = "Cluster 4")) 

# All clusters 


hcpc.vtest.filtered<-hcpc.vtest
hcpc.vtest.filtered<-hcpc.vtest.filtered[which(hcpc.vtest.filtered$p.value<0.0001),]
hcpc.vtest.filtered<-hcpc.vtest.filtered[which(hcpc.vtest.filtered$v.test>quantile(hcpc.vtest.filtered$v.test)[4]
                                               | hcpc.vtest.filtered$v.test<quantile(hcpc.vtest.filtered$v.test)[2]),]

hcpc.vtest.filtered$var=c("Pneumonia (mild)",
                          "Pneumonia pripheral pattern",
                          "Pnumonia (severe)",
                          "Pneumonia (moderate)",
                          "Pneumonia diffuse pattern",
                          "Pattern.1.central.2.peripheral.3.diffuse=0",
                          "Pneumonia.1.MILD.2.MODERATE.3.SEVERE=0",
                          "Pneumonia diffuse pattern",
                          "Pneumonia diffuse pattern",
                          "Pneumonia (moderate)",
                          "No PTCA",
                          "No cardiopathy",
                          "Pnumonia (severe)",
                          "Cardiopathy",
                          "PTCA",
                          "Pneumonia diffuse pattern",
                          "Pneumonia (mild)",
                          "Cardiopathy",
                          "PTCA",
                          "Coronaropathy",
                          "CKD",
                          "Hypertension",
                          "Peripheral vasculopathy",
                          "No diabetes",
                          "No BPAC",
                          "No peripheral vasculopathy",
                          "No hypertension",
                          "No CKD",
                          "No coronaropathy",
                          "No PTCA",
                          "No cardiopathy",
                          "Pulmonary breathing volume",
                          "PaO2 worst",
                          "SatO2 at hospitalization",
                          "Age",
                          "LDH at hospitalization",
                          "Number of lobes",
                          "Number of lobes",
                          "Number of lobes",
                          "LDH at hospitalization",
                          "SatO2 at hospitalization",
                          "PaO2 worst",
                          "TA calcium volume",
                          "CSV LAD",
                          "CSV TC",
                          "CSV CDX",
                          "CSV total volume",
                          "Pulmonary breathing volume",
                          "CSV total volume",
                          "CSV CDX",
                          "CSV TC",
                          "CSV LAD",
                          "TA calcium volume",
                          "CSV LCX",
                          "Creatinine at hospitalization",
                          "Age",
                          "AV calcium volume",
                          "EF basal"
)

write.csv2(hcpc.vtest.filtered, file=paste0(outdir, "/hcpc.vtest.filtered.csv"), row.names = F)

#hcpc.vtest.filtered<-hcpc.vtest.filtered[-selected,]

hcpc.vtest.filtered$var<-gsub("=Yes", "", hcpc.vtest.filtered$var)
hcpc.vtest.filtered$var<-gsub("..mm.", "", hcpc.vtest.filtered$var)
hcpc.vtest.filtered$var<-gsub("[ms.]", "", hcpc.vtest.filtered$var)
hcpc.vtest.filtered$var<-gsub("[.]", " ", hcpc.vtest.filtered$var)
hcpc.vtest.filtered$var<-gsub("_", " ", hcpc.vtest.filtered$var)

# All clusters
jpeg(filename = paste0(outdir, "/clusters.variance.jpg"), height=4, width=8, unit="in", res=600)

p<-ggplot(data=hcpc.vtest.filtered, aes(reorder(var, v.test), v.test)) + #aes(var, v.test)
  geom_bar(stat="identity", aes(fill = v.test > 0 ),  width = 0.8) + 
  facet_wrap(~CLUSTER, labeller = cluster.names) +
  xlab("Variable") +
  ylab("Variance contribution") +
  theme_pubr() +
  theme(text = element_text(size=7),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1),
        legend.position='none',
        panel.grid.major.y = element_line()
  ) +
  scale_fill_manual(values=c("black", "gray40"))
# Horizontal bar plot
p +  ggtitle("Variance contribution within clusters") # coord_flip() 

dev.off()

# cluster 1
cluster1<-hcpc.vtest[which(hcpc.vtest$CLUSTER==1),]

jpeg(filename = paste0(outdir, "/clusters1.variance.jpg"), height=8, width=4, unit="in", res=600)

p<-ggplot(data=cluster1, aes(reorder(var, v.test), v.test)) +
  geom_bar(stat="identity") +
  #facet_wrap(~CLUSTER, labeller = cluster.names) +
  xlab("Variable") +
  ylab("Variance") +
  theme(text = element_text(size=6), 
        axis.text.x = element_text(angle=90, hjust=1))
p
# Horizontal bar plot
p + coord_flip() + ggtitle("Cluster 1")

dev.off()

# cluster 2
cluster2<-hcpc.vtest[which(hcpc.vtest$CLUSTER==2),]

jpeg(filename = paste0(outdir, "/clusters2.variance.jpg"), height=8, width=4, unit="in", res=600)

p<-ggplot(data=cluster2, aes(reorder(var, v.test), v.test)) +
  geom_bar(stat="identity") +
  #facet_wrap(~CLUSTER, labeller = cluster.names) +
  xlab("Variable") +
  ylab("Variance") +
  theme(text = element_text(size=6), 
        axis.text.x = element_text(angle=90, hjust=1))
p
# Horizontal bar plot
p + coord_flip()+ ggtitle("Cluster 2")

dev.off()

# cluster 3
cluster3<-hcpc.vtest[which(hcpc.vtest$CLUSTER==3),]

jpeg(filename = paste0(outdir, "/clusters3.variance.jpg"), height=8, width=4, unit="in", res=600)

p<-ggplot(data=cluster3, aes(reorder(var, v.test), v.test)) +
  geom_bar(stat="identity") +
  #facet_wrap(~CLUSTER, labeller = cluster.names) +
  xlab("Variable") +
  ylab("Variance") +
  theme(text = element_text(size=6), 
        axis.text.x = element_text(angle=90, hjust=1))
p
# Horizontal bar plot
p + coord_flip()+ ggtitle("Cluster 3")

dev.off()

hcpc$desc.var$quanti
hcpc$desc.var$category


# ANOVA approach
summary(aov(dset.ana$Age~dset.ana$cluster))
#pairwise.t.test(dset.ana$Age, dset.ana$cluster, p.adj = "bonf")
TukeyHSD(aov(dset.ana$Age~factor(dset.ana$cluster)))

# 




