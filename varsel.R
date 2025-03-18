# VARSEL v 1.0

source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

library(cluster)
library(ClustOfVar)
library(dplyr)
library(ggplot2)
library(readr)
library(Rtsne)
library(factoextra)

#library(devtools)
#install_github("DataSlingers/ExclusiveLasso")

# PATH ----

wd<-setwd("~/Analisi/Manfrini/PHENO")
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out/varsel")

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

# CLUSTVAR ----

dset.varclust<-dset.ana[, nm[which(nm$noas==1), "id"]]

X.quanti <- dset.varclust[,getNumbersFromDataset(dset.varclust)]
X.quali <- dset.varclust[,getFactorsFromDataset(dset.varclust)]
tree <- hclustvar(X.quanti, X.quali)
plot(tree)

stab<-stability(tree,B=20)
plot(stab,main="Stability of the partitions")

part<-cutreevar(tree,27) #cut of the tree

table(part$cluster)

# GRPLASSO ----

library(gglasso)

outcome<-ifelse(dset.ana$Outcome==5,1,0)

gr<-gglasso(data.matrix(dset.varclust), outcome, lambda = 0.2,
              group = part$cluster, loss="ls",
              intercept = F)

# cross validation

gr_cv<-cv.gglasso(x=data.matrix(dset.varclust), y=outcome, group=part$cluster, 
                    pred.loss="L2", 
                    intercept = F, nfolds=5)

plot(gr_cv)
paste(gr_cv$lambda.min, gr_cv$lambda.1se)

# selected model

gr<-gglasso(data.matrix(dset.varclust), factor(outcome), lambda = gr_cv$lambda.1se+0.1,
              group = part$cluster, loss="ls",
              intercept = F)

write.csv2(cbind(rownames(gr$beta), gr$beta), file=paste0(outdir, "/gglasso.csv"), row.names = F)

dset.clust<-dset.ana[, c(2,
                         3,
                         16:20,
                         23,
                         29,
                         36,
                         37,
                         40,
                         42,
                         44,
                         45,
                         46)]

# ------------------------------------------------------------------------------------------------

# NBCLUST ----



library(NbClust)

gower_dist <- daisy(dset.clust, metric = "gower")
gower_mat <- as.matrix(gower_dist)

#nb.clust<-NbClust(data = dset.clust, diss = gower_mat, distance = NULL, min.nc = 2, max.nc = 5,
#        method = "ward.D", index = "silhouette", alphaBeale = 0.1)

nb.clust<-NbClust(data = dset.clust, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 4, 
                    method = "ward.D", index = "all", alphaBeale = 0.1)

print(nb.clust$Best.nc)

# PAM CLUSTERING ----
# Compute Gower distance
gower_dist <- daisy(dset.clust, metric = "gower")
gower_mat <- as.matrix(gower_dist)

# Compute euclidean distance
euc_dist <- daisy(dset.clust, metric = "euclidean")
euc_mat <- as.matrix(euc_dist)

# SERACH FOR CLUSTER NUMBER WITH SILHOUETTE METHOD
sil_width <- c(NA)
for(i in 2:8){  
  pam_fit <- pam(gower_dist, diss = TRUE, k = i)  
  sil_width[i] <- pam_fit$silinfo$avg.width  
}

plot(1:8, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:8, sil_width)

# CLUSTERS SUMMARY (BEST K=2)
k <- 3
pam_fit <- cluster::pam(euc_dist, k, diss = T)
#plot(pam_fit)
dset.ana$cluster<-pam_fit$clustering

write.csv2(dset.ana, file=paste0(outdir, "/dset.ana.cluster.last.csv"), na="", row.names = F)

pam_results <- dset.ana %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary

# VISUALIZATION
tsne_obj <- Rtsne(euc_dist, is_distance = TRUE)
tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering))
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

pam_fit$data<-dset.clust

fviz_cluster(pam_fit, geom = "point", ellipse.type = "norm", stand = F)

#fviz_nbclust(gower_mat, pam, method="wss")

#fviz_nbclust(gower_mat, pam, method="silhouette")

# CLUSTER DESCRIPTION ----

fviz_dend(pam_fit,
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



