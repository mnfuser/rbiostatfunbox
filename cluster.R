# CLUSTER TEMPLATE v1

source("D:/codice/rbiostatfunbox/funbox.R")

# PATH ----

wd<-setwd("H:/Drive condivisi/Lavori in corso/Manfrini/PHENO")
setwd(wd)
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out/cluster.testpipe")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana <- read.csv2(paste0(indir,"/dset.ana.csv"),
                        stringsAsFactors=FALSE) #, na.strings = ""
}

# Load dataset descriptor file
nm<-read.csv2(paste0(indir,"/nm.dset.ana.filtered.last.csv"), 
              stringsAsFactors=T)

selected<-nm[which(nm$rg==1 | nm$oc==1), "id"]
dset.ana<-dset.ana[,c(selected)]

# SET VAR TYPE ----

setVarType(paste0(indir, "/nm.dset.ana.filtered.last.csv"))

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


# CLUSTERING ----
library(cluster)
library(dplyr)
library(ggplot2)
library(readr)
library(Rtsne)
library(factoextra)

dset.cluster<-dset.ana[,c(1:15, 17, 18, 19, 20, 23, 25, 27, 29, 31, 33:77)]
#dset.cluster<-dset.ana

# Compute Gower distance
gower_dist <- daisy(dset.cluster, metric = "gower")
gower_mat <- as.matrix(gower_dist)

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
k <- 2
pam_fit <- pam(gower_dist, diss = TRUE, k)
dset.cluster$cluster<-pam_fit$clustering
dset.ana$cluster<-pam_fit$clustering

write.csv2(dset.ana, file=paste0(indir, "/dset.ana.cluster.last.csv"), na="", row.names = F)

pam_results <- dset.cluster %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary

# VISUALIZATION

# tsne

tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering))
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))

# fviz

# PAM clustering

# Visualize pam clustering
fviz_cluster(pam_fit, geom = "point", ellipse.type = "norm")

# Hierarchical clustering

# Use hcut() which compute hclust and cut the tree
hc.cut <- hcut(dset.cluster, k = 2, hc_method = "complete")
# Visualize dendrogram
fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE)
# Visualize cluster
fviz_cluster(hc.cut, ellipse.type = "convex")



# HIERARCHICAL CLUSTERING 
# heatmap(as.matrix(gower_dist), symm = F,
#         distfun = function(x) as.dist(x))

# Variable importance ----

pam_results <- dset.cluster %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary[[1]]
pam_results$the_summary[[2]]

dset.cluster1<-dset.cluster %>%
  filter(cluster == 1)

dset.cluster2<-dset.cluster %>%
  filter(cluster == 2)



#numericList<-getNumbersFromDataset(dset.cluster1)

#for(i in numericList){

factorList<-getFactorsFromDataset(dset.cluster)

for(i in factorList){
  dset.cluster[,i]<-as.numeric(as.character(dset.cluster[, i]))
}

factorList<-getFactorsFromDataset(dset.cluster1)

for(i in factorList){
  dset.cluster1[,i]<-as.numeric(as.character(dset.cluster1[, i]))
}

factorList<-getFactorsFromDataset(dset.cluster2)

for(i in factorList){
  dset.cluster2[,i]<-as.numeric(as.character(dset.cluster2[, i]))
}

vt1<-data.frame()

for(i in 1:dim(dset.cluster1)[2]){
  
  print(colnames(dset.cluster1[i]))
  print(var.test(dset.cluster1[, i], dset.cluster[, i]))
  f<-var.test(dset.cluster1[, i], dset.cluster[, i])
  v<-var(dset.cluster1[, i], na.rm = T)
  vt1<-rbind(
    vt1, data.frame(
      Name=colnames(dset.cluster1[i]),
      F=f$statistic,
      p.value=f$p.value
    )
  )
}

barplot(vt1$F, 
        names.arg = vt1$Name, 
        las=2,
        #space = 5,
        width = 10,
        cex.names = 0.6,
        main = "Variable contribution to cluster 1")
abline(h=1, col="red")

vt2<-data.frame()

for(i in 1:dim(dset.cluster2)[2]){
  
  print(colnames(dset.cluster2[i]))
  print(var.test(dset.cluster2[, i], dset.cluster[, i]))
  f<-var.test(dset.cluster2[, i], dset.cluster[, i])
  v<-var(dset.cluster2[, i], na.rm = T)
  vt2<-rbind(
    vt2, data.frame(
      Name=colnames(dset.cluster2[i]),
      F=f$statistic,
      p.value=f$p.value
    )
  )
}

barplot(vt2$F, 
        names.arg = vt2$Name, 
        las=2,
        #space = 5,
        width = 10,
        cex.names = 0.6,
        main = "Variable contribution to cluster 2")
abline(h=1, col="red")

write.csv2(vt1, file=paste0(outdir, "/vt1.csv"))

write.csv2(vt2, file=paste0(outdir, "/vt2.csv"))

  
  
