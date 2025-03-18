# TREES TEMPLATE v1

source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

library(caret)

# PATH ----

wd<-setwd("~/Analisi/Manfrini/PHENO")
indir=paste0(wd, "/out/treeexp")
outdir=paste0(wd,"/out/treeexp")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana.train <- read.csv2(paste0(indir,"/dset.ana.train.csv"),
                        stringsAsFactors=FALSE) #, na.strings = ""
}

if(RAW==F){
  dset.ana.test <- read.csv2(paste0(indir,"/dset.ana.test.csv"),
                              stringsAsFactors=FALSE) #, na.strings = ""
}

# if(RAW==F){
#   dset.ana <- read.csv2(paste0(indir,"/dset.ana.cluster.csv"),
#                              stringsAsFactors=FALSE) #, na.strings = ""
# }

# Load dataset descriptor file
nm<-read.csv2(paste0(indir,"/nm.dset.ana.filtered.cluster.csv"), 
              stringsAsFactors=T)


# SET VAR TYPE ----

# setVarType(paste0(indir, "/nm.dset.ana.filtered.tree.1.csv"))
# 
# if(length(numericList)>0){
#   dset.ana<- applyNumeric(dset.ana, numericList)  
# } else {
#   print("No numeric vars")
# }
# if(length(integerList)>0){
#   dset.ana <- applyInteger(dset.ana, integerList)
# } else {
#   print("No integer vars")
# }
# if(length(factorList)>0){
#   dset.ana <- applyFactors(dset.ana, factorList)
# } else {
#   print("No factor vars")
# }
# if(length(orderedList)>0){
#   dset.ana <- applyOrdered(dset.ana, orderedList)
# } else {
#   print("No ordered vars")
# }
# if(length(logicalList)>0){
#   dset.ana <- applyLogical(dset.ana, logicalList)
# } else {
#   print("No logical vars")
# }
# if(length(dates)>0){
#   dset.ana <- applyDate(dset.ana, dates, "%d/%m/%Y")
# } else {
#   print("No date vars")
# }

# SET VAR TYPE TRAIN ----

setVarType(paste0(indir, "/nm.dset.ana.filtered.cluster.csv"))

if(length(numericList)>0){
  dset.ana.train <- applyNumeric(dset.ana.train, numericList)  
} else {
  print("No numeric vars")
}
if(length(integerList)>0){
  dset.ana.train <- applyInteger(dset.ana.train, integerList)
} else {
  print("No integer vars")
}
if(length(factorList)>0){
  dset.ana.train <- applyFactors(dset.ana.train, factorList)
} else {
  print("No factor vars")
}
if(length(orderedList)>0){
  dset.ana.train <- applyOrdered(dset.ana.train, orderedList)
} else {
  print("No ordered vars")
}
if(length(logicalList)>0){
  dset.ana.train <- applyLogical(dset.ana.train, logicalList)
} else {
  print("No logical vars")
}
if(length(dates)>0){
  dset.ana.train <- applyDate(dset.ana.train, dates, "%d/%m/%Y")
} else {
  print("No date vars")
}

# SET VAR TYPE TEST ----

setVarType(paste0(indir, "/nm.dset.ana.filtered.cluster.csv"))

if(length(numericList)>0){
  dset.ana.test <- applyNumeric(dset.ana.test, numericList)  
} else {
  print("No numeric vars")
}
if(length(integerList)>0){
  dset.ana.test <- applyInteger(dset.ana.test, integerList)
} else {
  print("No integer vars")
}
if(length(factorList)>0){
  dset.ana.test <- applyFactors(dset.ana.test, factorList)
} else {
  print("No factor vars")
}
if(length(orderedList)>0){
  dset.ana.test <- applyOrdered(dset.ana.test, orderedList)
} else {
  print("No ordered vars")
}
if(length(logicalList)>0){
  dset.ana.test <- applyLogical(dset.ana.test, logicalList)
} else {
  print("No logical vars")
}
if(length(dates)>0){
  dset.ana.test <- applyDate(dset.ana.test, dates, "%d/%m/%Y")
} else {
  print("No date vars")
}

# VISUALIZE ANALYSIS DATASET
vis_dat_wrapper(dset.ana.train)
vis_dat_wrapper(dset.ana.test)
# BART TRAINIG SET ----

# selected<-nm[which(nm$rg==1), "id"]
# 
# library(BayesTree)
# 
# trainset<-dset.ana.train[,selected]
# testset<-dset.ana.test[,selected]
# outcome.train<-ifelse(dset.ana.train$Outcome==5,1,0)
# 
# bart.fit<-bart(x.train=trainset, y.train=outcome.train, x.test=testset, ntree=200)
# 
# plot(bart.fit)
# 
# bart.fit$yhat.train
# 
# library("BART")
# 
# bart.fit<-pbart(x.train=trainset, y.train=outcome.train, x.test=testset)


# RPART ----

set.seed(998)

library(rpart)
library(rattle)
library(rpart.plot)

selected<-nm[which(nm$tree==1), "id"]

trainset<-dset.ana.train[,selected]

testset<-dset.ana.test[,selected]

# testset$n.lobes<-as.numeric(as.character(testset$n.lobes))
# 
# trainset$n.lobes<-as.numeric(as.character(trainset$n.lobes))

# colnames(trainset)<-c("Age", 
#                       "Cardiopathy", 
#                       "Coronaropathy",
#                       "PTCA",
#                       "CKD",
#                       "Creatinine",
#                       "FE",
#                       "CSV",
#                       "TC",
#                       "LAD",
#                       "LCX",
#                       "CDX",
#                       "VA",
#                       "Pneumonia",
#                       "Lobes",
#                       "Pattern",
#                       "cluster"
# )
# 
# colnames(testset)<-c("Age", 
#                       "Cardiopathy", 
#                       "Coronaropathy",
#                       "PTCA",
#                       "CKD",
#                       "Creatinine",
#                       "FE",
#                       "CSV",
#                       "TC",
#                       "LAD",
#                       "LCX",
#                       "CDX",
#                       "VA",
#                       "Pneumonia",
#                       "Lobes",
#                       "Pattern",
#                       "cluster"
# )


control <- trainControl(method="repeatedcv",
                        number=20,
                        repeats=20)

metric <- "Accuracy"

rpart.fit <- train(cluster ~ ., data = trainset, 
                 method = "rpart", 
                 trControl = control,
                 # Now specify the exact models 
                 # to evaluate:
                 #tuneLength = 10
                 )

# Plot model accuracy vs different values of
# cp (complexity parameter)

plot(rpart.fit) 
plot(rpart.fit, metric = "Kappa")

# Print the best tuning parameter cp that
# maximizes the model accuracy

rpart.fit$bestTune

# Plot the final tree model

jpeg(filename = paste0(outdir, "/rpart.jpg"), height=4, width=8, unit="in", res=300)
  rpart.plot(rpart.fit$finalModel, type = 5, extra = 0, clip.right.labs = T, cex = 0.8, fallen.leaves=F, clip.facs=F)
dev.off()

# Decision rules in the model

rpart.fit$finalModel

# PERF

testset$pred<-predict(rpart.fit, newdata = testset)

confusionMatrix(data = testset$pred, reference = testset$cluster, mode = "everything") # prec_recall

# C50 ----

set.seed(998)

library(C50)

selected<-nm[which(nm$tree==1), "id"]

trainset<-dset.ana.train[,selected]
testset<-dset.ana.test[,selected]

# testset$n.lobes<-as.numeric(as.character(testset$n.lobes))
# 
# trainset$n.lobes<-as.numeric(as.character(trainset$n.lobes))
# 
# colnames(trainset)<-c("Age", 
#                       "Cardiopathy", 
#                       "Coronaropathy",
#                       "PTCA",
#                       "CKD",
#                       "Creatinine",
#                       "FE",
#                       "TC",
#                       "LAD",
#                       "LCX",
#                       "CDX",
#                       "VA",
#                       "Pneumonia",
#                       "Lobes",
#                       "Pattern",
#                       "cluster"
# )
# 
# colnames(testset)<-c("Age", 
#                      "Cardiopathy", 
#                      "Coronaropathy",
#                      "PTCA",
#                      "CKD",
#                      "Creatinine",
#                      "FE",
#                      "TC",
#                      "LAD",
#                      "LCX",
#                      "CDX",
#                      "VA",
#                      "Pneumonia",
#                      "Lobes",
#                      "Pattern",
#                      "cluster"
# )

control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats=10,
                        returnResamp="all")

metric <- "Accuracy"

tunegrid <- expand.grid( .winnow = c(TRUE,FALSE), .trials=c(1,5,10,15,20), .model="tree" )

c50.fit <- train(cluster ~., 
                 data=trainset, 
                 method="C5.0", 
                 metric=metric, 
                 tuneGrid=tunegrid, 
                 trControl=control)

print(c50.fit)
c50.fit$finalModel
c50.fit$bestTune
summary(c50.fit)

testset$pred<-predict(c50.fit, newdata = testset)

confusionMatrix(data = testset$pred, reference = testset$cluster, mode = "prec_recall") #, mode = "prec_recall"

# Plot the final tree model

summary(c50.fit$finalModel)

jpeg(filename = paste0(outdir, "/c50.last.jpg"), height=8, width=10, unit="in", res=300)
  plot(C5.0(trainset[,c(1:17)], trainset[, 18], trials = 10, winnow = T), cex = 0.6)
dev.off()
# RANDOM FOREST ----

set.seed(998)

library(randomForest)

selected<-nm[which(nm$tree==1), "id"]

 trainset<-dset.ana.train[,selected]
 testset<-dset.ana.test[,selected]

control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats=3)

metric <- "Accuracy"

tunegrid <- expand.grid(.mtry=c(1:5))

rf.fit <- train(cluster ~., 
                data=trainset, 
                method="rf", 
                metric=metric, 
                tuneGrid=tunegrid, 
                trControl=control)

print(rf.fit)
rf.fit$finalModel
plot(rf.fit)

rf.fit$bestTune

testset$pred<-predict(rf.fit, newdata = testset)

confusionMatrix(data = testset$pred, reference = testset$cluster, mode = "prec_recall") #, mode = "prec_recall"

varImpPlot(rf.fit$finalModel, sort=TRUE,
           type=NULL, class=NULL, scale=TRUE)

# PLOT TREE FROM RANDOMFOREST PACKAGE ----

library(dplyr)
library(ggraph)
library(igraph)

#getTree(rf.fit$finalModel, 1, labelVar = T)

tree_func <- function(final_model, 
                      tree_num) {
  
  # get tree by index
  tree <- randomForest::getTree(final_model, 
                                k = tree_num, 
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
  
  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))
  
  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  
  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$split <- as.character(round(tree$`split point`, digits = 2))
  
  # plot
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 8))
  
  print(plot)
}

#tree_num <- which(rf.fit$finalModel$forest$ndbigtree == min(rf.fit$finalModel$forest$ndbigtree))

jpeg(filename = paste0(outdir, "/rf.jpg"), height=20, width=20, unit="in", res=600)

tree_func(final_model = rf.fit$finalModel, 1)

dev.off()





