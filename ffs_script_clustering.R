
#libraries
library(caret)
library(CAST)
library(randomForest)
library(tidyverse)
library(doParallel)
library(parallel)
source('utils.R')

#read in data generated from pre-processing

data34 <- read_csv("data/data34.csv") %>% mutate(across(where(is.character), factor))

#set.seed for splitting data
set.seed(1234)

train <- data34
 
#create a list of clusters by sequential centers
#use 'train' data
# how many iterations is up to the user, below is an example with 5:9 but we 
# did way more!

clusters <- list()

for (i in 5:9){
  val <- i
  cluster_val <- nrow(train)/val
  
  Mycluster <- kmeans(train[,c(17,18)], cluster_val) 
  
  name <- paste("cluster",val)
  
  Mycluster <- Mycluster$cluster
  
  Mycluster <- data.frame(Mycluster) 
  
  colnames(Mycluster) <- name
  
  Myclusterl <- list(Mycluster)
  
  clusters <- append(clusters, Myclusterl)
}

#now take the list and create a data.frame

flat <- clusters %>% flatten() %>% as.data.frame()


#now create list of spatial indices by clustering amounts (e.g. data.frame)

index_list <- list()
for (i in 5:9) {
  
  set.seed(1234)
indices <- CAST::CreateSpacetimeFolds(flat, spacevar = paste("cluster.", i, sep = ""), k = 10)
name <- paste("cluster.", i, sep = "")
indices <- list(indices)
names(indices) <- name
index_list <- append(index_list, indices)

}

#now add the clusters groups back to the 'train' data

train <- train %>% cbind(flat)


#now iterate over the list length (e.g clusters) to get ffs results per cluster stratification

#Now with GBM

ffs_results_gbm <- data.frame() 

for (i in seq_along(index_list)){
  
  predictorsRF <- train[,c("twi30agg","tpi30agg","accum30","vv30agg", "vvsd30agg","npol30agg","NDVI_med", "nppmmid30agg","deficitRS","cad30RS","decid30RS","cpg30precip","Green_med","NIR_med")]
  
  responseRF <- train$stream
  
  
  #bring in cores
  cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
  registerDoParallel(cluster)
  
  
  #set up control
  ctrl <- trainControl(method="repeatedcv",
                       repeats = 5,number = 10,
                       allowParallel = TRUE,
                       returnResamp = "all", 
                       verbose = FALSE, 
                       classProbs = TRUE,
                       summaryFunction = fourStats,
                       index = index_list[[i]][[1]],
                       savePredictions = 'all')
  #run
  set.seed(1234)
  ffshw <- ffs(predictorsRF,responseRF, 
               metric = "Accuracy",
               method = "gbm",
               trControl = ctrl)
  
  vars <- ffshw$selectedvars
  
  metrics <- ffshw$results
  
  metrics <- metrics %>% mutate(vars = list(vars), cluster = paste(names(index_list[i])),
                                model = "gbm")
  
  
  ffs_results_gbm <- plyr::rbind.fill(ffs_results_gbm, metrics)
  #now remove cluster. consistency...
  stopCluster(cluster)
  registerDoSEQ()
  
  
}

library(dplyr)

#this removes `vars` as list column and makes character
ffs_results <- ffs_results %>%
  mutate(char_vars = str_remove_all(vars, c("^c" = "", "\\(|\\)" = "", "[\"]" = ""))) %>% 
  select(-vars)

write_csv(ffs_results, 'all_results_final.csv')

