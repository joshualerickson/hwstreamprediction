
# MSCV Methods ------------------------------------------------------------

#libraries
library(caret)
library(CAST)
library(tidyverse)
library(doParallel)
library(parallel)
library(gbm)
source('utils.R')

### We need to get global statistic metrics from 10-100 clusters
## we determined features selected by the most common occurrences in ffs_script_clustering.R results

#read in data

data34 <- read_csv("data/data34.csv") %>% mutate(across(where(is.character), factor))

#set.seed for splitting data
set.seed(1234)

train <- data34

#create a list of clusters by sequential centers
#use 'train' data

clusters <- list()

for (i in seq(10, 100, 10)){
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
for (i in seq(10, 100, 10)) {
  
  set.seed(1234)
  indices <- CAST::CreateSpacetimeFolds(flat, spacevar = paste("cluster.", i, sep = ""), k = 10)
  name <- paste("cluster.", i, sep = "")
  indices <- list(indices)
  names(indices) <- name
  index_list <- append(index_list, indices)
  
}

# GBM Topoclimatic --------------------------------------------------------


#now add the clusters groups back to the 'train' data

train <- train %>% cbind(flat)

gbm_topoclim_results <- data.frame() 
gbm_topoclim_varimp <- data.frame()

#now run for topolim model
for (i in seq_along(index_list)){
  
  #different names for the recipe
  topo_clim_names <- c("accum30", "tpi30agg", "deficitRS", "nppmmid30agg")
  
  #bring in cores
  cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
  registerDoParallel(cluster)
  
  gbm_rec <-  recipe(stream ~ ., data = train) %>% 
    update_role(-stream,-all_of(topo_clim_names), new_role = "bring along")
  
  gbmGrid <- expand.grid(interaction.depth = c(2,4,6,8), n.trees = c(500,1000,1500, 2000), shrinkage = c(0.001,.01, 0.05), n.minobsinnode = c(10,15,25,50))
  
  cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
  registerDoParallel(cluster)
  
  library(gbm)
  set.seed(1234)
  gbm_topoclim <- train(gbm_rec, data = train,
                     method = 'gbm',
                     distribution = 'bernoulli',
                     metric = "Accuracy",
                     tuneGrid = gbmGrid,
                     trControl = trainControl(method = "repeatedcv",
                                              repeats = 5,
                                              number = 10,
                                              classProbs = TRUE,
                                              returnResamp = "all",
                                              savePredictions = 'all',
                                              index = index_list[[i]][[1]],
                                              allowParallel = TRUE,
                                              verbose = FALSE, 
                                              summaryFunction = fourStats))

  bt <- gbm_topoclim$bestTune
  
  metrics <- gbm_topoclim$results %>% filter(n.trees %in% bt$n.trees, 
                               interaction.depth %in% bt$interaction.depth,
                               shrinkage %in% bt$shrinkage,
                               n.minobsinnode %in% bt$n.minobsinnode)
  
  metrics <- metrics %>% mutate(cluster = paste(names(index_list[i])),
                                model = "gbm_topoclim")
  
  vimp <- summary.gbm(gbm_topoclim$finalModel, method = permutation.test.gbm)
  rownames(vimp) <- NULL
  vimp <- vimp %>% mutate(cluster = paste(names(index_list[i])),
                          model = "gbm_topoclim")
  
  gbm_topoclim_results <- plyr::rbind.fill(gbm_topoclim_results, metrics)
  gbm_topoclim_varimp <- plyr::rbind.fill(gbm_topoclim_varimp, vimp)
  #now remove cluster. consistency...
  stopCluster(cluster)
  registerDoSEQ()
  
  
}


write_csv(gbm_topoclim_results, 'data/gbm_topoclim_results.csv')
write_csv(gbm_topoclim_varimp, 'data/gbm_topoclim_varimp.csv')


# GBM Topography only -----------------------------------------------------


# now do topo-only model but also add nhdplus to the different hold outs.

#add common identifier to training set
train$rowIndex <- seq(1,8802,1)

#now give the 0's"X0" and the values "X1"
train$NHDPlus_char <- ifelse(train$NHDPlus == 0, paste("X1"), paste("X2"))

train$NHDPlus_fact <- factor(train$NHDPlus_char)


gbm_topo_results <- data.frame() 
gbm_topo_varimp <- data.frame()
gbm_nhd_results <- data.frame()

#now run for topo-only model
for (i in seq_along(index_list)){
  
  #different names for the recipe
  
  topo_names <- c("twi30agg", "tpi30agg", "accum30")
  
  #bring in cores
  cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
  registerDoParallel(cluster)
  
  gbm_rec <-  recipe(stream ~ ., data = train) %>% 
    update_role(-stream,-all_of(topo_names), new_role = "bring along")
  
  gbmGrid <- expand.grid(interaction.depth = c(2,4,6,8), n.trees = c(500,1000,1500, 2000), shrinkage = c(0.001,.01, 0.05), n.minobsinnode = c(10,15,25,50))
  
  cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
  registerDoParallel(cluster)
  
  library(gbm)
  set.seed(1234)
  gbm_topo <- train(gbm_rec, data = train,
                        method = 'gbm',
                        distribution = 'bernoulli',
                        metric = "Accuracy",
                        tuneGrid = gbmGrid,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 5,
                                                 number = 10,
                                                 classProbs = TRUE,
                                                 returnResamp = "all",
                                                 savePredictions = 'all',
                                                 index = index_list[[i]][[1]],
                                                 allowParallel = TRUE,
                                                 verbose = FALSE, 
                                                 summaryFunction = fourStats))
  
  bt <- gbm_topo$bestTune
  
  metrics <- gbm_topo$results %>% filter(n.trees %in% bt$n.trees, 
                                             interaction.depth %in% bt$interaction.depth,
                                             shrinkage %in% bt$shrinkage,
                                             n.minobsinnode %in% bt$n.minobsinnode)
  
  metrics <- metrics %>% mutate(cluster = paste(names(index_list[i])),
                                model = "gbm_topo")
  
  vimp <- summary.gbm(gbm_topo$finalModel, method = permutation.test.gbm)
  rownames(vimp) <- NULL
  vimp <- vimp %>% mutate(cluster = paste(names(index_list[i])),
                          model = "gbm_topo")
  
  gbm_topo_results <- plyr::rbind.fill(gbm_topo_results, metrics)
  gbm_topo_varimp <- plyr::rbind.fill(gbm_topo_varimp, vimp)
  
  
  nhd_metrics <- gbm_topo$pred %>% filter(n.trees %in% bt$n.trees, 
                                              interaction.depth %in% bt$interaction.depth,
                                              shrinkage %in% bt$shrinkage,
                                              n.minobsinnode %in% bt$n.minobsinnode)
  
  
  #now join for train and resamps
  join_resamp <- left_join(nhd_metrics,train, by = "rowIndex" )
  
  nhdtune <- join_resamp %>% select(obs, NHDPlus_fact, Resample, rowIndex)
  
  blah <- nhdtune %>% split(.$Resample) %>% map(fourStatsNhdtune)
  
  blah_names <- data.frame(Resample = names(blah))
  
  blah <- blah %>% map_df(unlist)
  
  nhd <- cbind(blah_names, blah)
  
  nhd <- nhd %>% mutate(model = 'NHDPlus')
  
  
  nhd <- nhd %>% summarise(across(everything(), ~mean(., na.rm = T)))
  
  gbm_nhd_results <- plyr::rbind.fill(gbm_nhd_results, nhd)
  
  
  #now remove cluster. consistency...
  stopCluster(cluster)
  registerDoSEQ()
  
  
}


write_csv(gbm_nhd_results, 'data/gbm_nhd_results.csv')
write_csv(gbm_topo_results, 'data/gbm_topo_results.csv')
write_csv(gbm_topo_varimp, 'data/gbm_topo_varimp.csv')

##### need to save a final model for pdp's and plotting.

# Single Models for Plotting ----------------------------------------------

#topography only

topo_names <- c("twi30agg", "tpi30agg", "accum30")

#bring in cores
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

gbm_rec <-  recipe(stream ~ ., data = train) %>% 
  update_role(-stream,-all_of(topo_names), new_role = "bring along")

gbmGrid <- expand.grid(interaction.depth = c(2,4,6,8), n.trees = c(500,1000,1500, 2000), shrinkage = c(0.001,.01, 0.05), n.minobsinnode = c(10,15,25,50))

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

library(gbm)
set.seed(1234)
gbm_topo <- train(gbm_rec, data = train,
                  method = 'gbm',
                  distribution = 'bernoulli',
                  metric = "Accuracy",
                  tuneGrid = gbmGrid,
                  trControl = trainControl(method = "repeatedcv",
                                           repeats = 5,
                                           number = 10,
                                           classProbs = TRUE,
                                           returnResamp = "all",
                                           savePredictions = 'all',
                                           index = index_list[[i]][[1]],
                                           allowParallel = TRUE,
                                           verbose = FALSE, 
                                           summaryFunction = fourStats))

