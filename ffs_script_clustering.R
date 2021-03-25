
#libraries
library(caret)
library(CAST)
library(randomForest)
library(tidyverse)
library(doParallel)
library(parallel)

#read in data

data34 <- read_csv("data/data34.csv") %>% mutate(across(where(is.character), factor))

#set.seed for splitting data
set.seed(1234)

train <- data34
 
#create a list of clusters by sequential centers
#use 'train' data

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

ffs_results <- data.frame() 

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
              method = "rf",
              trControl = ctrl)

vars <- ffshw$selectedvars

metrics <- ffshw$results

metrics <- metrics %>% mutate(vars = list(vars), cluster = paste(names(index_list[i])),
                              model = "rf")


ffs_results <- plyr::rbind.fill(ffs_results, metrics)
#now remove cluster. consistency...
stopCluster(cluster)
registerDoSEQ()


}


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
  
  
  ffs_results <- plyr::rbind.fill(ffs_results, metrics)
  #now remove cluster. consistency...
  stopCluster(cluster)
  registerDoSEQ()
  
  
}

library(tidyverse)
ffs_results <- read_csv("D:/sending/ffs_results_full_52_clusters.csv")
#this removes `vars` as list column and makes character
ffs_results <- ffs_results %>%
  mutate(char_vars = str_remove_all(vars, c("^c" = "", "\\(|\\)" = "", "[\"]" = ""))) %>% 
  select(-vars)

wied_ex %>% count(value, sort = T)
cor_wise <- wied_ex %>% pairwise_cor(value, group, use = "pairwise.complete.obs", sort = T)

library(ggraph)
library(igraph)
library(tidygraph)
cor_wise %>% count(value, group, sort = TRUE)
cor_wise %>% head(62) %>% 
  as_tbl_graph() %>% 
  ggraph() + 
  geom_node_point() +
  geom_edge_link(aes(edge_colour = correlation), edge_width = 1.5) +
  geom_node_label(aes(label = name), repel = T) +
  theme_void()




write_csv(ffs_results, "ffs_results.csv")

# count the number of times two variables appear together

na_count <- na_count[!duplicated(t(apply(na_count, 1, sort))),]
counts_cor <- counts_cor[!duplicated(t(apply(counts_cor, 1, sort))),]
counts_cor <- counts_cor %>% mutate(comb = str_c(item1, item2, sep = " - "),
                                    correlation = round(correlation,2))

counts_cor %>% head(20) %>% mutate(comb = fct_reorder(comb, correlation)) %>% 
  ggplot(aes(comb, correlation)) + geom_col(alpha = .25) + geom_point() + geom_text(aes(label = correlation), nudge_y = .03)+
  labs(x = "variables", y = "correlation", title = "Top 20 Pairwise Variable Correlations") +
  coord_flip() + theme_light()

counts_vars %>% filter(is.na(item1)) %>% mutate(item2 = fct_reorder(item2, n)) %>% 
  ggplot(aes(item2, n)) + geom_col(alpha = .25) + geom_point() + geom_text(aes(label = n), nudge_y = 3)+
  labs(x = "variables", y = "count", title = "Variables used in each cluster model") +
  coord_flip()+ theme_light()


ffs_results %>% ggplot()

library(dplyr)
library(gapminder)
library(widyr)

gapminder %>%
  pairwise_cor(country, year, lifeExp)

gapminder %>%
  pairwise_cor(country, year, lifeExp, sort = TRUE)

