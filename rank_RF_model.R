require(switchBox)
require(randomForest)
library(dplyr)
library(magrittr)
library(Boruta)
library(tibble)
library(tidyr)
library(hash)
library(pROC)
library(PRROC)
require(svMisc)
library(matrixStats)
library(foreach)
library(doParallel)
library(caret)

#_____________________ RANK VER ________________________

class_def <- function(matrix){
  
  class <- c()
  for (x in 1:(dim(matrix)[1])){
    if (grepl("tum", rownames(matrix)[x]) == TRUE){
      class <- c(class, 1)
    }
    else if(grepl("adj", rownames(matrix)[x]) == TRUE){
      class <- c(class, 2)
    }
    else if(grepl("Normal", rownames(matrix)[x]) == TRUE){
      class <- c(class, 2)
    }
#    if (grepl("LumA", rownames(matrix)[x]) == TRUE){
#      class <- c(class, 1)
#    }
#    else if(grepl("LumB", rownames(matrix)[x]) == TRUE){
#      class <- c(class, 2)
#    }
#    else if(grepl("Her2", rownames(matrix)[x]) == TRUE){
#      class <- c(class, 2)
#    }
#    else if(grepl("Basal", rownames(matrix)[x]) == TRUE){
#      class <- c(class, 2)
#    }
    else{
      print(rownames(matrix)[x])
    }
  }
  class <- factor(class)
  
  return(class)
}

#_________________________ SIMPLE RANKING & INDICATOR_______________________________

ranking_matrix <- function(mat){
  rank_matrix <- matrix(nrow = dim(mat)[1], ncol = 0)
  
  rank_matrix = foreach(i=1:dim(mat)[2], .combine=cbind) %dopar% {
    one_col <- matrix(nrow = dim(mat)[1], ncol = 1)
    one_col[,1] <- rank(mat[,i])
    as.matrix(one_col)
  }
  colnames(rank_matrix) <- colnames(mat)
  rownames(rank_matrix) <- rownames(mat)
  
  return(rank_matrix)
}

#_________________________ BORUTA Feature Selection _______________________________

Boruta_Features <- function(train){
  #Boruta feature selection

  x <- train[, colnames(train) != 'cl']
  y <- factor(as.character(train$cl))
  
  boruta_rank <- Boruta(x=x, y=y, doTrace = 2, ntree = 345)
  
  #decision for tentative attributes
  boruta_rank_ten <- TentativeRoughFix(boruta_rank)
  features <- getSelectedAttributes(boruta_rank_ten)
  
  return (features)
}

#_________________________ RANDOM FOREST ________________________________

Quality <- function(pred, test){
  
  index_class0 <- test$cl == 1
  index_class1 <- test$cl == 2
  
  roc <- roc.curve(pred[index_class1, 2], pred[index_class0, 2], curve = TRUE)
  #png("ROC.png")
  #plot(roc)
  #dev.off()
  
  pr <- pr.curve(pred[index_class1, 2], pred[index_class0, 2], curve = TRUE)
  #png("PR.png")
  ### plot(pr)
  #dev.off()
  
  print(paste0("ROC-AUC : ", roc$auc))
  print(paste0("PR-AUC : ", pr$auc.integral))
  
  return (list(roc$auc, pr$auc.integral))
}


RF_function <- function(train_mat, test_mat){
  
  x <- train_mat[, colnames(train_mat) != 'cl']
  y <- factor(as.character(train_mat$cl))
  
  #building rf model
  rf <- foreach(ntree=rep(30, 15), .combine=randomForest::combine,
                .multicombine=TRUE, .packages='randomForest') %dopar% {
                  randomForest(x=x, y=y, importance=TRUE, ntree=ntree)}  
  #prediction
  pred <- predict(rf, test_mat, type="prob")
  pred_cl <- predict(rf, test_mat, type="class")
  
  #quality
  t_pred<-confusionMatrix(pred_cl, test_mat$cl)
  print(t_pred)
  
  #feat_imp <- randomForest::importance(rf)
  #table_imp_feat <- feat_imp[feat,]
  q <- Quality(pred, test_mat)
  
  feat<-Boruta_Features(train_mat)
  
  return (feat)
  #return (list(t_pred$byClass['Balanced Accuracy'], q[1], q[2]))
}

#______________________ DATA ________________________________________
set.seed(123)

train_ov <- read.table(file = "TRAIN_data.csv", sep = "\t", header = TRUE)
test_ov <- read.table(file = "TEST_data.csv", sep = "\t", header = TRUE)

rownames(train_ov) <- train_ov$GeneSymbol
rownames(test_ov) <- test_ov$Gene

train_ov <- train_ov[,-which(names(train_ov) %in% c("Gene"))]
test_ov <- test_ov[,-which(names(test_ov) %in% c("Gene"))]

cols <- intersect(rownames(train_ov), rownames(test_ov))
train_ov <- train_ov[cols,]
test_ov <- test_ov[cols,]

PUFA <- read.table("Pufa_cascade.tsv", sep="\t", header = TRUE)

#________________________ MAIN FUNCTION _______________

RF <- function(train, test, genes){
  
  tmp <- sample(1:dim(train)[2], size = dim(train)[2])
  train <- train[,tmp]
  
  train <- ranking_matrix(as.matrix(train)) 
  test <- ranking_matrix(as.matrix(test)) 
  
  inter <- intersect(rownames(train), genes)
  train <- train[inter,]
  test <- test[inter,]
  
  train <- as.data.frame(t(train)); 
  test <- as.data.frame(t(test))
  
  train$cl <- class_def(train)
  test$cl <- class_def(test)
  
  output <- RF_function(train, test)
  
  return(output)
}

main_select_genes <- function(train, test, genes){
  for(s in seq(1, 2, by=1)){
    file_name <- paste("Boruta/ImpGenes", s, ".csv", sep="_");
    features_tmp <- RF(train, test, genes);
    write.csv(features_tmp, file=file_name)
  } 
}