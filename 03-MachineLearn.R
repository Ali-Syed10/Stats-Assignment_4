#!/usr/bin/env Rscript

library(randomForest)
library(caret)
library(tidyverse)
library(e1071)
library(caTools)
library(caret)
library(kernlab)
library(plyr)
library(gglasso)


random_forest_test <- function(df, rng = 1212) {
  # Runs random forest test on the dataframe and returns the performance.

  set.seed(rng)

  # Split into test and train
  df.tst <- df %>%
    rownames_to_column('ID') %>%
    group_by(RACE) %>%
    slice_sample(prop = 0.25) %>%
    column_to_rownames('ID')

  df.trg <- df[-which(row.names(df) %in% row.names(df.tst)),] %>%
    rownames_to_column('ID') %>%
    group_by(RACE) %>%
    column_to_rownames('ID')

  # Change to factor
  df.tst$RACE <- as.factor(df.tst$RACE)
  df.trg$RACE <- as.factor(df.trg$RACE)
  df$RACE <- as.factor(df$RACE)

  ## Train

  # Control
  model.control <- trainControl(method = "cv", number = 10)

  model.bagged <- train(
    RACE ~ .,
    data = df.trg,
    method = "treebag",
    trControl = model.control,
    importance = TRUE
  )

  ## Test

  # Predict with testing df
  pred <- predict(model.bagged, df.tst)
  actual <- df.tst$RACE

  df.pred <- data.frame(pred = pred, actual = actual)

  ## Calculate stats

  # Confusion matrix
  cm <- table(actual, pred)

  # Accuracy
  accuracy <- sum(diag(cm)) / sum(cm)

  # Precision
  precision <- posPredValue(pred, actual, positive = "1")

  # Recall
  recall <- sensitivity(pred, actual, positive = "1")

  # F-score
  fscore <- (2 * precision * recall) / (precision + recall)

  # Summary DF
  df.summary <- data.frame(accuracy = accuracy, precision = precision, recall = recall, fscore = fscore, method = "RF")

  # Variable Importance
  var.imp <- varImp(model.bagged)$importance %>% rownames_to_column()
  names(var.imp) <- c("snp", "importance")
  
  return(list(pred = df.pred, vimp = var.imp, summary = df.summary))

}

svm_test <- function(df, rng = 1212) {
  # Same but for svm

  set.seed(rng)

  # Split into test and train
  df.tst <- df %>%
    rownames_to_column('ID') %>%
    group_by(RACE) %>%
    slice_sample(prop = 0.25) %>%
    column_to_rownames('ID')

  df.trg <- df[-which(row.names(df) %in% row.names(df.tst)),] %>%
    rownames_to_column('ID') %>%
    group_by(RACE) %>%
    column_to_rownames('ID')

  # Make factors
  # Change 1 and 0 to y/n
  df.tst$RACE <- as.factor(df.tst$RACE)
  df.tst$RACE <- as.factor(c('no', 'yes')[df.tst$RACE])
  df.trg$RACE <- as.factor(df.trg$RACE)
  df.trg$RACE <- as.factor(c('no', 'yes')[df.trg$RACE])

  # X and Y
  x.trg <- df.trg[, -1]
  y.trg <- df.trg[, 1]

  x.tst <- df.tst[, -1]
  y.tst <- df.tst[, 1]

  ## Train
  model.control <- trainControl(method = "cv",
                                number = 2,
                                summaryFunction = twoClassSummary,
                                classProbs = TRUE)

  # Grid search
  grid <- expand.grid(sigma = c(.01, .015, 0.2),
                      C = c(0.75, 0.9, 1, 1.1, 1.25)
  )

  # Scale data ???? TODO: Since the data is all 1 or 0 does it stil need to be scaled?
  model.tunned <- train(x = x.trg,
                        y = y.trg,
                        method = "svmRadial",
                        metric = "ROC",
                        tuneGrid = grid,
                        trControl = model.control)

  ## Test

  # Predict with testing df
  pred <- model.tunned %>% predict(x.tst)
  actual <- y.tst

  df.pred <- data.frame(pred = pred, actual = actual)

  ## Calculate stats

  # Confusion matrix
  cm <- table(actual, pred)

  # Accuracy
  accuracy <- sum(diag(cm)) / sum(cm)

  # Precision
  precision <- posPredValue(pred, actual, positive = "yes")

  # Recall
  recall <- sensitivity(pred, actual, positive = "yes")

  # F-score
  fscore <- (2 * precision * recall) / (precision + recall)

  # Summary DF
  df.summary <- data.frame(accuracy = accuracy, precision = precision, recall = recall, fscore = fscore, method = "SVM")

  # Variable Importance
  var.imp <- varImp(model.tunned)$importance %>% rownames_to_column()
  names(var.imp) <- c("snp", "importance", "importance_2")
  
  var.imp <- var.imp %>%
    select(snp, importance) %>%
    merge_snps()
  
  return(list(pred = df.pred, vimp = var.imp, summary = df.summary))
  
}

lasso_test <- function(df, rng = 1212) {
  # Same but with lasso

  set.seed(rng)

  # Split into test and train
  df.tst <- df %>%
    rownames_to_column('ID') %>%
    group_by(RACE) %>%
    slice_sample(prop = 0.25) %>%
    column_to_rownames('ID')

  df.trg <- df[-which(row.names(df) %in% row.names(df.tst)),] %>%
    rownames_to_column('ID') %>%
    group_by(RACE) %>%
    column_to_rownames('ID')

  # Making the ten folds for the training set. Again, I want even representation of the two populations, so I make the ten folds for each separately.
  lsFold0 <- rep(1:10, each = round_any(sum(df.trg$RACE == 0), 10, f = ceiling) / 10) %>%
    sample()
  lsFold0 <- lsFold0[-(1:(round_any(sum(df.trg$RACE == 0), 10, f = ceiling) - sum(df.trg$RACE == 0)))]

  lsFold1 <- rep(1:10, each = round_any(sum(df.trg$RACE == 1), 10, f = ceiling) / 10) %>%
    sample()
  lsFold1 <- lsFold1[-(1:(round_any(sum(df.trg$RACE == 1), 10, f = ceiling) - sum(df.trg$RACE == 1)))]

  # Adding the fold IDs to the training data frame.
  df.trg$FoldID <- c(lsFold0, lsFold1)
  df.trg <- df.trg %>%
    relocate(RACE, FoldID)

  # Making the vector that explains which predictors are part of the same group. So in this case, all SNPs with "SNP_1" will be part of group 1.
  lsSNPgrp <- colnames(df)[2:ncol(df)] %>%
    str_remove_all(pattern = "SNP_") %>%
    str_remove_all(pattern = "_[0-2].[0-2]") %>%
    as.integer()

  # Training the model.
  model <- cv.gglasso(x = as.matrix(df.trg[, 3:ncol(df.trg)]), y = df.trg$RACE, group = lsSNPgrp, foldid = df.trg$FoldID, pred.loss = "misclass")

  # Plotting the model. Tbh I don't know what's going on at the top. It should be how many predictors there are, so I don't know how there are negative numbers.

  # Letting you know how many predictors are actually used in the model.
  sum(coef(model, s = "lambda.min") != 0)

  # Sorting the coefficients by importance (they're sorted by absolute value, but their sign isn't actually lost).
  trg.coef <- coef(model, s = "lambda.min")[order(-abs(coef(model, s = "lambda.min"))),, drop = FALSE]
  trg.coef <- trg.coef[1:sum(coef(model, s = "lambda.min") != 0),, drop = FALSE]

  lsVIPSNPs <- row.names(trg.coef)[-1] %>%
    str_remove_all(pattern = "_[0-2].[0-2]$") %>%
    unique()

  # Making the prediction on the test set. Since this is still technically regression, the predictions are just values from ~0-~1, so I'm rounding them to the nearest integer.
  prediction <- predict(model, newx = df.tst[, 2:ncol(df.tst)], s = "lambda.min") %>%
    round()

  # Making the full model to visualize coefficients dropping out. I don't know why they plotted the same thing, in his script they're mirror images of each other (glmnet.R).
  model.full <- gglasso(x = as.matrix(df[, 2:ncol(df)]), y = df$RACE, group = lsSNPgrp, lambda = model$lambda)

  # Switching the values to factors, and then inputting them to the confusion matrix. It did pretty well!
  actual <- as.factor(df.tst$RACE)
  pred <- as.factor(prediction)

  df.pred <- data.frame(pred = pred, actual = actual)

  ## Calculate stats

  # Confusion matrix
  cm <- table(actual, pred)

  # Accuracy
  accuracy <- sum(diag(cm)) / sum(cm)

  # Precision
  precision <- posPredValue(pred, actual, positive = "1")

  # Recall
  recall <- sensitivity(pred, actual, positive = "1")

  # F-score
  fscore <- (2 * precision * recall) / (precision + recall)

  # Summary DF
  df.summary <- data.frame(accuracy = accuracy, precision = precision, recall = recall, fscore = fscore, method = "LASSO")

  # Variable Importance
  var.imp <- bind_cols(
    row.names(trg.coef),
    trg.coef %>% as_tibble()
  )[-1,]
  names(var.imp) <- c("snp", "importance")
  
  var.imp <- var.imp %>%
    merge_snps()
  
  return(list(pred = df.pred, vimp = var.imp, summary = df.summary))

}