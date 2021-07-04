#!/usr/bin/env Rscript 
library(tidyverse)
library(stringr)
library(rlist)

source("src/00-Util.R")
source("src/03-MachineLearn.R")

# Load the built and encoded data frame
df.ohe <- read_csv("./data/data_onehot_encoded.csv") %>%
  select(-ID)
df.label <- read_csv("./data/data_label_encoded.csv") %>%
  select(-ID)


# The number of replicates to use. Each algorithm will be run n times and the final summary of performance will be n*number of methods used long
n <- 50

# Seed 
set.seed(1212)

# To store performance summaries of each run
list.perf <- list()

for (i in 1:n) {
  
  # Set random seed for each run
  rng <- sample(1:43674, 1)
  
  # Generate a random seed to use, this will be what the splitting and training will be based off
  rf.perf <- random_forest_test(df = df.label, rng = rng)
  svm.perf <- svm_test(df = df.ohe, rng = rng)
  lasso.perf <- lasso_test(df = df.ohe, rng = rng)

  # Add the method to each of the objects
  rf.perf <- c(rf.perf, list(method = "RF"))
  svm.perf <- c(svm.perf, list(method = "SVM"))
  lasso.perf <- c(lasso.perf, list(method = "LASSO"))
  
  # Append to list
  list.perf <- c(list.perf, list(rf.perf), list(svm.perf), list(lasso.perf))
  
}

# Concatenate into a single performance data frame
df.perf <- data.frame()
for (i in list.perf) {
  df.perf <- rbind(df.perf, i$summary)
}

# Save to list to yaml
if (file.exists("./data/perfourmance_summary.yaml")) {
  print("skipping")
} else {
  list.save(list.perf, "./data/perfourmance_summary.yaml")
}

# Write to CSV.
if (file.exists("./data/perfourmance_summary.csv")) {
  print("skipping")
} else{
  write.csv(df.perf, "./data/perfourmance_summary.csv", row.names = FALSE)
}
