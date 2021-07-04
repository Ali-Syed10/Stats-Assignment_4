#!/usr/bin/env Rscript

# This file contains functions and stuff to use

library(tidyverse)
library(stringr)

one_hot_encode <- function(df, categories = c("0:1", "0:0", "1:1", "0:2", "1:2", "2:2")) {
  # This function takes a data frame (df), and returns a one hot encoded version.
  # Categories is a vector of all categories present in the data and will be appended to the column names.

  # Define categories
  cat.matches <- paste0(categories, collapse = "|")

  # Column names - Each is a SNP
  ncols <- colnames(df)[1:length(colnames(df))]

  # Appending to list is faster. Convert to dataframe at end.
  data <- list()

  # Iterate through each column/SNP
  for (i in 1:length(ncols)) {
    snp <- ncols[i]

    d <- df %>%
      select(snp) %>%
      mutate(n = 1) %>%
      mutate(ID = 1:nrow(df)) %>%
      distinct() %>%
      pivot_wider(id_cols = ID, names_from = snp, values_from = n) %>%
      mutate_at(vars(matches(cat.matches)), replace_na, 0) %>%
      select(-ID) # dont need ID anymore

    # Rename columns. Append category name.
    for (nt in categories) {
      colnames(d)[colnames(d) == nt] <- paste0(snp, "_", nt)
    }

    data[[i]] <- d
  }

  # Finally return a dataframe made from the list
  return(do.call(cbind, data))
}

arrange_gene <- function(sz) {
  # This function takes a string '1|0' and sorts it and returns the outputted string free of regex characters

  x <- unlist(str_split(str_replace_all(sz, "[^[:alnum:]]", "_"), "_"))

  paste(sort(x), collapse = ":")
}


merge_snps <- function(df) {
  # This function will merge SNPs in back from ohe by adding values across all

  # Remove ohe labels
  df$snp <- df %>%
    select(snp) %>%
    unlist() %>%
    str_remove_all(pattern = "_[0-2].[0-2]$")

  # Unique snps
  snps <- df %>%
    as_tibble() %>%
    select(snp) %>%
    unique() %>%
    unlist() %>%
    as.character()

  vals.list <- list()

  for (s in snps) {
    import <- df %>%
      filter(snp == s) %>%
      select(importance) %>%
      sum() %>%
      as.numeric()
    d <- data.frame(snp = s, importance = import)

    vals.list <- c(vals.list, list(d))
  }

  return(do.call(rbind, vals.list))
}


scree_plot <- function(pca, num.pc = 7) {
  #### Scree Plot Function ####
  # Takes a pca object from prcomp and plots the scree using ggplot

  # find PC explained variance
  df <- data.frame(comp = paste0("PC", 1:dim(pca$x)[2]), vare = 100 * (pca$sdev)^2 / sum((pca$sdev)^2))
  ggplot(data = df[1:num.pc, ], aes(x = comp, y = vare, group = 1)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    xlab("Component") +
    ylab("% of Variance Explained") +
    geom_text(aes(label = round(vare, digits = 2)), vjust = -0.25, color = "black", size = 3.5) +
    expand_limits(y = max(df$vare * 1.15)) # prevent numbers from being cut off
}


mean_sd <- function(x, n = 4) {
  # Returns the means ± sd as a string
  m = as.character(round(mean(x), n))
  d = as.character(round(sd(x), n))
  
  return(paste(m, "±", d))
}
