#!/usr/bin/env Rscript

library(vcfR)
library(tidyverse)
library(stringr)

# These snps are essential and should not be removed in any filtering even if it matches threshold
KEY.SNPS <- c("61214790", "61217542")

# Loop through all the VCF files
data.dir <- "./data"
vcf.list <- list.files(path = data.dir, pattern = "*.vcf.gz", full.names = TRUE, recursive = FALSE)

# Appending to list is faster. Convert to df AFTER.
list.data <- list()

# Iterate through every file
for (file in vcf.list) {

  # Skip EUR so it doesnt happen twice
  if (str_detect(file, "EUR")) {
    next
  }

  # Filerting Threshold
  THRESHOLD <- 0.001

  # Gene and file names
  file.basename <- unlist(str_split(file, "/"))[length(unlist(str_split(file, "/")))]
  dir <- "./data/"
  gene <- unlist(str_split(file.basename, "_"))[1]
  ext <- ".vcf.gz"
  populations <- c("EUR", "EAS")

  # Filtering list
  rmlist <- list()
  j <- 1

  # Inner data 
  df.data <- data.frame()

  # Do for each population. Find filtering. TODO: make this faster but whatever
  for (pop in populations) {
    f <- paste0(dir, gene, "_", pop, ext)

    vcf <- read.vcfR(f)
    mf <- maf(vcf, element = 2)

    # Filter step
    for (i in 1:nrow(mf)) {
      if (mf[i, "Frequency"] < THRESHOLD) {
        # Only remove if NOT esential
        if (str_detect(row.names(mf)[i], paste(KEY.SNPS, collapse = "|"), negate = TRUE)) {
          rmlist[j] <- row.names(mf)[i]
          j <- j + 1
        }
      }
    }

  }

  # Now do the SNP and filter
  for (pop in populations) {

    f <- paste0(dir, gene, "_", pop, ext)

    vcf <- read.vcfR(f)
    mf <- maf(vcf, element = 2)
    gt <- t(extract.gt(vcf))

    # Remove list of filtered snps
    gt <- data.frame(gt[, - which(colnames(gt) %in% rmlist)])

    gt <- gt %>%
      rownames_to_column(var = "rowname") %>%
      relocate(rowname)

    gt <- apply(gt, c(1, 2), arrange_gene)

    # Add Race information
    gt <- as_tibble(gt) %>%
      mutate(race = if_else(pop == "EUR", 0, 1)) %>%
      mutate(gene = gene) %>%
      relocate(gene, race)

    df.data <- rbind(df.data, gt)
  }

  list.data <- c(list.data, list(df.data))

  # Remove large files from memory
  rm(gt, gene, file.basename, rmlist, j, mf, vcf, pop, i, f)
}

# Merge Into 1 data frame
df1 <- list.data[[1]] %>% select(-gene)
df2 <- list.data[[2]] %>% select(-gene)
df3 <- list.data[[3]] %>% select(-gene)

x <- merge(df1, df2, by = c("rowname", "race"))
df.all <- merge(x, df3, by = c("rowname", "race"))

# Change from snp name to snp number, saving the snp name
snp.names <- colnames(df.all %>% select(-rowname, -race))

names(df.all) <- c("ID", "RACE", paste0("SNP_", 1:(length(df.all) - 2)))

# Additionally, Create a mapping of SNP number to gene and SNP name
df <- list.data[[1]]
  
df.mappings <- data.frame(
  SNP_NUM = paste0("SNP_", 1:(length(df.all) - 2)),
  SNP_NAME = snp.names
)

# Add in gene annotation
df.genes <- bind_rows(
  data.frame(GENE = unique(list.data[[1]]$gene),
                  SNP_NAME = colnames(list.data[[1]])),
  data.frame(GENE = unique(list.data[[2]]$gene),
             SNP_NAME = colnames(list.data[[2]])),
  data.frame(GENE = unique(list.data[[3]]$gene),
             SNP_NAME = colnames(list.data[[3]]))
)
# Merge with the SNP numbers.
df.mappings <- merge(df.mappings, df.genes, by = c("SNP_NAME")) %>%
  relocate(SNP_NUM)

# Save dataframes to CSV
write.csv(df.all, "./data/raw_data_all.csv", row.names = FALSE)
write.csv(df.mappings, "./data/snp_mappings.csv", row.names = FALSE)


rm(df1, df2, df3, x, data.dir, dir, ext, file, vcf.list, THRESHOLD)

