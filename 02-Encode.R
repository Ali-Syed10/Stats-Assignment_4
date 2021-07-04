#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)

# Load data
df.all <- read_csv("./data/raw_data_all.csv")

## Encode - One hot
df.ohe <- df.all %>%
  select(-ID, - RACE) %>%
  one_hot_encode() %>%
  mutate(RACE = df.all$RACE) %>%
  mutate(ID = df.all$ID) %>%
  relocate(ID, RACE)

## Encode Label
df.label <- df.all %>%
  select(-ID, - RACE) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_if(is.factor, as.numeric) %>%
  mutate(RACE = df.all$RACE) %>%
  mutate(ID = df.all$ID) %>%
  relocate(ID, RACE)


## Save each df to file
write.csv(df.label, "./data/data_label_encoded.csv", row.names = FALSE)
write.csv(df.ohe, "./data/data_onehot_encoded.csv", row.names = FALSE)
