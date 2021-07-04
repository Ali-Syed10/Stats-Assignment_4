#!/usr/bin/env Rscript 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(rlist)

source("./src/00-Util.R")

# Read in the summary data
list.perf <- list.load("./data/perfourmance_summary.yaml")
df.perf <- read_csv("./data/perfourmance_summary.csv")

# Gene/snp mapping
df.mappings <- read_csv("./data/snp_mappings.csv")
names(df.mappings) <- c("snp", "name", "gene")

## Variable importance plots
var.imp.rf <- list.perf[[1]]$vimp %>% as_tibble()
var.imp.svm <- list.perf[[2]]$vimp %>% as_tibble()

# Add the rest of the SNPs in as zero for lasso
var.imp.lasso <- full_join(list.perf[[3]]$vimp %>% as_tibble(),
                           var.imp.svm %>%
                             select(snp),
                           by=c('snp')) %>%
  replace(is.na(.), 0)

# Add gene annotations/SNP annotations, replace SNP_1 with the actual name of the SNP
var.imp.rf <- full_join(var.imp.rf, df.mappings, by="snp")
var.imp.svm <- full_join(var.imp.svm, df.mappings, by="snp")
var.imp.lasso <- full_join(var.imp.lasso, df.mappings, by="snp")

# Random Forest
p.vimp.rf <- var.imp.rf %>%
  arrange(importance, snp) %>%
  mutate(snp = factor(snp, levels = snp)) %>%
  ggplot() +
  geom_segment(aes(x = snp, y = 0, xend = snp, yend = importance), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = snp, y = importance, col = snp), 
             size = 4, show.legend = F) +
  coord_flip() +
  xlab("SNP") +
  ylab("Relative Importance (%)") + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(breaks = c(0,25,50,75,100), labels=c(0,25,50,75,100)) +
  ggtitle("A. Random Forest")

# SVM
p.vimp.svm <- var.imp.svm %>%
  arrange(importance, snp) %>%
  mutate(snp = factor(snp, levels = snp)) %>%
  ggplot() +
  geom_segment(aes(x = snp, y = 0, xend = snp, yend = importance), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = snp, y = importance, col = snp), 
             size = 4, show.legend = F) +
  coord_flip()+
  xlab("SNP") +
  ylab("Relative Importance (%)") + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(breaks = c(0,50,100,150,200), labels=c(0,25,50,75,100)) +
  ggtitle("B. Support Vector Machine")

# LASSO
p.vimp.lasso <- var.imp.lasso %>%
  arrange(importance, snp) %>%
  mutate(snp = factor(snp, levels = snp)) %>%
  ggplot() +
  geom_segment(aes(x = snp, y = 0, xend = snp, yend = importance), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = snp, y = importance, col = snp), 
             size = 4, show.legend = F) +
  coord_flip() +
  xlab("SNP") +
  ylab("Relative Importance (%)") + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(breaks = c(0,
                                max(var.imp.lasso[2]*0.25),
                                max(var.imp.lasso[2]*0.5),
                                max(var.imp.lasso[2]*0.75),
                                max(var.imp.lasso[2]*1)),
                     labels=c(0,25,50,75,100)) +
  ggtitle("C. LASSO")

## Same as above but with only the top x genes, along with their labels

# Random Forest
p.vimp.rf.names <- var.imp.rf %>%
  mutate(gene.name = paste(gene,name)) %>%
  arrange(importance, gene.name) %>%
  mutate(snp = factor(gene.name, levels = gene.name)) %>%
  tail(25) %>%
  ggplot() +
  geom_segment(aes(x = snp, y = 0, xend = snp, yend = importance), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = snp, y = importance, col = snp), 
             size = 4, show.legend = F) +
  coord_flip() +
  xlab("SNP") +
  ylab("Relative Importance (%)") + 
  scale_y_continuous(breaks = c(0,25,50,75,100), labels=c(0,25,50,75,100)) +
  ggtitle("A. Random Forest")

# SVM
p.vimp.svm.names <- var.imp.svm %>%
  mutate(gene.name = paste(gene,name)) %>%
  arrange(importance, gene.name) %>%
  mutate(snp = factor(gene.name, levels = gene.name)) %>%
  tail(25) %>%
  ggplot() +
  geom_segment(aes(x = snp, y = 0, xend = snp, yend = importance), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = snp, y = importance, col = snp), 
             size = 4, show.legend = F) +
  coord_flip()+
  xlab("SNP") +
  ylab("Relative Importance (%)") + 
  scale_y_continuous(breaks = c(0,50,100,150,200), labels=c(0,25,50,75,100)) +
  ggtitle("B. Support Vector Machine")

# LASSO
p.vimp.lasso.names <- var.imp.lasso %>%
  mutate(gene.name = paste(gene,name)) %>%
  arrange(importance, gene.name) %>%
  mutate(snp = factor(gene.name, levels = gene.name)) %>%
  tail(25) %>%
  ggplot() +
  geom_segment(aes(x = snp, y = 0, xend = snp, yend = importance), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = snp, y = importance, col = snp), 
             size = 4, show.legend = F) +
  coord_flip() +
  xlab("SNP") +
  ylab("Relative Importance (%)") + 
  scale_y_continuous(breaks = c(0,
                                max(var.imp.lasso[2]*0.25),
                                max(var.imp.lasso[2]*0.5),
                                max(var.imp.lasso[2]*0.75),
                                max(var.imp.lasso[2]*1)),
                     labels=c(0,25,50,75,100)) +
  ggtitle("C. LASSO")


## Performance Box plots - With error bars of paired T tests between

# Accuracy
p.acc <- ggboxplot(df.perf, x = "method", y = "accuracy", fill = "method") +
  stat_pvalue_manual(df.perf %>%
                       t_test(accuracy ~ method, paired = TRUE) %>%
                       adjust_pvalue(method = "bonferroni") %>%
                       add_significance("p.adj") %>%
                       add_xy_position(x = "method"), 
                     label = "p.adj.signif", tip.length = 0.025,
                     y.position = c(1,1.005,1.01), 
                     bracket.shorten = 0.05) + 
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_discrete(name = "Method", labels = c("LASSO", "Random Forest", "Support Vector Machine"))

# Precision
p.pre <- ggboxplot(df.perf, x = "method", y = "precision", fill = "method") +
  stat_pvalue_manual(df.perf %>%
                       t_test(precision ~ method, paired = TRUE) %>%
                       adjust_pvalue(method = "bonferroni") %>%
                       add_significance("p.adj") %>%
                       add_xy_position(x = "method"), 
                     label = "p.adj.signif", tip.length = 0.025,
                     y.position = c(1,1.005,1.01), 
                     bracket.shorten = 0.05) + 
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_discrete(name = "Method", labels = c("LASSO", "Random Forest", "Support Vector Machine"))

# Recall
p.rec <- ggboxplot(df.perf, x = "method", y = "recall", fill = "method") +
  stat_pvalue_manual(df.perf %>%
                       t_test(recall ~ method, paired = TRUE) %>%
                       adjust_pvalue(method = "bonferroni") %>%
                       add_significance("p.adj") %>%
                       add_xy_position(x = "method"), 
                     label = "p.adj.signif", tip.length = 0.025,
                     y.position = c(1,1.005,1.01), 
                     bracket.shorten = 0.05) + 
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_discrete(name = "Method", labels = c("LASSO", "Random Forest", "Support Vector Machine"))

# F-score
p.fsc <- ggboxplot(df.perf, x = "method", y = "fscore", fill = "method") +
  stat_pvalue_manual(df.perf %>%
                       t_test(fscore ~ method, paired = TRUE) %>%
                       adjust_pvalue(method = "bonferroni") %>%
                       add_significance("p.adj") %>%
                       add_xy_position(x = "method"), 
                     label = "p.adj.signif", tip.length = 0.025,
                     y.position = c(1,1.005,1.01), 
                     bracket.shorten = 0.05) + 
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_discrete(name = "Method", labels = c("LASSO", "Random Forest", "Support Vector Machine"))


## Arrange
ggarrange(p.acc, p.pre, p.rec, p.fsc, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")

ggarrange(p.vimp.rf,
          p.vimp.svm + theme(axis.title.y = element_blank()),
          p.vimp.lasso + theme(axis.title.y = element_blank()),
          nrow = 1, ncol = 3)

ggarrange(p.vimp.rf.names + theme(axis.title.x = element_blank()),
          p.vimp.svm.names + theme(axis.title.x = element_blank()),
          p.vimp.lasso.names,
          nrow = 3, ncol = 1)

p.vimp.rf
p.vimp.svm
p.vimp.lasso
p.vimp.rf.names
p.vimp.svm.names
p.vimp.lasso.names


## Summary Tables

# Taken from above plots
sig.label <- list(accuracy = c("a","a","a"),
                  precision = c("a","b","b"),
                  recall = c("a","b","b"),
                  fscore = c("a","a","ab"))
# Accuracy
df.perf %>%
  group_by(method) %>%
  summarise(min = min(accuracy),
            max = max(accuracy),
            mean = mean_sd(accuracy)) %>%
  mutate(significance = sig.label$accuracy) %>%
  mutate(method = c("LASSO", "Random Forest", "SVM")) %>%
  knitr::kable(., col.names = c(
    "(A) Method Accuracy", 
    "Minimum",
    "Maximum",
    "Mean ± Standard Deviation",
    "Significance Group"))
# Precision
df.perf %>%
  group_by(method) %>%
  summarise(min = min(precision),
            max = max(precision),
            mean = mean_sd(precision)) %>%
  mutate(significance = sig.label$precision) %>%
  mutate(method = c("LASSO", "Random Forest", "SVM")) %>% 
  knitr::kable(., col.names = c(
    "(B) Method Precision", 
    "Minimum",
    "Maximum",
    "Mean ± Standard Deviation",
    "Significance Group"))
# Recall
df.perf %>%
  group_by(method) %>%
  summarise(min = min(recall),
            max = max(recall),
            mean = mean_sd(recall)) %>%
  mutate(significance = sig.label$recall) %>%
  mutate(method = c("LASSO", "Random Forest", "SVM")) %>% 
  knitr::kable(., col.names = c(
    "(C) Method Recall", 
    "Minimum",
    "Maximum",
    "Mean ± Standard Deviation",
    "Significance Group"))
# F-score
df.perf %>%
  group_by(method) %>%
  summarise(min = min(fscore),
            max = max(fscore),
            mean = mean_sd(fscore)) %>%
  mutate(significance = sig.label$fscore) %>%
  mutate(method = c("LASSO", "Random Forest", "SVM")) %>%
  knitr::kable(., col.names = c(
    "(D) Method F-Score", 
    "Minimum",
    "Maximum",
    "Mean ± Standard Deviation",
    "Significance Group"))

