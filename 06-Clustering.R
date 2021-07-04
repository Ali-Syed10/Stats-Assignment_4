#### Loading Packages ####
library(ggplot2)
library(tidyverse)
library(ggfortify)
library(cluster)
library(clValid)
library(factoextra)
library(bpca)
library(pca3d)
library(GGally)
library(ggpubr)

source("./src/00-Util.R")

#### Loading the data and Preparing a custom data frame ####

# Load the built and encoded data frame
df.ohe <- read_csv("./data/data_onehot_encoded.csv") %>%
  select(-ID) %>%
  as.data.frame()

#### Loading the data and Preparing a custom data frame ####
df.ohe$RACE = as.factor(df.ohe$RACE) # Transforming the race as factors otherwise in plots it is read as a scale rather than individuals. 

df.ohe$RACE = fct_recode(df.ohe$RACE, Asian = "1", European = "0") # The factors 0 and 1 are changed to their designated names as Asian and European 


#### Creation of Data-Frame For Gene by Gene Analysis

MCIR = df.ohe[,448:491] # The dataset has a combination of all the genes SNP's thus they have been extracted by genes to asses population diviisions and if genes are able to produce disrcrimination between populations and also to what degree is their genetic diversity. The gene here is MCIR 
MCIR$Race = df.ohe[,1] # The data.frame extracted for a certain gene does not contain the column that defines each row by RACE thus that has been included as well. 

EOGT = df.ohe[,2:447] # Extraction of EOGT gene SNP data 
EOGT$Race = df.ohe[,1]

PUS10 = df.ohe[,492:817] # Extraction of PUS10 gene SNP data 
PUS10$Race = df.ohe[,1]

#### Prinicpal Component Computation and PCA Plot generation for MCIR gene ####


dat <- as.data.frame(sapply(MCIR[,-45], as.numeric)) # sapply is here to convert the data enteries to numeric from character as saved in object dat
MICR.pr = prcomp(dat) # computation of the principal components 
autoplot(MICR.pr, data = MCIR, colour = "Race") # generation of autoplot for prinicipal component analysis 

#### Prinicpal Component Computation and PCA Plot generation for EOGT gene ####

dat.1 <- as.data.frame(sapply(EOGT[,-447], as.numeric)) 
EOGT.pr = prcomp(dat.1)
autoplot(EOGT.pr, data = EOGT, colour = "Race")

#### Prinicpal Component Computation and PCA Plot generation for PUS10 gene ####

dat.2 <- as.data.frame(sapply(PUS10[,-327], as.numeric)) 
PUS.pr <- prcomp(dat.2)
autoplot(PUS.pr , data = PUS10, colour = "Race")

#### Scree plot generation for PCA ####

scree_plot(MICR.pr) # application of the scree plot function 
scree_plot(EOGT.pr)
scree_plot(PUS.pr)

#### Computing the Distance Matrix for Clusting Analysis ####

MICR.clust = dist(dat) # dist() function is used to compute the distance matrix needed for clustering 
EOGT.clust = dist(dat.1)
PUS.clust = dist(dat.2)

#### Visualizing the Distance Matrix ####

V.M <- fviz_dist(MICR.clust, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) #Fviz_dist() function is required for visualization of the distance matrix. 
V.E <- fviz_dist(EOGT.clust, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
V.P <- fviz_dist(PUS.clust, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#### Converting the dist class objects to numeric matrixes for cluster validation to identify appropriate number of clusters ####

valid.M = as.matrix(MICR.clust) # The conversion of dist class objects to matrix object which is required for cluster validation. 
valid.E = as.matrix(EOGT.clust)
valid.P = as.matrix(PUS.clust)

intern.M <- clValid(valid.M, 2:4, clMethods=c("kmeans", "pam"),
                  validation="internal") # clValid() function is used to determine optimum clusters. Pam and K-means method have been applied in the validation procedure. 
intern.E <- clValid(valid.E , 2:4, clMethods=c("kmeans", "pam"),
                  validation="internal")
intern.P <- clValid(valid.P, 2:4, clMethods=c("kmeans", "pam"),
                  validation="internal")

summary(intern.M) # View the summary to get an idea of the optimum clusters and dunn and silhouette clusters, MCIR has 3 medoids clusters optimum according to silhoutte index
summary(intern.E) # EOGT has 4 k clusters optimum according to silhoutte index
summary(intern.P) #PUS10 has 2 k clusters optimum according to silhoutte index


######****** Clustering for EOGT Gene #######*****

#Clustering for MCIR Gene

#+++++_____******++++++#############################


MICR.2 <- pam(MICR.clust, 2) # Determination of pam clusters a more robust version of k-means
MICR.3 <- pam(MICR.clust, 3)
MICR.4 <- pam(MICR.clust, 4)

cl2 <- kmeans(MICR.clust, 2) # Determination of k clusters a more robust version of k-means
cl3 <- kmeans(MICR.clust, 3)
cl4 <- kmeans(MICR.clust, 4)

MICR.2$medoids <- cl2$centers # this is done to make it possible to plot kmeans
MICR.2$clustering <- cl2$cluster

MICR.3$medoids <- cl3$centers # this is done to make it possible to plot kmeans
MICR.3$clustering <- cl3$cluster

MICR.4$medoids <- cl4$centers # this is done to make it possible to plot kmeans
MICR.4$clustering <- cl4$cluster

clusplot(MICR.2, labels = 0, main = "2 K-means Clusters of SNPs MCIR Gene for Asian and European Race") #clustplot function to visualize the k clusters 

clusplot(MICR.3, labels = 0, main = "3 K- means Clusters of SNPs MCIR Gene for Asian and European Race")

clusplot(MICR.4, labels = 0, main = "4-Kmeans Clusters of SNPs MCIR Gene for Asian and European Race")


MICR.22 <- pam(MICR.clust, 2)  # Determination of pam clusters a more robust version of k-means
MICR.33 <- pam(MICR.clust, 3)
MICR.44 <- pam(MICR.clust, 4)

clusplot(MICR.22, labels = 0, main = "2-Medoids Clusters of SNPs MCIR Gene for Asian and European Race") #clustplot function to visualize the medoids clusters 

clusplot(MICR.33, labels = 0, main = "3-Medoids Clusters of SNPs MCIR Gene for Asian and European Race", color = "TRUE")

clusplot(MICR.44, labels = 0, main = "4-Medoids Clusters of SNPs MCIR Gene for Asian and European Race")

#### Unfortunately you cannot visualize the clusters observations by race in the clusplot function thus the below operation is conducted to view the clusters by race in the legend #### Function retrieved from https://rpkgs.datanovia.com/ggpubr/reference/ggscatter.html

# Coordinates of individuals
M.coord <- as.data.frame(get_pca_ind(MICR.pr)$coord)
# Add clusters obtained using the K-means algorithm
M.coord$cluster <- factor(MICR.33$clustering) # for visulizaton of 3 medoids-clusters
# Add Racefrom the original data sett
M.coord$Race <- df.ohe$RACE
eigenvalue <- round(get_eigenvalue(MICR.pr), 1)
variance.percent <- eigenvalue$variance.percent
ggscatter(
  M.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Race", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" ) , title = "3 Medoids Clusters of SNPs MCIR Gene for Asian and European Race"
) + stat_mean(aes(color = cluster), size = 4)


######****** Clustering for EOGT Gene #######*****

#Clustering for EOGT Gene

#+++++_____******++++++#############################


EOGT.2 <- pam(EOGT.clust, 2)
EOGT.3 <- pam(EOGT.clust, 3)
EOGT.4 <- pam(EOGT.clust, 4)

cl2.1 <- kmeans(EOGT.clust, 2)
cl3.1 <- kmeans(EOGT.clust, 3)
cl4.1 <- kmeans(EOGT.clust, 4)

EOGT.2$medoids <- cl2.1$centers # this is done to make it possible to plot kmeans
EOGT.2$clustering <- cl2.1$cluster

EOGT.3$medoids <- cl3.1$centers # this is done to make it possible to plot kmeans
EOGT.3$clustering <- cl3.1$cluster

EOGT.4$medoids <- cl4.1$centers # this is done to make it possible to plot kmeans
EOGT.4$clustering <- cl4.1$cluster

clusplot(EOGT.2, labels = 0, main = "2 K-means Clusters of SNPs EOGT Gene for Asian and European Race")

clusplot(EOGT.3, labels = 0, main = "3 K-means Clusters of SNPs EOGT Gene for Asian and European Race")

clusplot(EOGT.4, labels = 0, main = "4 K-means Clusters of SNPs EOGT Gene for Asian and European Race")


EOGT.22 <- pam(EOGT.clust, 2)
EOGT.33 <- pam(EOGT.clust, 3)
EOGT.44 <- pam(EOGT.clust, 4)

clusplot(EOGT.22, labels = 0, main = "2 Medoids Clusters of SNPs EOGT Gene for Asian and European Race")

clusplot(EOGT.33, labels = 0, main = "3 Medoids Clusters of SNPs EOGT Gene for Asian and European Race")

clusplot(EOGT.44, labels = 0, main = "4 Medoids Clusters of SNPs EOGT Gene for Asian and European Race")

# Coordinates of individuals
E.coord <- as.data.frame(get_pca_ind(EOGT.pr)$coord)
# Add clusters obtained using the K-means algorithm
E.coord$cluster <- factor(EOGT.4$clustering) # for visulizaton of 4 k-clusters
# Add Racefrom the original data sett
E.coord$Race <- df.ohe$RACE
eigenvalue <- round(get_eigenvalue(EOGT.pr), 1)
variance.percent <- eigenvalue$variance.percent
ggscatter(
  E.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Race", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )"), title = "4 K-means Clusters of SNPs EOGT Gene for Asian and European Race"
) + stat_mean(aes(color = cluster), size = 4) 


######****** Clustering for PUS10 Gene #######*****

#Clustering for PUS10 Gene

#+++++_____******++++++#############################

PUS.2 <- pam(PUS.clust, 2)
PUS.3 <- pam(PUS.clust, 3)
PUS.4 <- pam(PUS.clust, 4)

cl2.2 <- kmeans(PUS.clust, 2)
cl3.2 <- kmeans(PUS.clust, 3)
cl4.2 <- kmeans(PUS.clust, 4)

PUS.2$medoids <- cl2.2$centers # this is done to make it possible to plot kmeans
PUS.2$clustering <- cl2.2$cluster

PUS.3$medoids <- cl3.2$centers # this is done to make it possible to plot kmeans
PUS.3$clustering <- cl3.2$cluster

PUS.4$medoids <- cl4.2$centers # this is done to make it possible to plot kmeans
PUS.4$clustering <- cl4.2$cluster

clusplot(PUS.2, labels = 0 , main = "2 K-means Clusters of SNPs PUS10 Gene for Asian and European Race")

clusplot(PUS.3, labels = 0, main = "3 K-means Clusters of SNPs PUS10 Gene for Asian and European Race")

clusplot(PUS.4, labels = 0, main = "4 K-means Clusters of SNPs PUS10 Gene for Asian and European Race")

PUS.22 <- pam(PUS.clust, 2)
PUS.33 <- pam(PUS.clust, 3)
PUS.44 <- pam(PUS.clust, 4)

clusplot(PUS.22, labels = 0, main = "2 Medoids Clusters of SNPs PUS10 Gene for Asian and European Race")

clusplot(PUS.33, labels = 0, main = "3 Medoids Clusters of SNPs PUS10 Gene for Asian and European Race")

clusplot(PUS.44, labels = 0, main = "4 Medoids Clusters of SNPs PUS10 Gene for Asian and European Race")

#### Creation of 2-k clusters with points labelled by race and this has been done because previously we determined 2 clusters to be optimum. ####

# Coordinates of individuals
P.coord <- as.data.frame(get_pca_ind(PUS.pr)$coord)
# Add clusters obtained using the K-means algorithm
P.coord$cluster <- factor(PUS.2$clustering) # for visulizaton of 4 k-clusters
# Add Racefrom the original data sett
P.coord$Race <- df.ohe$RACE
eigenvalue <- round(get_eigenvalue(PUS.pr), 1)
variance.percent <- eigenvalue$variance.percent
ggscatter(
  P.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Race", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" ) , title = "2 K-means Clusters of SNPs PUS10 Gene for Asian and European Race"
) + stat_mean(aes(color = cluster), size = 4)


#### Creation of Biplots###

pca2d(MICR.pr, biplot = TRUE, group = MCIR[,45], legend = "left") # Creaton of PCA biplot woth loadings shown as the SNP's
pca2d(EOGT.pr, biplot = TRUE, group = EOGT[,447], legend = "left")
pca2d(PUS.pr, biplot = TRUE, group = PUS10[,327], legend = "left")
pca3d(MICR.pr, legend = "left", group = MCIR[,45]) # creation of a 3d biplot 
pca3d(EOGT.pr, legend = "left", group = EOGT[,447])
pca3d(PUS.pr, legend = "left", group = PUS10[,327])

#### Creation of Heatmaps ####
heatmap(valid.M) # Heat maps allows to visualize the areas that seem the most similar. 
heatmap(valid.E)
heatmap(valid.P)

# Production of PC1 vs PC1 plot for EOGT gene and PUS10 gene. This has been done to capture the genetic diversity these two genes produce when combined together. 

PUS10.PC1 = PUS.pr$x[,1]
EOGT.PC1 = EOGT.pr$x[,1]

Race = as.character(df.ohe$RACE)

final  = cbind(PUS10.PC1,EOGT.PC1,Race)

final = as.data.frame(final)

ggplot(final, aes(x = PUS10.PC1, y = EOGT.PC1)) +
  geom_point(aes(color = Race)) + theme(axis.text.x=element_blank(), axis.text.y=element_blank()) + ggtitle("EOGT PC1 vs PUS10 PC1")
