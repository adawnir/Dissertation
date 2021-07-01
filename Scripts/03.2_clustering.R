## Clustering (Pooled LUX/GS)
## Rin Wada 1 July

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(cluster)
library(dbscan)
library(focus)
library(summarise)
source("function.R")

# Load data sets
mat = readRDS("../Processed/Pooled2/Chemical_compound_matrix_subset_trans_imp.rds")
chem = readRDS("../Processed/Pooled2/Chemical_compound_info_subset.rds")
covar = readRDS("../Processed/Pooled2/Participant_covariate_info_subset.rds")

# Change rownames to family sibling id
rownames(mat) = covar$Family.ID

# Standardisation
mat = scale(mat)

# Dissimilarity matrix
dist = dist(mat, method = "euclidean")

## Hierarchical clustering using Complete Linkage
hc = hclust(dist, method = "complete")

### Focus: Stability selection based HC----
out = Clustering(mat, K = 100, tau = 0.5, seed = 290621)

# Save outputs
saveRDS(out, "../Results/Pooled2/Stab_clustering_output.rds")

pdf("../Figures/Pooled2/Stab_clustering_calib.pdf")
par(mar=c(7, 5, 7, 6))
CalibrationPlot(out)
dev.off()

lambda = GetArgmax(out)[,1]
lambda

covar$hc_focus = cutree(hc, k=lambda)

# Plot the obtained dendrogram
pdf("../Figures/Pooled2/HC_complete_focus_dendrogram.pdf")
plot(hc, cex = 0.6, hang = -1)
rect.hclust(hc , k = lambda, border = 2:6)
dev.off()

# Families reconstructed
covar %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_focus))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar$Family.ID)!="Isolated")

### Summarise: HC with F-test based cluster score----
fs = NULL
for(i in seq(2,nrow(mat)-1)){
  member = cutree(hc, k=i)
  fs = c(fs,ClusteringScore(mat, member)$score)
}

pdf("../Figures/Pooled2/HC_complete_fs.pdf")
plot(y = fs, x = seq(2,nrow(mat)-1),
     col = "navy", type = "b",
     ylab = "F-test based clustering score",
     xlab = "Number of clusters")
dev.off()

# set k
lambda = seq(2,nrow(mat)-1)[which.max(fs)]

# cut trees
covar$hc_fs = cutree(hc, k=lambda)

# Plot the obtained dendrogram
pdf("../Figures/Pooled2/HC_complete_fs_dendrogram.pdf")
plot(hc, cex = 0.6, hang = -1)
rect.hclust(hc , k = lambda, border = 2:6)
dev.off()

# Families reconstructed
covar %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_fs))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar$Family.ID)!="Isolated")

### HC with silhouette score----
# Compute and plot average ss
ss = NULL
for(i in seq(2,nrow(mat)-1)){
  ss = c(ss,avg_sil_hc(i, hc, dist))
}

pdf("../Figures/Pooled2/HC_complete_ss.pdf")
plot(y = ss, x = seq(2,nrow(mat)-1),
     col = "navy", type = "b",
     ylab = "Average silhouette score",
     xlab = "Number of clusters")
dev.off()

# set k
lambda = seq(2,nrow(mat)-1)[which.max(ss)]
# cut trees
covar$hc_ss = cutree(hc, k=lambda)

# Plot the obtained dendrogram
pdf("../Figures/Pooled2/HC_complete_ss_dendrogram.pdf")
plot(hc, cex = 0.6, hang = -1)
rect.hclust(hc , k = lambda, border = 2:6)
dev.off()

# Families reconstructed
covar %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_ss))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar$Family.ID)!="Isolated")

### HC with total within-cluster sum of squares----
# Elbow (within-cluster sum of squares)
set.seed(123)
wss_hc = map_dbl(seq(2, nrow(mat)-1), wss_hc, hc, mat)

pdf("../Figures/Pooled2/HC_complete_wss.pdf")
plot(seq(2, nrow(mat)-1), wss_hc, col = "navy",
     type="b", xlab="Number of clusters",
     ylab="Total within-clusters sum of squares")
dev.off()

# Plot the obtained dendrogram
pdf("../Figures/Pooled2/HC_complete_wss_dendrogram.pdf")
plot(hc, cex = 0.6, hang = -1)
rect.hclust(hc , k = 20, border = 2:6)
dev.off()

# cut trees
covar$hc_wss = cutree(hc, k=20)

# Families reconstructed
covar %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_wss))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar$Family.ID)!="Isolated")

### Save data sets ----
saveRDS(covar, "../Results/Pooled2/Clustering_membership.rds")

