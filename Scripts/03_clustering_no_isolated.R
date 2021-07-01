## Clustering (excluding isolated children)
## Rin Wada 30 June

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(cluster)
library(focus)
library(summarise)
source("function.R")

# Load data sets
# mat_lux = readRDS("../Processed/Luxembourg/Chemical_compound_matrix_subset_trans_imp.rds")
# chem_lux = readRDS("../Processed/Luxembourg/Chemical_compound_info_subset.rds")
# covar_lux = readRDS("../Processed/Luxembourg/Participant_covariate_info_subset.rds")

mat_fra = readRDS("../Processed/France/Chemical_compound_matrix_subset_trans_imp.rds")
chem_fra = readRDS("../Processed/France/Chemical_compound_info_subset.rds")
covar_fra = readRDS("../Processed/France/Participant_covariate_info_subset.rds")

mat_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_matrix_subset_trans.rds")
chem_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_info_subset.rds")
covar_gs = readRDS("../Processed/GrandeSynthe/Participant_covariate_info_subset.rds")

# Exclude isolated children
# covar_lux = covar_lux %>%
#   filter(Family.ID != "Isolated")
covar_fra = covar_fra %>%
  filter(Family.ID != "Isolated")
covar_gs = covar_gs %>%
  filter(Family.ID != "Isolated")

# mat_lux = mat_lux[covar_lux$Indiv.ID,]
mat_fra = mat_fra[covar_fra$Indiv.ID,]
mat_gs = mat_gs[covar_gs$Indiv.ID,]

# Change rownames to family sibling id
# rownames(mat_lux) = covar_lux$Family.ID
rownames(mat_fra) = covar_fra$Family.ID
rownames(mat_gs) = covar_gs$Family.ID

# Standardisation
# mat_lux = scale(mat_lux)
mat_fra = scale(mat_fra)
mat_gs = scale(mat_gs)

# Dissimilarity matrix
# dist_lux = dist(mat_lux, method = "euclidean")
dist_fra = dist(mat_fra, method = "euclidean")
dist_gs = dist(mat_gs, method = "euclidean")

## Hierarchical clustering using Complete Linkage
# hc_lux = hclust(dist_lux, method = "complete")
hc_fra = hclust(dist_fra, method = "complete")
hc_gs = hclust(dist_gs, method = "complete")

### Focus: Stability selection based HC----
# out_lux = Clustering(mat_lux, K = 100, tau = 0.5, seed = 290621)
out_fra = Clustering(mat_fra, K = 100, tau = 0.5, seed = 290621)
out_gs = Clustering(mat_gs, K = 100, tau = 0.5, seed = 290621)

# Save outputs
# saveRDS(out_lux, "../Results/Luxembourg/Stab_clustering_output_no_isolated.rds")
saveRDS(out_fra, "../Results/France/Stab_clustering_output_no_isolated.rds")
saveRDS(out_gs, "../Results/GrandeSynthe/Stab_clustering_output_no_isolated.rds")

# pdf("../Figures/Luxembourg/Stab_clustering_calib_no_isolated.pdf")
# par(mar=c(7, 5, 7, 6))
# CalibrationPlot(out_lux)
# dev.off()
# lambda_lux = GetArgmax(out_lux)[,1]
# lambda_lux

pdf("../Figures/France/Stab_clustering_calib_no_isolated.pdf")
par(mar=c(7, 5, 7, 6))
CalibrationPlot(out_fra)
dev.off()
lambda_fra = GetArgmax(out_fra)[,1]
lambda_fra

pdf("../Figures/GrandeSynthe/Stab_clustering_calib_no_isolated.pdf")
par(mar=c(7, 5, 7, 6))
CalibrationPlot(out_gs)
dev.off()
lambda_gs = GetArgmax(out_gs)[,1]
lambda_gs

# covar_lux$hc_focus = cutree(hc_lux, k=lambda_lux)
covar_fra$hc_focus = cutree(hc_fra, k=lambda_fra)
covar_gs$hc_focus = cutree(hc_gs, k=lambda_gs)

# Plot the obtained dendrogram
# pdf("../Figures/Luxembourg/HC_complete_focus_dendrogram_no_isolated.pdf")
# plot(hc_lux, cex = 0.6, hang = -1)
# rect.hclust(hc_lux , k = lambda_lux, border = 2:6)
# dev.off()

pdf("../Figures/France/HC_complete_focus_dendrogram_no_isolated.pdf")
plot(hc_fra, cex = 0.6, hang = -1)
rect.hclust(hc_fra , k = lambda_fra, border = 2:6)
dev.off()

pdf("../Figures/GrandeSynthe/HC_complete_focus_dendrogram_no_isolated.pdf")
plot(hc_gs, cex = 0.6, hang = -1)
rect.hclust(hc_gs , k = lambda_gs, border = 2:6)
dev.off()

# Families reconstructed
# covar_lux %>%
#   filter(Family.ID != "Isolated") %>%
#   group_by(Family.ID) %>%
#   summarise(count = length(unique(hc_focus))) %>%
#   filter(count==1) %>%
#   nrow(.)/sum(unique(covar_lux$Family.ID)!="Isolated")

covar_fra %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_focus))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_fra$Family.ID)!="Isolated")

covar_gs %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_focus))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_gs$Family.ID)!="Isolated")

### Summarise: HC with F-test based cluster score----
# fs_lux = NULL
# for(i in seq(2,nrow(mat_lux)-1)){
#   member = cutree(hc_lux, k=i)
#   fs_lux = c(fs_lux,ClusteringScore(mat_lux, member)$score)
# }

fs_fra = NULL
for(i in seq(2,nrow(mat_fra)-1)){
  member = cutree(hc_fra, k=i)
  fs_fra = c(fs_fra,ClusteringScore(mat_fra, member)$score)
}

fs_gs = NULL
for(i in seq(2,nrow(mat_gs)-1)){
  member = cutree(hc_gs, k=i)
  fs_gs = c(fs_gs,ClusteringScore(mat_gs, member)$score)
}

# pdf("../Figures/Luxembourg/HC_complete_fs_no_isolated.pdf")
# plot(y = fs_lux, x = seq(2,nrow(mat_lux)-1),
#      col = "navy", type = "b",
#      ylab = "F-test based clustering score",
#      xlab = "Number of clusters")
# dev.off()

pdf("../Figures/France/HC_complete_fs_no_isolated.pdf")
plot(y = fs_fra, x = seq(2,nrow(mat_fra)-1),
     col = "navy", type = "b",
     ylab = "F-test based clustering score",
     xlab = "Number of clusters")
dev.off()

pdf("../Figures/GrandeSynthe/HC_complete_fs_no_isolated.pdf")
plot(y = fs_gs, x = seq(2,nrow(mat_gs)-1),
     col = "navy", type = "b",
     ylab = "F-test based clustering score",
     xlab = "Number of clusters")
dev.off()

# set k
# lambda_lux = seq(2,nrow(mat_lux)-1)[which.max(fs_lux)]
lambda_fra = seq(2,nrow(mat_fra)-1)[which.max(fs_fra)]
lambda_fra
lambda_gs = seq(2,nrow(mat_gs)-1)[which.max(fs_gs)]
lambda_gs

# cut trees
# covar_lux$hc_fs = cutree(hc_lux, k=lambda_lux)
covar_fra$hc_fs = cutree(hc_fra, k=lambda_fra)
covar_gs$hc_fs = cutree(hc_gs, k=lambda_gs)

# Plot the obtained dendrogram
# pdf("../Figures/Luxembourg/HC_complete_fs_dendrogram_no_isolated.pdf")
# plot(hc_lux, cex = 0.6, hang = -1)
# rect.hclust(hc_lux , k = lambda_lux, border = 2:6)
# dev.off()

pdf("../Figures/France/HC_complete_fs_dendrogram_no_isolated.pdf")
plot(hc_fra, cex = 0.6, hang = -1)
rect.hclust(hc_fra , k = lambda_fra, border = 2:6)
dev.off()

pdf("../Figures/GrandeSynthe/HC_complete_fs_dendrogram_no_isolated.pdf")
plot(hc_gs, cex = 0.6, hang = -1)
rect.hclust(hc_gs , k = lambda_gs, border = 2:6)
dev.off()

# Families reconstructed
# covar_lux %>%
#   filter(Family.ID != "Isolated") %>%
#   group_by(Family.ID) %>%
#   summarise(count = length(unique(hc_fs))) %>%
#   filter(count==1) %>%
#   nrow(.)/sum(unique(covar_lux$Family.ID)!="Isolated")

covar_fra %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_fs))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_fra$Family.ID)!="Isolated")

covar_gs %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_fs))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_gs$Family.ID)!="Isolated")

### HC with silhouette score----
# Compute and plot average ss
# ss_lux = NULL
# for(i in seq(2,nrow(mat_lux)-1)){
#   ss_lux = c(ss_lux,avg_sil_hc(i, hc_lux, dist_lux))
# }

ss_fra = NULL
for(i in seq(2,nrow(mat_fra)-1)){
  ss_fra = c(ss_fra,avg_sil_hc(i, hc_fra, dist_fra))
}

ss_gs = NULL
for(i in seq(2,nrow(mat_gs)-1)){
  ss_gs = c(ss_gs,avg_sil_hc(i, hc_gs, dist_gs))
}

# pdf("../Figures/Luxembourg/HC_complete_ss_no_isolated.pdf")
# plot(y = ss_lux, x = seq(2,nrow(mat_lux)-1),
#      col = "navy", type = "b",
#      ylab = "Average silhouette score",
#      xlab = "Number of clusters")
# dev.off()

pdf("../Figures/France/HC_complete_ss_no_isolated.pdf")
plot(y = ss_fra, x = seq(2,nrow(mat_fra)-1),
     col = "navy", type = "b",
     ylab = "Average silhouette score",
     xlab = "Number of clusters")
dev.off()

pdf("../Figures/GrandeSynthe/HC_complete_ss_no_isolated.pdf")
plot(y = ss_gs, x = seq(2,nrow(mat_gs)-1),
     col = "navy", type = "b",
     ylab = "Average silhouette score",
     xlab = "Number of clusters")
dev.off()

# set k
# lambda_lux = seq(2,nrow(mat_lux)-1)[which.max(ss_lux)]
lambda_fra = seq(2,nrow(mat_fra)-1)[which.max(ss_fra)]
lambda_fra
lambda_gs = seq(2,nrow(mat_gs)-1)[which.max(ss_gs)]
lambda_gs

# cut trees
# covar_lux$hc_ss = cutree(hc_lux, k=lambda_lux)
covar_fra$hc_ss = cutree(hc_fra, k=lambda_fra)
covar_gs$hc_ss = cutree(hc_gs, k=lambda_gs)

# Plot the obtained dendrogram
# pdf("../Figures/Luxembourg/HC_complete_ss_dendrogram_no_isolated.pdf")
# plot(hc_lux, cex = 0.6, hang = -1)
# rect.hclust(hc_lux , k = lambda_lux, border = 2:6)
# dev.off()

pdf("../Figures/France/HC_complete_ss_dendrogram_no_isolated.pdf")
plot(hc_fra, cex = 0.6, hang = -1)
rect.hclust(hc_fra , k = lambda_fra, border = 2:6)
dev.off()

pdf("../Figures/GrandeSynthe/HC_complete_ss_dendrogram_no_isolated.pdf")
plot(hc_gs, cex = 0.6, hang = -1)
rect.hclust(hc_gs , k = lambda_gs, border = 2:6)
dev.off()

# Families reconstructed
# covar_lux %>%
#   filter(Family.ID != "Isolated") %>%
#   group_by(Family.ID) %>%
#   summarise(count = length(unique(hc_ss))) %>%
#   filter(count==1) %>%
#   nrow(.)/sum(unique(covar_lux$Family.ID)!="Isolated")

covar_fra %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_ss))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_fra$Family.ID)!="Isolated")

covar_gs %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_ss))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_gs$Family.ID)!="Isolated")

### HC with total within-cluster sum of squares----
# Elbow (within-cluster sum of squares)
set.seed(123)
# wss_hc_lux = map_dbl(seq(2, nrow(mat_lux)-1), wss_hc, hc_lux, mat_lux)
wss_hc_fra = map_dbl(seq(2, nrow(mat_fra)-1), wss_hc, hc_fra, mat_fra)
wss_hc_gs = map_dbl(seq(2, nrow(mat_gs)-1), wss_hc, hc_gs, mat_gs)

# pdf("../Figures/Luxembourg/HC_complete_wss_no_isolated.pdf")
# plot(seq(2, nrow(mat_lux)-1), wss_hc_lux, col = "navy",
#      type="b", xlab="Number of clusters",
#      ylab="Total within-clusters sum of squares")
# dev.off()

pdf("../Figures/France/HC_complete_wss_no_isolated.pdf")
plot(seq(2, nrow(mat_fra)-1), wss_hc_fra, col = "navy",
     type="b", xlab="Number of clusters",
     ylab="Total within-clusters sum of squares")
dev.off()

pdf("../Figures/GrandeSynthe/HC_complete_wss_no_isolated.pdf")
plot(seq(2, nrow(mat_gs)-1), wss_hc_gs, col = "navy",
     type="b", xlab="Number of clusters",
     ylab="Total within-clusters sum of squares")
dev.off()

# Plot the obtained dendrogram
# pdf("../Figures/Luxembourg/HC_complete_wss_dendrogram_no_isolated.pdf")
# plot(hc_lux, cex = 0.6, hang = -1)
# rect.hclust(hc_lux , k = 20, border = 2:6)
# dev.off()

pdf("../Figures/France/HC_complete_wss_dendrogram_no_isolated.pdf")
plot(hc_fra, cex = 0.6, hang = -1)
rect.hclust(hc_fra , k = 13, border = 2:6)
dev.off()

pdf("../Figures/GrandeSynthe/HC_complete_wss_dendrogram_no_isolated.pdf")
plot(hc_gs, cex = 0.6, hang = -1)
rect.hclust(hc_gs , k = 7, border = 2:6)
dev.off()

# cut trees
# covar_lux$hc_wss = cutree(hc_lux, k=20)
covar_fra$hc_wss = cutree(hc_fra, k=13)
covar_gs$hc_wss = cutree(hc_gs, k=7)

# Families reconstructed
# covar_lux %>%
#   filter(Family.ID != "Isolated") %>%
#   group_by(Family.ID) %>%
#   summarise(count = length(unique(hc_wss))) %>%
#   filter(count==1) %>%
#   nrow(.)/sum(unique(covar_lux$Family.ID)!="Isolated")

covar_fra %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_wss))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_fra$Family.ID)!="Isolated")

covar_gs %>% as.data.frame %>%
  filter(Family.ID != "Isolated") %>%
  group_by(Family.ID) %>%
  summarise(count = length(unique(hc_wss))) %>%
  filter(count==1) %>%
  nrow(.)/sum(unique(covar_gs$Family.ID)!="Isolated")

### Save data sets ----
# saveRDS(covar_lux, "../Results/Luxembourg/Clustering_membership_no_isolated.rds")
saveRDS(covar_fra, "../Results/France/Clustering_membership_no_isolated.rds")
saveRDS(covar_gs, "../Results/GrandeSynthe/Clustering_membership_no_isolated.rds")

# saveRDS(mat_lux, "../Results/Luxembourg/Chemical_compound_matrix_no_isolated.rds")
saveRDS(mat_fra, "../Results/France/Chemical_compound_matrix_no_isolated.rds")
saveRDS(mat_gs, "../Results/GrandeSynthe/Chemical_compound_matrix_no_isolated.rds")

