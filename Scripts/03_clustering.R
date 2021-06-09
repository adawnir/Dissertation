## Clustering
## Rin Wada 24 May

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(cluster)

# Load data
mat=readRDS("../Processed/Chemical_component_matrix_subset_transformed_imputed.rds")
covar=readRDS("../Processed/Participant_covariate_info.rds")
chem=readRDS("../Processed/Chemical_component_info_subset.rds")

# Rename rows with unique id
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
sibling.number = substrRight(covar$SIBLING, 1)
sibling.number[-which(sibling.number%in%c(1:3))]=""
family.number = gsub("amily ","",covar$FAMILY.ID)
myid = paste0(family.number,".",sibling.number)
rownames(mat) = myid

# Correlation matrix
chem$Group = str_to_title(chem$Group)
mat_col = data.frame(Family = chem$Group)
rownames(mat_col) = chem$Component
mat_colors = list(Family = brewer.pal(length(unique(chem$Group)), "Set3"))
names(mat_colors$Family) = unique(chem$Group)
pdf("../Figures/Corr_mat.pdf", width = 8, height = 6)
pheatmap(cor(mat), cellwidth = 10, cellheight = 10,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

# Dissimilarity matrix
dist = dist(mat, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc = hclust(dist, method = "complete")

# Plot the obtained dendrogram
pdf("../Figures/HC_dendogram.pdf")
plot(hc, cex = 0.6, hang = -1)
dev.off()

avg_sil = function(k) {
  res = cutree(hc, k=k)
  ss = silhouette(res, dist)
  mean(ss[,3])
}

# Compute and plot wss for k = 2 to k = 19
avg_ss = NULL
for(i in 2:19){
  avg_ss = c(avg_ss,avg_sil(i))
}

pdf("../Figures/HC_ss_complete.pdf")
plot(avg_ss, col = "navy", type = "b", ylab = "Average silhouette score",
     xlab = "Number of clusters")
dev.off()

res_complete = cutree(hc, k=2)

# Hierarchical clustering using Ward Linkage
hc = hclust(dist, method = "ward.D2")

# Compute and plot wss for k = 2 to k = 19
avg_ss = NULL
for(i in 2:19){
  avg_ss = c(avg_ss,avg_sil(i))
}

pdf("../Figures/HC_ss_ward.pdf")
plot(avg_ss, col = "navy", type = "b", ylab = "Average silhouette score",
     xlab = "Number of clusters")
dev.off()

res_ward = cutree(hc, k=2)

# Plot the obtained dendrogram
pdf("../Figures/HC_dendogram_ward.pdf")
plot(hc, cex = 0.6, hang = -1)
dev.off()

# Kmeans
avg_sil = function(k) {
  km.res = kmeans(mat, centers = k, nstart = 25)
  ss = silhouette(km.res$cluster, dist(mat))
  mean(ss[, 3])
}

pdf("../Figures/Kmeans_ss.pdf")
plot(avg_ss, col = "navy", type = "b", ylab = "Average silhouette score",
     xlab = "Number of clusters")
dev.off()

# Compute and plot wss for k = 2 to k = 19
avg_ss = NULL
for(i in 2:19){
  avg_ss = c(avg_ss,avg_sil(i))
}

res_kmeans = kmeans(mat, centers = 2, nstart = 25)

# PCA
mypca=prcomp(mat)
S = mypca$x
ev = summary(mypca)$importance[2,]

unique = unique(family.number)
tab = as.data.frame(table(family.number)) %>%
  slice(match(unique,family.number))
rep = tab[,2]
colour = rep(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(family.number))),
             rep)
compare = sapply(family.number, function(x) x==family.number)
start = which(compare, arr.ind=TRUE)[,1]
end = which(compare, arr.ind=TRUE)[,2]

# Creating the score plots along the first 2 PCs
pdf("../Figures/PCA_hc_complete.pdf")
par(mar = c(5,5,1,1))
plot(S[,1:2], pch=ifelse(res_complete==1, 19, 17), las=1, col=colour, cex = 1.5,
     xlab=substitute(PC[1]*" ("*a*"% e.v.)", list(a=round(ev[1]*100, digits=2))), 
     ylab=substitute(PC[2]*" ("*a*"% e.v.)", list(a=round(ev[2]*100, digits=2))), 
     cex.lab=1.5)
text(S[,1:2], labels=myid, pos=3, offset=0.5, cex = 0.5)
segments(S[start,1],S[start,2],S[end,1],S[end,2])
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()

# Creating the score plots along the first 2 PCs
pdf("../Figures/PCA_hc_ward.pdf")
par(mar = c(5,5,1,1))
plot(S[,1:2], pch=ifelse(res_ward==1, 19, 17), las=1, col=colour, cex = 1.5,
     xlab=substitute(PC[1]*" ("*a*"% e.v.)", list(a=round(ev[1]*100, digits=2))), 
     ylab=substitute(PC[2]*" ("*a*"% e.v.)", list(a=round(ev[2]*100, digits=2))), 
     cex.lab=1.5)
text(S[,1:2], labels=myid, pos=3, offset=0.5, cex = 0.5)
segments(S[start,1],S[start,2],S[end,1],S[end,2])
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()

pdf("../Figures/PCA_kmeans.pdf")
par(mar = c(5,5,1,1))
plot(S[,1:2], pch=ifelse(res_kmeans$cluster==1, 19, 17), las=1, col=colour, cex = 1.5,
     xlab=substitute(PC[1]*" ("*a*"% e.v.)", list(a=round(ev[1]*100, digits=2))), 
     ylab=substitute(PC[2]*" ("*a*"% e.v.)", list(a=round(ev[2]*100, digits=2))), 
     cex.lab=1.5)
text(S[,1:2], labels=myid, pos=3, offset=0.5, cex = 0.5)
segments(S[start,1],S[start,2],S[end,1],S[end,2])
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()

all(res_complete==res_ward)
all(res_complete==res_kmeans$cluster)

covar$FAMILY.ID
covar$Cluster.k2 = res_complete

# Concentration matrix with cluster annotations
mat_row = data.frame(Cluster = ifelse(covar$Cluster.k2==1,"A","B"))
rownames(mat_row) = rownames(mat)

chem$Group = str_to_title(chem$Group)
mat_col = data.frame(Family = chem$Group)
rownames(mat_col) = chem$Component

mat_colors = list(Cluster = c("forestgreen","orange"),
                  Family = brewer.pal(length(unique(chem$Group)), "Set3"))
names(mat_colors$Cluster) = unique(mat_row$Cluster)
names(mat_colors$Family) = unique(mat_col$Family)

pdf("../Figures/Concentration_mat.pdf")
pheatmap(mat,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_row,annotation_col = mat_col,annotation_colors = mat_colors)
dev.off()
