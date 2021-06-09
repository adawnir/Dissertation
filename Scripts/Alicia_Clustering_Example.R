# set path to include my libraries
.libPaths(c("/rds/general/user/cls1017/home/R/x86_64-redhat-linux-gnu-library/3.6",
            "/usr/lib64/R/library",
            "/usr/share/R/library"))

rm(list=ls())

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(FactoMineR)

train_df<- readRDS("/rds/general/user/cls1017/home/ML_project_mentalhealth/Data/imp_train_data.rds")
test_df<- readRDS("/rds/general/user/cls1017/home/ML_project_mentalhealth//Data/imp_test_data.rds")

df<- rbind(train_df, test_df)

df_X = df[,colnames(df) != "HDoct6"]
df_y = df[,colnames(df) == "HDoct6"]

# Agglomerative Hierarchical Clustering -----------------------------------

d <- daisy(df_X, metric = c("gower"))

# hc1 <- hclust(d, method = "complete")
# 
# plot(hc1, cex=0.6, hang=-1)
# 
# ## Compute with agnes
# hc2 <- agnes(X_train, method = "complete")
# ## Agglomerative coefficient
# hc2$ac
# 
# ## methods to assess
# m <- c( "average", "single", "complete", "ward")
# names(m) <- c( "average", "single", "complete", "ward")
# 
# # function to compute coefficient
# ac <- function(x) {
#   agnes(X_train, method = x)$ac
# }
# 
# map_dbl(m, ac)
# 
# hc3 <- agnes(X_train, method = "ward")
# pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 

# # Ward's method
hc5 <- hclust(d, method = "ward.D2" )

# Enhanced hierarchical clustering
res.hc <- eclust(df_X, "hclust", k=2,
                 hc_method = "ward.D2", graph = FALSE) 
print(head(res.hc$cluster, 15))

# Dendrogram
pdf("/rds/general/user/cls1017/home/ML_project_mentalhealth/Figures/cluster_dendro.pdf")
fviz_dend(hc5, rect = TRUE, show_labels = FALSE) 
dev.off()

# Silhouette analysis
pdf("/rds/general/user/cls1017/home/ML_project_mentalhealth/Figures/cluster_sil_plot.pdf")
fviz_silhouette(hc5)
dev.off()

# Silhouette information
silinfo <- res.hc$silinfo
print(names(silinfo))
# Average silhouette width of each cluster
print(silinfo$clus.avg.widths)

# Check scores
table(df$HDoct6, res.hc$cluster)

library("fpc")
# Compute cluster stats
case <- as.numeric(df$HDoct6)
clust_stats <- cluster.stats(d = dist(df_X), 
                             case, res.hc$cluster)
# Corrected Rand index
clust_stats$corrected.rand



# Determining Optimal Clusters --------------------------------------------
pdf("/rds/general/user/cls1017/home/ML_project_mentalhealth/Figures/cluster_elbow.pdf")
fviz_nbclust(df_X, FUN = hcut, method = "wss",
             k.max= 10) +
  ggtitle("Elbow method")
dev.off()
pdf("/rds/general/user/cls1017/home/ML_project_mentalhealth/Figures/cluster_silhouette.pdf")
fviz_nbclust(df_X, FUN = hcut, method = "silhouette",
             k.max= 10)+
  ggtitle("Silouette method")# Silhouette 
dev.off()

# # Ward's method
# hc5 <- hclust(d, method = "ward.D2" )
# Enhanced hierarchical clustering
res.hc <- eclust(df_X, "hclust", k=2,
                 hc_method = "ward.D2", graph = FALSE) 
print(head(res.hc$cluster, 15))

# Dendrogram
pdf("/rds/general/user/cls1017/home/ML_project_mentalhealth/Figures/cluster_dendro.pdf")
fviz_dend(res.hc, rect = TRUE, show_labels = FALSE) 
dev.off()

# Silhouette analysis
pdf("/rds/general/user/cls1017/home/ML_project_mentalhealth/Figures/cluster_sil_plot.pdf")
fviz_silhouette(res.hc)
dev.off()

# Silhouette information
silinfo <- res.hc$silinfo
print(names(silinfo))
# Average silhouette width of each cluster
print(silinfo$clus.avg.widths)

# Check scores
table(df$HDoct6, res.hc$cluster)

library("fpc")
# Compute cluster stats
case <- as.numeric(df$HDoct6)
clust_stats <- cluster.stats(d = dist(df_X), 
                             case, res.hc$cluster)
# Corrected Rand index
clust_stats$corrected.rand

