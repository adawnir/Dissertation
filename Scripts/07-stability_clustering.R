## Clustering membership
## Rin Wada 22 July

# Load packages
library(ape)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(focus)

### Cluster membership ascertainment ----
# Initialisation
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

## Parameters
args=commandArgs(trailingOnly=TRUE)
m=as.numeric(args[1])

suffix = c("lux","fra","gs","pooled3","pooled2")
# Load data
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))

# Standardisation
expo = scale(expo)

# Dissimilarity matrix
d=dist(expo)

# Hierarchical clustering using Complete Linkage
h=hclust(d, method = "complete")
print(all(covars$Indiv.ID==h$labels))
h$labels=paste0(covars$Family.ID, "-",h$labels)

## Focus: Stability selection based HC
out = Clustering(expo, K = 100, tau = 0.5, seed = 290621,
                 implementation = HierarchicalClustering)

# Save outputs
saveRDS(out, paste0("../Results/",filepaths[m],"/Stability_clustering_output.rds"))

pdf(paste0("../Figures/",filepaths[m],"/Stability_clustering_output.pdf"))
par(mar=c(7, 5, 7, 6))
CalibrationPlot(out, clustering = TRUE)
dev.off()

# Save memberships
saveRDS(Clusters(out), paste0("../Results/",filepaths[m],"/Cluster_memberships.rds"))


# ### Plotting ----
# # Initialisation
# rm(list=ls())
# path="~/Dissertation/Scripts"
# setwd(path)
# 
# # Load custom
# source("functions.R")
# source("graph_param.R")
# 
# suffix = c("lux","fra","gs","pooled3","pooled2")
# 
# for (m in 1:length(bathces)){
#   # Load data
#   out = readRDS(paste0("../Results/",filepaths[m],"/Stability_clustering_output.rds"))
#   sc = readRDS(paste0("../Results/",filepaths[m],"/Cluster_memberships.rds"))
#   
#   # Load data
#   expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
#   covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
#   print(all(rownames(expo)==rownames(covars)))
#   
#   expo = scale(expo)
#   d=dist(expo)
#   h=hclust(d, method = "complete")
#   print(all(covars$Indiv.ID==h$labels))
#   h$labels=covars$Family.ID
#   
#   myphylo=as.phylo(h)
#   graph_edges = myphylo$edge
#   graph_net=graph.edgelist(graph_edges, directed=FALSE)
#   set.seed(1)
#   graph_layout = layout_with_kk(graph_net)
#   cluster_colour = rainbow(length(unique(sc)), alpha = 0.6)[membership(lc)]
#   names(cluster_colour) = membership(lc)
#   cluster_colour
#   {pdf(paste0("../Figures/",filepaths[m],"/Stab_clustering_graph.pdf"))
#     par(mar=c(0,5,0,5))
#     plot(graph_net, layout = graph_layout,
#          edge.width = 1,
#          vertex.size = 9,
#          vertex.label.cex = 0.7,
#          vertex.frame.color = darken(cluster_colour,0.5),
#          mark.groups = communities(lc),
#          mark.col = cluster_colour[as.character(1:length(communities(lc)))],
#          mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
#          vertex.color = mynode_colour)
#     legend("right",pch=19, inset = c(-0.15,0),
#            col = unique(cluster_colour[order(names(cluster_colour))]),
#            legend = 1:length(unique(cluster_colour)), cex = 0.6,
#            ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
#     legend("top",pch=19, inset = c(0,-0.15),
#            col = depart.colours[levels(covars$Department)],
#            legend = levels(covars$Department), cex = 0.6,
#            ncol = ceiling(length(levels(covars$Department))/3), bty = "n")
#     dev.off() 
#     dev.off()
#   }
# }
