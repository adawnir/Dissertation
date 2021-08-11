# Load packages
library(ape)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(focus)
library(colorspace)

### Plotting ----
# Initialisation
rm(list=ls())
path="~/HDAML/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")

for (m in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))

  out = readRDS(paste0("../Results/",filepaths[m],"/Stability_clustering_output.rds"))
  sc = readRDS(paste0("../Results/",filepaths[m],"/Cluster_memberships.rds"))
  sc = sc$stab_cluster
  names(sc) = rownames(covars)

  expo = scale(expo)
  d=dist(expo)
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))

  myphylo=as.phylo(h)
  mygraph = as.igraph(myphylo,directed=FALSE)
  
  mynode_colour = family.colours[as.character(covars$Family.ID)]
  names(mynode_colour) = rownames(covars)

  V(mygraph)$size = ifelse(grepl("Node",V(mygraph)$name), 0, 9)
  V(mygraph)$label = as.character(covars[V(mygraph)$name,"Family.ID"])
  V(mygraph)$color = mynode_colour[V(mygraph)$name]
  V(mygraph)$frame.color = darken(cluster_colour,0.5)[sc[V(mygraph)$name]]
  V(mygraph)$shape = ifelse(grepl("Node",V(mygraph)$name),"none","circle")
  V(mygraph)$label.family = "sans"
  V(mygraph)$label.cex = 0.7
  V(mygraph)$label.color = "grey20"
  
  # Formatting edges
  E(mygraph)$color = "grey60"
  E(mygraph)$width = 1
  
  cluster_colour = rainbow(length(unique(sc)), alpha = 0.6)[sc]
  names(cluster_colour) = sc
  
  communities = lapply(unique(sc), function(i) rownames(covars)[sc==i])
  names(communities) = unique(sc)
  
  {pdf(paste0("../Figures/",filepaths[m],"/Stab_clustering_graph.pdf"))
    par(mar=c(0,0,0,5))
    set.seed(1)
    plot(mygraph,
         mark.groups = communities,
         mark.col = cluster_colour[as.character(1:length(communities))],
         mark.border = darken(cluster_colour[as.character(1:length(communities))], 0.5))
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/20), bty = "n")
    # legend("top",pch=19, inset = c(0,-0.15),
    #        col = depart.colours[levels(covars$Department)],
    #        legend = levels(covars$Department), cex = 0.6,
    #        ncol = ceiling(length(levels(covars$Department))/3), bty = "n")
    dev.off()
  }
}
