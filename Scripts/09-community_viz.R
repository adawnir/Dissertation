# Community detection Visulisation
# Rin 4 Aug

## Load packages
library(tidyverse)
library(igraph)
library(colorspace)
library(focus)

# Initialise
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

suffix = c("lux","fra","gs","pooled3","pooled2")

# Load custom
source("functions.R")
source("graph_param.R")

for (m in 1:length(batches)){
  A = readRDS(paste0("../Results/",filepaths[m],"/Graphical_network_children_adjacency_matrix.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  
  mynode_colour = family.colours[as.character(covars$Family.ID)]
  names(mynode_colour) = rownames(covars)
  
  # Getting igraph object
  mygraph=Graph(adjacency=A,
                node_label = as.character(covars$Family.ID),
                node_colour=mynode_colour,
                satellites = TRUE)
  
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children.pdf"))
  # par(mar=c(0,0,5,0), xpd = TRUE)
  par(mar=c(0,0,0,0))
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       vertex.label.cex = 0.7,
       edge.width = 1,
       vertex.size = 9)
  # legend("top",pch=19, inset = c(0,-0.15),
  #        col = family.colours[levels(covars$Family.ID)],
  #        legend = levels(covars$Family.ID), cex = 0.6,
  #        ncol = ceiling(length(levels(covars$Family.ID))/5), bty = "n")
  dev.off()
  
  ### Community detection
  # Louvain
  lc = cluster_louvain(mygraph)
  
  # Set cluster colour
  cluster_colour = rainbow(length(unique(lc$membership)), alpha = 0.6)[membership(lc)]
  names(cluster_colour) = membership(lc)
  
  E(mygraph)$color = ifelse(igraph::crossing(lc, mygraph),
                            alpha('grey',0.6),
                            darken(cluster_colour[as.character(membership(lc)[ends(mygraph, E(mygraph))[,1]])],0.5))  
  
  # Plot clusters in graph
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community.pdf"))
  par(mar=c(0,0,0,5), xpd = TRUE)
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       edge.width = 1,
       vertex.size = 9,
       vertex.label.cex = 0.7,
       vertex.frame.color = darken(cluster_colour,0.5),
       mark.groups = communities(lc),
       mark.col = cluster_colour[as.character(1:length(communities(lc)))],
       mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
       vertex.color = mynode_colour)
  legend("right",pch=19, inset = c(-0.15,0),
         col = unique(cluster_colour[order(names(cluster_colour))]),
         legend = 1:length(unique(cluster_colour)), cex = 0.6,
         ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
  dev.off()
  
  mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
  # Plot clusters in graph separately
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_split.pdf"))
  par(mar=c(0,0,0,5), xpd = TRUE)
  set.seed(1)
  plot(mygraph_split, layout=layout_with_fr(mygraph_split),
       edge.width = 1,
       vertex.size = 9,
       vertex.label.cex = 0.7,
       vertex.frame.color = darken(cluster_colour,0.5),
       mark.groups = communities(lc),
       mark.col = cluster_colour[as.character(1:length(communities(lc)))],
       mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
       vertex.color = mynode_colour)
  legend("right",pch=19, inset = c(-0.15,0),
         col = unique(cluster_colour[order(names(cluster_colour))]),
         legend = 1:length(unique(cluster_colour)), cex = 0.6,
         ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
  dev.off()
  
  if(m %in% c(4,5)){
    mynode_colour = batch.colours[as.character(covars$Batch)]
    names(mynode_colour) = rownames(covars)
    
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_batch_community.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph, layout=layout_with_fr(mygraph),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = batch.colours[levels(covars$Batch)],
           legend = levels(covars$Batch), cex = 0.6,
           horiz = TRUE, bty = "n")
    dev.off()
    
    mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
    # Plot clusters in graph separately
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_batch_split.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph_split, layout=layout_with_fr(mygraph_split),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = batch.colours[covars$Batch])
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = batch.colours[levels(covars$Batch)],
           legend = levels(covars$Batch), cex = 0.6,
           horiz = TRUE, bty = "n")
    dev.off()  
  }
  
  if(m %in% c(2,4)){
    mynode_colour = region.colours[as.character(covars$Region)]
    names(mynode_colour) = rownames(covars)
    
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_region_community.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph, layout=layout_with_fr(mygraph),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = region.colours[levels(covars$Region)],
           legend = levels(covars$Region), cex = 0.6,
           ncol = ceiling(length(levels(covars$Region))/3), bty = "n")
    dev.off()
    
    mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
    
    # Plot clusters in graph separately
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_region_split.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph_split, layout=layout_with_fr(mygraph_split),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = region.colours[levels(covars$Region)],
           legend = levels(covars$Region), cex = 0.6,
           ncol = ceiling(length(levels(covars$Region))/3), bty = "n")
    dev.off()  
    
    # Department
    mynode_colour = depart.colours[as.character(covars$Department)]
    names(mynode_colour) = rownames(covars)
    
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_depart_community.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph, layout=layout_with_fr(mygraph),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = depart.colours[levels(covars$Department)],
           legend = levels(covars$Department), cex = 0.6,
           ncol = ceiling(length(levels(covars$Department))/3), bty = "n")
    dev.off()
    
    mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
    
    # Plot clusters in graph separately
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_depart_split.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph_split, layout=layout_with_fr(mygraph_split),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = depart.colours[levels(covars$Department)],
           legend = levels(covars$Department), cex = 0.6,
           ncol = ceiling(length(levels(covars$Department))/3), bty = "n")
    dev.off()  
  }
}