## Graphical Lasso Network Visualisation
## 14 July 

## Load packages
library(tidyverse)
library(igraph)
library(colorspace)
library(focus)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
suffix = c("lux","fra","gs","pooled3","pooled2")

# Load custom
source("functions.R")
source("graph_param.R")

for (m in 1:length(batches)){
  A = readRDS(paste0("../Results/",filepaths[m],"/Graphical_network_exposures_adjacency_matrix.rds"))
  expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  
  annot_sub = annot[colnames(expo)]
  mynode_colour = annot.colours[annot_sub]
  names(mynode_colour) = colnames(expo)
  
  # Getting igraph object
  mygraph=Graph(adjacency=A,node_colour=mynode_colour)
  
  # # Increase vertex size
  # V(mygraph)$size = V(mygraph)$size + 3
  # 
  # pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_exposures.pdf"))
  # par(mar=c(0,0,5,0), xpd = TRUE)
  # set.seed(1)
  # plot(mygraph, layout=layout_with_fr(mygraph),
  #      vertex.label.cex = 0.65)
  # legend("top",pch=19,
  #        inset=c(0,-0.15),
  #        col = annot.colours[unique(annot_sub)],
  #        legend = unique(annot_sub), cex = 0.6,
  #        ncol = ceiling(length(unique(annot_sub)))/5,
  #        bty = "n")
  # dev.off()
  
  ### Community detection
  # Louvain
  lc = cluster_louvain(mygraph)
  
  # Set cluster colour
  cluster_colour = rainbow(length(unique(lc$membership)), alpha = 0.1)[membership(lc)]
  names(cluster_colour) = membership(lc)
  
  E(mygraph)$color = ifelse(igraph::crossing(lc, mygraph),
                            alpha('grey',0.6),
                            darken(cluster_colour[as.character(membership(lc)[ends(mygraph, E(mygraph))[,1]])],0.5))  
  
  # Increase vertex size
  V(mygraph)$size = V(mygraph)$size + 5
  
  # Plot clusters in graph
  pdf(paste0("../Figures/Section1/Graphical_network_exposures_community_",suffix[m],".pdf"))
  par(mar=rep(1,4), xpd = TRUE)
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       vertex.label.cex = 0.8,
       mark.groups = communities(lc),
       mark.col = unique(cluster_colour[order(names(cluster_colour))]),
       mark.border = unique(cluster_colour[order(names(cluster_colour))]))
  dev.off()
}

{pdf("../Figures/Section1/Graphical_network_exposures_legend.pdf", height = 5, width = 4)
  grid.newpage()
  for(e in 1:length(annot.colours)) {
    grid.circle(
      2,
      e*0.55,
      0.25,
      default.units = "cm",
      gp=gpar( 
        fill=rev(annot.colours)[e],
        lwd=0
      )
    )
    grid.text(
      label = rev(names(annot.colours))[e],
      x = 2.5,
      y = e*0.55,
      default.units = "cm",
      hjust = 0
    )
  }
  dev.off()
}
