# Community detection
# Rin 4 Aug

# Load packages
library(tidyverse)
library(igraph)
library(colorspace)
library(focus)

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

### PFER calibration ----
suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  X = t(expo)
  
  out = GraphicalModel(xdata = X, verbose = TRUE)
  
  zero = NULL
  for (n in c(Inf, seq(1000,200,-100), seq(100,10,-10))){
    out = GraphicalModel(xdata = X, PFER_thr = n, verbose = FALSE)
    zero = c(zero, sum(rowSums(Adjacency(out))==0))
  }
  
  pdf(paste0("../Figures/",filepaths[i],"/Graphical_network_children_PFER_thres.pdf"))
  par(mar = c(5, 5, 1, 1))
  plot(1:length(zero), zero,
       ylab = "Number of nodes with zero edges",
       xlab = "PFER threshold",
       xaxt = "n",
       cex.lab = 1.5,
       pch = 19, type = "b", col = batch.colours[i])
  axis(side = 1, at = 1:length(zero), labels = c(Inf, seq(1000,200,-100), seq(100,10,-10)))
  dev.off()
  
  if(sum(zero != min(zero))!=0){
    start = c(Inf, seq(1000,200,-100), seq(100,10,-10))[which(zero != min(zero))[1]-1]
    stop = c(Inf, seq(1000,200,-100), seq(100,10,-10))[which(zero != min(zero))[1]]
    
    if(is.infinite(start)){
      for (n in c(Inf, seq(10000,1000,-1000))){
        out = GraphicalModel(xdata = X, PFER_thr = n, verbose = FALSE)
        zero = c(zero, sum(rowSums(Adjacency(out))==0))
        
        pdf(paste0("../Figures/",filepaths[i],"/Graphical_network_children_PFER_thres_broad.pdf"))
        par(mar = c(5, 5, 1, 1))
        plot(1:length(zero), zero,
             ylab = "Number of nodes with zero edges",
             xlab = "PFER threshold",
             xaxt = "n",
             cex.lab = 1.5,
             pch = 19, type = "b", col = batch.colours[i])
        axis(side = 1, at = 1:length(zero), labels = c(Inf, seq(10000,1000,-1000)))
        dev.off()
      }
    } else {
      int = ifelse(start > 100, -10, -1)
      
      zero = NULL
      for (n in seq(start,stop, int)){
        out = GraphicalModel(xdata = X, PFER_thr = n, verbose = FALSE)
        zero = c(zero, sum(rowSums(Adjacency(out))==0))
      }
      
      pdf(paste0("../Figures/",filepaths[i],"/Graphical_network_children_PFER_thres_fine.pdf"))
      par(mar = c(5, 5, 1, 1))
      plot(1:length(zero), zero,
           ylab = "Number of nodes with zero edges",
           xlab = "PFER threshold",
           xaxt = "n",
           cex.lab = 1.5,
           pch = 19, type = "b", col = batch.colours[i])
      axis(side = 1, at = 1:length(zero), labels = seq(start,stop, int))
      dev.off()
    } 
  } else {
    print("Select lowest PFER threshold")
  }
}

### Network ----
suffix = c("lux","fra","gs","pooled3","pooled2")
for(i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  X = t(expo)
  
  n = c(110,920,79)[i]
  
  out = GraphicalModel(xdata = X, PFER_thr = n, verbose = FALSE)
  
  pdf(paste0("../Figures/",filepaths[i],"/Graphical_network_children_output.pdf"))
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(out)
  dev.off()
  
  # Adjacency matrix of the calibrated network
  A=Adjacency(out)
  saveRDS(A, paste0("../Results/",filepaths[i],"/Graphical_network_children_adjacency_matrix.rds"))
  
  mynode_colour = family.colours[as.character(covars$Family.ID)]
  names(mynode_colour) = rownames(expo)
  
  # Getting igraph object
  mygraph=Graph(adjacency=A,
                node_label = as.character(covars$Family.ID),
                node_colour=mynode_colour,
                satellites = TRUE)
  
  pdf(paste0("../Figures/",filepaths[i],"/Graphical_network_children.pdf"))
  par(mar=rep(0,4))
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph))
  legend("topleft",pch=19,
         col = family.colours[levels(covars$Family.ID)],
         legend = levels(covars$Family.ID), cex = 0.6,
         ncol = ceiling(length(levels(covars$Family.ID))/10))
  dev.off()
  
  ### Community detection
  # Louvain
  lc <- cluster_louvain(mygraph)
  membership(lc)
  communities(lc)
  
  # Set cluster colour
  cluster_colour = rainbow(length(unique(lc$membership)), alpha = 0.6)[membership(lc)]
  names(cluster_colour) = membership(lc)
  
  E(mygraph)$color = ifelse(crossing(lc, mygraph),
                            alpha('grey',0.6),
                            darken(cluster_colour[as.character(membership(lc)[ends(mygraph, E(mygraph))[,1]])],0.5))  
  
  # Plot clusters in graph
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       edge.width = 1,
       vertex.size = 15,
       vertex.label.cex = 0.8,
       vertex.frame.color = darken(cluster_colour,0.5),
       mark.groups = communities(lc),
       mark.col = unique(cluster_colour[order(names(cluster_colour))]),
       mark.border = unique(cluster_colour[order(names(cluster_colour))]),
       vertex.color = mynode_colour)
}
