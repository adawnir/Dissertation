## Graphical Lasso Network
## Rin on 14 July 

## Load packages
library(igraph)
library(colorspace)
library(focus)

# Initialise
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

## Parameters
args=commandArgs(trailingOnly=TRUE)
m=as.numeric(args[1])

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (m in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  # Running stability selection
  t0=Sys.time()
  out=GraphicalModel(Xdata=expo, PFER_thr=20)
  t1=Sys.time()
  print(t1-t0)
  ifelse(dir.exists(paste0("../Results/",filepaths[m])),"",dir.create(paste0("../Results/",filepaths[m])))
  saveRDS(out, paste0("../Results/",filepaths[m],"/Graphical_network_exposures_output.rds"))

  # out = readRDS(paste0("../Results/",filepaths[m],"/Graphical_network_exposures_output.rds"))
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_exposures_output.pdf"))
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(out)
  dev.off()
  
  # Adjacency matrix of the calibrated network
  A=Adjacency(out)
  saveRDS(A, paste0("../Results/",filepaths[m],"/Graphical_network_exposures_adjacency_matrix.rds"))
  
  annot_sub = annot[colnames(expo)]
  mynode_colour = annot.colours[annot_sub]
  names(mynode_colour) = colnames(expo)
  
  # Getting igraph object
  mygraph=Graph(adjacency=A,node_color=mynode_colour)
  
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_exposures.pdf"))
  par(mar=rep(0,4))
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       vertex.label.cex = 0.7)
  legend("topleft",pch=19,
         col = annot.colours[unique(annot_sub)],
         legend = unique(annot_sub), cex = 0.6)
  dev.off()
  
  ### Community detection
  # Louvain
  lc <- cluster_louvain(mygraph)
  
  # Set cluster colour
  cluster_colour = rainbow(length(unique(lc$membership)), alpha = 0.6)[membership(lc)]
  names(cluster_colour) = membership(lc)
  
  E(mygraph)$color = ifelse(crossing(lc, mygraph),
                            alpha('grey',0.6),
                            darken(cluster_colour[as.character(membership(lc)[ends(mygraph, E(mygraph))[,1]])],0.5))  
  
  # Plot clusters in graph
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_exposures_community.pdf"))
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       vertex.label.cex = 0.7,
       mark.groups = communities(lc),
       mark.col = unique(cluster_colour[order(names(cluster_colour))]),
       mark.border = unique(cluster_colour[order(names(cluster_colour))]))
  legend("topleft",pch=19,
         col = annot.colours[unique(annot_sub)],
         legend = unique(annot_sub), cex = 0.6)
  dev.off()
}
