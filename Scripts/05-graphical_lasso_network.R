## Graphical Lasso Network
## Rin on 14 July 

## Load packages and source functions
LoadPackages=function(packages){
  for (i in 1:length(packages)){
    suppressPackageStartupMessages(library(packages[i], character.only=TRUE))
  }
}

LoadPackages(c("pheatmap","corpcor","abind","parallel",
               "RColorBrewer","igraph","ppcor","mvtnorm",
               "pROC","glasso","stabs","huge","pulsar",
               "QUIC","glassoFast","colorspace","glmnet","tidyverse"))
source("penalisation_functions.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  # Running stability selection
  t0=Sys.time()
  out=CalibrateNetwork(data=expo, K=100, tau=0.5, PFER_thr=20, refine_calib_grid=FALSE, verbose=FALSE)
  t1=Sys.time()
  print(t1-t0)
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(out, paste0("../Results/",filepaths[i],"/Graphical_lasso_output.rds"))

  # out = readRDS(paste0("../Results/",filepaths[i],"/Graphical_lasso_output.rds"))
  pdf(paste0("../Figures/",filepaths[i],"/Graphical_lasso_output.pdf"))
  CalibrationPlot(out)
  dev.off()
  
  # Adjacency matrix of the calibrated network
  A=CalibratedAdjacency(out)
  saveRDS(A, paste0("../Results/",filepaths[i],"/Graphical_lasso_adjacency_matrix.rds"))
  
  annot_sub = annot[colnames(expo)]
  mynode_colour = annot.colours[annot_sub]
  names(mynode_colour) = colnames(expo)
  
  # Getting igraph object
  mygraph=GetGraph(adjacency=A,node_color=mynode_colour)
  
  pdf(paste0("../Figures/",filepaths[i],"/Graphical_lasso_network.pdf"))
  par(mar=rep(0,4))
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph))
  legend("topleft",pch=19,
         col = annot.colours[unique(annot_sub)],
         legend = unique(annot_sub), cex = 0.6)
  dev.off()
  }
