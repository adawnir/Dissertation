## Fixed number clustering
## Rin on 10 Aug

# Load packages
library(tidyverse)
library(RColorBrewer)

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  # Standardisation
  expo = scale(expo)
  
  # Dissimilarity matrix
  d=dist(expo)
  
  # Hierarchical clustering using Complete Linkage
  h=hclust(d, method = "complete")
  
  #Fixed: k = # of families
  fc = cutree(h, k=length(unique(covars$Family.ID)))
  
  # Save memberships
  saveRDS(fc, paste0("../Results/",filepaths[i],"/Fixed_clusters.rds"))
}
