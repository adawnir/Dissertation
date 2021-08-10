## Clustering membership visualisation
## Rin Wada 22 July

# Load packages
library(ape)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(focus)
library(colorspace)

# Initialisation
rm(list=ls())
path="~/HDAML/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")

summary = NULL
for (m in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))

  sc = readRDS(paste0("../Results/",filepaths[m],"/Cluster_memberships.rds"))
  summary = c(summary, max(sc))

}
