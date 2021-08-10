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

## Focus: Stability selection based HC
out = Clustering(expo, K = 100, tau = 0.5, seed = 290621)
# Save outputs
saveRDS(out, paste0("../Results/",filepaths[m],"/Stability_clustering_output.rds"))
# Save memberships
saveRDS(Clusters(out), paste0("../Results/",filepaths[m],"/Stable_clusters.rds"))
