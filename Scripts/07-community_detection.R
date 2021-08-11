# Community detection
# Rin 4 Aug

# Load packages
library(tidyverse)
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

### Network ----
suffix = c("lux","fra","gs","pooled3","pooled2")
# Load data
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))

# Standardisation
expo = scale(expo)

X = t(expo)

n = c(120,240,110,n,800)[m]
start.type = ifelse(m %in% c(4,5),"cold","warm")
q
out = GraphicalModel(xdata = X, PFER_thr = n, start = start.type)

pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_output.pdf"))
par(mar = c(7, 5, 7, 6))
CalibrationPlot(out)
dev.off()

# Adjacency matrix of the calibrated network
A=Adjacency(out)
saveRDS(A, paste0("../Results/",filepaths[m],"/Graphical_network_children_adjacency_matrix.rds"))