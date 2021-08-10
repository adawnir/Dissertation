## Graphical Lasso Network
## Rin on 14 July 

## Load packages
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

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

# Load data
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh.rds"))
print(all(rownames(expo)==rownames(covars)))

# Running stability selection
t0=Sys.time()
out=GraphicalModel(xdata=expo, PFER_thr=20)
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

