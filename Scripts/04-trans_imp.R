## Log transformation, imputation and scaling
## Rin Wada 8 July

# Load packages
library(tidyverse)
library(imputeLCMD)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
for (i in 1:length(batches)){
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh.rds"))
  # Transformation
  expo = log10(expo)
  saveRDS(expo, paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log.rds"))
  # Missing rate per participant
  covars$NA_prop = apply(expo, 1, function(x) sum(is.na(x))/ncol(expo))
  print(max(covars$NA_prop) < 0.8)
  # impute.QRILC
  set.seed(7)
  expo_imp=impute.QRILC(t(expo))
  expo_imp=t(expo_imp[[1]])
  rownames(expo_imp) = rownames(expo)
  colnames(expo_imp) = colnames(expo)
  # Save data sets
  saveRDS(expo_imp, paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  saveRDS(covars, paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_naprop.rds"))
}
