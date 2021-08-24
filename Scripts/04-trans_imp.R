## Log transformation and imputation
## 10 July

# Load packages
library(tidyverse)
library(imputeLCMD)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom parameters
source("graph_param.R")

for (i in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh.rds"))
  # Transformationlod
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

# for (i in 4:5){
#   # Load data
#   covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
#   expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh.rds"))
#   # Transformation
#   expo = log10(expo)
#   saveRDS(expo, paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log.rds"))
#   # Missing rate per participant
#   covars$NA_prop = apply(expo, 1, function(x) sum(is.na(x))/ncol(expo))
#   print(max(covars$NA_prop) < 0.8)
#   
#   # impute.QRILC
#   expo_imp = NULL
#   for (k in levels(covars$Batch)){
#     set.seed(7)
#     tmp = impute.QRILC(t(expo[covars$Batch == k,]))
#     tmp=t(tmp[[1]])
#     expo_imp = rbind(expo_imp, tmp)
#   }
#   all(rownames(expo_imp) == rownames(expo))
#   all(colnames(expo_imp) == colnames(expo))
#   
#   # Save data sets
#   saveRDS(expo_imp, paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
#   saveRDS(covars, paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_naprop.rds"))
# }
