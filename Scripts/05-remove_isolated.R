## Removing isolated children for family association analysis
## 29 July

# Load packages
library(tidyverse)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("graph_param.R")

for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  ## Remove isolated children (no siblings)
  covars = covars %>%
    filter(Family.ID != "Isolated") %>%
    mutate_if(is.factor, droplevels)
  expo = expo[rownames(covars),]
  saveRDS(covars, paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  saveRDS(expo, paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
}