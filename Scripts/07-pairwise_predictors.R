## Pairwise: covariate mean and difference/delta exposure
## Rin Wada 15 July

# Load packages
library(tidyverse)

# Initialisation
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

#### Fixed clustering ####
for (m in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  expo  = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  ### Pair-wise covariates and delta exposures ----
  compare_family = sapply(covars$Family.ID, function(x) x==covars$Family.ID)
  rownames(compare_family)=colnames(compare_family)=rownames(covars)
  compare_family = compare_family[lower.tri(compare_family)]
  
  pairs = sapply(rownames(covars), function(x) paste0(x,"--",rownames(covars)))
  pairs = pairs[lower.tri(pairs)]
  
  mygrep = which(!compare_family)
  pairs = pairs[mygrep]
  
  tmp = combn(nrow(covars), 2,
              function(x) ifelse(sum(is.na(covars$Age[x[1:2]]))==0,
                                 mean(covars$Age[x[1]], covars$Age[x[2]]),NA))
  X=matrix(tmp[mygrep], nrow = length(pairs))
  colnames(X)[ncol(X)] = "age_mu"
  rownames(X) = pairs
  
  X=cbind(as.data.frame(X), as.numeric(dist(covars$Age))[mygrep])
  colnames(X)[ncol(X)] = "age_diff"
  
  if (m==1){
    tmp = combn(nrow(covars), 2, function(x) mean(covars$Weight[x[1]], covars$Weight[x[2]]))
    X=cbind(X, tmp[mygrep])
    colnames(X)[ncol(X)] = "weight_mu"
    
    X=cbind(X, as.numeric(dist(covars$Weight))[mygrep])
    colnames(X)[ncol(X)] = "weight_diff"
    
    tmp = combn(nrow(covars), 2, function(x) mean(covars$Length[x[1]], covars$Length[x[2]]))
    X=cbind(X, tmp[mygrep])
    colnames(X)[ncol(X)] = "length_mu"
    
    X=cbind(X, as.numeric(dist(covars$Length))[mygrep])
    colnames(X)[ncol(X)] = "length_diff"
  }
  
  gender = combn(nrow(covars), 2, function(x) as.character(covars$Gender)[x])
  tmp = apply(gender, 2,
              function(x) ifelse(sum(is.na(x))>0, NA,
                                 ifelse(length(unique(x))>1,"Different gender",
                                        ifelse(unique(x)=="Male", "All male", "All female"))))
  X=cbind(X, tmp[mygrep])
  colnames(X)[ncol(X)] = "gender_diff"
  
  ## Delta exposure
  for (k in 1:ncol(expo)){
    X=cbind(X, as.numeric(dist(expo[,k]))[mygrep])
    colnames(X)[ncol(X)] = colnames(expo)[k]
  }
  
  # Save covariates and exposures
  saveRDS(X, paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))
  
}
