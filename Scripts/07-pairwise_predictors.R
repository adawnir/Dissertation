## Pairwise: covariate mean and difference/delta exposure
## 15 July

# Load packages
library(tidyverse)

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

for (m in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  expo  = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  ### Family-wise covariates and delta exposures ----
  families=levels(covars$Family.ID)
  
  X=matrix(NA, nrow = length(families))
  for (f in 1:length(families)){
    X[f]=mean(covars$Age[covars$Family.ID==families[f]])
  }
  colnames(X)[ncol(X)] = "age_mu"
  rownames(X) = families
  
  tmp = NULL
  for (f in 1:length(families)){
    tmp = c(tmp, abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2]))
  }
  X=cbind(as.data.frame(X), tmp)
  colnames(X)[ncol(X)] = "age_diff"
  
  if (m==1){
    tmp = NULL
    for (f in 1:length(families)){
      tmp = c(tmp, mean(covars$Weight[covars$Family.ID==families[f]]))
    }
    X=cbind(X, tmp)
    colnames(X)[ncol(X)] = "weight_mu"
    
    tmp = NULL
    for (f in 1:length(families)){
      
      tmp = c(tmp, abs(covars$Weight[covars$Family.ID==families[f]][1]-covars$Weight[covars$Family.ID==families[f]][2]))
    }
    X=cbind(X, tmp)
    colnames(X)[ncol(X)] = "weight_diff"
    
    tmp = NULL
    for (f in 1:length(families)){
      tmp = c(tmp, mean(covars$Length[covars$Family.ID==families[f]]))
    }
    X=cbind(X, tmp)
    colnames(X)[ncol(X)] = "length_mu"
    
    tmp = NULL
    for (f in 1:length(families)){
      
      tmp = c(tmp, abs(covars$Length[covars$Family.ID==families[f]][1]-covars$Length[covars$Family.ID==families[f]][2]))
    }
    X=cbind(X, tmp)
    colnames(X)[ncol(X)] = "length_diff"
  }
  
  tmp = NULL
  for (f in 1:length(families)){
    gender=covars$Gender[covars$Family.ID==families[f]]
    if (sum(is.na(gender))>0) {
      tmp[f] = NA
    } else if (length(unique(gender))==1){
      if (unique(gender)=="Male"){
        tmp[f]="All male"
      }
      if (unique(gender)=="Female"){
        tmp[f]="All female"
      }
    } else {
      tmp[f]="Different gender"
    }
  }
  X=cbind(X, tmp)
  colnames(X)[ncol(X)] = "gender_diff"
  
  if (m %in% c(2,4)){
    tmp = NULL
    for (f in 1:length(families)){
      region=as.character(covars$Region)[covars$Family.ID==families[f]]
      if (unique(region)=="Île-de-France"){
        tmp[f] = "Yes"
      } else {tmp[f] = "No"}
    }
    X=cbind(X, tmp)
    colnames(X)[ncol(X)] = "region"
  }
  
  ## Delta exposure
  for (k in 1:ncol(expo)){
    for (f in 1:length(families)){
      delta = expo[covars$Family.ID==families[f],k]
      tmp[f] = abs(delta[1] - delta [2])
    }
    X=cbind(X, as.numeric(tmp))
    colnames(X)[ncol(X)] = colnames(expo)[k]
  }
  
  # Save covariates and exposures
  saveRDS(X, paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  
  
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
  
  if (m %in% c(2,4,5)){
  region = combn(nrow(covars), 2, function(x) as.character(covars$Region)[x])
  tmp = apply(region, 2,
              function(x) ifelse(sum(is.na(x))>0, NA,
                                 ifelse(length(unique(x))>1,"No", "Yes")))
  X=cbind(X, tmp[mygrep])
  colnames(X)[ncol(X)] = "region"
  }
  
  ## Delta exposure
  for (k in 1:ncol(expo)){
    X=cbind(X, as.numeric(dist(expo[,k]))[mygrep])
    colnames(X)[ncol(X)] = colnames(expo)[k]
  }
  
  # Save covariates and exposures
  saveRDS(X, paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))
  
}
