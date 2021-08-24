## Multivariate clustering - Association with shortest path
## 7 Aug

# Load packages
library(focus)
library(colorspace)
library(RColorBrewer)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# ## Parameters
# args=commandArgs(trailingOnly=TRUE)
# m=as.numeric(args[1])

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

for (m in 1:length(batches)){
  X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  fsp = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))
  
  betas = pvals = NULL
  f1='Y ~ X[,k]'
  t0=Sys.time()
  for (k in 1:ncol(X)){
    Y = fsp
    model1=lm(as.formula(f1))
    if(colnames(X)[k] == "gender_diff"){
      betas=c(betas, coefficients(model1)[2:length(coefficients(model1))])
      pvals=c(pvals, summary(model1)$coefficients[2:nrow(summary(model1)$coefficients),4])
    } else {
      betas=c(betas, coefficients(model1)[2])
      pvals=c(pvals, summary(model1)$coefficients[2,4])
    }
  }
  t1=Sys.time()
  print(t1-t0)
  
  head(colnames(X),10)
  
  mylabels = c(colnames(X)[1:which(colnames(X)=="gender_diff")-1],
               "All male (ref. All female)","Different gender (ref. All female)",
               colnames(X)[(which(colnames(X)=="gender_diff")+1):ncol(X)])
  
  length(betas)==length(mylabels)
  
  mylabels = gsub("age_mu","Mean age",mylabels)
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("region","Île-de-France (Y/N)",mylabels)
  
  names(pvals)=names(betas)=mylabels
  
  saveRDS(pvals, paste0("../Results/",filepaths[m],"/Shortest_path_univar_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[m],"/Shortest_path_univar_betas.rds"))
  
  ### Multivariate analysis using stability selection sPLS regression ----
  Y = fsp[complete.cases(X)]
  X = X[complete.cases(X),]
  X = model.matrix(~., X)[,-1]
  
  X = scale(X)
  
  stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)
  
  saveRDS(stab, paste0("../Results/",filepaths[m],"/Shortest_path_multivar_output.rds"))
  
}


for (m in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  expo  = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  fsp = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  families=unique(covars$Family.ID)
  
  Y = rep(fsp,each = 2)

  if (m==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else if(m %in% c(2,4)) {
    X = covars %>% select(Age, Gender, Region) %>%
      mutate(Region = ifelse(Region=="Île-de-France","Yes","No"))
  } else {X = covars %>% select(Age, Gender)}
  
  X = cbind(X, expo)
  
  ### Univariate analysis ----
  betas = pvals = NULL
  f1='Y ~ X[,k]'
  t0=Sys.time()
  for (k in 1:ncol(X)){
    model1=lm(as.formula(f1))
    betas=c(betas, coefficients(model1)[2])
    pvals=c(pvals, summary(model1)$coefficients[2,4])
  }
  t1=Sys.time()
  print(t1-t0)
  
  mylabels = colnames(X)
  
  mylabels=gsub("Region","Île-de-France (Y/N)", mylabels)
  
  length(betas)==length(mylabels)
  
  names(pvals)=names(betas)=mylabels
  
  saveRDS(pvals, paste0("../Results/",filepaths[m],"/Shortest_path_individual_univar_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[m],"/Shortest_path_individual_univar_betas.rds"))
  
  ### Multivariate analysis using stability selection sPLS regression ----
  
  Y = Y[complete.cases(X)]
  X = X[complete.cases(X),]
  X = model.matrix(~., X)[,-1]
  
  X = scale(X)
  
  stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)
  
  saveRDS(stab, paste0("../Results/",filepaths[m],"/Shortest_path_individual_multivar_output.rds"))
}
