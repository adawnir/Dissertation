## Univariate regression (including isolated children)
## Rin 12 July

# Load packages
library(tidyverse)
library(RColorBrewer)

### Exposure ~ Covariate ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  if (i==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else {X = covars %>% select(Age, Gender)}
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
  t0=Sys.time()
  for (k in 1:ncol(expo)){
    for (j in 1:ncol(X)){
      model1=lm(as.formula(f1))
      betas=c(betas, coefficients(model1)[2])
      pvals=c(pvals, summary(model1)$coefficients[2,4])
    }
  }
  t1=Sys.time()
  print(t1-t0)
  betas = matrix(betas, nrow = ncol(expo), ncol = ncol(X))
  pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
  rownames(pvals)=rownames(betas)=colnames(expo)
  colnames(pvals)=colnames(betas)=colnames(X)
  
  annot_sub = annot[rownames(betas)]
  
  ## Volcano plot
  pdf(paste0("../Figures/",filepaths[i],"/Univariate_exposure_covariate.pdf"), width = 10, height = 5*ncol(betas)/2)
  par(mar=c(5,5,1,1), mfrow = c(ncol(betas)/2,2))
  for (n in 1:ncol(betas)){
    plot(betas[,n], -log10(pvals[,n]), pch=19,
         col=ifelse(pvals[,n] < 0.05/length(betas[,n]), annot.colours[annot_sub], "grey"),
         cex.lab=1, cex = 0.7, xlim = c(-max(abs(betas[,n])), max(abs(betas[,n]))),
         ylab=expression(-log[10](p-value)), 
         xlab=substitute(paste(beta,"(",tmp,")"), list(tmp = colnames(betas)[n])))
    text(betas[,n]+sign(betas[,n])*0.15, -log10(pvals[,n]),
         labels = ifelse(pvals[,n] < 0.05/length(betas[,n]), rownames(betas), ""),
         col = annot.colours[annot_sub])
    abline(h = -log10(0.05/length(betas[,n])), lty = 2)
  }
  dev.off()
  
  ifelse(dir.exists("../Results"),"",dir.create("../Results"))
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas.rds"))
}


### Detection ~ Covariate ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_raw.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  if (i==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else {X = covars %>% select(Age, Gender)}
  
  # Convert matrix into binary
  expo = ifelse(expo=="nd",1,0)
  # Exclude chemicals with 0% non-detects
  expo = expo[,which(colSums(expo,na.rm = TRUE)!=0)]
  saveRDS(expo, paste0("../Processed/",filepaths[i],"/Exposure_matrix_nd.rds"))
  
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
  t0=Sys.time()
  for (k in 1:ncol(expo)){
    for (j in 1:ncol(X)){
      model1=glm(as.formula(f1), family = "binomial")
      betas=c(betas, coefficients(model1)[2])
      pvals=c(pvals, summary(model1)$coefficients[2,4])
    }
  }
  t1=Sys.time()
  print(t1-t0)
  betas = matrix(betas, nrow = ncol(expo), ncol = ncol(X))
  pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
  rownames(pvals)=rownames(betas)=colnames(expo)
  colnames(pvals)=colnames(betas)=colnames(X)
  
  annot_sub = annot[rownames(betas)]
  
  ## Volcano plot
  pdf(paste0("../Figures/",filepaths[i],"/Univariate_detection_covariate.pdf"), width = 10, height = 5*ncol(betas)/2)
  par(mar=c(5,5,1,1), mfrow = c(ncol(betas)/2,2))
  for (n in 1:ncol(betas)){
    plot(betas[,n], -log10(pvals[,n]), pch=19,
         col=ifelse(pvals[,n] < 0.05/length(betas[,n]), annot.colours[annot_sub], "grey"),
         cex.lab=1, cex = 0.7, xlim = c(-max(abs(betas[,n])), max(abs(betas[,n]))),
         ylab=expression(-log[10](p-value)), 
         xlab=substitute(paste(beta,"(",tmp,")"), list(tmp = colnames(betas)[n])))
    text(betas[,n]+sign(betas[,n])*5.5, -log10(pvals[,n]),
         labels = ifelse(pvals[,n] < 0.05/length(betas[,n]), rownames(betas), ""),
         col = annot.colours[annot_sub])
    abline(h = -log10(0.05/length(betas[,n])), lty = 2)
  }
  dev.off()
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_detection_covariate_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_detection_covariate_betas.rds"))
  # assign(paste0("betas_",suffix[i]),betas)
  # assign(paste0("pvals_",suffix[i]),pvals)
}
