## Univariate regression modelling association with exposure
## Rin 29 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(sgPLS)
library(focus)

### Exposure ~ Covariates ----
# Initialisation
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in c(2,4)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  ### Region (Ile-de-France Y/N) ----
  Y = ifelse(covars$Region=="Île-de-France",1,0)
  X = expo
  
  stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS, family = "binomial")
  
  saveRDS(stab, paste0("../Results/",filepaths[i],"/Multivariate_exposure_region_output.rds"))
  
  pdf(paste0("../Figures/",filepaths[i],"/Multivariate_exposure_region_output.pdf"))
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-((length(a)-1):length(a))]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-((length(b)-1):length(b))]
  pdf(paste0("../Figures/",filepaths[i],"/Multivariate_exposure_region_beta_consistency.pdf"))
  par(mar=c(5,5,1,1))
  plot(a/(a+b), selprop, las=1, pch=19, col="navy", cex.lab=1.5,
       xlab="Proportion of positive beta among non-zero betas", 
       ylab="Selection Proportion") # Ideally no point around 0.5 with high selection proportion
  dev.off()
  
  # Extract beta
  beta = stab$Beta[ArgmaxId(stab)[1],,]
  beta_mu = rowMeans(beta[-nrow(beta),] %*% diag(beta[nrow(beta),]))
  names(beta_mu) = names(selprop)
  
  selprop_ranked <- sort(selprop, decreasing = TRUE)
  beta_mu_ranked <- beta_mu[order(-selprop)]
  print(all(names(selprop_ranked)==names(beta_mu_ranked)))
  
  pdf(paste0("../Figures/",filepaths[i],"/Multivariate_exposure_region_selprop.pdf"), width = 14)
  par(mar=c(10, 5, 1, 1))
  plot(selprop_ranked,
       type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
       col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[i], no = "grey"),
       xaxt = "n", xlab = "", ylab = "Selection proportions"
  )
  abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
  axis(side = 1, at = 1:length(selprop_ranked), las = 2, labels = names(selprop_ranked))
  dev.off()
}

