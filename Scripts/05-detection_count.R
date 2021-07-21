## Detection count Poisson regression (including isolated children)
## Rin 12 July

# Load packages
library(tidyverse)
library(RColorBrewer)

### Detection frequency per participant ----
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
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_raw_thresh.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  # Count non-detects for each participant
  covars$nd_count = apply(expo, 1, function(x) sum(x=="nd", na.rm = T))
  # Association with gender
  mean(covars$nd_count)
  var(covars$nd_count)
  summary(covars$nd_count)
  if (i==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else {X = covars %>% select(Age, Gender)}
  pdf(paste0("../Figures/",filepaths[i],"/Detection_count_covariate.pdf"), width=10, height=5*ceiling(ncol(X)/2))
  par(mar=c(5,5,1,1), mfrow = c(ceiling(ncol(X)/2),2))
  for (n in 1:ncol(X)){
    model = glm(nd_count ~ X[,n], poisson, covars)
    if (is.factor(X[,n])){
      boxplot(nd_count ~ X[,n], data = covars, col = batch.colours[i],
              xlab = colnames(X)[n], ylab = "Number of non-detects", main = NULL,
              cex.lab = 1.5) 
    } else {
      plot(nd_count ~ X[,n], data = covars, pch = 19, col = batch.colours[i],
           xlab = colnames(X)[n], ylab = "Number of non-detects", main = NULL,
           cex.lab = 1.5)
    }
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
  }
  dev.off()
  
  if (i %in% 4:5){
    if (i==4){X = covars %>% select(Batch, Region, Department)}
    if (i==5){X = covars %>% select(Batch)}
    pdf(paste0("../Figures/",filepaths[i],"/Detection_count_geo.pdf"), width=5*(ceiling(ncol(X)/ceiling(ncol(X)/3))), height=5*ceiling(ncol(X)/3))
    if (i==4){par(mar=c(5,13,1,1), mfrow = c(ceiling(ncol(X)/3),ceiling(ncol(X)/ceiling(ncol(X)/3))))}
    if (i==5){par(mar=c(5,5,1,1), mfrow = c(ceiling(ncol(X)/3),ceiling(ncol(X)/ceiling(ncol(X)/3))))}
    for (n in 1:ncol(X)){
      model1=glm(nd_count ~ X[,n], poisson, covars)
      model0=glm(nd_count ~ 1, poisson, covars)
      pval = anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2]
      pval = ifelse(pval==0,
                    paste0("<",formatC(.Machine$double.xmin, format="e", digits=2)),
                    formatC(pval, format="e", digits=2))
      if (colnames(X)[n]=="Batch") {mycolours = batch.colours}
      if (colnames(X)[n]=="Region") {mycolours = region.colours}
      if (colnames(X)[n]=="Department") {mycolours = depart.colours}
      if (colnames(X)[n]=="Batch"){
        boxplot(nd_count ~ X[,n], data = covars, col = mycolours, las = 1,
                xlab = "Number of non-detects", ylab = colnames(X)[n], main = NULL,
                horizontal = TRUE) 
      } else{
        boxplot(nd_count ~ X[,n], data = covars, col = mycolours, las = 1,
                xlab = "Number of non-detects", ylab = "", main = NULL,
                cex.axis = 0.7, horizontal = TRUE) 
      }
      legend("topright", bty="n", cex=1.5,legend=paste0("p=", pval))
    }
    dev.off()
    }
}
