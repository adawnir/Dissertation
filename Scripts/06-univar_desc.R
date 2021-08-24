## Univariate regression: Detection vs Family ID
## Rin 29 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)

options(warn = 0)

### Detection ~ Family ID ----
# Initialise
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
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))

  # Convert matrix into binary
  expo = ifelse(expo=="nd",1,0)
  expo = expo[rownames(covars),]
  print(all(rownames(expo)==rownames(covars)))
  
  # Exclude chemicals with 0% non-detects
  expo = expo[,which(colSums(expo,na.rm = TRUE)!=0)]
  saveRDS(expo, paste0("../Processed/",filepaths[i],"/Exposure_matrix_nd_thresh_no_isolated.rds"))
  
  ## Detection vs family ID by compound
  # Luxembourg
  pvals = NULL
  f1='expo[,k] ~ covars$Family.ID'
  f0='expo[,k] ~ 1'
  t0=Sys.time()
  for (k in 1:ncol(expo)){
    model1=glm(as.formula(f1), family = "binomial")
    model0=glm(as.formula(f0), family = "binomial")
    pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
  }
  t1=Sys.time()
  print(t1-t0)
  names(pvals) = colnames(expo)
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_detection_family_pvals.rds"))
  assign(paste0("pvals_",suffix[i]),pvals)
}

### Check significant pvals for each data set ----
for (i in 1:length(batches)){
  pvals = eval(parse(text = paste0("pvals_",suffix[i])))
  bonf = 0.05/length(pvals)
  
  list = names(pvals)[pvals < bonf]
  assign(paste0("list_",suffix[i]),list)
}
list_lux
list_fra
list_gs
list_pooled2
list_pooled3

### Plotting ----
expo = readRDS(paste0("../Processed/",filepaths[4],"/Exposure_matrix_nd_thresh_no_isolated.rds"))

nd = round(colMeans(expo, na.rm = T),2)

values = -log10(pvals_pooled3)
annot_sub = annot[names(values)]

# Scatter plot
bonf = 0.05/length(pvals_pooled3)
bonf = -log10(bonf)

ifelse(dir.exists("../Figures/Section1"), "", dir.create("../Figures/Section1"))
{pdf("../Figures/Section1/Univariate_detection_family_pooled3.pdf", width=7, height=5)
  par(mar=c(5,5,1,10))
  plot(nd, values,
       col=annot.colours[annot_sub], cex = 3,
       xlab="Proportion of non-detects", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="p", pch = 19)
  abline(h = -log10(0.05), lty = 2, col = "grey")
  text(x = 0, y = -log10(0.05), labels="Nominal threshold", adj = c(0,1),col = "grey", cex = 0.8)
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  text(x = 0, y = bonf,labels="Bonferroni threshold", adj = c(0,1),col = batch.colours[4], cex = 0.8)
  abline(h = bonf, lty = 2, col = darken(batch.colours[4], 0.5))
  for (k in 1:length(values)){
    text(nd[k],values[k],
         labels = ifelse(values[k] > bonf, which(names(values)[k] == names(values)[values > bonf]), ""),
         cex = 0.8)
  }
  par(xpd = TRUE)
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(values>bonf),". ", names(values)[values>bonf]),
         text.col = darken(annot.colours[annot_sub[values>bonf]],0.5),
         bty = "n",  x.intersp = 0)
  dev.off()
}

for (i in c(1:3,5)){
  pvals = eval(parse(text = paste0("pvals_",suffix[i])))
  
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_nd_thresh_no_isolated.rds"))
  
  nd = round(colMeans(expo, na.rm = T),2)
  
  values = -log10(pvals)
  annot_sub = annot[names(values)]
  
  # Scatter plot
  bonf = 0.05/length(pvals)
  bonf = -log10(bonf)
  
  xseq = seq(1, length(values))
  
  {pdf(paste0("../Figures/Supplementary/Section1/Univariate_detection_family_",suffix[i],".pdf"), width=7, height=5)
    par(mar=c(5,5,1,10))
    plot(nd, values,
         col=annot.colours[annot_sub], cex = 3,
         xlab="Proportion of non-detects", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
         type="p", pch = 19)
    abline(h = -log10(0.05), lty = 2, col = "grey")
    text(x = 0, y = -log10(0.05), labels="Nominal threshold", adj = c(0,1),col = "grey", cex = 0.8)
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    text(x = 0, y = bonf,labels="Bonferroni threshold", adj = c(0,1),col = batch.colours[i], cex = 0.8)
    abline(h = bonf, lty = 2, col = darken(batch.colours[i], 0.5))
    for (k in 1:length(values)){
      text(nd[k],values[k],
           labels = ifelse(values[k] > bonf, which(names(values)[k] == names(values)[values > bonf]), ""),
           cex = 0.8)
    }
    par(xpd = TRUE)
    coord = par("usr")
    legend(x = coord[2]*1.05, y = coord[4],
           legend = paste0(1:sum(values>bonf),". ", names(values)[values>bonf]),
           text.col = darken(annot.colours[annot_sub[values>bonf]],0.5),
           bty = "n",  x.intersp = 0)
    dev.off()
  }
}


### Detection ~ Sample length/weight ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

# expo = readRDS(paste0("../Processed/",filepaths[1],"/Exposure_matrix_raw_thresh.rds"))
# covars = readRDS(paste0("../Processed/",filepaths[1],"/Participant_covariate_info_thresh.rds"))
# 
# # Convert matrix into binary
# expo = ifelse(expo=="nd",1,0)
# expo = expo[rownames(covars),]
# print(all(rownames(expo)==rownames(covars)))
# 
# # Exclude chemicals with 0% non-detects
# expo = expo[,which(colSums(expo,na.rm = TRUE)!=0)]
# saveRDS(expo, paste0("../Processed/",filepaths[1],"/Exposure_matrix_nd_thresh.rds"))

# Load data
expo = readRDS(paste0("../Processed/",filepaths[1],"/Exposure_matrix_nd_thresh.rds"))
covars = readRDS(paste0("../Processed/",filepaths[1],"/Participant_covariate_info_thresh.rds"))

## Detection vs length by compound
# Luxembourg
pvals = NULL
f1='expo[,k] ~ covars$Length'
t0=Sys.time()
for (k in 1:ncol(expo)){
  model1=glm(as.formula(f1), family = "binomial")
  pvals=c(pvals,summary(model1)$coefficients[2,4])
}
t1=Sys.time()
print(t1-t0)
names(pvals) = colnames(expo)
saveRDS(pvals, paste0("../Results/",filepaths[1],"/Univariate_detection_length_pvals.rds"))

bonf = 0.05/length(pvals)
names(pvals)[pvals < bonf]

nd = round(colMeans(expo, na.rm = T),2)

values = -log10(pvals)
annot_sub = annot[names(values)]

bonf = -log10(bonf)

xseq = seq(1, length(values))
{pdf("../Figures/Section1/Univariate_detection_length_lux.pdf", width=5, height=5)
  par(mar=c(5,5,1,1))
  plot(nd, values,
       col=annot.colours[annot_sub], cex = 3,
       xlab="Proportion of non-detects", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="p", pch = 19, ylim = c(min(values),bonf))
  abline(h = -log10(0.05), lty = 2, col = "grey")
  text(x = 0, y = -log10(0.05), labels="Nominal threshold", adj = c(0,1),col = "grey", cex = 0.8)
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  text(x = 0, y = bonf,labels="Bonferroni threshold", adj = c(0,1),col = batch.colours[1], cex = 0.8)
  abline(h = bonf, lty = 2, col = darken(batch.colours[1], 0.5))
  # text(nd,values,
  #      labels = ifelse(values > bonf, which(names(values)[values > bonf] %in% names(values)), ""),
  #      cex = 0.8)
  # par(xpd = TRUE)
  # coord = par("usr")
  # legend(x = coord[2]*1.05, y = coord[4],
  #        legend = paste0(1:sum(values>bonf),". ", names(values)[values>bonf]),
  #        text.col = darken(annot.colours[annot_sub[values>bonf]],0.5),
  #        bty = "n",  x.intersp = 0)
  dev.off()
}

## Detection vs weight by compound
# Luxembourg
pvals = NULL
f1='expo[,k] ~ covars$Weight'
t0=Sys.time()
for (k in 1:ncol(expo)){
  model1=glm(as.formula(f1), family = "binomial")
  pvals=c(pvals,summary(model1)$coefficients[2,4])
}
t1=Sys.time()
print(t1-t0)
names(pvals) = colnames(expo)
saveRDS(pvals, paste0("../Results/",filepaths[1],"/Univariate_detection_weight_pvals.rds"))

bonf = 0.05/length(pvals)
names(pvals)[pvals < bonf]

nd = round(colMeans(expo, na.rm = T),2)

values = -log10(pvals)
annot_sub = annot[names(values)]

bonf = -log10(bonf)

xseq = seq(1, length(values))
{pdf("../Figures/Section1/Univariate_detection_weight_lux.pdf", width=5, height=5)
  par(mar=c(5,5,1,1))
  plot(nd, values,
       col=annot.colours[annot_sub], cex = 3,
       xlab="Proportion of non-detects", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="p", pch = 19, ylim = c(min(values),bonf))
  abline(h = -log10(0.05), lty = 2, col = "grey")
  text(x = 0, y = -log10(0.05), labels="Nominal threshold", adj = c(0,1),col = "grey", cex = 0.8)
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  text(x = 0, y = bonf,labels="Bonferroni threshold", adj = c(0,1),col = batch.colours[1], cex = 0.8)
  abline(h = bonf, lty = 2, col = darken(batch.colours[1], 0.5))
  # text(nd,values,
  #      labels = ifelse(values > bonf, which(names(values)[values > bonf] %in% names(values)), ""),
  #      cex = 0.8)
  # par(xpd = TRUE)
  # coord = par("usr")
  # legend(x = coord[2]*1.05, y = coord[4],
  #        legend = paste0(1:sum(values>bonf),". ", names(values)[values>bonf]),
  #        text.col = darken(annot.colours[annot_sub[values>bonf]],0.5),
  #        bty = "n",  x.intersp = 0)
  dev.off()
}


## Univariate regression modelling association with exposure
## Rin 29 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)

### Exposure ~ Covariates ----
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
  X = covars %>% select(Age, Gender)
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
  f0='expo[,k] ~ 1'
  t0=Sys.time()
  for (j in 1:ncol(X)){
    for (k in 1:ncol(expo)){
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
  
  ifelse(dir.exists("../Results"),"",dir.create("../Results"))
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas.rds"))
}
### Check significant pvals for each data set ----
for (i in 1:length(batches)){
  pvals = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals.rds"))
  bonf = 0.05/nrow(pvals)
  
  list = apply(pvals, 2, function(x) x < bonf)
  assign(paste0("list_",suffix[i]),list)
}
colSums(list_lux)
colSums(list_fra)
colSums(list_gs)
colSums(list_pooled2)
colSums(list_pooled3)
rownames(list_lux)[list_lux[,"Age"]]
rownames(list_fra)[list_fra[,"Age"]]
rownames(list_gs)[list_gs[,"Age"]]
rownames(list_pooled2)[list_pooled2[,"Age"]]
rownames(list_pooled3)[list_pooled3[,"Age"]]

rownames(list_lux)[list_lux[,"Gender"]]
rownames(list_fra)[list_fra[,"Gender"]]
rownames(list_gs)[list_gs[,"Gender"]]
rownames(list_pooled2)[list_pooled2[,"Gender"]]
rownames(list_pooled3)[list_pooled3[,"Gender"]]

### Volcano plots ----
# Age
for (i in 1:length(batches)){
  pvals = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas.rds"))
  annot_sub = annot[rownames(betas)]
  if(i ==4){
    pdf(paste0("../Figures/Section1/Univariate_exposure_age_",suffix[i],".pdf"), width = 5, height = 5)
  } else{pdf(paste0("../Figures/Supplementary/Section1/Univariate_exposure_age_",suffix[i],".pdf"), width = 5, height = 5)
  }
  par(mar=c(5,5,1,1))
  plot(betas[,"Age"], -log10(pvals[,"Age"]), pch=19,
       col=ifelse(pvals[,"Age"] < 0.05/length(betas[,"Age"]), annot.colours[annot_sub], "grey"),
       cex.lab=1.5, cex = 1.5,
       ylim = c(0, max(-log10(pvals[,"Age"]))+0.25),
       xlim = c(-max(abs(betas[,"Age"]))-0.2, max(abs(betas[,"Age"]))+0.2),
       ylab=expression(-log[10](italic(p))),
       xlab=expression(beta))
  text(betas[,"Age"]+sign(betas[,"Age"])*0.05, -log10(pvals[,"Age"])+0.1,
       labels = ifelse(pvals[,"Age"] < 0.05/length(betas[,"Age"]), rownames(betas), ""),
       col = annot.colours[annot_sub], cex = 0.8)
  abline(h = -log10(0.05/length(betas[,"Age"])), lty = 2, col = batch.colours[i])
  text(-max(abs(betas[,"Age"]))-0.2, -log10(0.05/length(betas[,"Age"])),  labels = "Bonferroni threshold", adj=c(0,1), col = batch.colours[i], cex = 0.8)
  dev.off()
}

### Gender
for (i in 1:length(batches)){
  pvals = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas.rds"))
  annot_sub = annot[rownames(betas)]
  ## Volcano plots
  if(i == 4){
    pdf(paste0("../Figures/Section1/Univariate_exposure_gender_",suffix[i],".pdf"), width = 5, height = 5)
  } else{pdf(paste0("../Figures/Supplementary/Section1/Univariate_exposure_gender_",suffix[i],".pdf"), width = 5, height = 5)
  }
  par(mar=c(5,5,1,1))
  plot(betas[,"Gender"], -log10(pvals[,"Gender"]), pch=19,
       col=ifelse(pvals[,"Gender"] < 0.05/length(betas[,"Gender"]), annot.colours[annot_sub], "grey"),
       cex.lab=1.5, cex = 1.5,
       ylim = c(0, max(-log10(pvals[,"Gender"]))+0.25),
       xlim = c(-max(abs(betas[,"Gender"]))-0.2, max(abs(betas[,"Gender"]))+0.2),
       ylab=expression(-log[10](italic(p))),
       xlab=expression(beta))
  text(betas[,"Gender"]+sign(betas[,"Gender"])*0.05, -log10(pvals[,"Gender"])+0.3,
       labels = ifelse(pvals[,"Gender"] < 0.05/length(betas[,"Gender"]), rownames(betas), ""),
       col = annot.colours[annot_sub], cex = 0.8)
  abline(h = -log10(0.05/length(betas[,"Gender"])), lty = 2, col = batch.colours[i])
  text(-max(abs(betas[,"Gender"]))-0.2, -log10(0.05/length(betas[,"Gender"])), labels = "Bonferroni threshold", adj=c(0,1), col = batch.colours[i], cex = 0.8)
  dev.off()
}
