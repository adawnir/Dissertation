## Misclassified analysis
## Rin Wada 15 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)

### Stability-based ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

for (i in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Results/",filepaths[i],"/Cluster_memberships.rds"))
  expo  = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  families=unique(covars$Family.ID)
  
  misclass=rep(NA, length(families))
  for (f in 1:length(families)){
    tmp=covars$stab_cluster[covars$Family.ID==families[f]]
    if (length(unique(tmp))==1){
      misclass[f] = 0
    } else {
      misclass[f] = 1
      }
  }
  covars$stab_outcome = rep(misclass, as.vector(table(covars$Family.ID)))
  
  ncol(expo)
  betas = pvals = NULL
  f1='covars$stab_outcome ~ expo[,k]'
  t0=Sys.time()
  for (k in 1:ncol(expo)){
    model1=lm(as.formula(f1))
    betas=c(betas, coefficients(model1)[2])
    pvals=c(pvals, summary(model1)$coefficients[2,4])
  }
  t1=Sys.time()
  print(t1-t0)
  names(pvals)=names(betas)=colnames(expo)
  
  annot_sub = annot[names(betas)]
  
  ## Volcano plot
  pdf(paste0("../Figures/",filepaths[i],"/Univariate_stab_misclassified_exposure.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(betas, -log10(pvals), pch=19,
         col=ifelse(pvals < 0.05/length(betas), annot.colours[annot_sub], "grey"),
         cex.lab=1, cex = 0.7,
         ylim = c(0, max(-log10(pvals))+0.25),
         xlim = c(-max(abs(betas))-0.4, max(abs(betas))+0.4),
         ylab=expression(-log[10](p)), 
         xlab=expression(beta))
  text(betas+sign(betas)*0.1, -log10(pvals)+0.25,
       labels = ifelse(pvals < 0.05/length(betas), names(betas), ""),
       col = annot.colours[annot_sub])
  abline(h = -log10(0.05/length(betas)), lty = 2)
  dev.off()
  
  ## Delta exposure
  delta_mat=matrix(NA, nrow=length(families), ncol=ncol(expo))
  for (f in 1:length(families)){
    for (k in 1:ncol(expo)){
      tmp = expo[covars$Family.ID==families[f],k]
      delta_mat[f,k] = abs(tmp[1] - tmp [2])
    }
  }
  rownames(delta_mat)=families
  colnames(delta_mat)=colnames(expo)
  
  betas = pvals = NULL
  f1='misclass ~ delta_mat[,k]'
  t0=Sys.time()
  for (k in 1:ncol(delta_mat)){
    model1=lm(as.formula(f1))
    betas=c(betas, coefficients(model1)[2])
    pvals=c(pvals, summary(model1)$coefficients[2,4])
  }
  t1=Sys.time()
  print(t1-t0)
  names(pvals)=names(betas)=colnames(delta_mat)
  
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_vs_delta_expo.pdf"))
    par(mar=c(5,5,1,1))
    plot(betas, -log10(pvals), pch=19,
         col=ifelse(pvals < 0.05/length(betas), annot.colours[annot_sub], "grey"),
         cex.lab=1, cex = 0.7,
         ylim = c(0, max(-log10(pvals))+0.25),
         xlim = c(-max(abs(betas))-0.4, max(abs(betas))+0.4),
         ylab=expression(-log[10](p)), 
         xlab=expression(beta))
    text(betas+sign(betas)*0.1, -log10(pvals)+0.25,
         labels = ifelse(pvals < 0.05/length(betas), names(betas), ""),
         col = annot.colours[annot_sub])
    abline(h = -log10(0.05/length(betas)), lty = 2)
    dev.off()
  }
  
  ## Relationship with age
  age_mu=rep(NA, length(families))
  for (f in 1:length(families)){
    age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
  }
  names(age_mu)=families
  
  model=glm(misclass~age_mu, family = "binomial")
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_vs_age_mean.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_mu, misclass, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(age_mu)],
         xlab="Within-family age mean", ylab="Misclassified", ylim = c(0,1.05))
    text(age_mu, misclass, labels=names(age_mu), pos=3, col=darken(family.colours[names(age_mu)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("right", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
    dev.off()
  }
  
  age_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    age_diff[f]=abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2])
  }
  names(age_diff)=families
  
  model=glm(misclass~age_diff, family = "binomial")
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_vs_age_diff.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_diff, misclass, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(age_diff)],
         xlab="Within-family age absolute difference", ylab="Misclassified", ylim = c(0,1.05))
    text(age_diff, misclass, labels=names(age_diff), pos=3, col=darken(family.colours[names(age_diff)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("right", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
    dev.off()
  }
  
  if (i==1){
    weight_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_mu[f]=mean(covars$Weight[covars$Family.ID==families[f]])
    }
    names(weight_mu)=families
    
    model=glm(misclass~weight_mu, family = "binomial")
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_vs_weight_mean.pdf"))
      par(mar=c(5,5,1,1))
      plot(weight_mu, misclass, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(weight_mu)],
           xlab="Within-family weight mean", ylab="Misclassified", ylim = c(0,1.05))
      text(weight_mu, misclass, labels=names(weight_mu), pos=3, col=darken(family.colours[names(weight_mu)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("left", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    weight_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_diff[f]=abs(covars$Weight[covars$Family.ID==families[f]][1]-covars$Weight[covars$Family.ID==families[f]][2])
    }
    names(weight_diff)=families
    
    model=glm(misclass~weight_diff, family = "binomial")
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_vs_weight_diff.pdf"))
      par(mar=c(5,5,1,1))
      plot(weight_diff, misclass, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(weight_diff)],
           xlab="Within-family weight absolute difference", ylab="Misclassified", ylim = c(0,1.05))
      text(weight_diff, misclass, labels=names(weight_diff), pos=3, col=darken(family.colours[names(weight_diff)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("right", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    length_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      length_mu[f]=mean(covars$Length[covars$Family.ID==families[f]])
    }
    names(length_mu)=families
    
    model=glm(misclass~length_mu, family = "binomial")
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_vs_length_mean.pdf"))
      par(mar=c(5,5,1,1))
      plot(length_mu, misclass, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(length_mu)],
           xlab="Within-family length mean", ylab="Misclassified", ylim = c(0,1.05))
      text(length_mu, misclass, labels=names(length_mu), pos=3, col=darken(family.colours[names(length_mu)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("right", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    length_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      length_diff[f]=abs(covars$Length[covars$Family.ID==families[f]][1]-covars$Length[covars$Family.ID==families[f]][2])
    }
    names(length_diff)=families
    
    model=glm(misclass~length_diff, family = "binomial")
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_vs_length_diff.pdf"))
      par(mar=c(5,5,1,1))
      plot(length_diff, misclass, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(length_diff)],
           xlab="Within-family length absolute difference", ylab="Misclassified", ylim = c(0,1.05))
      text(length_diff, misclass, labels=names(length_diff), pos=3, col=darken(family.colours[names(length_diff)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("right", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
  }
  
  gender_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    tmp=covars$Gender[covars$Family.ID==families[f]]
    if (length(unique(tmp))==1){
      if (unique(tmp)=="Male"){
        gender_diff[f]=1
      }
      if (unique(tmp)=="Female"){
        gender_diff[f]=2
      }
    } else {
      gender_diff[f]=3
    } 
  }
  names(gender_diff)=families
  
  model1 = glm(misclass ~ as.factor(gender_diff), family = "binomial")
  model0 = glm(misclass ~ 1, family = "binomial")
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_by_gender.pdf"))
    par(mar=c(5,5,1,1))
    plot(misclass, pch=19, cex=1, cex.lab=1.5, las=1, 
         col=c("skyblue", "pink", "tan")[gender_diff],
         xlab="Family ID", ylab="Misclassified", ylim = c(0,1.05))
    text(misclass, labels=names(gender_diff), pos=3, col=darken(family.colours[names(gender_diff)], amount=0.5))
    legend("left", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
           legend=c("All Male", "All Female", "Different genders"))
    legend("right", bty="n", cex=1.5,
           legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
    dev.off()
  }
  if (i %in% 4:5){
    Batch=rep(NA, length(families))
    for (f in 1:length(families)){
      tmp=covars$Batch[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Batch[f] = unique(tmp)
      }
    }
    names(Batch)=families
    
    model1 = glm(misclass ~ as.factor(Batch), family = "binomial")
    model0 = glm(misclass ~ 1, family = "binomial")
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_by_Batch.pdf"))
      par(mar=c(5,5,1,1))
      plot(misclass, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=batch.colours[Batch],
           xlab="Family ID", ylab="Misclassified", ylim = c(0,1.05))
      text(misclass, labels=names(Batch), pos=3, col=darken(family.colours[names(Batch)], amount=0.5))
      legend("left", pch=19, col=batch.colours[unique(Batch)],
             legend=levels(covars$Batch))
      legend("right", bty="n", cex=1.5,
             legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
      dev.off()
    }
  }
  if (i == 4){
    Region=rep(NA, length(families))
    Department=rep(NA, length(families))
    for (f in 1:length(families)){
      tmp=covars$Region[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Region[f] = unique(tmp)
      }
      tmp=covars$Department[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Department[f] = unique(tmp)
      }
    }
    names(Region)=families
    names(Department)=families
    
    region.colours=brewer.pal(n=12,name='Paired')
    region.colours=colorRampPalette(region.colours)(length(levels(covars$Region)))
    names(region.colours)=levels(covars$Region)
    
    depart.colours=brewer.pal(n=12,name='Paired')
    depart.colours=colorRampPalette(depart.colours)(length(levels(covars$Department)))
    names(depart.colours)=levels(covars$Department)
    
    model1 = glm(misclass ~ as.factor(Region), family = "binomial")
    model0 = glm(misclass ~ 1, family = "binomial")
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_by_Region.pdf"))
      par(mar=c(5,5,1,1))
      plot(misclass, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=region.colours[Region],
           xlab="Family ID", ylab="Misclassified", ylim = c(0,1.05))
      text(misclass, labels=names(Region), pos=3, col=darken(family.colours[names(Region)], amount=0.5))
      legend("left", pch=19, col=region.colours[unique(Region)],
             legend=levels(covars$Region), cex = 0.6, ncol = 2)
      legend("right", bty="n", cex=1.5,
             legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
      dev.off()
    }
    
    model1 = glm(misclass ~ as.factor(Department), family = "binomial")
    model0 = glm(misclass ~ 1, family = "binomial")
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_misclassified_by_Department.pdf"))
      par(mar=c(5,5,1,1))
      plot(misclass, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=depart.colours[Department],
           xlab="Family ID", ylab="Misclassified", ylim = c(0,1.05))
      text(misclass, labels=names(Department), pos=3, col=darken(family.colours[names(Department)], amount=0.5))
      legend("left", pch=19, col=depart.colours[unique(Department)],
             legend=levels(covars$Department), cex = 0.6, ncol = 2)
      legend("right", bty="n", cex=1.5,
             legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
      dev.off()
    }
  }
}
