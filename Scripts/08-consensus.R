## Stability-based clustering co-membership proportion
## Rin Wada 2 Aug

# Load packages
library(tidyverse)
library(RColorBrewer)
library(focus)
library(colorspace)

### Consensus probability association ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  out = readRDS(paste0("../Results/",filepaths[i],"/Stability_clustering_output.rds"))
  prop = SelectionProportions(out)
  
  families=levels(covars$Family.ID)
  
  fprop=NULL
  for (f in 1:length(families)){
    tmpmat=prop[covars$Family.ID==families[f],covars$Family.ID==families[f]]
    fprop=c(fprop,mean(tmpmat[upper.tri(tmpmat)]))
  }
  names(fprop)=families
  
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all.pdf"),width=14)
    par(mar=c(5,5,1,1))
    plot(fprop, pch=19, col=family.colours[families], xaxt="n", las=1,
         panel.first=abline(v=1:length(families), lty=3, col="grey"),
         xlab="Family ID", cex=2,
         ylab="Co-membership proportion between siblings", cex.lab=1.5)
    axis(side=1, at=1:length(families), labels=families, las = 2)
    dev.off()
  }
  
  ## Relationship with age
  age_mu=rep(NA, length(families))
  for (f in 1:length(families)){
    age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
  }
  names(age_mu)=families
  
  model=lm(fprop~age_mu)
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_vs_age_mean.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_mu, fprop, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(age_mu)],
         xlab="Within-family age mean", ylab="Co-membership proportion between siblings")
    text(age_mu, fprop, labels=names(age_mu), pos=3, col=darken(family.colours[names(age_mu)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
    dev.off()
  }
  
  age_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    age_diff[f]=abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2])
  }
  names(age_diff)=families
  
  model=lm(fprop~age_diff)
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_vs_age_diff.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_diff, fprop, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(age_diff)],
         xlab="Within-family age absolute difference", ylab="Co-membership proportion between siblings")
    text(age_diff, fprop, labels=names(age_diff), pos=3, col=darken(family.colours[names(age_diff)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
    dev.off()
  }
  
  if (i==1){
    weight_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_mu[f]=mean(covars$Weight[covars$Family.ID==families[f]])
    }
    names(weight_mu)=families
    
    model=lm(fprop~weight_mu)
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_vs_weight_mean.pdf"))
      par(mar=c(5,5,1,1))
      plot(weight_mu, fprop, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(weight_mu)],
           xlab="Within-family weight mean", ylab="Co-membership proportion between siblings")
      text(weight_mu, fprop, labels=names(weight_mu), pos=3, col=darken(family.colours[names(weight_mu)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topleft", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    weight_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_diff[f]=abs(covars$Weight[covars$Family.ID==families[f]][1]-covars$Weight[covars$Family.ID==families[f]][2])
    }
    names(weight_diff)=families
    
    model=lm(fprop~weight_diff)
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_vs_weight_diff.pdf"))
      par(mar=c(5,5,1,1))
      plot(weight_diff, fprop, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(weight_diff)],
           xlab="Within-family weight absolute difference", ylab="Co-membership proportion between siblings")
      text(weight_diff, fprop, labels=names(weight_diff), pos=3, col=darken(family.colours[names(weight_diff)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    length_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      length_mu[f]=mean(covars$Length[covars$Family.ID==families[f]])
    }
    names(length_mu)=families
    
    model=lm(fprop~length_mu)
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_vs_length_mean.pdf"))
      par(mar=c(5,5,1,1))
      plot(length_mu, fprop, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(length_mu)],
           xlab="Within-family length mean", ylab="Co-membership proportion between siblings")
      text(length_mu, fprop, labels=names(length_mu), pos=3, col=darken(family.colours[names(length_mu)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    length_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      length_diff[f]=abs(covars$Length[covars$Family.ID==families[f]][1]-covars$Length[covars$Family.ID==families[f]][2])
    }
    names(length_diff)=families
    
    model=lm(fprop~length_diff)
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_vs_length_diff.pdf"))
      par(mar=c(5,5,1,1))
      plot(length_diff, fprop, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(length_diff)],
           xlab="Within-family length absolute difference", ylab="Co-membership proportion between siblings")
      text(length_diff, fprop, labels=names(length_diff), pos=3, col=darken(family.colours[names(length_diff)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
  }
  
  gender_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    tmp=covars$Gender[covars$Family.ID==families[f]]
    if (sum(is.na(tmp))>0) {
      gender_diff[f] = NA
    } else if (length(unique(tmp))==1){
      if (unique(tmp)=="Male"){
        gender_diff[f]="All male"
      }
      if (unique(tmp)=="Female"){
        gender_diff[f]="All female"
      }
    } else {
      gender_diff[f]="Different genders"
    }
  }
  names(gender_diff)=families
  
  gender_col = c("skyblue", "pink", "tan")
  names(gender_col) = c("All male", "All female", "Different genders")
  
  model1 = lm(fprop[!is.na(gender_diff)] ~ gender_diff[!is.na(gender_diff)])
  print(summary(model1))
  model0 = lm(fprop[!is.na(gender_diff)] ~ 1)
  {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_by_gender.pdf"))
    par(mar=c(5,5,1,1))
    plot(fprop, pch=19, cex=1, cex.lab=1.5, las=1, 
         col=gender_col[gender_diff],
         xlab="Family", ylab="Co-membership proportion between siblings")
    text(fprop, labels=ifelse(!is.na(gender_diff),names(gender_diff),""), pos=3,
         col=darken(family.colours[names(gender_diff)], amount=0.5))
    legend("topleft", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
           legend=c("All Male", "All Female", "Different genders"))
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
    dev.off()
  }
  if (i %in% c(2,4)){
    Region=rep(NA, length(families))
    Department=rep(NA, length(families))
    for (f in 1:length(families)){
      tmp=as.character(covars$Region)[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Region[f] = unique(tmp)
      }
      tmp=as.character(covars$Department)[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Department[f] = unique(tmp)
      }
    }
    names(Region)=families
    names(Department)=families
    
    Department = ifelse(Department=="Paris",1,0)
    Region = ifelse(Region=="Île-de-France",1,0)
    
    model = lm(fprop ~ as.factor(Region))
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_by_Region.pdf"))
      par(mar=c(5,5,1,1))
      plot(fprop, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=ifelse(Region==1,region.colours["Île-de-France"],"grey"),
           xlab="Family", ylab="Co-membership proportion between siblings")
      text(fprop, labels=names(Region), pos=3, col=darken(family.colours[names(Region)], amount=0.5))
      legend("top", pch=19, col=c(region.colours["Île-de-France"],"grey"),
             legend=c("Île-de-France","Other region"))
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    model = lm(fprop ~ as.factor(Department))
    {pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_cont_all_by_Department.pdf"))
      par(mar=c(5,5,1,1))
      plot(fprop, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=ifelse(Department==1,depart.colours["Paris"],"grey"),
           xlab="Family", ylab="Co-membership proportion between siblings")
      text(fprop, labels=names(Department), pos=3, col=darken(family.colours[names(Department)], amount=0.5))
      legend("top", pch=19, col=c(depart.colours["Paris"],"grey"),
             legend=c("Paris","Other department"))
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
  }
  
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
  
  # Multivariate analysis using stability selection LASSO regression
  if (i == 1){
    X = cbind(age_mu, age_diff, gender_diff,
              length_mu, length_diff, weight_mu, weight_diff,
              as.data.frame(delta_mat))
  }
  if (i %in% c(2,4)){
    X = cbind(age_mu, age_diff, gender_diff,
              Region,
              as.data.frame(delta_mat))
  }
  if (i %in% c(3,5)){
    X = cbind(age_mu, age_diff, gender_diff,
              as.data.frame(delta_mat))
  }
  Y = fprop[complete.cases(X)]
  X = model.matrix(~., X)[,-1]

  stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)
  
  pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_multivariate_output.pdf"))
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_multivariate_beta_consistency.pdf"))
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
  
  mylabels = gsub("`","",names(selprop_ranked))
  mylabels = gsub("age_mu","Mean age",mylabels)
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("gender_diff","",mylabels)
  mylabels = gsub("All male","All male (ref. All female)",mylabels)
  mylabels = gsub("Different genders","Different genders (ref. All female)",mylabels)
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("Region","Île-de-France (Y/N)",mylabels)

  pdf(paste0("../Figures/",filepaths[i],"/Stab_prop_multivariate_selprop.pdf"), width = 14)
  par(mar=c(15, 5, 1, 1))
  plot(selprop_ranked,
       type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
       col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[i], no = "grey"),
       xaxt = "n", xlab = "", ylab = "Selection proportions"
  )
  abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
  axis(side = 1, at = 1:length(selprop_ranked), labels = NA)
  for (k in 1:length(selprop_ranked)){
    if (beta_mu_ranked[k] > 0){
      col = "red"
    } else if (beta_mu_ranked[k] < 0){
      col = "blue"
    } else {
      col = "black"
    }
    if (mylabels[k] %in% colnames(delta_mat)){
      label = substitute(Delta~tmp, list(tmp = mylabels[k]))
    } else {
      label = mylabels[k]
    }
    axis(side = 1, at = k, las = 2, labels = label,
         col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
  }
  dev.off()
}
