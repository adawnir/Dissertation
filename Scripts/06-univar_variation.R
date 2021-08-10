## Within family variation (excluding isolated children)
## Rin 30 July

# Load packages
library(tidyverse)
library(RColorBrewer)

### Within-family variation ----
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
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  families=levels(covars$Family.ID)
  family_sd=matrix(NA, nrow=length(families), ncol=ncol(expo))
  for (p in 1:ncol(expo)){
    for (f in 1:length(families)){
      family_sd[f,p]=sd(expo[covars$Family.ID==families[f],p])
    }
  }
  colnames(family_sd)=colnames(expo)
  rownames(family_sd)=families
  
  overall_sd=apply(expo,2,sd)

  x=as.vector(rownames(family_sd))
  y=as.vector(family_sd)
  z=as.vector(col(family_sd))
  
  annot_sub=annot[colnames(family_sd)]
  
  # One-way ANOVA
  pvals = NULL
  for (k in 1:ncol(expo)){
    pvals = c(pvals, summary(aov(expo[,k] ~ covars$Family.ID))[[1]]$`Pr(>F)`[1])
  }
  names(pvals) = colnames(expo)
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_sd_family_overall_pvals.rds"))
  bonf = 0.05/ncol(expo)
  
  {pdf(paste0("../Figures/",filepaths[i],"/Univariate_sd_family_overall_cont.pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(z,y, pch=19, cex=0.5, col=family.colours[x], las=1, xaxt="n",
           xlab="", ylab="Standard deviation", cex.lab=1.5,
           panel.first=abline(v=levels(z),lty=3,col="grey"),
           ylim=range(c(family_sd,overall_sd)))
      points(overall_sd, pch=15)
      for (k in 1:ncol(family_sd)){
        axis(side=1, at=k, labels=colnames(family_sd)[k], cex.axis=0.8, las=2,
             col.axis = ifelse(pvals[k] < bonf, darken(annot.colours[annot_sub[k]],0.5), "black"))
      }
      xgroup=c(which(!duplicated(annot_sub))-0.5, length(colnames(family_sd))+0.5)
      axis(side=1, line=8, at=xgroup, labels=NA)
      tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
      for (k in 1:length(unique(annot_sub))){
        axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
             col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
      }
      # legend("top", col=c(family.colours[families], "black"),
      #        pch=c(rep(19,length(families)),15),
      #        pt.cex=c(rep(0.5,length(families)),1),
      #        legend=c(as.character(families), "Overall"), ncol=10, bg="white")
      dev.off()
  }
  if (i %in% 4:5){
    print(all(rownames(expo)==rownames(covars)))
    
    batch_sd=matrix(NA, nrow=length(levels(covars$Batch)), ncol=ncol(expo))
    for (p in 1:ncol(expo)){
      for (b in 1:length(levels(covars$Batch))){
        batch_sd[b,p]=sd(expo[covars$Batch==levels(covars$Batch)[b],p])
      }
    }
    colnames(batch_sd)=colnames(expo)
    rownames(batch_sd)=levels(covars$Batch)
    
    x=as.vector(rownames(batch_sd))
    y=as.vector(batch_sd)
    z=as.vector(col(batch_sd))
    
    # One-way ANOVA
    pvals = NULL
    for (k in 1:ncol(expo)){
      pvals = c(pvals, summary(aov(expo[,k] ~ covars$Batch))[[1]]$`Pr(>F)`[1])
    }
    names(pvals) = colnames(expo)
    saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_sd_batch_overall_pvals.rds"))
    
    {pdf(paste0("../Figures/",filepaths[i],"/Univariate_sd_batch_overall_cont.pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(z,y, pch=19, cex=0.5, col=batch.colours[x], las=1, xaxt="n",
           xlab="", ylab="Standard deviation", cex.lab=1.5,
           panel.first=abline(v=unique(z),lty=3,col="grey"),
           ylim=range(c(batch_sd,overall_sd)))
      points(overall_sd, pch=15)
      for (k in 1:ncol(batch_sd)){
        axis(side=1, at=k, labels=colnames(batch_sd)[k], cex.axis=0.8, las=2,
             col.axis = ifelse(pvals[k] < bonf, darken(annot.colours[annot_sub[k]], amount=0.5), "black"))
      }
      xgroup=c(which(!duplicated(annot_sub))-0.5, length(colnames(batch_sd))+0.5)
      axis(side=1, line=8, at=xgroup, labels=NA)
      tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
      for (k in 1:length(unique(annot_sub))){
        axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
             col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
      }
      legend("top", col=c(batch.colours[x], "black"),
             pch=c(rep(19,length(x)),15),
             pt.cex=c(rep(0.5,length(x)),1),
             legend=c(batches[1:3], "Overall"),
             ncol=length(x)+1, bg="white", cex = 0.8)
      dev.off()
    }
  }
  if (i %in% c(2,4)){
    covars$Region = ifelse(as.character(covars$Region) %in% c("Brittany","Normandy","Pays de la Loire"),
                           "Brittany/Normandy/Pays de la Loire", as.character(covars$Region))
    covars$Region = factor.order(covars$Region)
    table(covars$Region)
    names(region.colours)[names(region.colours)=="Normandy"] = "Brittany/Normandy/Pays de la Loire"
    region_sd=matrix(NA, nrow=length(levels(covars$Region)), ncol=ncol(expo))
    for (p in 1:ncol(expo)){
      for (r in 1:length(levels(covars$Region))){
        region_sd[r,p]=sd(expo[covars$Region==levels(covars$Region)[r],p])
      }
    }
    colnames(region_sd)=colnames(expo)
    rownames(region_sd)=levels(covars$Region)
    
    x=as.vector(rownames(region_sd))
    y=as.vector(region_sd)
    z=as.vector(col(region_sd))
    
    # One-way ANOVA
    pvals = NULL
    for (k in 1:ncol(expo)){
      pvals = c(pvals, summary(aov(expo[,k] ~ covars$Region))[[1]]$`Pr(>F)`[1])
    }
    names(pvals) = colnames(expo)
    saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_sd_region_overall_pvals.rds"))
    
    {pdf(paste0("../Figures/",filepaths[i],"/Univariate_sd_region_overall_cont.pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(z,y, pch=19, cex=0.5, col=region.colours[x], las=1, xaxt="n",
           xlab="", ylab="Standard deviation", cex.lab=1.5,
           panel.first=abline(v=unique(z),lty=3,col="grey"),
           ylim=range(c(region_sd,overall_sd)))
      points(overall_sd, pch=15)
      for (k in 1:ncol(region_sd)){
        axis(side=1, at=k, labels=colnames(region_sd)[k], cex.axis=0.8, las=2,
             col.axis = ifelse(pvals[k] < bonf, darken(annot.colours[annot_sub[k]], amount=0.5), "black"))
      }
      xgroup=c(which(!duplicated(annot_sub))-0.5, length(colnames(region_sd))+0.5)
      axis(side=1, line=8, at=xgroup, labels=NA)
      tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
      for (k in 1:length(unique(annot_sub))){
        axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
             col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
      }
      legend("top", col=c(region.colours[x], "black"),
             pch=c(rep(19,length(x)),15),
             pt.cex=c(rep(0.5,length(x)),1),
             legend=c(x, "Overall"), cex = 0.8,
             ncol=ceiling(length(x)/2), bg="white")
      dev.off()
    }
  }
}

## Family vs batch (Pooled3)
family = readRDS(paste0("../Results/",filepaths[4],"/Univariate_sd_family_overall_pvals.rds"))
batch = readRDS(paste0("../Results/",filepaths[4],"/Univariate_sd_batch_overall_pvals.rds"))

annot_sub = annot[names(family)]

# Scatter plot
{pdf(paste0("../Figures/",filepaths[4],"/Univariate_sd_family_batch_pvals.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(-log10(family), -log10(batch), pch=19,
       col=ifelse(family < 0.05/length(family) & batch < 0.05/length(batch), annot.colours[annot_sub],
                  ifelse(family < 0.05/length(family)|batch < 0.05/length(batch), alpha(annot.colours[annot_sub],0.5),"grey")),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, max(-log10(batch))+2),
       xlim = c(0, max(-log10(family))+2),
       ylab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-batch")), 
       xlab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-family")))
  text(-log10(family), -log10(batch)+1, cex = 0.7,
       labels = ifelse(family < 0.05/length(family)&batch < 0.05/length(batch), names(family), ""),
       col = annot.colours[annot_sub])
  abline(v = -log10(0.05/length(family)), lty = 2, col = batch.colours[4])
  abline(h = -log10(0.05/length(batch)), lty = 2, col = batch.colours[4])    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

## Family vs batch (Pooled2)
family = readRDS(paste0("../Results/",filepaths[5],"/Univariate_sd_family_overall_pvals.rds"))
batch = readRDS(paste0("../Results/",filepaths[5],"/Univariate_sd_batch_overall_pvals.rds"))

annot_sub = annot[names(family)]

# Scatter plot
{pdf(paste0("../Figures/",filepaths[5],"/Univariate_sd_family_batch_pvals.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(-log10(family), -log10(batch), pch=19,
       col=ifelse(family < 0.05/length(family) & batch < 0.05/length(batch), annot.colours[annot_sub],
                  ifelse(family < 0.05/length(family)|batch < 0.05/length(batch), alpha(annot.colours[annot_sub],0.5),"grey")),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, max(-log10(batch))+2),
       xlim = c(0, max(-log10(family))+2),
       ylab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-batch")), 
       xlab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-family")))
  text(-log10(family), -log10(batch)+1, cex = 0.7,
       labels = ifelse(family < 0.05/length(family)&batch < 0.05/length(batch), names(family), ""),
       col = annot.colours[annot_sub])
  abline(v = -log10(0.05/length(family)), lty = 2, col = batch.colours[5])
  abline(h = -log10(0.05/length(batch)), lty = 2, col = batch.colours[5])    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

## Family vs Region (Pooled3)
family = readRDS(paste0("../Results/",filepaths[4],"/Univariate_sd_family_overall_pvals.rds"))
region = readRDS(paste0("../Results/",filepaths[4],"/Univariate_sd_region_overall_pvals.rds"))

annot_sub = annot[names(family)]

# Scatter plot
{pdf(paste0("../Figures/",filepaths[4],"/Univariate_sd_family_region_pvals.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(-log10(family), -log10(region), pch=19,
       col=ifelse(family < 0.05/length(family) & region < 0.05/length(region), annot.colours[annot_sub],
                  ifelse(family < 0.05/length(family)|region < 0.05/length(region), alpha(annot.colours[annot_sub],0.5),"grey")),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, max(-log10(region))+2),
       xlim = c(0, max(-log10(family))+2),
       ylab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-region")), 
       xlab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-family")))
  text(-log10(family), -log10(region)+1, cex = 0.7,
       labels = ifelse(family < 0.05/length(family)&region < 0.05/length(region), names(family), ""),
       col = annot.colours[annot_sub])
  abline(v = -log10(0.05/length(family)), lty = 2, col = region.colours[4])
  abline(h = -log10(0.05/length(region)), lty = 2, col = region.colours[4])    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

## Family vs Region (France)
family = readRDS(paste0("../Results/",filepaths[2],"/Univariate_sd_family_overall_pvals.rds"))
region = readRDS(paste0("../Results/",filepaths[2],"/Univariate_sd_region_overall_pvals.rds"))

annot_sub = annot[names(family)]

# Scatter plot
{pdf(paste0("../Figures/",filepaths[2],"/Univariate_sd_family_region_pvals.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(-log10(family), -log10(region), pch=19,
       col=ifelse(family < 0.05/length(family) & region < 0.05/length(region), annot.colours[annot_sub],
                  ifelse(family < 0.05/length(family)|region < 0.05/length(region), alpha(annot.colours[annot_sub],0.5),"grey")),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, max(-log10(region))+2),
       xlim = c(0, max(-log10(family))+2),
       ylab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-region")), 
       xlab=substitute(paste(-log[10](italic(p))," (",tmp,")"), list(tmp = "Within-family")))
  text(-log10(family), -log10(region)+1, cex = 0.7,
       labels = ifelse(family < 0.05/length(family)&region < 0.05/length(region), names(family), ""),
       col = annot.colours[annot_sub])
  abline(v = -log10(0.05/length(family)), lty = 2, col = region.colours[2])
  abline(h = -log10(0.05/length(region)), lty = 2, col = region.colours[2])    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}


# ### Within-family absolute difference ----
# # Initialise
# rm(list=ls())
# path=dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(path)
# 
# # Load custom
# source("functions.R")
# source("graph_param.R")
# 
# # Load data sets
# annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
# 
# suffix = c("lux","fra","gs","pooled3","pooled2")
# for (i in 1:length(batches)){
#   # Load data
#   expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
#   covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
#   print(all(rownames(expo)==rownames(covars)))
#   
#   families=unique(covars$Family.ID)
#   family_diff=matrix(NA, nrow=length(families), ncol=ncol(expo))
#   for (p in 1:ncol(expo)){
#     for (f in 1:length(families)){
#       family_diff[f,p]=abs(expo[covars$Family.ID==families[f],p][1]-expo[covars$Family.ID==families[f],p][2])
#     }
#   }
#   colnames(family_diff)=colnames(expo)
#   rownames(family_diff)=families
#   
#   # Relationship wtihin age
#   age_mu=rep(NA, length(families))
#   for (f in 1:length(families)){
#     age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
#   }
#   names(age_mu)=families
#   age_diff=rep(NA, length(families))
#   for (f in 1:length(families)){
#     age_diff[f]=abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2])
#   }
#   names(age_diff)=families
#   gender_diff=rep(NA, length(families))
#   for (f in 1:length(families)){
#     tmp=covars$Gender[covars$Family.ID==families[f]]
#     if (sum(is.na(tmp))>0) {
#       gender_diff[f] = NA
#     } else if (length(unique(tmp))==1){
#       if (unique(tmp)=="Male"){
#         gender_diff[f]=1
#       }
#       if (unique(tmp)=="Female"){
#         gender_diff[f]=2
#       }
#     } else {
#       gender_diff[f]=3
#     }
#   }
#   names(gender_diff)=families
#   if (i==1){
#     weight_mu=rep(NA, length(families))
#     for (f in 1:length(families)){
#       weight_mu[f]=mean(covars$Weight[covars$Family.ID==families[f]])
#     }
#     names(weight_mu)=families
#     weight_diff=rep(NA, length(families))
#     for (f in 1:length(families)){
#       weight_diff[f]=abs(covars$Weight[covars$Family.ID==families[f]][1]-covars$Weight[covars$Family.ID==families[f]][2])
#     }
#     names(weight_diff)=families
#     length_mu=rep(NA, length(families))
#     for (f in 1:length(families)){
#       length_mu[f]=mean(covars$Length[covars$Family.ID==families[f]])
#     }
#     names(length_mu)=families
#     length_diff=rep(NA, length(families))
#     for (f in 1:length(families)){
#       length_diff[f]=abs(covars$Length[covars$Family.ID==families[f]][1]-covars$Length[covars$Family.ID==families[f]][2])
#     }
#     names(length_diff)=families
#   }
#   if (i %in% 4:5){
#     Batch=rep(NA, length(families))
#     for (f in 1:length(families)){
#       tmp=covars$Batch[covars$Family.ID==families[f]]
#       if (length(unique(tmp))==1){
#         Batch[f] = unique(tmp)
#       }
#     }
#     names(Batch)=families
#   }
#   if (i == 4){
#     Region=rep(NA, length(families))
#     Department=rep(NA, length(families))
#     for (f in 1:length(families)){
#       tmp=covars$Region[covars$Family.ID==families[f]]
#       if (length(unique(tmp))==1){
#         Region[f] = unique(tmp)
#       }
#       tmp=covars$Department[covars$Family.ID==families[f]]
#       if (length(unique(tmp))==1){
#         Department[f] = unique(tmp)
#       }
#     }
#     names(Region)=families
#     names(Department)=families
#   }
#   pvals = NULL
#   for (k in 1:ncol(family_diff)){
#     model = lm(family_diff[,k] ~ age_mu)
#     pvals = c(pvals, summary(model)$coefficients[2,4])
#   }
#   for (k in 1:ncol(family_diff)){
#     model = lm(family_diff[,k] ~ age_diff)
#     pvals = c(pvals, summary(model)$coefficients[2,4])
#   }
#   for (k in 1:ncol(family_diff)){
#     tmp = family_diff[!is.na(gender_diff),]
#     model1 = lm(tmp[,k] ~ as.factor(gender_diff[!is.na(gender_diff)]))
#     model0 = lm(tmp[,k] ~ 1)
#     pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
#   }
#   if (i==1){
#     for (k in 1:ncol(family_diff)){
#       model = lm(family_diff[,k] ~ weight_mu)
#       pvals = c(pvals, summary(model)$coefficients[2,4])
#     }
#     for (k in 1:ncol(family_diff)){
#       model = lm(family_diff[,k] ~ weight_diff)
#       pvals = c(pvals, summary(model)$coefficients[2,4])
#     }
#     for (k in 1:ncol(family_diff)){
#       model = lm(family_diff[,k] ~ length_mu)
#       pvals = c(pvals, summary(model)$coefficients[2,4])
#     }
#     for (k in 1:ncol(family_diff)){
#       model = lm(family_diff[,k] ~ length_diff)
#       pvals = c(pvals, summary(model)$coefficients[2,4])
#     }
#   }
#   if (i %in% 4:5){
#     for (k in 1:ncol(family_diff)){
#       model1 = lm(family_diff[,k] ~ as.factor(Batch))
#       model0 = lm(family_diff[,k] ~ 1)
#       pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
#     }
#   }
#   if (i==4){
#     for (k in 1:ncol(family_diff)){
#       model1 = lm(family_diff[,k] ~ as.factor(Region))
#       model0 = lm(family_diff[,k] ~ 1)
#       pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
#     }
#     for (k in 1:ncol(family_diff)){
#       model1 = lm(family_diff[,k] ~ as.factor(Department))
#       model0 = lm(family_diff[,k] ~ 1)
#       pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
#     }
#   }
#   X = c("age_mu","age_diff","gender_diff")
#   if (i==1){X = c(X, "weight_mu","weight_diff","length_mu","length_diff")}
#   if (i %in% 4:5){X = c(X, "Batch")}
#   if (i==4){X = c(X, "Region","Department")}
#   # If p-value = 0 replace wtihin smallest double in R
#   pvals = ifelse(pvals == 0, .Machine$double.xmin, pvals)
#   pvals = matrix(pvals, nrow = ncol(family_diff), ncol = length(X))
#   rownames(pvals)=colnames(family_diff)
#   colnames(pvals)= X
#   saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_family_diff_covariate_pvals.rds"))
#   annot_sub = annot[rownames(pvals)]
#   # Manhattan plot
#   xseq = seq(1, nrow(pvals))
#   for (j in colnames(pvals)){
#     {pdf(paste0("../Figures/",filepaths[i],"/Univariate_family_diff_",j,".pdf"), width=14, height=8)
#       par(mar=c(20,5,1,1))
#       plot(-log10(pvals[,j]),
#            pch = ifelse(i %in% 4:5, 17, 19), col=annot.colours[annot_sub],
#            xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
#            ylim = c(min(-log10(pvals[,j]), na.rm = TRUE),-log10(0.05/nrow(pvals))+0.1))
#       abline(h = -log10(0.05), lty = 2, col = "grey")
#       abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
#       abline(h = -log10(0.05/nrow(pvals)), lty = 2, col = darken(batch.colours[4], 0.5))
#       for(k in 1:length(xseq)){
#         axis(1, at=xseq[k], labels = rownames(pvals)[k], las=2, cex.axis = 0.8,
#              col.axis = ifelse(isTRUE(pvals[k,j] < 0.05/nrow(pvals)),
#                                darken(annot.colours[(annot_sub)[k]], amount=0.5),"black"))
#       }
#       xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals))+0.5)
#       axis(side=1, line=8, at=xgroup, labels=NA)
#       tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#       for (k in 1:length(unique(annot_sub))){
#         axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#              col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#       }
#       legend("top", lty = c(rep(2,2)),
#              col=c(darken(batch.colours[i], 0.5), "grey"),
#              legend = c("Bonferroni threshold","Nominal threshold"),
#              horiz = TRUE,
#              bg="white", cex = 0.8)
#       dev.off()
#     }
#   }
# }


