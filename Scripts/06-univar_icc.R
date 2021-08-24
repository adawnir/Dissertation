## Intra-class correlation (excluding isolated children)
## Rin 30 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(lme4)

### Intra-class correlation ----
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
  
  ### Univariate linear mixed models
  if (i %in% c(2,4,5)){
    icc_family = icc_region = NULL
    covars$Region = ifelse(as.character(covars$Region) %in% c("Brittany","Normandy","Pays de la Loire"),
                           "Brittany/Normandy/Pays de la Loire", as.character(covars$Region))
    covars$Region = factor.order(covars$Region)
    table(covars$Region)
    names(region.colours)[names(region.colours)=="Brittany"] = "Brittany/Normandy/Pays de la Loire"
    for (p in 1:ncol(expo)){
      x=expo[,p]
      model=lmer(x~(1|covars$Family.ID) + (1|covars$Region))
      vcov = as.data.frame(VarCorr(model))$vcov
      icc_family=c(icc_family,vcov[1]/sum(vcov))
      icc_region=c(icc_region,vcov[2]/sum(vcov))
    }
    names(icc_family) = colnames(expo)
    names(icc_region) = colnames(expo)
    assign(paste0("icc_family_",suffix[i]),icc_family)
    assign(paste0("icc_region_",suffix[i]),icc_region)
  } else {
    icc_family = NULL
    for (p in 1:ncol(expo)){
      x=expo[,p]
      model=lmer(x~(1|covars$Family.ID))
      vcov = as.data.frame(VarCorr(model))$vcov
      icc_family=c(icc_family,vcov[1]/sum(vcov))
    }
    names(icc_family) = colnames(expo)
    assign(paste0("icc_family_",suffix[i]),icc_family)
  }
}

### Plotting ----
## Family vs Region
# Scatter plot
annot_sub = annot[names(icc_family_pooled3)]
{pdf(paste0("../Figures/Section2/Intra_class_correlation_family_region_pooled3.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(icc_family_pooled3, icc_region_pooled3, pch=19,
       col=ifelse(icc_family_pooled3 > 0.7 | icc_region_pooled3 > 0.7, annot.colours[annot_sub], alpha(annot.colours[annot_sub],0.5)),
       cex.lab=1.5,
       ylim = c(0, 1),
       xlim = c(0, 1),
       ylab="Region ICC", 
       xlab="Family ICC", cex = 1.5)
  text(icc_family_pooled3, icc_region_pooled3+0.03,
       labels = ifelse(icc_family_pooled3 > 0.7 | icc_region_pooled3 > 0.7, names(icc_family_pooled3), ""),
       col = annot.colours[annot_sub])
  abline(0,1, lty = 2, col = "grey")    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

# Scatter plot
annot_sub = annot[names(icc_family_fra)]
{pdf(paste0("../Figures/Supplementary/Section2/Intra_class_correlation_family_region_fra.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(icc_family_fra, icc_region_fra, pch=19,
       col=ifelse(icc_family_fra > 0.7 | icc_region_fra > 0.7, annot.colours[annot_sub], alpha(annot.colours[annot_sub],0.5)),
       cex.lab=1.5,
       ylim = c(0, 1),
       xlim = c(0, 1),
       ylab="Region ICC", 
       xlab="Family ICC", cex = 1.5)
  text(icc_family_fra, icc_region_fra+0.03,
       labels = ifelse(icc_family_fra > 0.7 | icc_region_fra > 0.7, names(icc_family_fra), ""),
       col = annot.colours[annot_sub])
  abline(0,1, lty = 2, col = "grey")    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}


# Scatter plot
annot_sub = annot[names(icc_family_pooled2)]
{pdf(paste0("../Figures/Section2/Intra_class_correlation_family_batch_pooled2.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(icc_family_pooled2, icc_region_pooled2, pch=19,
       col=ifelse(icc_family_pooled2 > 0.7 | icc_region_pooled2 > 0.7, annot.colours[annot_sub], alpha(annot.colours[annot_sub],0.5)),
       cex.lab=1.5,
       ylim = c(0, 1),
       xlim = c(0, 1),
       ylab="Country ICC", 
       xlab="Family ICC", cex = 1.5)
  mygrep = icc_family_pooled2 > 0.7
  jitter = seq(0.7,0.9,length.out = sum(mygrep))
  names(jitter) = names(icc_family_pooled2)[mygrep][order(icc_family_pooled2[mygrep])]
  yjitter = seq(0.2,0.4,length.out = sum(mygrep))
  names(yjitter) = names(icc_family_pooled2)[mygrep][order(icc_family_pooled2[mygrep])]
  
  mygrep2 = icc_region_pooled2 > 0.7
  jitter2 = seq(0.7,0.9,length.out = sum(mygrep2))
  names(jitter2) = names(icc_region_pooled2)[mygrep2][order(icc_region_pooled2[mygrep2])]
  for(k in 1:length(icc_family_pooled2)){
    if(icc_family_pooled2[k] > 0.7){
        x = jitter[names(icc_family_pooled2)[k]]
        y = yjitter[names(icc_family_pooled2)[k]]
        segments(icc_family_pooled2[k], icc_region_pooled2[k], x,  y, col = annot.colours[annot_sub][k])
        text(x, y, labels = names(icc_family_pooled2)[k],col = annot.colours[annot_sub][k], adj = c(0.5,0))
    } else if (icc_region_pooled2[k] > 0.7){
      y = jitter2[names(icc_region_pooled2)[k]]
      x = 0.3
      segments(icc_family_pooled2[k], icc_region_pooled2[k], x,  y, col = annot.colours[annot_sub][k])
      text(x, y, labels = names(icc_family_pooled2)[k],col = annot.colours[annot_sub][k], adj = c(0,0.5))
    }
  }
  abline(0,1, lty = 2, col = "grey")    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

## Family (LUX)
annot_sub = annot[names(icc_family_lux)]
xseq = seq(1, length(icc_family_lux))

{pdf(paste0("../Figures/Supplementary/Section2/Intra_class_correlation_family_",suffix[1],".pdf"), width=13, height=6.5)
  par(mar=c(20,5,1,1))
  plot(icc_family_lux,
       col=annot.colours[annot_sub],
       xaxt="n", xlab="", ylab = "Family ICC", cex.lab=1.5,
       type="p", pch = 19, cex = 1.5,)
  abline(v = xseq, lty = 3, col = "grey")
  for (k in xseq){
    if(icc_family_lux[k] > 0.7){
      colour = darken(annot.colours[annot_sub][k],0.5)
    } else {colour = "black"}
    axis(1, at=k, labels = names(icc_family_lux)[k],col.axis = colour, las=2, cex.axis = 0.8)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(names(icc_family_lux))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  dev.off()
}

## Family (GS)
annot_sub = annot[names(icc_family_gs)]
xseq = seq(1, length(icc_family_gs))

{pdf(paste0("../Figures/Supplementary/Section2/Intra_class_correlation_family_",suffix[3],".pdf"), width=13, height=6.5)
  par(mar=c(20,5,1,1))
  plot(icc_family_gs,
       col=annot.colours[annot_sub],
       xaxt="n", xlab="", ylab = "Family ICC", cex.lab=1.5,
       type="p", pch = 19, cex = 1.5,)
  abline(v = xseq, lty = 3, col = "grey")
  for (k in xseq){
    if(icc_family_gs[k] > 0.7){
      colour = darken(annot.colours[annot_sub][k],0.5)
    } else {colour = "black"}
    axis(1, at=k, labels = names(icc_family_gs)[k],col.axis = colour, las=2, cex.axis = 0.8)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(names(icc_family_gs))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  dev.off()
}

