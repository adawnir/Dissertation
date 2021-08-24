## Within family variation (excluding isolated children)
## 30 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)

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
  
  ifelse(dir.exists("../Figures/Section2"), "", dir.create("../Figures/Section2"))
  if (i == 4){
    pdf(paste0("../Figures/Section2/Univariate_sd_family_overall_",suffix[i],".pdf"), width=14, height=8)
  } else {    pdf(paste0("../Figures/Supplementary/Section2/Univariate_sd_family_overall_",suffix[i],".pdf"), width=14, height=8)
    }
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
  
  if (i %in% c(2,4,5)){
    covars$Region = ifelse(as.character(covars$Region) %in% c("Brittany","Normandy","Pays de la Loire"),
                           "Brittany/Normandy/Pays de la Loire", as.character(covars$Region))
    covars$Region = factor.order(covars$Region)
    table(covars$Region)
    names(region.colours)[names(region.colours)=="Brittany"] = "Brittany/Normandy/Pays de la Loire"
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
    
    {pdf(paste0("../Figures/Supplementary/Section2/Univariate_sd_region_overall_",suffix[i],".pdf"), width=14, height=8)
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


## Family vs Region (Pooled3)
family = readRDS(paste0("../Results/",filepaths[4],"/Univariate_sd_family_overall_pvals.rds"))
region = readRDS(paste0("../Results/",filepaths[4],"/Univariate_sd_region_overall_pvals.rds"))

annot_sub = annot[names(family)]

# Scatter plot
{pdf(paste0("../Figures/Section2/Univariate_sd_family_region_pvals_pooled3.pdf"), width = 9, height = 5)
  par(mar=c(5,5,1,20))
  plot(-log10(family), -log10(region), pch=19,
       col=ifelse(-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region)),
                  annot.colours[annot_sub], "grey"),
       cex.lab=1.5, cex = 3,
       ylab=substitute(paste(-log[10](italic(p))," ",tmp), list(tmp = "Within-region")), 
       xlab=substitute(paste(-log[10](italic(p))," ",tmp), list(tmp = "Within-family")))
  text(x = max(-log10(family)), y = -log10(0.05/length(family)),labels="Bonferroni threshold", adj = c(1,1),col = batch.colours[4], cex = 0.8)
  abline(h = -log10(0.05/length(family)), lty = 2, col = batch.colours[4])
  text(x = -log10(0.05/length(region)), y = max(-log10(region)),labels="Bonferroni threshold", adj = c(0,0),col = batch.colours[4], cex = 0.8, srt = 270)
  abline(v = -log10(0.05/length(region)), lty = 2, col = batch.colours[4])
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  for (k in 1:length(family)){
    if(-log10(family)[k] > -log10(0.05/length(family))| -log10(region)[k] > -log10(0.05/length(region))){
      text(-log10(family)[k], -log10(region)[k],
           labels = which(names(family)[k] == names(family)[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]),
           cex = 0.8)
    }
  }
  par(xpd = TRUE)
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))),
                         ". ", names(family)[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]),
         text.col = darken(annot.colours[annot_sub[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]],0.5),
         bty = "n",  x.intersp = 0, ncol = 2)
  dev.off()
}

## Family vs Region (France)
family = readRDS(paste0("../Results/",filepaths[2],"/Univariate_sd_family_overall_pvals.rds"))
region = readRDS(paste0("../Results/",filepaths[2],"/Univariate_sd_region_overall_pvals.rds"))

annot_sub = annot[names(family)]

# Scatter plot
{pdf(paste0("../Figures/Supplementary/Section2/Univariate_sd_family_region_pvals_fra.pdf"), width = 9, height = 5)
  par(mar=c(5,5,1,20))
  plot(-log10(family), -log10(region), pch=19,
       col=ifelse(-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region)),
                  annot.colours[annot_sub], "grey"),
       cex.lab=1.5, cex = 3,
       ylab=substitute(paste(-log[10](italic(p))," ",tmp), list(tmp = "Within-region")), 
       xlab=substitute(paste(-log[10](italic(p))," ",tmp), list(tmp = "Within-family")))
  text(x = max(-log10(family)), y = -log10(0.05/length(family)),labels="Bonferroni threshold", adj = c(1,1),col = batch.colours[2], cex = 0.8)
  abline(h = -log10(0.05/length(family)), lty = 2, col = batch.colours[2])
  text(x = -log10(0.05/length(region)), y = max(-log10(region)),labels="Bonferroni threshold", adj = c(0,0),col = batch.colours[2], cex = 0.8, srt = 270)
  abline(v = -log10(0.05/length(region)), lty = 2, col = batch.colours[2])
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  for (k in 1:length(family)){
    if(-log10(family)[k] > -log10(0.05/length(family))| -log10(region)[k] > -log10(0.05/length(region))){
      text(-log10(family)[k], -log10(region)[k],
           labels = which(names(family)[k] == names(family)[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]),
           cex = 0.8)
    }
  }
  par(xpd = TRUE)
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))),
                         ". ", names(family)[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]),
         text.col = darken(annot.colours[annot_sub[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]],0.5),
         bty = "n",  x.intersp = 0, ncol = 2)
  dev.off()
}

## Family vs Region (Pooled2)
family = readRDS(paste0("../Results/",filepaths[5],"/Univariate_sd_family_overall_pvals.rds"))
region = readRDS(paste0("../Results/",filepaths[5],"/Univariate_sd_region_overall_pvals.rds"))

annot_sub = annot[names(family)]

# Scatter plot
{pdf(paste0("../Figures/Supplementary/Section2/Univariate_sd_family_region_pvals_pooled2.pdf"), width = 9, height = 5)
  par(mar=c(5,5,1,20))
  plot(-log10(family), -log10(region), pch=19,
       col=ifelse(-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region)),
                  annot.colours[annot_sub], "grey"),
       cex.lab=1.5, cex = 3,
       ylab=substitute(paste(-log[10](italic(p))," ",tmp), list(tmp = "Within-country")), 
       xlab=substitute(paste(-log[10](italic(p))," ",tmp), list(tmp = "Within-family")))
  text(x = max(-log10(family)), y = -log10(0.05/length(family)),labels="Bonferroni threshold", adj = c(1,1),col = batch.colours[5], cex = 0.8)
  abline(h = -log10(0.05/length(family)), lty = 2, col = batch.colours[5])
  text(x = -log10(0.05/length(region)), y = max(-log10(region)),labels="Bonferroni threshold", adj = c(0,0),col = batch.colours[5], cex = 0.8, srt = 270)
  abline(v = -log10(0.05/length(region)), lty = 2, col = batch.colours[5])
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  for (k in 1:length(family)){
    if(-log10(family)[k] > -log10(0.05/length(family))| -log10(region)[k] > -log10(0.05/length(region))){
      text(-log10(family)[k], -log10(region)[k],
           labels = which(names(family)[k] == names(family)[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]),
           cex = 0.8)
    }
  }
  par(xpd = TRUE)
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))),
                         ". ", names(family)[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]),
         text.col = darken(annot.colours[annot_sub[-log10(family) > -log10(0.05/length(family))| -log10(region) > -log10(0.05/length(region))]],0.5),
         bty = "n",  x.intersp = 0, ncol = 2)
  dev.off()
}



