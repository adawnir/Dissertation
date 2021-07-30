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
    names(region.colours)[names(region.colours)=="Normandy"] = "Brittany/Normandy/Pays de la Loire"
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

## Family vs Region
# Scatter plot
annot_sub = annot[names(icc_family_pooled3)]
{pdf(paste0("../Figures/",filepaths[4],"/Intra_class_correlation_family_region.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(icc_family_pooled3, icc_region_pooled3, pch=19,
       col=ifelse(icc_family_pooled3 > 0.5 | icc_region_pooled3 > 0.5, annot.colours[annot_sub], alpha(annot.colours[annot_sub],0.5)),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, 1),
       xlim = c(0, 1),
       ylab="Within-region ICC", 
       xlab="Within-family ICC")
  text(icc_family_pooled3, icc_region_pooled3+0.05, cex = 0.7,
       labels = ifelse(icc_family_pooled3 > 0.5 | icc_region_pooled3 > 0.5, names(icc_family_pooled3), ""),
       col = annot.colours[annot_sub])
  abline(0,1, lty = 2, col = batch.colours[4])    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

# Scatter plot
annot_sub = annot[names(icc_family_pooled2)]
{pdf(paste0("../Figures/",filepaths[5],"/Intra_class_correlation_family_batch.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(icc_family_pooled2, icc_region_pooled2, pch=19,
       col=ifelse(icc_family_pooled2 > 0.5 | icc_region_pooled2 > 0.5, annot.colours[annot_sub], alpha(annot.colours[annot_sub],0.5)),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, 1),
       xlim = c(0, 1),
       ylab="Within-batch ICC", 
       xlab="Within-family ICC")
  text(icc_family_pooled2, icc_region_pooled2+0.05, cex = 0.7,
       labels = ifelse(icc_family_pooled2 > 0.5 | icc_region_pooled2 > 0.5, names(icc_family_pooled2), ""),
       col = annot.colours[annot_sub])
  abline(0,1, lty = 2, col = batch.colours[5])    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

# Scatter plot
annot_sub = annot[names(icc_family_fra)]
{pdf(paste0("../Figures/",filepaths[2],"/Intra_class_correlation_family_region.pdf"), width = 5, height = 5)
  par(mar=c(5,5,1,1))
  plot(icc_family_fra, icc_region_fra, pch=19,
       col=ifelse(icc_family_fra > 0.5 | icc_region_fra > 0.5, annot.colours[annot_sub], alpha(annot.colours[annot_sub],0.5)),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, 1),
       xlim = c(0, 1),
       ylab="Within-region ICC", 
       xlab="Within-family ICC")
  text(icc_family_fra, icc_region_fra+0.05, cex = 0.7,
       labels = ifelse(icc_family_fra > 0.5 | icc_region_fra > 0.5, names(icc_family_fra), ""),
       col = annot.colours[annot_sub])
  abline(0,1, lty = 2, col = region.colours[2])    
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  dev.off()
}

## Family
icc = t(bind_rows(icc_family_lux, icc_family_fra, icc_family_gs, icc_family_pooled3,
                  icc_family_pooled2))
annot_sub = factor.order(annot)[rownames(icc)]
annot_sub = annot_sub[order(annot_sub)]
icc = icc[names(annot_sub),]
mylabels = rownames(icc)

myspacing = 1
xseq = 1:nrow(icc)

{pdf("../Figures/Intra_class_correlation_family_removed_region.pdf", width=14, height=12)
  par(oma=c(21,0,0,0),mar=c(0,5,2,1), mfrow = c(2,1))
  for (i in c(5,4)){
    plot(icc[,i],
         col=annot.colours[annot_sub],
         xaxt="n", ylab=paste0("Intra-class correlation"," (", batches[i] ,")"), xlab = "", cex.lab=1.2,
         type="p", pch=19, ylim = c(0,1.1),
         panel.first=abline(v=1:nrow(icc),lty=3,col="grey"))
  }
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.8)
  }
  xseqblack=c(xseq[!duplicated(annot_sub)]-myspacing/2, max(xseq)+myspacing/2)
  axis(side=1, line=9, at=xseqblack, labels=NA)
  tmp=apply(rbind(xseqblack[-length(xseqblack)],xseqblack[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=9, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  dev.off()
}



## Region
icc = t(bind_rows(icc_region_fra, icc_region_pooled3))
annot_sub = factor.order(annot)[rownames(icc)]
annot_sub = annot_sub[order(annot_sub)]
icc = icc[names(annot_sub),]
mylabels = rownames(icc)

myspacing = 1
xseq = 1:nrow(icc)

{pdf("../Figures/Intra_class_correlation_region_removed_family.pdf", width=14, height=11)
  par(oma=c(21,0,0,0),mar=c(0,5,2,1), mfrow = c(2,1))
  for (i in c(2,4)){
    plot(icc[,i/2],
         col=annot.colours[annot_sub],
         xaxt="n", ylab=paste0("Intra-class correlation"," (", batches[i] ,")"), xlab = "", cex.lab=1.2,
         type="p", pch=19, ylim = c(0,1.1),
         panel.first=abline(v=1:nrow(icc),lty=3,col="grey"))
  }
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.8)
  }
  xseqblack=c(xseq[!duplicated(annot_sub)]-myspacing/2, max(xseq)+myspacing/2)
  axis(side=1, line=9, at=xseqblack, labels=NA)
  tmp=apply(rbind(xseqblack[-length(xseqblack)],xseqblack[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=9, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  dev.off()
}

# # Manhattan plot (Family ID, Main)
# values = t(bind_rows(icc_family_lux, icc_family_fra, icc_family_gs, icc_family_pooled3, icc_family_pooled2))
# annot_sub = annot[rownames(values)]
# annot_sub = annot_sub[order(annot_sub)]
# values = values[names(annot_sub),]
# 
# xseq = seq(1, nrow(values))
# 
# {pdf(paste0("../Figures/Intra_class_correlation_univariate_expo_cont.pdf"), width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(values[,1], pch=19, las=1, xaxt="n", type = "n",
#        ylim = c(0,1),
#        xlab="", ylab="Intra-Class Correlation within families", cex.lab=1.2,
#        panel.first=abline(v=xseq,lty=3,col="grey"),
#        col=batch.colours[1])
#   points(values[,1], pch = 19, col = batch.colours[1], cex = 0.8)
#   points(values[,2], pch = 19, col = batch.colours[2], cex = 0.8)
#   points(values[,3], pch = 19, col = batch.colours[3], cex = 0.8)
#   points(values[,4], pch = 17, col = batch.colours[4], cex = 0.8)
#   for (k in 1:length(xseq)){
#     axis(side=1, at=xseq[k], labels=rownames(values)[k], cex.axis=0.8, las=2)
#   }
#   xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
#   axis(side=1, line=8, at=xgroup, labels=NA)
#   tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
#          cex.axis=0.9)
#   }
#   legend("top", pch=c(rep(19,3),17),
#          col=batch.colours[1:4],
#          legend = batches[1:4],
#          horiz = TRUE,
#          bg="white", cex = 0.8)
#   dev.off()
# }
# 
# {pdf(paste0("../Figures/Intra_class_correlation_univariate_expo_cont_Pooled2.pdf"), width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(values[,1], pch=19, las=1, xaxt="n", type = "n",
#        ylim = c(0,1),
#        xlab="", ylab="Intra-Class Correlation within families", cex.lab=1.2,
#        panel.first=abline(v=xseq,lty=3,col="grey"),
#        col=batch.colours[1])
#   points(values[,1], pch = 19, col = batch.colours[1], cex = 0.8)
#   points(values[,3], pch = 19, col = batch.colours[3], cex = 0.8)
#   points(values[,5], pch = 17, col = batch.colours[5], cex = 0.8)
#   for (k in 1:length(xseq)){
#     axis(side=1, at=xseq[k], labels=rownames(values)[k], cex.axis=0.8, las=2)
#   }
#   xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
#   axis(side=1, line=8, at=xgroup, labels=NA)
#   tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
#          cex.axis=0.9)
#   }
#   legend("top", pch=c(rep(19,2),17),
#          col=batch.colours[c(1,3,5)],
#          legend = batches[c(1,3,5)],
#          horiz = TRUE,
#          bg="white", cex = 0.8)
#   dev.off()
# }

# annot_sub=annot[names(icc_batch_pooled2)]
# {pdf("../Figures/Pooled2/Intra_class_correlation_batch_univariate_expo_cont.pdf", width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(icc_batch_pooled2, pch=19, cex=1, las=1, xaxt="n",
#        xlab="", ylab="Intra-Class Correlation wtihin batch", cex.lab=1.2,
#        panel.first=abline(v=1:length(icc_batch_pooled2),lty=3,col="grey"),
#        col=annot.colours[annot_sub])
#   for (k in 1:length(icc_batch_pooled2)){
#     axis(side=1, at=k, labels=names(icc_batch_pooled2)[k], cex.axis=0.8, las=2)
#   }
#   xgroup=c(which(!duplicated(annot_sub))-0.5, length(icc_batch_pooled2)+0.5)
#   axis(side=1, line=8, at=xgroup, labels=NA)
#   tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
#          cex.axis = 0.9)
#   }
#   dev.off()
# }
# 
# annot_sub=annot[names(icc_batch_pooled3)]
# {pdf("../Figures/Pooled3/Intra_class_correlation_batch_univariate_expo_cont.pdf", width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(icc_batch_pooled3, pch=19, cex=1, las=1, xaxt="n",
#        xlab="", ylab="Intra-Class Correlation wtihin batch", cex.lab=1.2,
#        panel.first=abline(v=1:length(icc_batch_pooled3),lty=3,col="grey"),
#        col=annot.colours[annot_sub])
#   for (k in 1:length(icc_batch_pooled3)){
#     axis(side=1, at=k, labels=names(icc_batch_pooled3)[k], cex.axis=0.8, las=3)
#   }
#   xgroup=c(which(!duplicated(annot_sub))-0.5, length(icc_batch_pooled3)+0.5)
#   axis(side=1, line=8, at=xgroup, labels=NA)
#   tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
#          cex.axis = 0.9)
#   }
#   dev.off()
# }
