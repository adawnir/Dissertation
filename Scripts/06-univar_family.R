## Univariate regression and within family variation (excluding isolated children)
## Rin 12 July

# Load packages
library(tidyverse)
library(RColorBrewer)

### Exposure ~ covariate ----
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
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  ## Remove isolated children (no siblings)
  covars = covars %>%
    filter(Family.ID != "Isolated") %>%
    mutate_if(is.factor, droplevels)
  expo = expo[rownames(covars),]
  saveRDS(covars, paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  saveRDS(expo, paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))

  ## Exposure vs family ID by compound
  # Luxembourg
  pvals = NULL
  f1='expo[,k] ~ covars$Family.ID'
  f0='expo[,k] ~ 1'
  t0=Sys.time()
  for (k in 1:ncol(expo)){
    model1=lm(as.formula(f1))
    model0=lm(as.formula(f0))
    pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
  }
  t1=Sys.time()
  print(t1-t0)
  names(pvals) = colnames(expo)
  # If p-value = 0 replace wtihin smallest double in R
  pvals = ifelse(pvals ==0, .Machine$double.xmin, pvals)
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_family_pvals.rds"))
  assign(paste0("pvals_",suffix[i]),pvals)
}

values = -log10(t(bind_rows(pvals_lux, pvals_fra, pvals_gs, pvals_pooled3, pvals_pooled2)))
annot_sub = annot[rownames(values)]
annot_sub = annot_sub[order(annot_sub)]
values = values[names(annot_sub),]

# Manhattan plot
bonf = c(0.05/length(pvals_lux), 0.05/length(pvals_fra), 0.05/length(pvals_gs),
         0.05/length(pvals_pooled3), 0.05/length(pvals_pooled2))
bonf = -log10(bonf)

xseq = seq(1, nrow(values))

options(scipen=999)
{pdf("../Figures/Univariate_exposure_family.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1],
       col=batch.colours[1],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="n", ylim = c(min(values, na.rm = TRUE), 500),
       log = "y")
  abline(h = -log10(0.05), lty = 2, col = "grey")
  abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
  for (i in 1:4){
    if (i == 4){
      pch = 17
    } else {pch = 19}
    points(values[,i], pch = pch, col = batch.colours[i], cex = 0.8)
    abline(h = bonf[i], lty = 2, col = darken(batch.colours[i], 0.5))
  }
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.8,
         col.axis = ifelse(isTRUE(values[i,4] > bonf[4]),
                           darken(annot.colours[annot_sub[i]], amount=0.5),"black"))
    }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("top", pch=c(rep(19,3),17, rep(NA,5)), lty = c(rep(NA,4),rep(2,5)),
         col=c(batch.colours[1:4],darken(batch.colours[1:4], 0.5), "grey"),
         legend = c(batches[1:4],
                    paste0("Bonferroni threshold (",c("LUX","FRA","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.8, ncol = 5)
  dev.off()
}

{pdf("../Figures/Univariate_exposure_family_Pooled2.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1],
       col=batch.colours[1],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="n", ylim = c(min(values, na.rm = TRUE), 500),
       log = "y")
  abline(h = -log10(0.05), lty = 2, col = "grey")
  abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
  for (i in c(1,3,5)){
    if (i == 5){
      pch = 17
    } else {pch = 19}
    points(values[,i], pch = pch, col = batch.colours[i], cex = 0.8)
    abline(h = bonf[i], lty = 2, col = darken(batch.colours[i], 0.5))
  }
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.8,
         col.axis = ifelse(isTRUE(values[i,5] > bonf[5]),
                           darken(annot.colours[annot_sub[i]], amount=0.5),"black"))
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("top", pch=c(rep(19,2),17, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
         col=c(batch.colours[c(1,3,5)],darken(batch.colours[c(1,3,5)], 0.5), "grey"),
         legend = c(batches[c(1,3,5)],
                    paste0("Bonferroni threshold (",c("LUX","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.8, ncol = 4)
  dev.off()
}
options(scipen=0)

### Detection ~ covariate ----
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
  nd = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_nd.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  ## Remove isolated children (no siblings)
  nd = nd[rownames(covars),]
  print(all(rownames(nd)==rownames(covars)))
  
  saveRDS(nd, paste0("../Processed/",filepaths[i],"/Exposure_matrix_nd_no_isolated.rds"))
  saveRDS(covars, paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))

  ## Detection vs family ID by compound
  # Luxembourg
  pvals = NULL
  f1='nd[,k] ~ covars$Family.ID'
  f0='nd[,k] ~ 1'
  t0=Sys.time()
  for (k in 1:ncol(nd)){
    model1=glm(as.formula(f1), family = "binomial")
    model0=glm(as.formula(f0), family = "binomial")
    pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
  }
  t1=Sys.time()
  print(t1-t0)
  names(pvals) = colnames(nd)
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_detection_family_pvals.rds"))
  assign(paste0("pvals_",suffix[i]),pvals)
}

values = -log10(t(bind_rows(pvals_lux, pvals_fra, pvals_gs, pvals_pooled3, pvals_pooled2)))
annot_sub = annot[rownames(values)]
annot_sub = annot_sub[order(annot_sub)]
values = values[names(annot_sub),]

# Manhattan plot
bonf = c(0.05/length(pvals_lux), 0.05/length(pvals_fra), 0.05/length(pvals_gs),
         0.05/length(pvals_pooled3), 0.05/length(pvals_pooled2))
bonf = -log10(bonf)

xseq = seq(1, nrow(values))


{pdf("../Figures/Univariate_detection_family.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1],
       col=batch.colours[1],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="n", ylim = c(min(values, na.rm = TRUE), 13))
  abline(h = -log10(0.05), lty = 2, col = "grey")
  abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
  for (i in 1:4){
    if (i == 4){
      pch = 17
    } else {pch = 19}
    points(values[,i], pch = pch, col = batch.colours[i], cex = 0.8)
    abline(h = bonf[i], lty = 2, col = darken(batch.colours[i], 0.5))
  }
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.8,
         col.axis = ifelse(isTRUE(values[i,4] > bonf[4]),
                           darken(annot.colours[annot_sub[i]], amount=0.5),"black"))
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("top", pch=c(rep(19,3),17, rep(NA,5)), lty = c(rep(NA,4),rep(2,5)),
         col=c(batch.colours[1:4],darken(batch.colours[1:4], 0.5), "grey"),
         legend = c(batches[1:4],
                    paste0("Bonferroni threshold (",c("LUX","FRA","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.8, ncol = 5)
  dev.off()
}

{pdf("../Figures/Univariate_detection_family_Pooled2.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1],
       col=batch.colours[1],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="n", ylim = c(min(values, na.rm = TRUE), 8))
  abline(h = -log10(0.05), lty = 2, col = "grey")
  abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
  for (i in c(1,3,5)){
    if (i == 5){
      pch = 17
    } else {pch = 19}
    points(values[,i], pch = pch, col = batch.colours[i], cex = 0.8)
    abline(h = bonf[i], lty = 2, col = darken(batch.colours[i], 0.5))
  }
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.8,
         col.axis = ifelse(isTRUE(values[i,5] > bonf[5]),
                           darken(annot.colours[annot_sub[i]], amount=0.5),"black"))
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("top", pch=c(rep(19,2),17, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
         col=c(batch.colours[c(1,3,5)],darken(batch.colours[c(1,3,5)], 0.5), "grey"),
         legend = c(batches[c(1,3,5)],
                    paste0("Bonferroni threshold (",c("LUX","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.8, ncol = 4)
  dev.off()
}

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
  
  # bonf = 0.05/ncol(expo)
  
  {pdf(paste0("../Figures/",filepaths[i],"/Univariate_sd_family_overall_cont.pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(z,y, pch=19, cex=0.5, col=family.colours[x], las=1, xaxt="n",
           xlab="", ylab="Standard deviation", cex.lab=1.5,
           panel.first=abline(v=levels(z),lty=3,col="grey"),
           ylim=range(c(family_sd,overall_sd)))
      points(overall_sd, pch=15)
      for (k in 1:ncol(family_sd)){
        axis(side=1, at=k, labels=colnames(family_sd)[k], cex.axis=0.8, las=2,
             col.axis = ifelse(pvals[k] < 0.05, darken(annot.colours[annot_sub[k]],0.5), "black"))
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
    
    {pdf(paste0("../Figures/",filepaths[i],"/Univariate_sd_batch_overall_cont.pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(z,y, pch=19, cex=0.5, col=batch.colours[x], las=1, xaxt="n",
           xlab="", ylab="Standard deviation", cex.lab=1.5,
           panel.first=abline(v=unique(z),lty=3,col="grey"),
           ylim=range(c(batch_sd,overall_sd)))
      points(overall_sd, pch=15)
      for (k in 1:ncol(batch_sd)){
        axis(side=1, at=k, labels=colnames(batch_sd)[k], cex.axis=0.8, las=2,
             col.axis = ifelse(pvals[k] < 0.05, darken(annot.colours[annot_sub[k]], amount=0.5), "black"))
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
  if (i == 4){
    region_sd=matrix(NA, nrow=length(levels(covars$Region)), ncol=ncol(expo))
    for (p in 1:ncol(expo)){
      for (r in 1:length(levels(covars$Region))){
        region_sd[r,p]=sd(expo[covars$Region==levels(covars$Region)[r],p])
      }
    }
    colnames(region_sd)=colnames(expo)
    rownames(region_sd)=levels(covars$Region)
    depart_sd=matrix(NA, nrow=length(levels(covars$Department)), ncol=ncol(expo))
    for (p in 1:ncol(expo)){
      for (d in 1:length(levels(covars$Department))){
        depart_sd[d,p]=sd(expo[covars$Department==levels(covars$Department)[d],p])
      }
    }
    colnames(depart_sd)=colnames(expo)
    rownames(depart_sd)=levels(covars$Department)
    
    x=as.vector(rownames(region_sd))
    y=as.vector(region_sd)
    z=as.vector(col(region_sd))
    
    # One-way ANOVA
    pvals = NULL
    for (k in 1:ncol(expo)){
      pvals = c(pvals, summary(aov(expo[,k] ~ covars$Region))[[1]]$`Pr(>F)`[1])
    }
    
    {pdf(paste0("../Figures/",filepaths[i],"/Univariate_sd_region_overall_cont.pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(z,y, pch=19, cex=0.5, col=region.colours[x], las=1, xaxt="n",
           xlab="", ylab="Standard deviation", cex.lab=1.5,
           panel.first=abline(v=unique(z),lty=3,col="grey"),
           ylim=range(c(region_sd,overall_sd)))
      points(overall_sd, pch=15)
      for (k in 1:ncol(region_sd)){
        axis(side=1, at=k, labels=colnames(region_sd)[k], cex.axis=0.8, las=2,
             col.axis = ifelse(pvals[k] < 0.05, darken(annot.colours[annot_sub[k]], amount=0.5), "black"))
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
             ncol=5, bg="white")
      dev.off()
    }
    
    x=as.vector(rownames(depart_sd))
    y=as.vector(depart_sd)
    z=as.vector(col(depart_sd))
    
    # One-way ANOVA
    pvals = NULL
    for (k in 1:ncol(expo)){
      pvals = c(pvals, summary(aov(expo[,k] ~ covars$Department))[[1]]$`Pr(>F)`[1])
    }
    
    {pdf(paste0("../Figures/",filepaths[i],"/Univariate_sd_depart_overall_cont.pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(z,y, pch=19, cex=0.5, col=depart.colours[x], las=1, xaxt="n",
           xlab="", ylab="Standard deviation", cex.lab=1.5,
           panel.first=abline(v=unique(z),lty=3,col="grey"),
           ylim=range(c(depart_sd,overall_sd)))
      points(overall_sd, pch=15)
      for (k in 1:ncol(depart_sd)){
        axis(side=1, at=k, labels=colnames(depart_sd)[k], cex.axis=0.8, las=2,
             col.axis = ifelse(pvals[k] < 0.05, darken(annot.colours[annot_sub[k]], amount=0.5), "black"))
      }
      xgroup=c(which(!duplicated(annot_sub))-0.5, length(colnames(depart_sd))+0.5)
      axis(side=1, line=8, at=xgroup, labels=NA)
      tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
      for (k in 1:length(unique(annot_sub))){
        axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
             col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
      }
      legend("top", col=c(depart.colours[x], "black"),
             pch=c(rep(19,length(x)),15),
             pt.cex=c(rep(0.5,length(x)),1),
             legend=c(x, "Overall"), cex = 0.8,
             ncol=5, bg="white")
      dev.off()
    }
  }
}

### Within-family absolute difference ----
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
  
  families=unique(covars$Family.ID)
  family_diff=matrix(NA, nrow=length(families), ncol=ncol(expo))
  for (p in 1:ncol(expo)){
    for (f in 1:length(families)){
      family_diff[f,p]=abs(expo[covars$Family.ID==families[f],p][1]-expo[covars$Family.ID==families[f],p][2])
    }
  }
  colnames(family_diff)=colnames(expo)
  rownames(family_diff)=families
  
  # Relationship wtihin age
  age_mu=rep(NA, length(families))
  for (f in 1:length(families)){
    age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
  }
  names(age_mu)=families
  age_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    age_diff[f]=abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2])
  }
  names(age_diff)=families
  gender_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    tmp=covars$Gender[covars$Family.ID==families[f]]
    if (sum(is.na(tmp))>0) {
      gender_diff[f] = NA
    } else if (length(unique(tmp))==1){
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
  if (i==1){
    weight_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_mu[f]=mean(covars$Weight[covars$Family.ID==families[f]])
    }
    names(weight_mu)=families
    weight_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_diff[f]=abs(covars$Weight[covars$Family.ID==families[f]][1]-covars$Weight[covars$Family.ID==families[f]][2])
    }
    names(weight_diff)=families
    length_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      length_mu[f]=mean(covars$Length[covars$Family.ID==families[f]])
    }
    names(length_mu)=families
    length_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      length_diff[f]=abs(covars$Length[covars$Family.ID==families[f]][1]-covars$Length[covars$Family.ID==families[f]][2])
    }
    names(length_diff)=families
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
  }
  pvals = NULL
  for (k in 1:ncol(family_diff)){
    model = lm(family_diff[,k] ~ age_mu)
    pvals = c(pvals, summary(model)$coefficients[2,4])
  }
  for (k in 1:ncol(family_diff)){
    model = lm(family_diff[,k] ~ age_diff)
    pvals = c(pvals, summary(model)$coefficients[2,4])
  }
  for (k in 1:ncol(family_diff)){
    tmp = family_diff[!is.na(gender_diff),]
    model1 = lm(tmp[,k] ~ as.factor(gender_diff[!is.na(gender_diff)]))
    model0 = lm(tmp[,k] ~ 1)
    pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
  }
  if (i==1){
    for (k in 1:ncol(family_diff)){
      model = lm(family_diff[,k] ~ weight_mu)
      pvals = c(pvals, summary(model)$coefficients[2,4])
    }
    for (k in 1:ncol(family_diff)){
      model = lm(family_diff[,k] ~ weight_diff)
      pvals = c(pvals, summary(model)$coefficients[2,4])
    }
    for (k in 1:ncol(family_diff)){
      model = lm(family_diff[,k] ~ length_mu)
      pvals = c(pvals, summary(model)$coefficients[2,4])
    }
    for (k in 1:ncol(family_diff)){
      model = lm(family_diff[,k] ~ length_diff)
      pvals = c(pvals, summary(model)$coefficients[2,4])
    }
  }
  if (i %in% 4:5){
    for (k in 1:ncol(family_diff)){
      model1 = lm(family_diff[,k] ~ as.factor(Batch))
      model0 = lm(family_diff[,k] ~ 1)
      pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
    }
  }
  if (i==4){
    for (k in 1:ncol(family_diff)){
      model1 = lm(family_diff[,k] ~ as.factor(Region))
      model0 = lm(family_diff[,k] ~ 1)
      pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
    }
    for (k in 1:ncol(family_diff)){
      model1 = lm(family_diff[,k] ~ as.factor(Department))
      model0 = lm(family_diff[,k] ~ 1)
      pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
    }
  }
  X = c("age_mu","age_diff","gender_diff")
  if (i==1){X = c(X, "weight_mu","weight_diff","length_mu","length_diff")}
  if (i %in% 4:5){X = c(X, "Batch")}
  if (i==4){X = c(X, "Region","Department")}
  # If p-value = 0 replace wtihin smallest double in R
  pvals = ifelse(pvals == 0, .Machine$double.xmin, pvals)
  pvals = matrix(pvals, nrow = ncol(family_diff), ncol = length(X))
  rownames(pvals)=colnames(family_diff)
  colnames(pvals)= X
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_family_diff_covariate_pvals.rds"))
  annot_sub = annot[rownames(pvals)]
  # Manhattan plot
  xseq = seq(1, nrow(pvals))
  for (j in colnames(pvals)){
    {pdf(paste0("../Figures/",filepaths[i],"/Univariate_family_diff_",j,".pdf"), width=14, height=8)
      par(mar=c(20,5,1,1))
      plot(-log10(pvals[,j]),
           pch = ifelse(i %in% 4:5, 17, 19), col=annot.colours[annot_sub],
           xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
           ylim = c(min(-log10(pvals[,j]), na.rm = TRUE),-log10(0.05/nrow(pvals))+0.1))
      abline(h = -log10(0.05), lty = 2, col = "grey")
      abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
      abline(h = -log10(0.05/nrow(pvals)), lty = 2, col = darken(batch.colours[4], 0.5))
      for(k in 1:length(xseq)){
        axis(1, at=xseq[k], labels = rownames(pvals)[k], las=2, cex.axis = 0.8,
             col.axis = ifelse(isTRUE(pvals[k,j] < 0.05/nrow(pvals)),
                               darken(annot.colours[(annot_sub)[k]], amount=0.5),"black"))
      }
      xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals))+0.5)
      axis(side=1, line=8, at=xgroup, labels=NA)
      tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
      for (k in 1:length(unique(annot_sub))){
        axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
             col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
      }
      legend("top", lty = c(rep(2,2)),
             col=c(darken(batch.colours[i], 0.5), "grey"),
             legend = c("Bonferroni threshold","Nominal threshold"),
             horiz = TRUE,
             bg="white", cex = 0.8)
      dev.off()
    }
  }
}

### Intra-class correlation ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load package
library(lme4)

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
  icc_family=NULL
  
  
  
  if (i %in% c(2,4,5)){
    for (p in 1:ncol(expo)){
      x=expo[,p]
      model=lmer(x~(1|covars$Family.ID) + (1|covars$Department) + (1|covars$Region))
      vcov = as.data.frame(VarCorr(model))$vcov
      icc_family=c(icc_family,vcov[1]/sum(vcov))
    }
    names(icc_family) = colnames(expo)
    assign(paste0("icc_family_",suffix[i]),icc_family)
  } else {
    for (p in 1:ncol(expo)){
      x=expo[,p]
      model=lmer(x~(1|covars$Family.ID))
      vcov = as.data.frame(VarCorr(model))$vcov
      icc_family=c(icc_family,vcov[1]/sum(vcov))
    }
    names(icc_family) = colnames(expo)
    assign(paste0("icc_family_",suffix[i]),icc_family)
  }
  
  if (i %in% 4:5){
    icc_batch=NULL
    for (p in 1:ncol(expo)){
      x=expo[,p]
      model=lmer(x~(1|covars$Batch))
      vcov = as.data.frame(VarCorr(model))$vcov
      icc_batch=c(icc_batch,vcov[1]/sum(vcov))
    }
    names(icc_batch) = colnames(expo)
    assign(paste0("icc_batch_",suffix[i]),icc_batch)
  }
  
  if (i==4){
    icc_region=NULL
    for (p in 1:ncol(expo)){
      x=expo[,p]
      model=lmer(x~(1|covars$Region))
      vcov = as.data.frame(VarCorr(model))$vcov
      icc_region=c(icc_region,vcov[1]/sum(vcov))
    }
    names(icc_region) = colnames(expo)
    
    icc_depart=NULL
    for (p in 1:ncol(expo)){
      x=expo[,p]
      model=lmer(x~(1|covars$Department))
      vcov = as.data.frame(VarCorr(model))$vcov
      icc_depart=c(icc_depart,vcov[1]/sum(vcov))
    }
    names(icc_depart) = colnames(expo)
  }
}

# Manhattan plot (Family ID, Main)
values = t(bind_rows(icc_family_lux, icc_family_fra, icc_family_gs, icc_family_pooled3, icc_family_pooled2))
annot_sub = annot[rownames(values)]
annot_sub = annot_sub[order(annot_sub)]
values = values[names(annot_sub),]

xseq = seq(1, nrow(values))

{pdf(paste0("../Figures/Intra_class_correlation_univariate_expo_cont.pdf"), width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1], pch=19, las=1, xaxt="n", type = "n",
       ylim = c(0,1),
       xlab="", ylab="Intra-Class Correlation within families", cex.lab=1.2,
       panel.first=abline(v=xseq,lty=3,col="grey"),
       col=batch.colours[1])
  points(values[,1], pch = 19, col = batch.colours[1], cex = 0.8)
  points(values[,2], pch = 19, col = batch.colours[2], cex = 0.8)
  points(values[,3], pch = 19, col = batch.colours[3], cex = 0.8)
  points(values[,4], pch = 17, col = batch.colours[4], cex = 0.8)
  for (k in 1:length(xseq)){
    axis(side=1, at=xseq[k], labels=rownames(values)[k], cex.axis=0.8, las=2)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis=0.9)
  }
  legend("top", pch=c(rep(19,3),17),
         col=batch.colours[1:4],
         legend = batches[1:4],
         horiz = TRUE,
         bg="white", cex = 0.8)
  dev.off()
}

{pdf(paste0("../Figures/Intra_class_correlation_univariate_expo_cont_Pooled2.pdf"), width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1], pch=19, las=1, xaxt="n", type = "n",
       ylim = c(0,1),
       xlab="", ylab="Intra-Class Correlation within families", cex.lab=1.2,
       panel.first=abline(v=xseq,lty=3,col="grey"),
       col=batch.colours[1])
  points(values[,1], pch = 19, col = batch.colours[1], cex = 0.8)
  points(values[,3], pch = 19, col = batch.colours[3], cex = 0.8)
  points(values[,5], pch = 17, col = batch.colours[5], cex = 0.8)
  for (k in 1:length(xseq)){
    axis(side=1, at=xseq[k], labels=rownames(values)[k], cex.axis=0.8, las=2)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis=0.9)
  }
  legend("top", pch=c(rep(19,2),17),
         col=batch.colours[c(1,3,5)],
         legend = batches[c(1,3,5)],
         horiz = TRUE,
         bg="white", cex = 0.8)
  dev.off()
}

annot_sub=annot[names(icc_batch_pooled2)]
{pdf("../Figures/Pooled2/Intra_class_correlation_batch_univariate_expo_cont.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(icc_batch_pooled2, pch=19, cex=1, las=1, xaxt="n",
       xlab="", ylab="Intra-Class Correlation wtihin batch", cex.lab=1.2,
       panel.first=abline(v=1:length(icc_batch_pooled2),lty=3,col="grey"),
       col=annot.colours[annot_sub])
  for (k in 1:length(icc_batch_pooled2)){
    axis(side=1, at=k, labels=names(icc_batch_pooled2)[k], cex.axis=0.8, las=2)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(icc_batch_pooled2)+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis = 0.9)
  }
  dev.off()
}

annot_sub=annot[names(icc_batch_pooled3)]
{pdf("../Figures/Pooled3/Intra_class_correlation_batch_univariate_expo_cont.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(icc_batch_pooled3, pch=19, cex=1, las=1, xaxt="n",
       xlab="", ylab="Intra-Class Correlation wtihin batch", cex.lab=1.2,
       panel.first=abline(v=1:length(icc_batch_pooled3),lty=3,col="grey"),
       col=annot.colours[annot_sub])
  for (k in 1:length(icc_batch_pooled3)){
    axis(side=1, at=k, labels=names(icc_batch_pooled3)[k], cex.axis=0.8, las=3)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(icc_batch_pooled3)+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis = 0.9)
  }
  dev.off()
}

annot_sub=annot[names(icc_region)]
{pdf("../Figures/Pooled3/Intra_class_correlation_region_univariate_expo_cont.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(icc_region, pch=19, cex=1, las=1, xaxt="n",
       xlab="", ylab="Intra-Class Correlation wtihin region", cex.lab=1.2,
       panel.first=abline(v=1:length(icc_region),lty=3,col="grey"),
       col=annot.colours[annot_sub])
  for (k in 1:length(icc_region)){
    axis(side=1, at=k, labels=names(icc_region)[k], cex.axis=0.8, las=3)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(icc_region)+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis = 0.9)
  }
  dev.off()
}

annot_sub=annot[names(icc_depart)]
{pdf("../Figures/Pooled3/Intra_class_correlation_depart_univariate_expo_cont.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(icc_depart, pch=19, cex=1, las=1, xaxt="n",
       xlab="", ylab="Intra-Class Correlation wtihin department", cex.lab=1.2,
       panel.first=abline(v=1:length(icc_depart),lty=3,col="grey"),
       col=annot.colours[annot_sub])
  for (k in 1:length(icc_depart)){
    axis(side=1, at=k, labels=names(icc_depart)[k], cex.axis=0.8, las=3)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(icc_depart)+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis = 0.9)
  }
  dev.off()
}


