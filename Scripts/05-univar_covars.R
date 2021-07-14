## Univariate regression (including isolated children)
## Rin 12 July

# Load packages
library(tidyverse)
library(RColorBrewer)

### Exposure ~ Covariate ----
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
  if (i==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else {X = covars %>% select(Age, Gender)}
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
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
  
  annot_sub = annot[rownames(betas)]
  
  ## Volcano plot
  pdf(paste0("../Figures/",filepaths[i],"/Univariate_exposure_covariate.pdf"), width = 10, height = 5*ceiling(ncol(betas)/2))
  par(mar=c(5,5,1,1), mfrow = c(ceiling(ncol(betas)/2),2))
  for (n in 1:ncol(betas)){
    plot(betas[,n], -log10(pvals[,n]), pch=19,
         col=ifelse(pvals[,n] < 0.05/length(betas[,n]), annot.colours[annot_sub], "grey"),
         cex.lab=1, cex = 0.7,
         ylim = c(0, max(-log10(pvals[,n]))+0.25),
         xlim = c(-max(abs(betas[,n]))-0.4, max(abs(betas[,n]))+0.4),
         ylab=expression(-log[10](p-value)), 
         xlab=substitute(paste(beta,"(",tmp,")"), list(tmp = colnames(betas)[n])))
    text(betas[,n]+sign(betas[,n])*0.1, -log10(pvals[,n])+0.25,
         labels = ifelse(pvals[,n] < 0.05/length(betas[,n]), rownames(betas), ""),
         col = annot.colours[annot_sub])
    abline(h = -log10(0.05/length(betas[,n])), lty = 2)
  }
  dev.off()
  
  ifelse(dir.exists("../Results"),"",dir.create("../Results"))
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas.rds"))
  
  # Pooled analysis: by Batch, Region and Department
  pvals = NULL
  if (i %in% 4:5){
    if (i==4){X = covars %>% select(Batch, Region, Department)}
    if (i==5){X = covars %>% select(Batch)}
    ncol(expo) * ncol(X) # number of tests
    betas = pvals = NULL
    f1='expo[,k] ~ X[,j]'
    f0='expo[,k] ~ 1'
    t0=Sys.time()
    for (j in 1:ncol(X)){
      for (k in 1:ncol(expo)){
        model1=lm(as.formula(f1))
        model0=lm(as.formula(f0))
        pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
      }
    }
    t1=Sys.time()
    print(t1-t0)
    pvals = ifelse(pvals ==0, .Machine$double.xmin, pvals)
    pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
    rownames(pvals)=colnames(expo)
    colnames(pvals)=colnames(X)
    # If p-value = 0 replace with smallest double in R
    saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_geo_pvals.rds"))
    assign(paste0("pvals_",suffix[i]),pvals)
  }
}

annot_sub = annot[rownames(pvals_pooled3)]

# Manhattan plot
xseq = seq(1, nrow(pvals_pooled3))

for (j in colnames(pvals_pooled3)){
  {pdf(paste0("../Figures/",filepaths[4],"/Univariate_exposure_",j,".pdf"), width=14, height=8)
    par(mar=c(20,5,1,1))
    plot(-log10(pvals_pooled3[,j]),
         pch = 17, col=annot.colours[annot_sub],
         xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
         ylim = c(min(-log10(pvals_pooled3[,j]), na.rm = TRUE), max(-log10(pvals_pooled3[,j]), na.rm = TRUE)))
    abline(h = -log10(0.05), lty = 2, col = "grey")
    abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
    abline(h = -log10(0.05/nrow(pvals_pooled3)), lty = 2, col = darken(batch.colours[4], 0.5))
    for(i in 1:length(xseq)){
      axis(1, at=xseq[i], labels = rownames(pvals_pooled3)[i], las=2, cex.axis = 0.8,
           col.axis = ifelse(isTRUE(pvals_pooled3[i,j] < 0.05/nrow(pvals_pooled3)),
                             darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
    }
    xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals_pooled3))+0.5)
    axis(side=1, line=8, at=xgroup, labels=NA)
    tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
    for (k in 1:length(unique(annot_sub))){
      axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
           col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
    }
    legend("topleft", lty = c(rep(2,2)),
           col=c(darken(batch.colours[4], 0.5), "grey"),
           legend = c("Bonferroni threshold","Nominal threshold"),
           bg="white", cex = 0.7)
    dev.off()
  }
}

annot_sub = annot[rownames(pvals_pooled2)]

# Manhattan plot
xseq = seq(1, nrow(pvals_pooled2))

for (j in colnames(pvals_pooled2)){
  {pdf(paste0("../Figures/",filepaths[5],"/Univariate_exposure_",j,".pdf"), width=14, height=8)
    par(mar=c(20,5,1,1))
    plot(-log10(pvals_pooled2[,j]),
         pch = 17, col=annot.colours[annot_sub],
         xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
         ylim = c(min(-log10(pvals_pooled2[,j]), na.rm = TRUE), max(-log10(pvals_pooled2[,j]), na.rm = TRUE)))
    abline(h = -log10(0.05), lty = 2, col = "grey")
    abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
    abline(h = -log10(0.05/nrow(pvals_pooled2)), lty = 2, col = darken(batch.colours[5], 0.5))
    for(i in 1:length(xseq)){
      axis(1, at=xseq[i], labels = rownames(pvals_pooled2)[i], las=2, cex.axis = 0.8,
           col.axis = ifelse(isTRUE(pvals_pooled2[i,j] < 0.05/nrow(pvals_pooled2)),
                             darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
    }
    xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals_pooled2))+0.5)
    axis(side=1, line=8, at=xgroup, labels=NA)
    tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
    for (k in 1:length(unique(annot_sub))){
      axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
           col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
    }
    legend("topleft", lty = c(rep(2,2)),
           col=c(darken(batch.colours[5], 0.5), "grey"),
           legend = c("Bonferroni threshold","Nominal threshold"),
           bg="white", cex = 0.7)
    dev.off()
  }
}

### Detection ~ Covariate ----
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
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_raw.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  if (i==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else {X = covars %>% select(Age, Gender)}
  
  # Convert matrix into binary
  expo = ifelse(expo=="nd",1,0)
  # Exclude chemicals with 0% non-detects
  expo = expo[,which(colSums(expo,na.rm = TRUE)!=0)]
  saveRDS(expo, paste0("../Processed/",filepaths[i],"/Exposure_matrix_nd.rds"))
  
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
  t0=Sys.time()
  for (j in 1:ncol(X)){
    for (k in 1:ncol(expo)){
      model1=glm(as.formula(f1), family = "binomial")
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
  
  annot_sub = annot[rownames(betas)]
  
  ## Volcano plot
  pdf(paste0("../Figures/",filepaths[i],"/Univariate_detection_covariate.pdf"), width = 10, height = 5*ceiling(ncol(betas)/2))
  par(mar=c(5,5,1,1), mfrow = c(ceiling(ncol(betas)/2),2))
  for (n in 1:ncol(betas)){
    plot(betas[,n], -log10(pvals[,n]), pch=19,
         col=ifelse(pvals[,n] < 0.05/length(betas[,n]), annot.colours[annot_sub], "grey"),
         cex.lab=1, cex = 0.7,
         ylim = c(0, max(-log10(pvals[,n])+0.25)),
         xlim = c(-max(abs(betas[,n])), max(abs(betas[,n]))),
         ylab=expression(-log[10](p-value)), 
         xlab=substitute(paste(beta,"(",tmp,")"), list(tmp = colnames(betas)[n])))
    text(betas[,n]+sign(betas[,n])*0.5, -log10(pvals[,n])+0.25,
         labels = ifelse(pvals[,n] < 0.05/length(betas[,n]), rownames(betas), ""),
         col = annot.colours[annot_sub])
    abline(h = -log10(0.05/length(betas[,n])), lty = 2)
  }
  dev.off()
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_detection_covariate_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_detection_covariate_betas.rds"))
  
  # Pooled analysis: by Batch, Region and Department
  pvals = NULL
  if (i %in% 4:5){
    if (i==4){X = covars %>% select(Batch, Region, Department)}
    if (i==5){X = covars %>% select(Batch)}
    ncol(expo) * ncol(X) # number of tests
    pvals = NULL
    f1='expo[,k] ~ X[,j]'
    f0='expo[,k] ~ 1'
    t0=Sys.time()
    for (j in 1:ncol(X)){
      for (k in 1:ncol(expo)){
        model1=lm(as.formula(f1))
        model0=lm(as.formula(f0))
        pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
      }
    }
    t1=Sys.time()
    print(t1-t0)
    pvals = ifelse(pvals ==0, .Machine$double.xmin, pvals)
    pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
    rownames(pvals)=colnames(expo)
    colnames(pvals)=colnames(X)
    # If p-value = 0 replace with smallest double in R
    saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_detection_geo_pvals.rds"))
    assign(paste0("pvals_",suffix[i]),pvals)
  }
}

annot_sub = annot[rownames(pvals_pooled3)]

# Manhattan plot
xseq = seq(1, nrow(pvals_pooled3))

for (j in colnames(pvals_pooled3)){
  {pdf(paste0("../Figures/",filepaths[4],"/Univariate_detection_",j,".pdf"), width=14, height=8)
    par(mar=c(20,5,1,1))
    plot(-log10(pvals_pooled3[,j]),
         pch = 17, col=annot.colours[annot_sub],
         xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
         ylim = c(min(-log10(pvals_pooled3[,j]), na.rm = TRUE), max(-log10(pvals_pooled3[,j]), na.rm = TRUE)))
    abline(h = -log10(0.05), lty = 2, col = "grey")
    abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
    abline(h = -log10(0.05/nrow(pvals_pooled3)), lty = 2, col = darken(batch.colours[4], 0.5))
    for(i in 1:length(xseq)){
      axis(1, at=xseq[i], labels = rownames(pvals_pooled3)[i], las=2, cex.axis = 0.7,
           col.axis = ifelse(isTRUE(pvals_pooled3[i,j] < 0.05/nrow(pvals_pooled3)),
                             darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
    }
    xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals_pooled3))+0.5)
    axis(side=1, line=8, at=xgroup, labels=NA)
    tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
    for (k in 1:length(unique(annot_sub))){
      axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
           col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
    }
    legend("topleft", lty = c(rep(2,2)),
           col=c(darken(batch.colours[4], 0.5), "grey"),
           legend = c("Bonferroni threshold","Nominal threshold"),
           bg="white", cex = 0.7)
    dev.off()
  }
}

annot_sub = annot[rownames(pvals_pooled2)]

# Manhattan plot
xseq = seq(1, nrow(pvals_pooled2))

for (j in colnames(pvals_pooled2)){
  {pdf(paste0("../Figures/",filepaths[5],"/Univariate_detection_",j,".pdf"), width=14, height=8)
    par(mar=c(20,5,1,1))
    plot(-log10(pvals_pooled2[,j]),
         pch = 17, col=annot.colours[annot_sub],
         xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
         ylim = c(min(-log10(pvals_pooled2[,j]), na.rm = TRUE), max(-log10(pvals_pooled2[,j]), na.rm = TRUE)))
    abline(h = -log10(0.05), lty = 2, col = "grey")
    abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
    abline(h = -log10(0.05/nrow(pvals_pooled2)), lty = 2, col = darken(batch.colours[5], 0.5))
    for(i in 1:length(xseq)){
      axis(1, at=xseq[i], labels = rownames(pvals_pooled2)[i], las=2, cex.axis = 0.7,
           col.axis = ifelse(isTRUE(pvals_pooled2[i,j] < 0.05/nrow(pvals_pooled2)),
                             darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
    }
    xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals_pooled2))+0.5)
    axis(side=1, line=8, at=xgroup, labels=NA)
    tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
    for (k in 1:length(unique(annot_sub))){
      axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
           col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
    }
    legend("topleft", lty = c(rep(2,2)),
           col=c(darken(batch.colours[5], 0.5), "grey"),
           legend = c("Bonferroni threshold","Nominal threshold"),
           bg="white", cex = 0.7)
    dev.off()
  }
}
