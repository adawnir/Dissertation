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
  
  ## Volcano plot - age and gender
  pdf(paste0("../Figures/",filepaths[i],"/Univariate_exposure_age_gender.pdf"), width = 10, height = 5)
  par(mar=c(5,5,1,1), mfrow = c(1,2))
  for (n in 1:2){
    plot(betas[,n], -log10(pvals[,n]), pch=19,
         col=ifelse(pvals[,n] < 0.05/length(betas[,n]), annot.colours[annot_sub], "grey"),
         cex.lab=1.5, cex = 0.7,
         ylim = c(0, max(-log10(pvals[,n]))+0.25),
         xlim = c(-max(abs(betas[,n]))-0.4, max(abs(betas[,n]))+0.4),
         ylab=expression(-log[10](italic(p))), 
         xlab=substitute(paste(beta,"(",tmp,")"), list(tmp = colnames(betas)[n])))
    text(betas[,n]+sign(betas[,n])*0.1, -log10(pvals[,n])+0.25,
         labels = ifelse(pvals[,n] < 0.05/length(betas[,n]), rownames(betas), ""),
         col = annot.colours[annot_sub])
    abline(h = -log10(0.05/length(betas[,n])), lty = 2)
  }
  dev.off()
  
  ## Volcano plot - length and weight
  if (i==1){
    pdf(paste0("../Figures/",filepaths[i],"/Univariate_exposure_length_weight.pdf"), width = 10, height = 5)
    par(mar=c(5,5,1,1), mfrow = c(1,2))
    for (n in 3:4){
      plot(betas[,n], -log10(pvals[,n]), pch=19,
           col=ifelse(pvals[,n] < 0.05/length(betas[,n]), annot.colours[annot_sub], "grey"),
           cex.lab=1.5, cex = 0.7,
           ylim = c(0, max(-log10(pvals[,n]))+0.25),
           xlim = c(-max(abs(betas[,n]))-0.4, max(abs(betas[,n]))+0.4),
           ylab=expression(-log[10](italic(p))), 
           xlab=substitute(paste(beta,"(",tmp,")"), list(tmp = colnames(betas)[n])))
      text(betas[,n]+sign(betas[,n])*0.1, -log10(pvals[,n])+0.25,
           labels = ifelse(pvals[,n] < 0.05/length(betas[,n]), rownames(betas), ""),
           col = annot.colours[annot_sub])
      abline(h = -log10(0.05/length(betas[,n])), lty = 2)
    }
    dev.off() 
  }
  
  ifelse(dir.exists("../Results"),"",dir.create("../Results"))
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas.rds"))
  
  # # Pooled analysis: by Batch, Region and Department
  # pvals = NULL
  # if (i %in% 4:5){
  #   if (i==4){X = covars %>% select(Batch, Region, Department)}
  #   if (i==5){X = covars %>% select(Batch)}
  #   ncol(expo) * ncol(X) # number of tests
  #   betas = pvals = NULL
  #   f1='expo[,k] ~ X[,j]'
  #   f0='expo[,k] ~ 1'
  #   t0=Sys.time()
  #   for (j in 1:ncol(X)){
  #     for (k in 1:ncol(expo)){
  #       model1=lm(as.formula(f1))
  #       model0=lm(as.formula(f0))
  #       pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
  #     }
  #   }
  #   t1=Sys.time()
  #   print(t1-t0)
  #   pvals = ifelse(pvals ==0, .Machine$double.xmin, pvals)
  #   pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
  #   rownames(pvals)=colnames(expo)
  #   colnames(pvals)=colnames(X)
  #   # If p-value = 0 replace with smallest double in R
  #   saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_geo_pvals.rds"))
  #   # assign(paste0("pvals_",suffix[i]),pvals)
  # }
}

# annot_sub = annot[rownames(pvals_pooled3)]
# 
# # Manhattan plot
# xseq = seq(1, nrow(pvals_pooled3))
# 
# for (j in colnames(pvals_pooled3)){
#   {pdf(paste0("../Figures/",filepaths[4],"/Univariate_exposure_",j,".pdf"), width=14, height=8)
#     par(mar=c(20,5,1,1))
#     plot(-log10(pvals_pooled3[,j]),
#          pch = 17, col=annot.colours[annot_sub],
#          xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
#          ylim = c(min(-log10(pvals_pooled3[,j]), na.rm = TRUE), max(-log10(pvals_pooled3[,j]), na.rm = TRUE)))
#     abline(h = -log10(0.05), lty = 2, col = "grey")
#     abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
#     abline(h = -log10(0.05/nrow(pvals_pooled3)), lty = 2, col = darken(batch.colours[4], 0.5))
#     for(i in 1:length(xseq)){
#       axis(1, at=xseq[i], labels = rownames(pvals_pooled3)[i], las=2, cex.axis = 0.8,
#            col.axis = ifelse(isTRUE(pvals_pooled3[i,j] < 0.05/nrow(pvals_pooled3)),
#                              darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
#     }
#     xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals_pooled3))+0.5)
#     axis(side=1, line=8, at=xgroup, labels=NA)
#     tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#     for (k in 1:length(unique(annot_sub))){
#       axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#            col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#     }
#     legend("top", lty = c(rep(2,2)),
#            col=c(darken(batch.colours[4], 0.5), "grey"),
#            legend = c("Bonferroni threshold","Nominal threshold"),
#            bg="white", cex = 0.8, horiz = TRUE)
#     dev.off()
#   }
# }

# annot_sub = annot[rownames(pvals_pooled2)]
# 
# # Manhattan plot
# xseq = seq(1, nrow(pvals_pooled2))
# 
# for (j in colnames(pvals_pooled2)){
#   {pdf(paste0("../Figures/",filepaths[5],"/Univariate_exposure_",j,".pdf"), width=14, height=8)
#     par(mar=c(20,5,1,1))
#     plot(-log10(pvals_pooled2[,j]),
#          pch = 17, col=annot.colours[annot_sub],
#          xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
#          ylim = c(min(-log10(pvals_pooled2[,j]), na.rm = TRUE), max(-log10(pvals_pooled2[,j]), na.rm = TRUE)))
#     abline(h = -log10(0.05), lty = 2, col = "grey")
#     abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
#     abline(h = -log10(0.05/nrow(pvals_pooled2)), lty = 2, col = darken(batch.colours[5], 0.5))
#     for(i in 1:length(xseq)){
#       axis(1, at=xseq[i], labels = rownames(pvals_pooled2)[i], las=2, cex.axis = 0.8,
#            col.axis = ifelse(isTRUE(pvals_pooled2[i,j] < 0.05/nrow(pvals_pooled2)),
#                              darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
#     }
#     xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(pvals_pooled2))+0.5)
#     axis(side=1, line=8, at=xgroup, labels=NA)
#     tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#     for (k in 1:length(unique(annot_sub))){
#       axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#            col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#     }
#     legend("top", lty = c(rep(2,2)),
#            col=c(darken(batch.colours[5], 0.5), "grey"),
#            legend = c("Bonferroni threshold","Nominal threshold"),
#            bg="white", cex = 0.8, horiz = TRUE)
#     dev.off()
#   }
# }

### Exposure ~ Family ID ----
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

values = -log10(pvals_pooled3)
annot_sub = annot[names(values)]

# Manhattan plot
bonf = 0.05/length(pvals_pooled3)
bonf = -log10(bonf)

xseq = seq(1, length(values))

{pdf("../Figures/Pooled3/Univariate_exposure_family.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values,
       col=annot.colours[annot_sub],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="p", ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)+1), pch = 19)
  abline(h = -log10(0.05), lty = 2, col = "grey")
  abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
  abline(h = bonf, lty = 2, col = darken(batch.colours[4], 0.5))
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = names(values)[i], las=2, cex.axis = 0.8,
         col.axis = ifelse(isTRUE(values[i] > bonf),
                           darken(annot.colours[annot_sub[i]], amount=0.5),"black"))
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(names(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("top", lty = c(rep(2,2)),
         col=c(darken(batch.colours[4], 0.5), "grey"),
         legend = c("Bonferroni threshold (Pooled)",
                    "Nominal threshold"),
         bg="white", cex = 0.8, horiz = TRUE)
  dev.off()
}

# values = -log10(t(bind_rows(pvals_lux, pvals_fra, pvals_gs, pvals_pooled3, pvals_pooled2)))
# annot_sub = annot[rownames(values)]
# annot_sub = annot_sub[order(annot_sub)]
# values = values[names(annot_sub),]
# 
# # Manhattan plot
# bonf = c(0.05/length(pvals_lux), 0.05/length(pvals_fra), 0.05/length(pvals_gs),
#          0.05/length(pvals_pooled3), 0.05/length(pvals_pooled2))
# bonf = -log10(bonf)
# 
# xseq = seq(1, nrow(values))
# 
# options(scipen=999)
# {pdf("../Figures/Univariate_exposure_family.pdf", width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(values[,1],
#        col=batch.colours[1],
#        xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
#        type="n", ylim = c(min(values, na.rm = TRUE), 500),
#        log = "y")
#   abline(h = -log10(0.05), lty = 2, col = "grey")
#   abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
#   for (i in 1:4){
#     if (i == 4){
#       pch = 17
#     } else {pch = 19}
#     points(values[,i], pch = pch, col = batch.colours[i], cex = 0.8)
#     abline(h = bonf[i], lty = 2, col = darken(batch.colours[i], 0.5))
#   }
#   for(i in 1:length(xseq)){
#     axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.8,
#          col.axis = ifelse(isTRUE(values[i,4] > bonf[4]),
#                            darken(annot.colours[annot_sub[i]], amount=0.5),"black"))
#   }
#   xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
#   axis(side=1, line=8, at=xgroup, labels=NA)
#   tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#   }
#   legend("top", pch=c(rep(19,3),17, rep(NA,5)), lty = c(rep(NA,4),rep(2,5)),
#          col=c(batch.colours[1:4],darken(batch.colours[1:4], 0.5), "grey"),
#          legend = c(batches[1:4],
#                     paste0("Bonferroni threshold (",c("LUX","FRA","GS","Pooled"),")"),
#                     "Nominal threshold"),
#          bg="white", cex = 0.8, ncol = 5)
#   dev.off()
# }
# 
# {pdf("../Figures/Univariate_exposure_family_Pooled2.pdf", width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(values[,1],
#        col=batch.colours[1],
#        xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
#        type="n", ylim = c(min(values, na.rm = TRUE), 500),
#        log = "y")
#   abline(h = -log10(0.05), lty = 2, col = "grey")
#   abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
#   for (i in c(1,3,5)){
#     if (i == 5){
#       pch = 17
#     } else {pch = 19}
#     points(values[,i], pch = pch, col = batch.colours[i], cex = 0.8)
#     abline(h = bonf[i], lty = 2, col = darken(batch.colours[i], 0.5))
#   }
#   for(i in 1:length(xseq)){
#     axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.8,
#          col.axis = ifelse(isTRUE(values[i,5] > bonf[5]),
#                            darken(annot.colours[annot_sub[i]], amount=0.5),"black"))
#   }
#   xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
#   axis(side=1, line=8, at=xgroup, labels=NA)
#   tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#   }
#   legend("top", pch=c(rep(19,2),17, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
#          col=c(batch.colours[c(1,3,5)],darken(batch.colours[c(1,3,5)], 0.5), "grey"),
#          legend = c(batches[c(1,3,5)],
#                     paste0("Bonferroni threshold (",c("LUX","GS","Pooled"),")"),
#                     "Nominal threshold"),
#          bg="white", cex = 0.8, ncol = 4)
#   dev.off()
# }
# options(scipen=0)

### Exposure ~ Paris/Ile-de-France ----
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
for (i in c(2,4)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  X = covars %>%
      select(Department, Region) %>%
      mutate(Paris = ifelse(.$Department=="Paris",1,0),
             `Ile-de-France` = ifelse(.$Region=="Île-de-France",1,0)) %>%
    select(-Department, -Region)
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
  
  ## Volcano plot - age and gender
  pdf(paste0("../Figures/",filepaths[i],"/Univariate_exposure_Paris_IDF.pdf"), width = 10, height = 5)
  par(mar=c(5,5,1,1), mfrow = c(1,2))
  for (n in 1:2){
    plot(betas[,n], -log10(pvals[,n]), pch=19,
         col=ifelse(pvals[,n] < 0.05/length(betas[,n]), annot.colours[annot_sub], "grey"),
         cex.lab=1.5, cex = 0.7,
         ylim = c(0, max(-log10(pvals[,n]))+0.25),
         xlim = c(-max(abs(betas[,n]))-0.5, max(abs(betas[,n]))+0.5),
         ylab=expression(-log[10](italic(p))), 
         xlab=substitute(paste(beta,"(",tmp,")"), list(tmp = colnames(betas)[n])))
    text(betas[,n]+sign(betas[,n])*0.1, -log10(pvals[,n])+0.25,
         labels = ifelse(pvals[,n] < 0.05/length(betas[,n]), rownames(betas), ""),
         col = annot.colours[annot_sub])
    abline(h = -log10(0.05/length(betas[,n])), lty = 2)
  }
  dev.off()
  
  ifelse(dir.exists("../Results"),"",dir.create("../Results"))
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_Paris_IDF_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_exposure_Paris_IDF_betas.rds"))
  assign(paste0("pvals_",suffix[i]),pvals)
}
