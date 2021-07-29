## Univariate regression and within family variation
## Rin 29 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)

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
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_detection_family_pvals.rds"))
  assign(paste0("pvals_",suffix[i]),pvals)
}
values = -log10(pvals_pooled3)
annot_sub = annot[names(values)]

# Manhattan plot
bonf = 0.05/length(pvals_pooled3)
bonf = -log10(bonf)

xseq = seq(1, length(values))

{pdf("../Figures/Pooled3/Univariate_detection_family.pdf", width=14, height=8)
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
# annot_sub = factor.order(annot)[rownames(values)]
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
# 
# {pdf("../Figures/Univariate_detection_family.pdf", width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(values[,1],
#        col=batch.colours[1],
#        xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
#        type="n", ylim = c(min(values, na.rm = TRUE), 13))
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

# {pdf("../Figures/Univariate_detection_family_Pooled2.pdf", width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(values[,1],
#        col=batch.colours[1],
#        xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
#        type="n", ylim = c(min(values, na.rm = TRUE), 8))
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
