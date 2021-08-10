## Pre-processed data overview
## Rin Wada 8 July

# Load packages
library(tidyverse)
library(colorspace)
library(RColorBrewer)

### Load data ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

covars_lux = readRDS(paste0("../Processed/",filepaths[1],"/Participant_covariate_info_thresh.rds"))
chem_lux = readRDS(paste0("../Processed/",filepaths[1],"/Chemical_compound_info.rds"))
expo_lux = readRDS(paste0("../Processed/",filepaths[1],"/Exposure_matrix_raw.rds"))

covars_fra = readRDS(paste0("../Processed/",filepaths[2],"/Participant_covariate_info_thresh.rds"))
chem_fra = readRDS(paste0("../Processed/",filepaths[2],"/Chemical_compound_info.rds"))
expo_fra = readRDS(paste0("../Processed/",filepaths[2],"/Exposure_matrix_raw.rds"))

covars_gs = readRDS(paste0("../Processed/",filepaths[3],"/Participant_covariate_info_thresh.rds"))
chem_gs = readRDS(paste0("../Processed/",filepaths[3],"/Chemical_compound_info.rds"))
expo_gs = readRDS(paste0("../Processed/",filepaths[3],"/Exposure_matrix_raw.rds"))

covars_pooled3 = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))
chem_pooled3 = readRDS(paste0("../Processed/",filepaths[4],"/Chemical_compound_info.rds"))
expo_pooled3 = readRDS(paste0("../Processed/",filepaths[4],"/Exposure_matrix_raw.rds"))

covars_pooled2 = readRDS(paste0("../Processed/",filepaths[5],"/Participant_covariate_info_thresh.rds"))
chem_pooled2 = readRDS(paste0("../Processed/",filepaths[5],"/Chemical_compound_info.rds"))
expo_pooled2 = readRDS(paste0("../Processed/",filepaths[5],"/Exposure_matrix_raw.rds"))

### Proportion detected ----
tmp1 = chem_lux$detect_rate
names(tmp1) = chem_lux$Compound
tmp2 = chem_fra$detect_rate
names(tmp2) = chem_fra$Compound
tmp3 = chem_gs$detect_rate
names(tmp3) = chem_gs$Compound

grep = intersect(chem_lux$Compound[which(chem_lux$detect_rate!=0)],
                 chem_fra$Compound[which(chem_fra$detect_rate!=0)]) %>%
  intersect(chem_gs$Compound[which(chem_gs$detect_rate!=0)])

tmp4 = chem_pooled3$detect_rate
names(tmp4) = chem_pooled3$Compound
tmp4 = tmp4[grep]

grep = intersect(chem_lux$Compound[which(chem_lux$detect_rate!=0)],
                 chem_gs$Compound[which(chem_gs$detect_rate!=0)])
tmp5 = chem_pooled2$detect_rate
names(tmp5) = chem_pooled2$Compound
tmp5 = tmp5[grep]

prop = list(tmp1, tmp2, tmp3, tmp4, tmp5)
sorted_prop = lapply(prop, sort)
N=list(NULL, NULL, NULL, NULL, NULL)
thrseq=seq(0,1,by=0.05)
for (i in 1:length(sorted_prop)){
  for (k in thrseq){
    N[[i]]=c(N[[i]],sum(sorted_prop[[i]]>k))
  }
}

ifelse(dir.exists("../Figures"), "", dir.create("../Figures"))
{pdf("../Figures/Thresholds_proportion_detection.pdf", width=10, height=7)
  par(mar=c(5,5,1,1))
  plot(thrseq, N[[1]], type="n", pch=19, col=batch.colours[1], las=1, cex.lab=1.5,
       xlab="Threshold in proportion of detection",
       ylab="Number of kept variables", ylim = c(0,152))
  pch = c(rep(19,3),rep(17,2))
  lty = c(rep(1,3),rep(2,2))
  for (k in 1:3){
    points(thrseq, N[[k]], type="b", pch=pch[k], lty = lty[k], col=batch.colours[k])
  }
  abline(v = 0.1, lty = 2)
  legend("top", lty=lty, pch = pch, lwd=2, col=batch.colours[1:3],
         legend = batches[1:3], bg="white", horiz = TRUE)
  dev.off()
}

### Overlay ----
prop = t(bind_rows(tmp1, tmp2, tmp3))
prop[which(is.na(prop))] = 0 # Replace NA with 0
prop = prop[order(factor.order(annot[rownames(prop)]),-prop[,1],-prop[,2],-prop[,3]),]
mylabels = rownames(prop)
prop = cbind(prop[,1],rep(NA,nrow(prop)),prop[,2],rep(NA,nrow(prop)),prop[,3],rep(NA,nrow(prop)),rep(NA,nrow(prop)))
prop = as.vector(t(prop))
prop=prop[c(length(prop),1:length(prop)-1)]

background = TRUE
myspacing = 7
xseq = seq(myspacing/2,length(mylabels)*myspacing, by=myspacing)
annot_sub = annot[mylabels]
{pdf("../Figures/Detection_prop.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(prop,
       col=c(NA, batch.colours[1],NA, batch.colours[2], NA, batch.colours[3], NA),
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1.1))
  xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
  if (background){
    for (k in seq(1,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(annot.colours[annot_sub[k]],0.95), border=NA)
    }
    for (k in seq(2,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(annot.colours[annot_sub[k]],0.99), border=NA)
    }
    box()
  }
  abline(v=xseqgreysep,lty=1,lwd=0.1,col="grey")
  par(new = TRUE)
  plot(prop,
       col=c(NA, batch.colours[1],NA, batch.colours[2], NA, batch.colours[3], NA),
       xaxt="n", ylab="Proportion detected", xlab = "", cex.lab=1.5,
       type="h", lwd=2, ylim = c(0,1.1))
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.6)
  }
  xseqblack=c(xseq[!duplicated(annot_sub)]-myspacing/2, max(xseq)+myspacing/2)
  abline(v=xseqblack,lty=3,col="black")
  axis(side=1, line=7, at=xseqblack, labels=NA)
  tmp=apply(rbind(xseqblack[-length(xseqblack)],xseqblack[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=7, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("top", lty=1, lwd=2, col=batch.colours[1:3], cex = 0.7,
         legend = batches[1:3], bg="white", horiz = TRUE)
  dev.off()
}

{pdf("../Figures/Detection_prop_thresh.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(prop,
       col=c(NA, batch.colours[1],NA, batch.colours[2], NA, batch.colours[3], NA),
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1.1))
  xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
  if (background){
    for (k in seq(1,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(annot.colours[annot_sub[k]],0.95), border=NA)
    }
    for (k in seq(2,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(annot.colours[annot_sub[k]],0.99), border=NA)
    }
    box()
  }
  abline(v=xseqgreysep,lty=1,lwd=0.1,col="grey")
  par(new = TRUE)
  plot(prop,
       col=ifelse(prop>0.1,
                  c(NA, batch.colours[1],NA, batch.colours[2], NA, batch.colours[3], NA),
                  alpha(c(NA, batch.colours[1],NA, batch.colours[2],NA, batch.colours[3], NA), 0.5)),
       xaxt="n", ylab="Proportion detected", xlab = "", cex.lab=1.5,
       type="h", lwd=2, ylim = c(0,1.1))
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.6,
         col.axis = ifelse(mylabels[i] %in% names(tmp4) & tmp4[mylabels[i]] > 0.1,
                           darken(annot.colours[(annot_sub)[i]], amount=0.5),
                           ifelse(mylabels[i] %in% names(tmp5) & tmp5[mylabels[i]] > 0.1,
                                  alpha(darken(annot.colours[(annot_sub)[i]], amount=0.5),0.5),
                                  "black")))
  }
  xseqblack=c(xseq[!duplicated(annot_sub)]-myspacing/2, max(xseq)+myspacing/2)
  abline(v=xseqblack,lty=3,col="black")
  abline(h = 0.1, lty = 2)
  axis(side=1, line=7, at=xseqblack, labels=NA)
  tmp=apply(rbind(xseqblack[-length(xseqblack)],xseqblack[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=7, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("top", lty=1, lwd=2, col=batch.colours[1:3], cex = 0.7,
         legend = batches[1:3], bg="white", horiz = TRUE)
  dev.off()
}

# ### Multi-panel ----
# prop = t(bind_rows(tmp1, tmp2, tmp3, tmp4))
# prop[which(is.na(prop))] = 0 # Replace NA with 0
# prop = prop[names(annot),]
# mylabels = rownames(prop)
# 
# myspacing = 1
# xseq = 1:nrow(prop)
# annot_sub = annot[mylabels]
# 
# {pdf("../Figures/Detection_prop.pdf", width=14, height=11)
#   par(oma=c(21,0,0,0),mar=c(0,5,1,1), mfrow = c(3,1))
#   for (i in 1:3){
#     plot(prop[,i],
#          col=annot.colours[annot_sub],
#          xaxt="n", ylab=paste0("Proportion detected"," (", batches[i] ,")"), xlab = "", cex.lab=1.2,
#          type="h", lwd=2, ylim = c(0,1.1))
#   }
#   for(i in 1:length(xseq)){
#     axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.8)
#   }
#   xseqblack=c(xseq[!duplicated(annot_sub)]-myspacing/2, max(xseq)+myspacing/2)
#   axis(side=1, line=9, at=xseqblack, labels=NA)
#   tmp=apply(rbind(xseqblack[-length(xseqblack)],xseqblack[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=9, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#   }
#   dev.off()
# }
# 
# {pdf("../Figures/Detection_prop_thresh.pdf", width=14, height=11)
#   par(oma=c(21,0,0,0),mar=c(0,5,1,1), mfrow = c(3,1))
#   for (i in 1:3){
#     plot(prop[,i],
#          col=ifelse(prop>0.1,annot.colours[annot_sub],alpha(annot.colours[annot_sub], 0.3)),
#          xaxt="n", ylab=paste0("Proportion detected"," (", batches[i] ,")"), xlab = "", cex.lab=1.2,
#          type="h", lwd=2, ylim = c(0,1.1))
#     abline(h = 0.1, lty = 2)
#   }
#   for(i in 1:length(xseq)){
#     axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.8,
#          col.axis = ifelse(mylabels[i] %in% names(tmp4) & tmp4[mylabels[i]] > 0.1,
#                            darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
#   }
#   xseqblack=c(xseq[!duplicated(annot_sub)]-myspacing/2, max(xseq)+myspacing/2)
#   axis(side=1, line=9, at=xseqblack, labels=NA)
#   tmp=apply(rbind(xseqblack[-length(xseqblack)],xseqblack[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=9, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#   }
#   dev.off()
# }
# 


# ### Age distribution ----
# {pdf("../Figures/Age_dist_thresh.pdf", width = 12, height = 8)
#   par(mfrow=c(2,3), mar = c(5,5,2,1))
#   hist(covars_lux$Age, main = batches[1], xlab = "Age (years)", col = batch.colours[1])
#   hist(covars_fra$Age, main = batches[2], xlab = "Age (years)", col = batch.colours[2])
#   hist(covars_gs$Age, main = batches[3], xlab = "Age (years)", col = batch.colours[3])
#   hist(covars_pooled3$Age, main = batches[4], xlab = "Age (years)", col = batch.colours[4])
#   hist(covars_pooled2$Age, main = batches[5], xlab = "Age (years)", col = batch.colours[5])
#   dev.off()
# }
# 
# ### Hair length and weight ----
# # Luxembourg only
# {pdf(paste0("../Figures/",filepaths[1],"/Sample_length_weight_dist_thresh.pdf"), width=8, height=4)
#   par(mfrow = c(1,2), mar = c(5,5,1,1))
#   hist(covars_lux$Length,
#        xlab = "Length of hair sample (cm)", main = NULL, col = batch.colours[1])
#   hist(covars_lux$Weight,
#        xlab = "Weight of hair sample (mg)", main = NULL, col = batch.colours[1])
#   dev.off()
# }
# 
# ### Proportion missing ----
# # Proportion missing
# tmp1 = chem_lux$NA_prop
# names(tmp1) = chem_lux$Compound
# tmp2 = chem_fra$NA_prop
# names(tmp2) = chem_fra$Compound
# tmp3 = chem_gs$NA_prop
# names(tmp3) = chem_gs$Compound
# tmp4 = chem_pooled3$NA_prop
# names(tmp4) = chem_pooled3$Compound
# tmp5 = chem_pooled2$NA_prop
# names(tmp5) = chem_pooled2$Compound
# 
# prop = t(bind_rows(tmp1, tmp2, tmp3, tmp4, tmp5))
# prop[which(is.na(prop))] = 1 # Replace NA with 1
# prop = prop[order(match(rownames(prop), names(annot))),]
# prop = cbind(prop,rep(NA,nrow(prop)),rep(NA,nrow(prop)))
# mylabels = rownames(prop)
# prop = as.vector(t(prop))
# prop=prop[c(length(prop),1:length(prop)-1)]
# 
# background = TRUE
# myspacing = 7
# xseq = seq(myspacing/2,length(mylabels)*myspacing, by=myspacing)
# annot_sub = annot[mylabels]
# {pdf("../Figures/Missing_prop.pdf", width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(prop,
#        col=c(NA, batch.colours, NA),
#        xaxt="n", ylab="", xlab = "",
#        type="n", lwd=2, ylim = c(0,1.1))
#   xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
#   if (background){
#     for (k in seq(1,length(xseqgreysep),by=2)){
#       polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
#               y=c(-10,-10,10,10), col=lighten(annot.colours[annot_sub[k]],0.95), border=NA)
#     }
#     for (k in seq(2,length(xseqgreysep),by=2)){
#       polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
#               y=c(-10,-10,10,10), col=lighten(annot.colours[annot_sub[k]],0.99), border=NA)
#     }
#     box()
#   }
#   abline(v=xseqgreysep,lty=1,lwd=0.1,col="grey")
#   par(new = TRUE)
#   plot(prop,
#        col=c(NA, batch.colours, NA),
#        xaxt="n", ylab="Proportion missing", xlab = "", cex.lab=1.5,
#        type="h", lwd=1, ylim = c(0,1.1))
#   for(i in 1:length(xseq)){
#     axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.6)
#   }
#   xseqblack=c(xseq[!duplicated(annot_sub)]-myspacing/2, max(xseq)+myspacing/2)
#   abline(v=xseqblack,lty=3,col="black")
#   axis(side=1, line=7, at=xseqblack, labels=NA)
#   tmp=apply(rbind(xseqblack[-length(xseqblack)],xseqblack[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=7, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
#   }
#   legend("top", lty=1, lwd=2, col=batch.colours, cex = 0.7,
#          legend = batches, bg="white", horiz = TRUE)
#   dev.off()
# }
# 