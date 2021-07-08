## Pre-processed data overview
## Rin Wada 8 July

# Load packages
library(openxlsx)
library(tidyverse)

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
chem_lux = readRDS(paste0("../Data/",filepaths[1],"/Chemical_compound_info.rds"))
expo_lux = readRDS(paste0("../Data/",filepaths[1],"/Exposure_matrix_raw.rds"))

covars_fra = readRDS(paste0("../Processed/",filepaths[2],"/Participant_covariate_info_thresh.rds"))
chem_fra = readRDS(paste0("../Data/",filepaths[2],"/Chemical_compound_info.rds"))
expo_fra = readRDS(paste0("../Data/",filepaths[2],"/Exposure_matrix_raw.rds"))

covars_gs = readRDS(paste0("../Processed/",filepaths[3],"/Participant_covariate_info_thresh.rds"))
chem_gs = readRDS(paste0("../Data/",filepaths[3],"/Chemical_compound_info.rds"))
expo_gs = readRDS(paste0("../Data/",filepaths[3],"/Exposure_matrix_raw.rds"))

{pdf("../Figures/Age_dist_thresh.pdf", width = 12, height = 4)
  par(mfrow=c(1,3), mar = c(5,5,2,1))
  hist(covars_lux$Age, main = batches[1], xlab = "Age (years)", col = batch.colours[1])
  hist(covars_fra$Age, main = batches[2], xlab = "Age (years)", col = batch.colours[2])
  hist(covars_gs$Age, main = batches[3], xlab = "Age (years)", col = batch.colours[3])
  dev.off()
}

{pdf(paste0("../Figures/",filepaths[1],"/Sample_length_weight_dist_thresh.pdf"), width=8, height=4)
  par(mfrow = c(1,2), mar = c(5,5,1,1))
  hist(covars_lux$Length,
       xlab = "Length of hair sample (cm)", main = NULL, col = batch.colours[1])
  hist(covars_lux$Weight,
       xlab = "Weight of hair sample (mg)", main = NULL, col = batch.colours[1])
  dev.off()
}

# Total number of families
sum(covars_lux$Family.ID=="Isolated") + sum(table(table(covars_lux$Family.ID[which(covars_lux$Family.ID!="Isolated")])))
sum(covars_fra$Family.ID=="Isolated") + sum(table(table(covars_fra$Family.ID[which(covars_fra$Family.ID!="Isolated")])))
sum(covars_gs$Family.ID=="Isolated") + sum(table(table(covars_gs$Family.ID[which(covars_gs$Family.ID!="Isolated")])))

# Number of each type of family
table(table(covars_lux$Family.ID))
table(table(covars_fra$Family.ID))
table(table(covars_gs$Family.ID))

# Age (Mean, SD)
cat(round(mean(covars_lux$Age, na.rm = T),2),
    " (", round(sd(covars_lux$Age, na.rm = T),2), ")\n", sep = "")
cat(round(mean(covars_fra$Age, na.rm = T),2),
    " (", round(sd(covars_fra$Age, na.rm = T),2), ")\n", sep = "")
cat(round(mean(covars_gs$Age, na.rm = T),2),
    " (", round(sd(covars_gs$Age, na.rm = T),2), ")\n", sep = "")

# Age (NA %)
cat(sum(is.na(covars_lux$Age)),
    " (", round(sum(is.na(covars_lux$Age))/length(covars_lux$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_fra$Age)),
    " (", round(sum(is.na(covars_fra$Age))/length(covars_fra$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_gs$Age)),
    " (", round(sum(is.na(covars_gs$Age))/length(covars_gs$Age)*100,2), "%)\n", sep = "")

# Gender (N, %)
table(covars_lux$Gender, useNA = "ifany")
round(prop.table(table(covars_lux$Gender, useNA = "ifany"))*100,2)

table(covars_fra$Gender, useNA = "ifany")
round(prop.table(table(covars_fra$Gender, useNA = "ifany"))*100,2)

table(covars_gs$Gender, useNA = "ifany")
round(prop.table(table(covars_gs$Gender, useNA = "ifany"))*100,2)

# Gender (NA %)
cat(sum(is.na(covars_lux$Gender)),
    " (", round(sum(is.na(covars_lux$Gender))/length(covars_lux$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_fra$Gender)),
    " (", round(sum(is.na(covars_fra$Gender))/length(covars_fra$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_gs$Gender)),
    " (", round(sum(is.na(covars_gs$Gender))/length(covars_gs$Gender)*100,2), "%)\n", sep = "")

tmp = chem_lux$NA_prop
names(tmp) = chem_lux$Compound
tmp2 = chem_fra$NA_prop
names(tmp2) = chem_fra$Compound
tmp3 = chem_gs$NA_prop
names(tmp3) = chem_gs$Compound

# Proportion observed (not missing)
prop = 1-t(bind_rows(tmp, tmp2, tmp3))
prop[which(is.na(prop))] = 0 # Replace NA with 0
prop = prop[order(annot,-prop[,1],-prop[,2],-prop[,3]),]
prop = cbind(prop,rep(NA,nrow(prop)),rep(NA,nrow(prop)))
mylabels = rownames(prop)
prop = as.vector(t(prop))
prop=prop[c(length(prop),1:length(prop)-1)]

background = TRUE
myspacing = 5
xseq = seq(myspacing/2,length(mylabels)*myspacing, by=myspacing)
annot_sub = annot[mylabels]
{pdf("../Figures/Observed_prop.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(prop,
       col=c(NA, batch.colours[1:3], NA),
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1))
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
       col=c(NA, batch.colours[1:3], NA),
       xaxt="n", ylab="Proportion observed", xlab = "", cex.lab=1.5,
       type="h", lwd=2, ylim = c(0,1))
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
  legend("right", lty=1, lwd=2, col=batch.colours[1:3], cex = 0.7,
         legend = batches[1:3], bg="white")
  dev.off()
}

# Proportion NOT detected
tmp = chem_lux$nd_prop + chem_lux$NA_prop
names(tmp) = chem_lux$Compound
tmp2 = chem_fra$nd_prop + chem_fra$NA_prop
names(tmp2) = chem_fra$Compound
tmp3 = chem_gs$nd_prop + chem_gs$NA_prop
names(tmp3) = chem_gs$Compound

prop = list(tmp, tmp2, tmp3)
sorted_prop = lapply(prop, sort)
N=list(NULL, NULL, NULL)
thrseq=seq(0,1,by=0.05)
for (i in 1:length(sorted_prop)){
  for (k in thrseq){
    N[[i]]=c(N[[i]],sum(sorted_prop[[i]]>k))
  }
}

{pdf("../Figures/Thresholds_proportion_detection.pdf", width=10, height=7)
  par(mar=c(5,5,1,1))
  plot(thrseq, N[[1]], type="b", pch=19, col=batch.colours[1], las=1, cex.lab=1.5,
       xlab="Threshold in proportion of detection",
       ylab="Number of kept variables", ylim = c(0,152))
  points(thrseq, N[[2]], type="b", pch=19, col=batch.colours[2])
  points(thrseq, N[[3]], type="b", pch=19, col=batch.colours[3])
  abline(v = 0.1, lty = 2)
  legend("top", lty=1, lwd=2, col=batch.colours[1:3], ncol = 3,
         legend = batches[1:3], cex=0.8, bg="white")
  dev.off()
}

# Proportion detected
prop = 1-t(bind_rows(tmp, tmp2, tmp3))
prop[which(is.na(prop))] = 0 # Replace NA with 0
prop = prop[order(annot,-prop[,1],-prop[,2],-prop[,3]),]
prop = cbind(prop,rep(NA,nrow(prop)),rep(NA,nrow(prop)))
mylabels = rownames(prop)
prop = as.vector(t(prop))
prop=prop[c(length(prop),1:length(prop)-1)]

background = TRUE
myspacing = 5
xseq = seq(myspacing/2,length(mylabels)*myspacing, by=myspacing)
annot_sub = annot[mylabels]
{pdf("../Figures/Detection_prop.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(prop,
       col=c(NA, batch.colours[1:3], NA),
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1))
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
       col=c(NA, batch.colours[1:3], NA),
       xaxt="n", ylab="Proportion detected", xlab = "", cex.lab=1.5,
       type="h", lwd=2, ylim = c(0,1))
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
  legend("topright", lty=1, lwd=2, col=batch.colours[1:3], cex = 0.7,
         legend = batches[1:3], bg="white")
  dev.off()
}

{pdf("../Figures/Detection_prop_thresh.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(prop,
       col=ifelse(prop > 0.1,
                  c(NA, batch.colours[1:3], NA),
                  alpha(c(NA, batch.colours[1:3], NA), 0.5)),
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1))
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
       col=ifelse(prop > 0.1,
                  c(NA, batch.colours[1:3], NA),
                  alpha(c(NA, batch.colours[1:3], NA), 0.5)),
       xaxt="n", ylab="Proportion detected", xlab = "", cex.lab=1.5,
       type="h", lwd=2, ylim = c(0,1))
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
  abline(h = 0.1, lty = 2)
  legend("topright", lty=1, lwd=2, col=batch.colours[1:3], cex = 0.7,
         legend = batches[1:3], bg="white")
  dev.off()
}
