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
  covars = covars[which(covars$Family.ID != "Isolated"),]
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


{pdf("../Figures/Univariate_exposure_family.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1],
       col=batch.colours[1],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="n", ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)),
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
                           darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
    }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("bottomright", pch=c(rep(19,3),17, rep(NA,5)), lty = c(rep(NA,4),rep(2,5)),
         col=c(batch.colours[1:4],darken(batch.colours[1:4], 0.5), "grey"),
         legend = c(batches[1:4],
                    paste0("Bonferroni threshold (",c("LUX","FRA","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.7, ncol = 2)
  dev.off()
}

{pdf("../Figures/Univariate_exposure_family_Pooled2.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1],
       col=batch.colours[1],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="n", ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)),
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
                           darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("bottomright", pch=c(rep(19,2),17, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
         col=c(batch.colours[c(1,3,5)],darken(batch.colours[c(1,3,5)], 0.5), "grey"),
         legend = c(batches[c(1,3,5)],
                    paste0("Bonferroni threshold (",c("LUX","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.7, ncol = 2)
  dev.off()
}

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
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(nd)==rownames(covars)))
  
  ## Remove isolated children (no siblings)
  covars = covars[which(covars$Family.ID != "Isolated"),]
  nd = nd[rownames(covars),]
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
       type="n", ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)))
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
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.6,
         col.axis = ifelse(isTRUE(values[i,4] > bonf[4]),
                           darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("topright", pch=c(rep(19,3),17, rep(NA,5)), lty = c(rep(NA,4),rep(2,5)),
         col=c(batch.colours[1:4],darken(batch.colours[1:4], 0.5), "grey"),
         legend = c(batches[1:4],
                    paste0("Bonferroni threshold (",c("LUX","FRA","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.7, ncol = 2)
  dev.off()
}

{pdf("../Figures/Univariate_detection_family_Pooled2.pdf", width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1],
       col=batch.colours[1],
       xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
       type="n", ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)))
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
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.6,
         col.axis = ifelse(isTRUE(values[i,5] > bonf[5]),
                           darken(annot.colours[(annot_sub)[i]], amount=0.5),"black"))
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("topright", pch=c(rep(19,2),17, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
         col=c(batch.colours[c(1,3,5)],darken(batch.colours[c(1,3,5)], 0.5), "grey"),
         legend = c(batches[c(1,3,5)],
                    paste0("Bonferroni threshold (",c("LUX","GS","Pooled"),")"),
                    "Nominal threshold"),
         bg="white", cex = 0.7, ncol = 2)
  dev.off()
}

### Within-family variation ----
families=unique(covars$Family.ID)
family_sd=matrix(NA, nrow=length(families), ncol=ncol(expo))
for (p in 1:ncol(expo)){
  for (f in 1:length(families)){
    family_sd[f,p]=sd(expo[covars$Family.ID==families[f],p])
  }
}
colnames(family_sd)=colnames(expo)
rownames(family_sd)=families
overall_sd=apply(expo,2,sd)

mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families

x=as.vector(row(family_sd))
y=as.vector(family_sd)
z=as.vector(col(family_sd))
mycolours = family.colours[unique(annot)]
annot_sub=annot[colnames(family_sd)]

{pdf(paste0("../Figures/",filepath,"/Univariate_sd_family_overall_cont.pdf"), width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(z,y, pch=19, cex=0.5, col=mycolours[x], las=1, xaxt="n",
       xlab="", ylab="Standard deviation", cex.lab=1.5,
       panel.first=abline(v=unique(z),lty=3,col="grey"),
       ylim=range(c(family_sd,overall_sd)))
  points(overall_sd, pch=15)
  for (k in 1:ncol(family_sd)){
    axis(side=1, at=k, labels=colnames(family_sd)[k], cex.axis=0.8, las=2)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(colnames(family_sd))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(mycolours[as.character(unique(annot_sub)[k])], amount=0.5), cex.axis = 0.9)
  }
  legend("top", col=c(mycolours, "black"),
         pch=c(rep(19,length(families)),15),
         pt.cex=c(rep(0.5,length(families)),1),
         legend=c(families, "Overall"), ncol=10, bg="white")
  dev.off()
}

# ### Intra-class correlation ----
# library(lme4)
# 
# expo = mat_lux
# covars = covar_lux
# 
# ### Univariate linear mixed models
# icc_lux=NULL
# for (p in 1:ncol(expo)){
#   x=expo[,p]
#   model=lmer(x~(1|covars$Family.ID))
#   vcov = as.data.frame(VarCorr(model))$vcov
#   icc_lux=c(icc_lux,vcov[1]/sum(vcov))
# }
# names(icc_lux) = colnames(expo)
# 
# expo = mat_fra
# covars = covar_fra
# 
# ### Univariate linear mixed models
# icc_fra=NULL
# for (p in 1:ncol(expo)){
#   x=expo[,p]
#   model=lmer(x~(1|covars$Family.ID))
#   vcov = as.data.frame(VarCorr(model))$vcov
#   icc_fra=c(icc_fra,vcov[1]/sum(vcov))
# }
# names(icc_fra) = colnames(expo)
# 
# expo = mat_gs
# covars = covar_gs
# 
# ### Univariate linear mixed models
# icc_gs=NULL
# for (p in 1:ncol(expo)){
#   x=expo[,p]
#   model=lmer(x~(1|covars$Family.ID))
#   vcov = as.data.frame(VarCorr(model))$vcov
#   icc_gs=c(icc_gs,vcov[1]/sum(vcov))
# }
# names(icc_gs) = colnames(expo)
# 
# # Manhattan plot
# mycolours = c("tomato","royalblue","forestgreen")
# values = merge(as.data.frame(icc_lux),as.data.frame(icc_fra), by=0, all=TRUE)
# rownames(values) = values[,1]
# values = merge(values[,-1], as.data.frame(icc_gs), by=0, all=TRUE)
# rownames(values) = values[,1]
# values = values[,-1]
# 
# # Sort rows
# values = values[intersect(names(annot), rownames(values)),]
# 
# mycolours=c("tomato","royalblue","forestgreen")
# annot_sub=annot[rownames(values)]
# xseq = seq(1, nrow(values))
# 
# {pdf(paste0("../Figures/Intra_class_correlation_univariate_expo_cont.pdf"), width=14, height=8)
#   par(mar=c(20,5,1,1))
#   plot(values[,1], pch=17, cex=0.7, las=1, xaxt="n", type = "n",
#        ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)),
#        xlab="", ylab="Intra-Class Correlation", cex.lab=1.5,
#        panel.first=abline(v=xseq,lty=3,col="grey"),
#        col=mycolours[1])
#   points(values[,1], pch = 17, col = mycolours[1], cex = 0.7)
#   points(values[,2], pch = 19, col = mycolours[2], cex = 0.7)
#   points(values[,3], pch = 15, col = mycolours[3], cex = 0.8)
#   for (k in 1:length(xseq)){
#     axis(side=1, at=xseq[k], labels=rownames(values)[k], cex.axis=0.8, las=2)
#   }
#   xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
#   axis(side=1, line=8, at=xgroup, labels=NA)
#   tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
#   for (k in 1:length(unique(annot_sub))){
#     axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
#          col.axis=darken(family.colours[unique(annot_sub)[k]], amount=0.5),
#          cex.axis=0.9)
#   }
#   legend("topright", pch=c(17,19,15),
#          col=mycolours,
#          legend = c("Luxembourg","France","Grande-Synthe"),
#          cex=0.7, bg="white")
#   dev.off()
# }
