## Univariate regression modelling association with exposure
## 29 July

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
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  if(i %in% c(2,4)){
    X = covars %>%
      select(Age, Gender, Family.ID, Region) %>%
      mutate(`Ile-de-France` = ifelse(.$Region=="Île-de-France",1,0)) %>%
      select(-Region)
  } else if(i == 5){
    X = covars %>% select(Age, Gender, Family.ID, Batch)
  } else {
    X = covars %>% select(Age, Gender, Family.ID)
  }
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
  f0='expo[,k] ~ 1'
  t0=Sys.time()
  for (j in 1:ncol(X)){
    for (k in 1:ncol(expo)){
      model1=lm(as.formula(f1))
      if(colnames(X)[j] %in% c("Family.ID")){
        model0=lm(as.formula(f0))
        pvals=c(pvals, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
    } else {
      betas=c(betas, coefficients(model1)[2])
      pvals=c(pvals, summary(model1)$coefficients[2,4])
    }
    }
    }
  t1=Sys.time()
  print(t1-t0)
  betas = matrix(betas, nrow = ncol(expo), ncol = sum(!colnames(X) %in% c("Family.ID","Region")))
  pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
  rownames(pvals)=rownames(betas)=colnames(expo)
  colnames(pvals)=colnames(X)
  colnames(betas)=colnames(X)[!colnames(X) %in% c("Family.ID")]
  
  ifelse(dir.exists("../Results"),"",dir.create("../Results"))
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas_no_isolated.rds"))
}
### Check significant pvals for each data set ----
for (i in 1:length(batches)){
  pvals = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
  bonf = 0.05/nrow(pvals)
  
  list = apply(pvals, 2, function(x) x < bonf)
  assign(paste0("list_",suffix[i]),list)
}
colSums(list_lux)
colSums(list_fra)
colSums(list_gs)
colSums(list_pooled2)
colSums(list_pooled3)
rownames(list_lux)[list_lux[,"Age"]]
rownames(list_fra)[list_fra[,"Age"]]
rownames(list_gs)[list_gs[,"Age"]]
rownames(list_pooled2)[list_pooled2[,"Age"]]
rownames(list_pooled3)[list_pooled3[,"Age"]]

rownames(list_lux)[list_lux[,"Gender"]]
rownames(list_fra)[list_fra[,"Gender"]]
rownames(list_gs)[list_gs[,"Gender"]]
rownames(list_pooled2)[list_pooled2[,"Gender"]]
rownames(list_pooled3)[list_pooled3[,"Gender"]]


rownames(list_lux)[list_lux[,"Family.ID"]]
rownames(list_fra)[list_fra[,"Family.ID"]]
rownames(list_gs)[list_gs[,"Family.ID"]]
rownames(list_pooled2)[list_pooled2[,"Family.ID"]]
rownames(list_pooled3)[list_pooled3[,"Family.ID"]]

rownames(list_fra)[list_fra[,"Ile-de-France"]]
rownames(list_pooled3)[list_pooled3[,"Ile-de-France"]]

# ### Volcano plots ----
# # Age
# for (i in 1:length(batches)){
#   pvals = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
#   betas = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas_no_isolated.rds"))
#   annot_sub = annot[rownames(betas)]
#   if(i ==4){
#     pdf(paste0("../Figures/Section2/Univariate_exposure_age_",suffix[i],".pdf"), width = 5, height = 5)
#   } else{pdf(paste0("../Figures/Supplementary/Section2/Univariate_exposure_age_",suffix[i],".pdf"), width = 5, height = 5)
#   }
#   par(mar=c(5,5,1,1))
#   plot(betas[,"Age"], -log10(pvals[,"Age"]), pch=19,
#        col=ifelse(pvals[,"Age"] < 0.05/length(betas[,"Age"]), annot.colours[annot_sub], "grey"),
#        cex.lab=1.5, cex = 0.7,
#        ylim = c(0, max(-log10(pvals[,"Age"]))+0.25),
#        xlim = c(-max(abs(betas[,"Age"]))-0.2, max(abs(betas[,"Age"]))+0.2),
#        ylab=expression(-log[10](italic(p))),
#        xlab=expression(beta))
#   text(betas[,"Age"]+sign(betas[,"Age"])*0.05, -log10(pvals[,"Age"])+0.1,
#        labels = ifelse(pvals[,"Age"] < 0.05/length(betas[,"Age"]), rownames(betas), ""),
#        col = annot.colours[annot_sub], cex = 0.8)
#   abline(h = -log10(0.05/length(betas[,"Age"])), lty = 2, col = batch.colours[i])
#   text(-max(abs(betas[,"Age"]))-0.2, -log10(0.05/length(betas[,"Age"])),  labels = "Bonferroni threshold", adj=c(0,1), col = batch.colours[i], cex = 0.8)
#   dev.off()
# }
# 
# ### Gender
# for (i in 1:length(batches)){
#   pvals = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
#   betas = readRDS(paste0("../Results/",filepaths[i],"/Univariate_exposure_covariate_betas_no_isolated.rds"))
#   annot_sub = annot[rownames(betas)]
#   ## Volcano plots
#   if(i == 4){
#     pdf(paste0("../Figures/Section2/Univariate_exposure_gender_",suffix[i],".pdf"), width = 5, height = 5)
#   } else{pdf(paste0("../Figures/Supplementary/Section2/Univariate_exposure_gender_",suffix[i],".pdf"), width = 5, height = 5)
#   }
#   par(mar=c(5,5,1,1))
#   plot(betas[,"Gender"], -log10(pvals[,"Gender"]), pch=19,
#        col=ifelse(pvals[,"Gender"] < 0.05/length(betas[,"Gender"]), annot.colours[annot_sub], "grey"),
#        cex.lab=1.5, cex = 0.7,
#        ylim = c(0, max(-log10(pvals[,"Gender"]))+0.25),
#        xlim = c(-max(abs(betas[,"Gender"]))-0.2, max(abs(betas[,"Gender"]))+0.2),
#        ylab=expression(-log[10](italic(p))),
#        xlab=expression(beta))
#   text(betas[,"Gender"]+sign(betas[,"Gender"])*0.05, -log10(pvals[,"Gender"])+0.2,
#        labels = ifelse(pvals[,"Gender"] < 0.05/length(betas[,"Gender"]), rownames(betas), ""),
#        col = annot.colours[annot_sub], cex = 0.8)
#   abline(h = -log10(0.05/length(betas[,"Gender"])), lty = 2, col = batch.colours[i])
#   text(-max(abs(betas[,"Gender"]))-0.2, -log10(0.05/length(betas[,"Gender"])), labels = "Bonferroni threshold", adj=c(0,1), col = batch.colours[i], cex = 0.8)
#   dev.off()
# }

### Comparison ----
pvals = readRDS(paste0("../Results/",filepaths[4],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
betas = readRDS(paste0("../Results/",filepaths[4],"/Univariate_exposure_covariate_betas_no_isolated.rds"))

annot_sub = annot[rownames(pvals)]
# Scatter plot
bonf = 0.05/length(pvals)
bonf = -log10(bonf)

x = -log10(pvals[,"Family.ID"])

y = -log10(pvals[,"Ile-de-France"])
y = sign(betas[,"Ile-de-France"])*y

{pdf("../Figures/Section2/Univariate_exposure_family_region_pooled3.pdf", width=9, height=5)
  par(mar=c(5,5,1,20))
  plot(x, y,
       col=ifelse(x > bonf | y > bonf | y < -bonf,
                  annot.colours[annot_sub], "grey"),
       cex = 3,
       xlab=substitute(paste(-log[10](italic(p))," ",x), list(x = "Family membership")),
       ylab=substitute(paste("Signed ",-log[10](italic(p))," ",x), list(x = "Île-de-France")),
       cex.lab=1.5,
       ylim = c(min(y), bonf+0.5),
       type="p", pch = 19)
  
  text(x = max(x), y = bonf,labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[4], cex = 0.8)
  abline(h = bonf, lty = 2, col = darken(batch.colours[4], 0.5))
  
  text(x = max(x), y = -bonf,labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[4], cex = 0.8)
  abline(h = -bonf, lty = 2, col = darken(batch.colours[4], 0.5))
  
  text(x = bonf, y = min(y),labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[4], cex = 0.8, srt = 270)
  abline(v = bonf, lty = 2, col = darken(batch.colours[4], 0.5))
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  for (k in 1:length(x)){
    if(x[k] > bonf| y[k] > bonf | y[k] < -bonf){
      text(x[k], y[k],
           labels = which(names(x)[k] == names(x)[x > bonf| y > bonf | y < -bonf]),
           cex = 0.8)
    }
  }
  par(xpd = TRUE)
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(x > bonf| y > bonf | y < -bonf),
                         ". ", names(x)[x > bonf| y > bonf | y < -bonf]),
         text.col = darken(annot.colours[annot_sub[names(x)[x > bonf| y > bonf | y < -bonf]]],0.5),
         bty = "n",  x.intersp = 0, ncol = 2)
  dev.off()
}

### Comparison 2----
pvals = readRDS(paste0("../Results/",filepaths[5],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
betas = readRDS(paste0("../Results/",filepaths[5],"/Univariate_exposure_covariate_betas_no_isolated.rds"))
annot_sub = annot[rownames(pvals)]
# Scatter plot
bonf = 0.05/length(pvals)
bonf = -log10(bonf)

x = -log10(pvals[,"Family.ID"])

y = -log10(pvals[,"Batch"])
y = sign(betas[,"Batch"])*y

{pdf("../Figures/Section2/Univariate_exposure_family_batch_pooled2.pdf", width=9, height=5)
  par(mar=c(5,5,1,20))
  plot(x, y,
       col=ifelse(x > bonf | y > bonf | y < -bonf,
                  annot.colours[annot_sub], "grey"), cex = 3,
       xlab=substitute(paste(-log[10](italic(p))," ",x), list(x = "Family membership")),
       ylab=substitute(paste("Signed ",-log[10](italic(p))," ",x), list(x = "LUX vs GS")),
       cex.lab=1.5,
       type="p", pch = 19)
  
  text(x = max(x), y = bonf,labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[5], cex = 0.8)
  abline(h = bonf, lty = 2, col = darken(batch.colours[5], 0.5))
  
  text(x = max(x), y = -bonf,labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[5], cex = 0.8)
  abline(h = -bonf, lty = 2, col = darken(batch.colours[5], 0.5))
  
  text(x = bonf, y = max(y),labels="Bonferroni threshold", adj = c(0,0),col = batch.colours[5], cex = 0.8, srt = 270)
  abline(v = bonf, lty = 2, col = darken(batch.colours[5], 0.5))
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  
  condition = x > bonf & y <bonf & y >-bonf & y <0
  jitter = seq(17,max(x),length.out = sum(condition))
  names(jitter) = names(x)[condition][order(x[condition])]
  
  condition2 = x > bonf & y <bonf & y >-bonf & y>0
  jitter2 = seq(20,max(x),length.out = sum(condition2))
  names(jitter2) = names(x)[condition2][order(x[condition2])]
  
  for(k in 1:length(x)){
  if(x[k] > bonf| y[k] > bonf | y[k] < -bonf){
    if(x[k] > bonf & y[k]<bonf & y[k]>-bonf){
      if (y[k] > 0){
        yjitter = 8
        xjitter = jitter2[names(x)[k]]
      } else{
        yjitter = -8
        xjitter = jitter[names(x)[k]]
      }
      segments(x[k], y[k], xjitter,  yjitter, col = annot.colours[annot_sub[names(x)][k]])
      points(xjitter, yjitter,
             col = annot.colours[annot_sub[names(x)][k]],
             pch = 19,
             cex = 3)
      text(xjitter, yjitter, labels = which(names(x)[k] == names(x)[x > bonf| y > bonf | y < -bonf]), cex=0.8)
    } else{
      text(x[k], y[k], labels = which(names(x)[k] == names(x)[x > bonf| y > bonf | y < -bonf]), cex=0.8)
    } 
  }
  }

  par(xpd = TRUE)
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(x > bonf| y > bonf | y < -bonf),
                         ". ", names(x)[x > bonf| y > bonf | y < -bonf]),
         text.col = darken(annot.colours[annot_sub[names(x)[x > bonf| y > bonf | y < -bonf]]],0.5),
         bty = "n",  x.intersp = 0, ncol = 2)
  dev.off()
}

### Comparison 3 ----
pvals = readRDS(paste0("../Results/",filepaths[2],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
betas = readRDS(paste0("../Results/",filepaths[2],"/Univariate_exposure_covariate_betas_no_isolated.rds"))
annot_sub = annot[rownames(pvals)]
# Scatter plot
bonf = 0.05/length(pvals)
bonf = -log10(bonf)

x = -log10(pvals[,"Family.ID"])

y = -log10(pvals[,"Ile-de-France"])
y = sign(betas[,"Ile-de-France"])*y

{pdf("../Figures/Supplementary/Section2/Univariate_exposure_family_region_fra.pdf", width=9, height=5)
  par(mar=c(5,5,1,20))
  plot(x, y,
       col=ifelse(x > bonf | y > bonf | y < -bonf,
                  annot.colours[annot_sub], "grey"),
       cex = 3,
       xlab=substitute(paste(-log[10](italic(p))," ",x), list(x = "Family membership")),
       ylab=substitute(paste("Signed ",-log[10](italic(p))," ",x), list(x = "Île-de-France")),
       cex.lab=1.5,
       type="p", pch = 19)
  
  text(x = max(x), y = bonf,labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[2], cex = 0.8)
  abline(h = bonf, lty = 2, col = darken(batch.colours[2], 0.5))
  
  text(x = max(x), y = -bonf,labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[2], cex = 0.8)
  abline(h = -bonf, lty = 2, col = darken(batch.colours[2], 0.5))
  
  text(x = bonf, y = min(y),labels="Bonferroni threshold", adj = c(1,0),col = batch.colours[2], cex = 0.8, srt = 270)
  abline(v = bonf, lty = 2, col = darken(batch.colours[2], 0.5))
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  for (k in 1:length(x)){
    if(x[k] > bonf| y[k] > bonf | y[k] < -bonf){
      text(x[k], y[k],
           labels = which(names(x)[k] == names(x)[x > bonf| y > bonf | y < -bonf]),
           cex = 0.8)
    }
  }
  par(xpd = TRUE)
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(x > bonf| y > bonf | y < -bonf),
                         ". ", names(x)[x > bonf| y > bonf | y < -bonf]),
         text.col = darken(annot.colours[annot_sub[names(x)[x > bonf| y > bonf | y < -bonf]]],0.5),
         bty = "n",  x.intersp = 0, ncol = 2)
  dev.off()
}



### Manhattan plot ----
for (m in 1:3){
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Univariate_exposure_covariate_pvals_no_isolated.rds"))
  values = -log10(pvals[,"Family.ID"])
  annot_sub = annot[names(values)]

  bonf = 0.05/length(pvals)
  bonf = -log10(bonf)
  
  xseq = seq(1, length(values))
  
  {pdf(paste0("../Figures/Supplementary/Section2/Univariate_exposure_family_",suffix[m],".pdf"), width=12, height=6.5)
    par(mar=c(20,5,1,1))
    plot(values,
         col=annot.colours[annot_sub],
         xaxt="n", xlab="", ylab = expression(-log[10](italic(p))), cex.lab=1.5,
         type="p", ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)+1), pch = 19)
    abline(h = bonf, lty = 2, col = batch.colours[m])
    abline(h = -log10(0.05), lty = 2, col = "grey")
    text(x = 0, y = bonf,labels="Bonferroni threshold", adj = c(0,1),col = batch.colours[m], cex = 0.8)
    text(x = 0, y = -log10(0.05),labels="Nominal threshold", adj = c(0,1),col = "grey", cex = 0.8)
    abline(v = xseq, lty = 3, col = "grey")
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
    dev.off()
  }
  }


