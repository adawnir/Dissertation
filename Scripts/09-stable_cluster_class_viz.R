## Characterisitng well/ill-classified families Visualisation
## Rin Wada 15 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(focus)

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

#### Consensus clustering ####
{
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Stable_cluster_class_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Stable_cluster_class_univar_betas.rds"))
  
  mylabels = names(betas)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]
  
  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
  
  {pdf(paste0("../Figures/",filepaths[m],"/Stable_cluster_class_univariate.pdf"))
    par(mar=c(5,5,1,1))
    plot(betas, -log10(pvals), pch=19,
         col=ifelse(pvals < 0.01, mycolours, "grey"),
         cex.lab=1.5, cex = 0.7,
         ylim = c(0, max(-log10(pvals))+0.25),
         xlim = c(-max(abs(betas))-0.1, max(abs(betas))+0.1),
         ylab=expression(-log[10](p)),
         xlab=expression(beta))
    for (k in 1:length(mylabels)){
      if(mylabels[k] %in% names(annot_sub)){
        if(pvals[k] < 0.01){
          text(betas[k], -log10(pvals[k])+0.1,
               labels = mylabels[k], col = mycolours[k], cex = 0.7)
        }
      } else {
        text(betas[k], -log10(pvals[k])+0.1,
             labels = mylabels[k], col = mycolours[k], cex = 0.7)
      }
    }
    text(min(betas)-0.1, -log10(0.05/length(betas)),
         paste0("Bonferroni threshold = ",formatC(0.05/length(betas), digits = 2, format = "e")),
         adj = c(0,0), col = "darkred", cex = 0.6)
    abline(h = -log10(0.05/length(betas)), lty = 2, col = "darkred")
    dev.off()
  }
}
{stab = readRDS(paste0("../Results/",filepaths[m],"/Stable_cluster_class_multivar_output.rds"))
  
  pdf(paste0("../Figures/",filepaths[m],"/Stable_cluster_class_multivariate_output.pdf"))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[m],"/Stable_cluster_class_multivariate_beta_consistency.pdf"))
  par(mar=c(5,5,1,1))
  plot(a/(a+b), selprop, las=1, pch=19, col="navy", cex.lab=1.5,
       xlab="Proportion of positive beta among non-zero betas",
       ylab="Selection Proportion") # Ideally no point around 0.5 with high selection proportion
  dev.off()
  
  # Extract beta
  beta = stab$Beta[ArgmaxId(stab)[1],,]
  beta_mu = rowMeans(beta[-nrow(beta),] %*% diag(beta[nrow(beta),]))
  names(beta_mu) = names(selprop)
  
  selprop_ranked <- sort(selprop, decreasing = TRUE)
  beta_mu_ranked <- beta_mu[order(-selprop)]
  print(all(names(selprop_ranked)==names(beta_mu_ranked)))
  
  pdf(paste0("../Figures/",filepaths[m],"/Stable_cluster_class_multivariate_selprop.pdf"), width = 14)
  par(mar=c(15, 5, 1, 1))
  plot(selprop_ranked,
       type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
       col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[m], no = "grey"),
       xaxt = "n", xlab = "", ylab = "Selection proportions"
  )
  text(length(selprop_ranked)-0.01, Argmax(stab)[2]+0.01, paste0("Selection proportion threshold = ",Argmax(stab)[2]), adj = c(1,0), col = "darkred")
  abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
  axis(side = 1, at = 1:length(selprop_ranked), labels = NA)
  for (k in 1:length(selprop_ranked)){
    if (beta_mu_ranked[k] > 0){
      col = "red"
    } else if (beta_mu_ranked[k] < 0){
      col = "blue"
    } else {
      col = "black"
    }
    axis(side = 1, at = k, las = 2, labels = mylabels[k],
         col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
  }
  dev.off()
}
