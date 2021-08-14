### Multivariate clusterings
### Rin on 10 Aug

# Load packages
library(focus)
library(colorspace)
library(RColorBrewer)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

## Univariate
{
  X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Community_family_class_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Community_family_class_univar_betas.rds"))
  
  mylabels = names(betas)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]

  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])

  {pdf(paste0("../Figures/",filepaths[m],"/Community_family_class_univariate.pdf"))
    par(mar=c(5,5,1,1))
    plot(betas, -log10(pvals), pch=19,
         col=ifelse(pvals < 0.01, mycolours, "grey"),
         cex.lab=1.5, cex = 0.7,
         ylim = c(0, max(-log10(pvals))+0.25),
         xlim = c(-max(abs(betas))-5, max(abs(betas))+5),
         ylab=expression(-log[10](p)),
         xlab=expression(beta))
    for (k in 1:length(mylabels)){
      if(mylabels[k] %in% names(annot_sub)){
        if(pvals[k] < 0.01){
          label = substitute(Delta~tmp, list(tmp=mylabels[k]))
          text(betas[k]+sign(betas[k])*1, -log10(pvals[k])+0.1,
               labels = label, col = mycolours[k], cex = 0.7)
        }
      } else {
        text(betas[k]+sign(betas[k])*1, -log10(pvals[k])+0.1,
             labels = mylabels[k], col = mycolours[k], cex = 0.7)
      }
    }
    text(max(abs(betas))+4.9, -log10(0.05/ncol(X))+0.01,
         paste0("Bonferroni threshold = ",formatC(0.05/ncol(X), digits = 2, format = "e")),
         adj = c(1,0), col = "darkred", cex = 0.6)
    abline(h = -log10(0.05/ncol(X)), lty = 2, col = "darkred")
    dev.off()
  }
}

## Multivariate
{stab = readRDS(paste0("../Results/",filepaths[m],"/Community_family_class_multivar_output.rds"))

  pdf(paste0("../Figures/",filepaths[m],"/Community_family_class_multivariate_output.pdf"))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[m],"/Community_family_class_multivariate_beta_consistency.pdf"))
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
  
  mylabels = gsub("`","",names(selprop_ranked))
  mylabels = gsub("age_mu","Mean age",mylabels)
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("gender_diff","",mylabels)
  mylabels = gsub("All male","All male (ref. All female)",mylabels)
  mylabels = gsub("Different gender","Different gender (ref. All female)",mylabels)
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("regionYes","Île-de-France (Y/N)",mylabels)
  
  pdf(paste0("../Figures/",filepaths[m],"/Community_family_class_multivariate_selprop.pdf"), width = 14)
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
    if (mylabels[k] %in% names(annot_sub)){
      label = substitute(Delta~tmp, list(tmp = mylabels[k]))
    } else {
      label = mylabels[k]
    }
    axis(side = 1, at = k, las = 2, labels = label,
         col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
  }
  dev.off()
}

### Stranger status: Is the stranger in the same cluster or not?
## Univariate
{
  X = readRDS(paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Community_stranger_class_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Community_stranger_class_univar_betas.rds"))
  
  if(sum(pvals==0)>0){
    warning(paste0("Replacing ", sum(pvals==0), " value(s) with minimum double"," (",.Machine$double.xmin,")"))
    pvals = ifelse(pvals==0, .Machine$double.xmin, pvals)
  }
  
  mylabels = names(betas)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]
  
  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
  
  {pdf(paste0("../Figures/",filepaths[m],"/Community_stranger_class_univariate.pdf"))
    par(mar=c(5,5,1,1))
    plot(betas, -log10(pvals), pch=19,
         col=ifelse(pvals < 0.05/ncol(X), mycolours, "grey"),
         cex.lab=1.5, cex = 0.7,
         ylim = c(0, max(-log10(pvals))+0.25),
         xlim = c(-max(abs(betas))-0.05, max(abs(betas))+0.05),
         ylab=expression(-log[10](p)),
         xlab=expression(beta))
    for (k in 1:length(mylabels)){
      if(mylabels[k] %in% names(annot_sub)){
        if(pvals[k] < 0.05/ncol(X)){
          label = substitute(Delta~tmp, list(tmp=mylabels[k]))
          text(betas[k]+sign(betas[k])*0.01, -log10(pvals[k])+0.1,
               labels = label, col = mycolours[k], cex = 0.7)
        }
      } else {
        text(betas[k]+sign(betas[k])*0.01, -log10(pvals[k])+0.1,
             labels = mylabels[k], col = mycolours[k], cex = 0.7)
      }
    }
    text(max(abs(betas))+0.04, -log10(0.05/ncol(X))+0.01,
         paste0("Bonferroni threshold = ",formatC(0.05/ncol(X), digits = 2, format = "e")),
         adj = c(1,0), col = "darkred", cex = 0.6)
    abline(h = -log10(0.05/ncol(X)), lty = 2, col = "darkred")
    dev.off()
  }
}

## Multivariate
{stab = readRDS(paste0("../Results/",filepaths[m],"/Community_stranger_class_multivar_output.rds"))
  
  pdf(paste0("../Figures/",filepaths[m],"/Community_stranger_class_multivariate_output.pdf"))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[m],"/Community_stranger_class_multivariate_beta_consistency.pdf"))
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
  
  mylabels = gsub("`","",names(selprop_ranked))
  mylabels = gsub("age_mu","Mean age",mylabels)
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("gender_diff","",mylabels)
  mylabels = gsub("All male","All male (ref. All female)",mylabels)
  mylabels = gsub("Different gender","Different gender (ref. All female)",mylabels)
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("region","Same region (Y/N)",mylabels)
  
  pdf(paste0("../Figures/",filepaths[m],"/Community_stranger_class_multivariate_selprop.pdf"), width = 14)
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
    if (mylabels[k] %in% names(annot_sub)){
      label = substitute(Delta~tmp, list(tmp = mylabels[k]))
    } else {
      label = mylabels[k]
    }
    axis(side = 1, at = k, las = 2, labels = label,
         col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
  }
  dev.off()
}

## Univariate vs Multivariate
{
  X = readRDS(paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Community_stranger_class_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Community_stranger_class_univar_betas.rds"))
  
  stab = readRDS(paste0("../Results/",filepaths[m],"/Community_stranger_class_multivar_output.rds"))
  selprop=SelectionProportions(stab)
  
  mylabels = gsub("`","",names(selprop))
  mylabels = gsub("age_mu","Mean age",mylabels)
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("gender_diff","",mylabels)
  mylabels = gsub("All male","All male (ref. All female)",mylabels)
  mylabels = gsub("Different gender","Different gender (ref. All female)",mylabels)
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("regionYes","Same region (Y/N)",mylabels)
  
  all(names(betas)==mylabels)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]
  
  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
  
  x = sign(betas)*-log10(pvals)
  
  {pdf(paste0("../Figures/",filepaths[m],"/Community_stranger_class_univar_vs_multivar.pdf"), width = 9, height = 7)
    par(mar=c(5,5,5,15), xpd = TRUE)
    plot(x, selprop, pch=19,
         col=ifelse(pvals < 0.05/ncol(X) | selprop > Argmax(stab)[2], mycolours, "grey"),
         cex.lab=1.5, cex = 3,
         ylim = c(0, 1),
         xlim = range(x),
         ylab="Selection Proportion",
         xlab=expression(Signed~-log[10](p)))
    # jitter = seq(10,-10,length.out = )
    for (k in 1:length(mylabels)){
      if(pvals[k] < 0.05/ncol(X)){
        if(min(abs(x[-k]-x[k]))<3 & min(abs(selprop[-k]-selprop[k]))<0.05 & selprop[k] < Argmax(stab)[2]){
          mygrep = pvals < 0.05/ncol(X) & selprop < Argmax(stab)[2]
          jitter = seq(min(selprop),max(selprop),length.out = sum(mygrep))
          names(jitter) = mylabels[mygrep][order(selprop[mygrep])]
          jitter = jitter[mylabels[k]]
          segments(x[k], selprop[k], max(x)+1, jitter, col = mycolours[k])
          points(max(x)+1, jitter,
                 col = mycolours[k],
                 pch = 19,
                 cex = 3)
          text(max(x)+1, jitter,
               labels = which(mylabels[k] == mylabels[pvals < 0.05/ncol(X)]),
               cex = 0.8)
        } else if(min(abs(selprop[-k]-selprop[k]))<0.05 & min(abs(x[-k]-x[k]))<3  & selprop[k] > Argmax(stab)[2]){
          mygrep = selprop > Argmax(stab)[2]
          jitter = seq(min(x),max(x),length.out = sum(mygrep))
          names(jitter) = mylabels[mygrep][order(x[mygrep])]
          jitter = jitter[mylabels[k]]
          segments(x[k], selprop[k], jitter, max(selprop)+0.1, col = mycolours[k])
          points(jitter, max(selprop)+0.1,
                 col = mycolours[k],
                 pch = 19,
                 cex = 3)
          text(jitter, max(selprop)+0.1,
               labels = which(mylabels[k] == mylabels[pvals < 0.05/ncol(X)]),
               cex = 0.8)
        } else {
          text(x[k],
               selprop[k],
               labels = which(mylabels[k] == mylabels[pvals < 0.05/ncol(X)]),
               cex = 0.8)
          
        }
      }
    }
    # abline(v= -log10(0.05/ncol(X)), lty = 2, col = "darkred")
    text(min(x), Argmax(stab)[2],
         paste0("Selection proportion threshold = ",formatC(Argmax(stab)[2], digits = 2, format = "f")),
         adj = c(0,1), col = "darkred", cex = 1)
    # text(-log10(0.05/ncol(X)), 1,
    #      paste0("Bonferroni threshold = ",formatC(0.05/ncol(X), digits = 2, format = "e")),
    #      adj = c(0,0), srt = 270, col = "darkred", cex = 0.8)
    coord = par("usr")
    legend(x = coord[2]*1.05, y = coord[4],
           legend = paste0(1:length(mylabels[pvals < 0.05/ncol(X)]),". ", mylabels[pvals < 0.05/ncol(X)]),
           text.col = darken(mycolours[pvals < 0.05/ncol(X)],0.5),
           bty = "n",  x.intersp = 0)
    par(xpd = FALSE)
    abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
    dev.off()
  }
}
