### Multivariate clusterinL Shortest path between siblings
### Rin on 10 Aug

# Load packages
library(focus)
library(colorspace)
library(RColorBrewer)

# Initialise
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

m = 1
m = 2
m = 3
m = 4
m = 5
{
  X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_univar_betas.rds"))
  
  mylabels = names(betas)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]

  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])

  {pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_univariate.pdf"))
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
{stab = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_multivar_output.rds"))

  pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_output.pdf"))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_beta_consistency.pdf"))
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
  mylabels = gsub("region","Île-de-France (Y/N)",mylabels)
  
  pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_selprop.pdf"), width = 14)
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
  X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_univar_betas.rds"))
  
  stab = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_multivar_output.rds"))
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
  mylabels = gsub("regionYes","Île-de-France (Y/N)",mylabels)
  
  all(names(betas)==mylabels)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]
  
  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
  
  {pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_univar_vs_multivar.pdf"))
    par(mar=c(5,5,1,1))
    plot(sign(betas)*-log10(pvals), selprop, pch=19,
         col=ifelse(pvals < 0.05/ncol(X) | selprop > Argmax(stab)[2], mycolours, "grey"),
         cex.lab=1.5, cex = 0.7,
         ylim = c(0, 1),
         xlim = c(-max(-log10(pvals))-1, max(-log10(pvals))+1),
         ylab="Selection Proportion",
         xlab=expression(Signed~-log[10](p)))
    for (k in 1:length(mylabels)){
      if(pvals[k] < 0.05/ncol(X)){
        if(mylabels[k] %in% names(annot_sub)){
          label = substitute(Delta~tmp, list(tmp=mylabels[k]))
          text(sign(betas[k])*-log10(pvals[k])*1, selprop[k]+0.05,
               labels = label, col = mycolours[k], cex = 1)
          } else {
            text(sign(betas[k])*-log10(pvals[k])*1, selprop[k]+0.05,
                 labels = mylabels[k], col = mycolours[k], cex = 1)
          }
      }
    }
    abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
    abline(v= -log10(0.05/ncol(X)), lty = 2, col = "darkred")
    text(-max(-log10(pvals))-1, Argmax(stab)[2],
         paste0("Selection proportion threshold = ",formatC(Argmax(stab)[2], digits = 2, format = "f")),
         adj = c(0,0), col = "darkred", cex = 0.8)
    text(-log10(0.05/ncol(X)), 1,
         paste0("Bonferroni threshold = ",formatC(0.05/ncol(X), digits = 2, format = "e")),
         adj = c(0,0), srt = 270, col = "darkred", cex = 0.8)
    
    dev.off()
  }
}

### Comparing clustering performance between Pooled, France and LUX/GS
fsp_pooled = readRDS(paste0("../Results/",filepaths[4],"/Shortest_path_cont_all.rds"))
fsp_luxgs = readRDS(paste0("../Results/",filepaths[5],"/Shortest_path_cont_all.rds"))
fsp_fra = readRDS(paste0("../Results/",filepaths[2],"/Shortest_path_cont_all.rds"))

y = c(fsp_luxgs[grepl("L",names(fsp_luxgs))], fsp_fra,fsp_luxgs[grepl("G",names(fsp_luxgs))])
{pdf("../Figures/Shortest_path_cont_all_compare.pdf")
  par(mar=c(5,5,1,1))
  plot(fsp_pooled, y, pch=19, col = family.colours[names(fsp_pooled)],
       xlim = c(2, max(c(fsp_pooled,y))),
       ylim = c(2, max(c(fsp_pooled,y))),
       cex.lab=1.5, cex = 3,
       ylab="Shortest path (FRA and LUX/GS combined)",
       xlab="Shortest path (Pooled)")
  abline(0,1, lty = 3)
  text(fsp_pooled, y, names(fsp_pooled), cex = 0.8)
  dev.off()
}
