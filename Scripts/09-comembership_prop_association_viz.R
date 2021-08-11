### Multivariate clustering: Shortest path between siblings
### Rin on 10 Aug

# Load packages
library(focus)
library(colorspace)
library(RColorBrewer)

# Initialise
rm(list=ls())
path="~/HDAML/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

#### Between siblings ####
## Co-membership proportion
{covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))

fprop = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop.rds"))
families=levels(covars$Family.ID)

{pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop.pdf"),width=14)
  par(mar=c(5,5,1,1))
  plot(fprop, pch=19, col=family.colours[families], xaxt="n", las=1,
       panel.first=abline(v=1:length(families), lty=3, col="grey"),
       xlab="Family ID", cex=2,
       ylab="Co-membership proportion between siblings", cex.lab=1.5)
  axis(side=1, at=1:length(families), labels=families, las = 2)
  dev.off()
}
}
## Univariate
{
  X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_univar_betas.rds"))
  
  mylabels = names(betas)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]

  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])

  {pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_univariate.pdf"))
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
{stab = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_multivar_output.rds"))

  pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_multivariate_output.pdf"))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_multivariate_beta_consistency.pdf"))
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
  
  pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_multivariate_selprop.pdf"), width = 14)
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


#### Between stranger ####
## Univariate
{
  X = readRDS(paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_betas.rds"))
  
  if(sum(pvals==0)>0){
    warning(paste0("Replacing ", sum(pvals==0), " value(s) with minimum double"," (",.Machine$double.xmin,")"))
    pvals = ifelse(pvals==0, .Machine$double.xmin, pvals)
  }
  
  mylabels = names(betas)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]
  
  mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
  
  {pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_stranger_univariate.pdf"))
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
{stab = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_multivar_output.rds"))
  
  pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_stranger_multivariate_output.pdf"))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_stranger_multivariate_beta_consistency.pdf"))
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
  mylabels = gsub("regionYes","Same region (Y/N)",mylabels)
  
  pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_stranger_multivariate_selprop.pdf"), width = 14)
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
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_betas.rds"))
  
  if(sum(pvals==0)>0){
    warning(paste0("Replacing ", sum(pvals==0), " value(s) with minimum double"," (",.Machine$double.xmin,")"))
    pvals = ifelse(pvals==0, .Machine$double.xmin, pvals)
  }
  
  stab = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_multivar_output.rds"))
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
  
  yjitter = rep(0.03,length.out = (length(selprop)))
  yjitter[(sign(betas)*selprop)==-1] = rep(c(c(1,-1)*rep(c(seq(0.03,0.12,0.03),0.09,0.06),each = 2)),
                                           times = ceiling(sum((sign(betas)*selprop)==-1)/4))[1:sum((sign(betas)*selprop)==-1)]
  
  xjitter = rep(c(-10,+10),length.out = (length(betas)))
  xjitter = xjitter[order((sign(betas)*-log10(pvals)), selprop)]
  
  {pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_stranger_univar_vs_multivar.pdf"))
    par(mar=c(5,5,1,1))
    plot(sign(betas)*-log10(pvals), selprop, pch=19,
         col=ifelse(pvals < 0.05/ncol(X) | selprop > Argmax(stab)[2], mycolours, "grey"),
         cex.lab=1.5, cex = 0.7,
         ylim = c(0, 1.15),
         xlim = c(-max(-log10(pvals)[sign(betas)<0])-20, max(-log10(pvals)[sign(betas)>0])+40),
         ylab="Selection Proportion",
         xlab=expression(Signed~-log[10](p)))
    for (k in 1:length(mylabels)){
      if(pvals[k] < 0.05/ncol(X)){
        xadj = ifelse(xjitter<0,1,0)
        yadj = ifelse(yjitter<0,1,0)
        if(mylabels[k] %in% names(annot_sub)){
          label = substitute(Delta~tmp, list(tmp=mylabels[k]))
          text(sign(betas[k])*-log10(pvals[k])+xjitter[k], selprop[k]+yjitter[k],
               labels = label, col = mycolours[k], cex = 1, adj = c(xadj,yadj))
        } else {
          text(sign(betas[k])*-log10(pvals[k])+xjitter[k], selprop[k]+yjitter[k],
               labels = mylabels[k], col = mycolours[k], cex = 1, adj = c(xadj,yadj))
        }
      }
    }
    abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
    # abline(v= -log10(0.05/ncol(X)), lty = 2, col = "darkred")
    text(max(-log10(pvals)[sign(betas)>0])+40, Argmax(stab)[2],
         paste0("Selection proportion threshold = ",formatC(Argmax(stab)[2], digits = 2, format = "f")),
         adj = c(1,1), col = "darkred", cex = 0.8)
    # text(-log10(0.05/ncol(X)), 1,
    #      paste0("Bonferroni threshold = ",formatC(0.05/ncol(X), digits = 2, format = "e")),
    #      adj = c(0,0), srt = 270, col = "darkred", cex = 0.8)
    dev.off()
  }
}

