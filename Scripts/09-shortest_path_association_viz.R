### Multivariate clustering Shortest path between siblings
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

## Individual characteristics: Univariate vs Multivariate ----
for (m in 1:length(batches)){
    pvals = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_individual_univar_pvals.rds"))
    betas = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_individual_univar_betas.rds"))
    
    stab = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_individual_multivar_output.rds"))
    selprop=SelectionProportions(stab)
    beta_stab=stab$Beta[ArgmaxId(stab)[1],,]
    beta_mu = rowMeans(beta_stab[-nrow(beta_stab),] %*% diag(beta_stab[nrow(beta_stab),]))
    names(beta_mu) = names(selprop)
    
    mylabels = gsub("`","",names(selprop))
    mylabels = gsub("Gender","",mylabels)
    mylabels = gsub("Male","Male (ref. Female)",mylabels)
    mylabels = gsub("RegionYes","Île-de-France (Y/N)",mylabels)
    
    annot_sub = annot[which(names(annot) %in% mylabels)]
    bonf = 0.05/length(betas) 
    mycolours = c(rep("magenta",length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
    
    x = sign(betas)*-log10(pvals)
   pdf(paste0("../Figures/Section3/Shortest_path_individual_univar_vs_multivar_",suffix[m],".pdf"), width = 7, height = 5.1)
   par(mar=c(5,5,2,10), xpd = TRUE)
   plot(x, selprop, pch=19,
        col=ifelse(pvals < bonf| selprop > Argmax(stab)[2], mycolours, "grey"),
        cex.lab=1.5, cex = 3,
        ylim = c(0, 1),
        xlim = c(min(log10(bonf),x), max(c(-log10(bonf), x))),
        ylab="Selection Proportion",
        xlab=expression(Signed~-log[10](p)))
   condition  = pvals < bonf| selprop > Argmax(stab)[2]
   for (k in 1:length(mylabels)){
     if(pvals[k] < bonf|selprop[k] > Argmax(stab)[2]){
       ## Check consistency of beta coefficients
       if(sign(betas[k]) != sign(beta_mu[k])){
         warning(paste0("Inconsistent direction of associaion for ", names(beta_mu)[k]))
       }
       xthres = sum(abs(c(range(x))))/60
       ythres = 0.05
       if(length(intersect(which(abs(x[-k]-x[k])<xthres),which(abs(selprop[-k]-selprop[k])<ythres))) > 0
          & selprop[k] < Argmax(stab)[2]){
         mygrep = pvals < bonf& selprop < Argmax(stab)[2]
         jitter = seq(min(selprop),max(selprop),length.out = sum(mygrep))
         names(jitter) = mylabels[mygrep][order(selprop[mygrep])]
         jitter = jitter[mylabels[k]]
         segments(x[k], selprop[k], max(x)+1, jitter, col = mycolours[k])
         points(max(x)+1, jitter,
                col = mycolours[k],
                pch = 19,
                cex = 3)
         text(max(x)+1, jitter,
              labels = which(mylabels[k] == mylabels[condition]),
              cex = 0.8)
       } else if(length(intersect(which(abs(x[-k]-x[k])<xthres),which(abs(selprop[-k]-selprop[k])<ythres))) > 0
         & selprop[k] > Argmax(stab)[2]){
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
              labels = which(mylabels[k] == mylabels[condition]),
              cex = 0.8)
       } else {
         text(x[k],
              selprop[k],
              labels = which(mylabels[k] == mylabels[condition]),
              cex = 0.8)
       }
     }
   }
   text(max(min(x),0), Argmax(stab)[2], "Selection proportion threshold",
        adj = c(0.5,1), col = batch.colours[m], cex = 1)
   text(-log10(bonf), 0, "Bonferroni threshold",
        adj = c(1,0), srt = 270, col = batch.colours[m])
   text(log10(bonf), 0, "Bonferroni threshold",
        adj = c(1,0), srt = 270, col = batch.colours[m])
   coord = par("usr")
   legend(x = coord[2]*1.05, y = coord[4],
          legend = paste0(1:sum(condition),". ", mylabels[condition]),
          text.col = darken(mycolours[condition],0.5),
          bty = "n",  x.intersp = 0)
   par(xpd = FALSE)
   abline(h = Argmax(stab)[2], lty = 2, col = batch.colours[m])
   abline(v= -log10(bonf), lty = 2, col = batch.colours[m])
   abline(v= log10(bonf), lty = 2, col = batch.colours[m])
   dev.off()
}

## Family characteristics: Univariate vs Multivariate ----
for (m in 1:length(batches)){
  X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_univar_betas.rds"))
  
  stab = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_multivar_output.rds"))
  selprop=SelectionProportions(stab)
  beta_stab=stab$Beta[ArgmaxId(stab)[1],,]
  beta_mu = rowMeans(beta_stab[-nrow(beta_stab),] %*% diag(beta_stab[nrow(beta_stab),]))
  names(beta_mu) = names(selprop)
  
  mylabels = gsub("`","",names(selprop))
  mylabels = gsub("age_mu","Mean age",mylabels)
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("gender_diff","",mylabels)
  mylabels = gsub("All male","All male (ref. All female)",mylabels)
  mylabels = gsub("Different gender","Different gender (ref. All female)",mylabels)
  mylabels = gsub("regionYes","Île-de-France (Y/N)",mylabels)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]
  bonf = 0.05/ncol(X)
  mycolours = c(rep("magenta",length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
  
  x = sign(betas)*-log10(pvals)
  pdf(paste0("../Figures/Section3/Shortest_path_univar_vs_multivar_",suffix[m],".pdf"), width = 7, height = 5.1)
  par(mar=c(5,5,2,10), xpd = TRUE)
  plot(x, selprop, pch=19,
       col=ifelse(pvals < bonf | selprop > Argmax(stab)[2], mycolours, "grey"),
       cex.lab=1.5, cex = 3,
       ylim = c(0, 1),
       xlim = c(min(log10(bonf),x), max(c(-log10(bonf), x))),
       ylab="Selection Proportion",
       xlab=expression(Signed~-log[10](p)))
  condition  = pvals < bonf | selprop > Argmax(stab)[2]
  for (k in 1:length(mylabels)){
    if(pvals[k] < bonf|selprop[k] > Argmax(stab)[2]){
      ## Check consistency of beta coefficients
      if(sign(betas[k]) != sign(beta_mu[k])){
        warning(paste0("Inconsistent direction of associaion for ", names(beta_mu)[k]))
      }
      xthres = sum(abs(range(x)))/60
      ythres = 0.05
      if(length(intersect(which(abs(x[-k]-x[k])<xthres),which(abs(selprop[-k]-selprop[k])<ythres))) > 0
         & selprop[k] < Argmax(stab)[2]){
        mygrep = pvals < bonf & selprop < Argmax(stab)[2]
        jitter = seq(min(selprop),max(selprop),length.out = sum(mygrep))
        names(jitter) = mylabels[mygrep][order(selprop[mygrep])]
        jitter = jitter[mylabels[k]]
        segments(x[k], selprop[k], max(x)+1, jitter, col = mycolours[k])
        points(max(x)+1, jitter,
               col = mycolours[k],
               pch = 19,
               cex = 3)
        text(max(x)+1, jitter,
             labels = which(mylabels[k] == mylabels[condition]),
             cex = 0.8)
      } else if(length(intersect(which(abs(x[-k]-x[k])<xthres),which(abs(selprop[-k]-selprop[k])<ythres))) > 0
                & selprop[k] > Argmax(stab)[2]){
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
             labels = which(mylabels[k] == mylabels[condition]),
             cex = 0.8)
      } else {
        text(x[k],
             selprop[k],
             labels = which(mylabels[k] == mylabels[condition]),
             cex = 0.8)
      }
    }
  }
  text(max(min(x),0), Argmax(stab)[2], "Selection proportion threshold",
       adj = c(0.5,1), col = batch.colours[m], cex = 1)
  text(-log10(bonf), 0, "Bonferroni threshold",
       adj = c(1,0), srt = 270, col = batch.colours[m])
  text(log10(bonf), 0, "Bonferroni threshold",
       adj = c(1,0), srt = 270, col = batch.colours[m])
  coord = par("usr")
  legend(x = coord[2]*1.05, y = coord[4],
         legend = paste0(1:sum(condition),". ", mylabels[condition]),
         text.col = darken(mycolours[condition],0.5),
         bty = "n",  x.intersp = 0)
  par(xpd = FALSE)
  abline(h = Argmax(stab)[2], lty = 2, col = batch.colours[m])
  abline(v= -log10(bonf), lty = 2, col = batch.colours[m])
  abline(v= log10(bonf), lty = 2, col = batch.colours[m])
  dev.off()
}

