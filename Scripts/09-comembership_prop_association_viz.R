### Multivariate clustering: Shortest path between siblings
### 10 Aug

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
suffix = c("lux","fra","gs","pooled3","pooled2")

#### Between siblings ####
## Co-membership proportion
### Compare ----
covars = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh_no_isolated.rds"))

param = NULL
for(m in 1:4){
  out = readRDS(paste0("../Results/",filepaths[m],"/Consensus_clustering_output.rds"))
  param = rbind(param, Argmax(out))
  fprop = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop.rds"))
  assign(paste0("fprop_",suffix[m]),fprop)
}

rownames(param) = batches[1:4]

{pdf(paste0("../Figures/Section4/Comembership_prop_compare.pdf"),width=,height=7)
  par(mar=c(5,7,1,3), mgp=c(4,1,0), pty="s")
  plot(fprop_pooled3, c(fprop_lux, fprop_fra, fprop_gs), pch=19,
       ylim = c(0,1),
       xlim = c(0,1),
       col=family.colours[names(fprop_pooled3)], las=1,
       panel.first=list(abline(v=axTicks(1), lty=3, col="grey"),
                        abline(h=axTicks(2), lty=3, col="grey"),
                        abline(0,1, lty=2, col="grey"),
                        abline(v = , lty=2, col="grey")),
       xlab="Co-membership proportion between siblings \n(Pooled)",
       ylab = "Co-membership proportion between siblings \n(By cohort)",
       cex=3, cex.lab=1.5)
  par(xpd = TRUE)
  condition = fprop_pooled3>0.9 & c(fprop_lux, fprop_fra, fprop_gs) >0.9
  jitter = seq(0,1,length.out = sum(condition))
  names(jitter) = names(fprop_pooled3)[condition][order(fprop_pooled3[condition])]
  
  for(k in 1:length(fprop_pooled3)){
    if(fprop_pooled3[k]>0.9 & c(fprop_lux, fprop_fra, fprop_gs)[k]>0.9){
      y = jitter[names(fprop_pooled3)[k]]
      x = 1.1
      segments(fprop_pooled3[k], c(fprop_lux, fprop_fra, fprop_gs)[k], x,  y, col = family.colours[names(fprop_pooled3)][k])
      points(x, y,
             col = family.colours[names(fprop_pooled3)][k],
             pch = 19,
             cex = 3)
      text(x, y, labels = names(fprop_pooled3)[k], cex=0.8)
    } else{
      text(fprop_pooled3[k], c(fprop_lux, fprop_fra, fprop_gs)[k], labels = names(fprop_pooled3)[k], cex=0.8)
    }
  }
  
  dev.off()
}

mycolours = region.colours[as.character(covars$Region[!duplicated(covars$Family.ID)])]
names(mycolours) = levels(covars$Family.ID)
{pdf(paste0("../Figures/Section4/Comembership_prop_compare_region.pdf"),width=9,height=7)
  par(mar=c(5,7,1,13), mgp=c(4,1,0), pty="s")
  plot(fprop_pooled3, c(fprop_lux, fprop_fra, fprop_gs), pch=19,
       ylim = c(0,1),
       xlim = c(0,1),
       col=mycolours[names(fprop_pooled3)], las=1,
       panel.first=list(abline(v=axTicks(1), lty=3, col="grey"),
                        abline(h=axTicks(2), lty=3, col="grey"),
                        abline(0,1, lty=2, col="grey")),
       xlab="Co-membership proportion between siblings \n(Pooled)",
       ylab = "Co-membership proportion between siblings \n(By cohort)",
       cex=3, cex.lab=1.5)
  par(xpd = TRUE)
  condition = fprop_pooled3>0.9 & c(fprop_lux, fprop_fra, fprop_gs) >0.9
  jitter = seq(0,1,length.out = sum(condition))
  names(jitter) = names(fprop_pooled3)[condition][order(fprop_pooled3[condition])]
  
  for(k in 1:length(fprop_pooled3)){
    if(fprop_pooled3[k]>0.9 & c(fprop_lux, fprop_fra, fprop_gs)[k]>0.9){
      y = jitter[names(fprop_pooled3)[k]]
      x = 1.1
      segments(fprop_pooled3[k], c(fprop_lux, fprop_fra, fprop_gs)[k], x,  y, col = mycolours[names(fprop_pooled3)][k])
      points(x, y,
             col = mycolours[names(fprop_pooled3)][k],
             pch = 19,
             cex = 3)
      text(x, y, labels = names(fprop_pooled3)[k], cex=0.8)
    } else{
      text(fprop_pooled3[k], c(fprop_lux, fprop_fra, fprop_gs)[k], labels = names(fprop_pooled3)[k], cex=0.8)
    }
  }
  coord = par("usr")
  legend(x = coord[2]+.1, y = coord[4],
         legend = levels(covars$Region),
         pch = 19,
         pt.cex = 2,
         col = region.colours[levels(covars$Region)],
         bty = "n")
  
  dev.off()
}

### Compare (Sensitivity)----
for(m in c(1,3,5)){
  fprop = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop.rds"))
  assign(paste0("fprop_",suffix[m]),fprop)
}

{pdf(paste0("../Figures/Supplementary/Comembership_prop_compare_pooled2.pdf"),width=7,height=7)
  par(mar=c(5,7,1,3), mgp=c(4,1,0), pty="s")
  plot(fprop_pooled2, c(fprop_lux, fprop_gs), pch=19,
       ylim = c(0,1),
       xlim = c(0,1),
       col=family.colours[names(fprop_pooled2)], las=1,
       panel.first=list(abline(v=axTicks(1), lty=3, col="grey"),
                        abline(h=axTicks(2), lty=3, col="grey"),
                        abline(0,1, lty=2, col="grey")),
       xlab="Co-membership proportion between siblings \n(Pooled)",
       ylab = "Co-membership proportion between siblings \n(By cohort)",
       cex=3, cex.lab=1.5)
  par(xpd = TRUE)
  condition = fprop_pooled2>0.9 & c(fprop_lux, fprop_gs) >0.9
  jitter = seq(0,1,length.out = sum(condition))
  names(jitter) = names(fprop_pooled2)[condition][order(fprop_pooled2[condition])]
  
  condition2 = fprop_pooled2>0.9 & c(fprop_lux, fprop_gs) <0.1
  jitter2 = seq(0,1,length.out = sum(condition2))
  names(jitter2) = names(fprop_pooled2)[condition2][order(c(fprop_lux, fprop_gs)[condition2])]
  
  for(k in 1:length(fprop_pooled2)){
    if(fprop_pooled2[k]>0.9 & c(fprop_lux, fprop_gs)[k]>0.9){
      y = 1.1
      x = jitter[names(fprop_pooled2)[k]]
      segments(fprop_pooled2[k], c(fprop_lux, fprop_gs)[k], x,  y, col = family.colours[names(fprop_pooled2)][k])
      points(x, y,
             col = family.colours[names(fprop_pooled2)][k],
             pch = 19,
             cex = 3)
      text(x, y, labels = names(fprop_pooled2)[k], cex=0.8)
    } else if(fprop_pooled2[k]>0.9 & c(fprop_lux, fprop_gs)[k]<0.1){
      y = jitter2[names(fprop_pooled2)[k]]
      x = 1.1
      segments(fprop_pooled2[k], c(fprop_lux, fprop_gs)[k], x,  y, col = family.colours[names(fprop_pooled2)][k])
      points(x, y,
             col = family.colours[names(fprop_pooled2)][k],
             pch = 19,
             cex = 3)
      text(x, y, labels = names(fprop_pooled2)[k], cex=0.8)
    } else {
      text(fprop_pooled2[k], c(fprop_lux, fprop_gs)[k], labels = names(fprop_pooled2)[k], cex=0.8)
    }
  }
  
  dev.off()
}

# {pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop.pdf"),width=14)
#   par(mar=c(5,5,1,1))
#   plot(fprop, pch=19, col=family.colours[families], xaxt="n", las=1,
#        panel.first=abline(v=1:length(families), lty=3, col="grey"),
#        xlab="Family ID", cex=2,
#        ylab="Co-membership proportion between siblings", cex.lab=1.5)
#   axis(side=1, at=1:length(families), labels=families, las = 2)
#   dev.off()
# }
# }
## Family characteristics: Univariate vs Multivariate ----
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
suffix = c("lux","fra","gs","pooled3","pooled2")
for (m in 1:length(batches)){
  X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_univar_betas.rds"))
  
  stab = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_multivar_output.rds"))
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
  pdf(paste0("../Figures/Section4/Comembership_prop_univar_vs_multivar_",suffix[m],".pdf"), width = 7, height = 5.1)
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
      if(length(intersect(which(abs(x[-k]-x[k])<xthres),which(abs(selprop[-k]-selprop[k])<ythres))) > 0 &
         selprop[k] < Argmax(stab)[2]){
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


#### Between stranger ####
## Univariate vs Multivariate
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
suffix = c("lux","fra","gs","pooled3","pooled2")
for (m in 1:length(batches)){
  X = readRDS(paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))
  pvals = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_pvals.rds"))
  betas = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_betas.rds"))
  
  stab = readRDS(paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_multivar_output.rds"))
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
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("regionYes","Same region (Y/N)",mylabels)
  
  annot_sub = annot[which(names(annot) %in% mylabels)]
  bonf = 0.05/ncol(X)
  mycolours = c(rep("magenta",length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
  
  if(sum(pvals==0)>0){
    warning(paste0("Replacing ", sum(pvals==0), " value(s) with minimum double"," (",.Machine$double.xmin,")"))
    mylabels = ifelse(pvals==0, paste0(mylabels,"*"), mylabels)
    pvals = ifelse(pvals==0, .Machine$double.xmin, pvals)
  }
  
  x = sign(betas)*-log10(pvals)
  
  if (m %in% c(2,4,5)){
    pdf(paste0("../Figures/Section4/Comembership_prop_stranger_univar_vs_multivar_",suffix[m],".pdf"), width = 7, height = 5.1)
    par(mar=c(5,5,2,10), xpd = TRUE)
    plot(x, selprop, pch=19,
         col=ifelse(pvals < bonf & selprop > Argmax(stab)[2], mycolours, "grey"),
         cex.lab=1.5, cex = 3,
         ylim = c(0, 1),
         xlim = c(min(log10(bonf),x), max(c(-log10(bonf), x))),
         ylab="Selection Proportion",
         xlab=expression(Signed~-log[10](p)))
    condition  = pvals < bonf & selprop > Argmax(stab)[2]
    for (k in 1:length(mylabels)){
      if(pvals[k] < bonf & selprop[k] > Argmax(stab)[2]){
        ## Check consistency of beta coefficients
        if(sign(betas[k]) != sign(beta_mu[k])){
          warning(paste0("Inconsistent direction of associaion for ", names(beta_mu)[k]))
        }
        xthres = sum(abs(range(x)))/60
        ythres = 0.05
        if(length(intersect(which(abs(x[-k]-x[k])<xthres),which(abs(selprop[-k]-selprop[k])<ythres))) > 0
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
  } else {pdf(paste0("../Figures/Section4/Comembership_prop_stranger_univar_vs_multivar_",suffix[m],".pdf"), width = 7, height = 5.1)
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
        if(length(intersect(which(abs(x[-k]-x[k])<xthres),which(abs(selprop[-k]-selprop[k])<ythres))) > 0 &
           selprop[k] < Argmax(stab)[2]){
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
}


