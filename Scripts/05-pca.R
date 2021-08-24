## Descriptive analysis: PCA
## Rin Wada 22 June

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(FactoMineR)
source("functions.R")
source("graph_param.R")

# Main ----
# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

expo = readRDS(paste0("../Processed/",filepaths[4],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
covars = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))
print(all(rownames(expo)==rownames(covars)))
print(ncol(expo))

mypca=PCA(expo, graph = FALSE)

ev=mypca$eig[,2]

CreateScorePlot(mypca=mypca, type=as.character(covars$Batch),
                mycolours=batch.colours[levels(covars$Batch)],
                filename=paste0("../Figures/Section1/PCA_score_plot_batch_pooled3.pdf"),
                segments = FALSE, ellipse = TRUE)
CreateScorePlot(mypca=mypca, type=as.character(covars$Region),
                mycolours=region.colours[levels(covars$Region)],
                filename=paste0("../Figures/Section1/PCA_score_plot_region_pooled3.pdf"),
                segments = FALSE, ellipse = TRUE, wide = TRUE)

mycor=cor(expo, mypca$ind$coord)

annot_sub = annot[colnames(expo)]
{pdf(paste0("../Figures/Section1/PCA_correlation_circle_pooled3.pdf"), width=18, height=5)
  par(mar = c(5,5,1,1))
  layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(2,2,2,1.2))
  comp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
  for (k in 1:nrow(comp)){
    xcomp=comp[k,1]
    ycomp=comp[k,2]
    plot(mycor[,c(xcomp,ycomp)], xlim=c(-1.25,1.25), ylim=c(-1.25,1.25), cex=0.1, pch=19, las=1, cex.lab=1.5,
         xlab=paste0("Comp ",xcomp," (", round(ev[xcomp], digits=2), "% e.v.)"),
         ylab=paste0("Comp ",ycomp," (", round(ev[ycomp], digits=2), "% e.v.)"),
         col=annot.colours[annot_sub])
    arrows(x0=rep(0, nrow(mycor)), y0=rep(0, nrow(mycor)),
           x1=mycor[,xcomp], y1=mycor[,ycomp], length=0.1, col=annot.colours[annot_sub])
    abline(h=0, lty=2)
    abline(v=0, lty=2)
    xseq=seq(-1,1,length.out=10000)
    lines(xseq, sqrt(1-xseq^2))
    lines(xseq, -sqrt(1-xseq^2))
    text(mycor[,xcomp]+sign(mycor[,xcomp])*0.25, mycor[,ycomp],
         labels=colnames(expo), cex=0.8,
         col=darken(annot.colours[annot_sub]),0.5)
  }
  par(mar = c(0,1,1,0))
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
  legend("left", col=annot.colours[unique(annot_sub)],
         legend=names(annot.colours[unique(annot_sub)]),
         ncol = ceiling(length(annot.colours[unique(annot_sub)])/25),
         lty = 1, lwd = 2, bty = "n", cex = 1.5)
  dev.off()
}
                                    

# Supplementary ----
suffix = c("lux","fra","gs","pooled3","pooled2")
ifelse(dir.exists("../Figures/Supplementary"),"",dir.create("../Figures/Supplementary"))
ifelse(dir.exists("../Figures/Supplementary/Section1"),"",dir.create("../Figures/Supplementary/Section1"))
for (i in 1:length(batches)){
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  print(ncol(expo))

  mypca=PCA(expo, graph = FALSE)

  
  {pdf(paste0("../Figures/Supplementary/Section1/PCA_scree_plot_",suffix[i],".pdf"))
    par(mar=c(5,5,1,1))
    plot(mypca$eig[,3], col="navy", pch=19, cex=0.5, cex.lab=1.5, las=1,
         xlab="Number of PCs", ylab="Cumulative explained variance", ylim=c(0,100))
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
    points(mypca$eig[,3], col=batch.colours[i], pch=19, cex=0.5)
    dev.off()
  }
  ev=mypca$eig[,2]
  if (i !=4 ){
  CreateScorePlot(mypca=mypca, type=as.character(covars$Family.ID),
                  mycolours=family.colours[levels(covars$Family.ID)],
                  filename=paste0("../Figures/Supplementary/Section1/PCA_score_plot_",suffix[i],".pdf"))
  }
    if (i==5){
    CreateScorePlot(mypca=mypca, type=as.character(covars$Batch),
                    mycolours=batch.colours[levels(covars$Batch)],
                    filename=paste0("../Figures/Supplementary/Section1/PCA_score_plot_batch_",suffix[i],".pdf"),
                    segments = FALSE, ellipse = TRUE)
  }
  if (i==2){
    CreateScorePlot(mypca=mypca, type=as.character(covars$Region),
                    mycolours=region.colours[levels(covars$Region)],
                    filename=paste0("../Figures/Supplementary/Section1/PCA_score_plot_region_",suffix[i],".pdf"),
                    segments = FALSE, ellipse = TRUE, wide = TRUE)
  }

  mycor=cor(expo, mypca$ind$coord)
  annot_sub = annot[colnames(expo)]
  if (i !=4 ){
  {pdf(paste0("../Figures/Supplementary/Section1/PCA_correlation_circle_",suffix[i],".pdf"), width=18, height=5)
    par(mar = c(5,5,1,1))
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(2,2,2,1.2))
    comp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
    for (k in 1:nrow(comp)){
      xcomp=comp[k,1]
      ycomp=comp[k,2]
      plot(mycor[,c(xcomp,ycomp)], xlim=c(-1.25,1.25), ylim=c(-1.25,1.25), cex=0.1, pch=19, las=1, cex.lab=1.5,
           xlab=paste0("Comp ",xcomp," (", round(ev[xcomp], digits=2), "% e.v.)"),
           ylab=paste0("Comp ",ycomp," (", round(ev[ycomp], digits=2), "% e.v.)"),
           col=annot.colours[annot_sub])
      arrows(x0=rep(0, nrow(mycor)), y0=rep(0, nrow(mycor)),
             x1=mycor[,xcomp], y1=mycor[,ycomp], length=0.1, col=annot.colours[annot_sub])
      abline(h=0, lty=2)
      abline(v=0, lty=2)
      xseq=seq(-1,1,length.out=10000)
      lines(xseq, sqrt(1-xseq^2))
      lines(xseq, -sqrt(1-xseq^2))
      text(mycor[,xcomp]+sign(mycor[,xcomp])*0.25, mycor[,ycomp],
           labels=colnames(expo), cex=0.8,
           col=darken(annot.colours[annot_sub]),0.5)
    }
    par(mar = c(0,1,1,0))
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("left", col=annot.colours[unique(annot_sub)],
           legend=names(annot.colours[unique(annot_sub)]),
           ncol = ceiling(length(annot.colours[unique(annot_sub)])/25),
           lty = 1, lwd = 2, bty = "n", cex = 1.5)
    dev.off()
  }
  }
}
