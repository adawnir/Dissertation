## Descriptive analysis
## Rin Wada 22 June

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(FactoMineR)
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  ### Correlation matrix ----
  annot_sub = annot[colnames(expo)]
  mat_col = data.frame(Family = annot_sub)
  rownames(mat_col) = colnames(expo)
  mat_colors = list(Family = annot.colours[unique(annot_sub)])
  names(mat_colors$Family) = unique(annot_sub)
  
  cor = cor(expo, method = "spearman")
  ifelse(dir.exists(paste0("../Figures/",filepaths[i])),"",dir.create(paste0("../Figures/",filepaths[i])))
  {pdf(paste0("../Figures/",filepaths[i],"/Spearman_correlation_matrix.pdf"), width = 14)
    pheatmap(cor,
             show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
             cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE, 
             annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors)
    dev.off()
  }
  {pdf(paste0("../Figures/",filepaths[i],"/Spearman_correlation_matrix_clustered.pdf"), width = 14)
    pheatmap(cor,
             show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
             treeheight_row = 0, treeheight_col = 0,
             annotation_names_row = FALSE, annotation_names_col = FALSE, 
             annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors)
    dev.off()
    }
  # assign(paste0("expo_",suffix[i]),expo)
 
  ### PCA for visualisation ----
  mypca=PCA(expo)
  
  {pdf(paste0("../Figures/",filepaths[i],"/PCA_scree_plot.pdf"))
    par(mar=c(5,5,1,1))
    plot(mypca$eig[,3], col="navy", pch=19, cex=0.5, cex.lab=1.5, las=1,
         xlab="Number of PCs", ylab="Cumulative explained variance", ylim=c(0,100))
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
    points(mypca$eig[,3], col=batch.colours[i], pch=19, cex=0.5)
    dev.off()
  }
  ev=mypca$eig[,2]
  families=unique(covars$Family.ID)
  families.number = as.numeric(gsub(".*[A-Z](\\d+).*", "\\1", covars$Family.ID))
  mycolours=brewer.pal(n=12,name='Paired')
  mycolours=colorRampPalette(mycolours)(length(families[families!="Isolated"]))
  names(mycolours)=families[families!="Isolated"]
  CreateScorePlot(mypca=mypca, type=families.number, mycolours=mycolours,
                  filename=paste0("../Figures/",filepaths[i],"/PCA_score_plot.pdf"))
  mycor=cor(expo, mypca$ind$coord)

  {pdf(paste0("../Figures/",filepaths[i],"/PCA_correlation_circle.pdf"), width=14, height=5)
    par(mar=c(5,5,1,1), mfrow=c(1,3))
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
           labels=colnames(expo), cex=0.75, 
           col=darken(annot.colours[annot_sub]),0.5)
    }
    dev.off()
  }
}





