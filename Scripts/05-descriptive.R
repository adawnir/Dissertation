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
  print(all(rownames(expo)==rownames(covars)))
  
  ### Correlation matrix ----
  annot_sub = annot[colnames(expo)]
  mat_col = data.frame(Family = annot_sub)
  rownames(mat_col) = colnames(expo)
  mat_colors = list(Family = annot.colours[unique(annot_sub)])
  names(mat_colors$Family) = unique(annot_sub)
  
  cor = cor(expo, method = "spearman")
  
  tmp = cor
  tmp[lower.tri(tmp, diag = T)] = NA
  high_corr = cbind(rownames(tmp)[which(abs(tmp) > 0.8, arr.ind = TRUE)[,1]],
                    colnames(tmp)[which(abs(tmp) > 0.8, arr.ind = TRUE)[,2]])
  ifelse(dir.exists(paste0("../Exports/",filepaths[i])),"",dir.create(paste0("../Exports/",filepaths[i])))
  write.csv(high_corr, file = paste0("../Exports/",filepaths[i],"/High_correlation_compound_pairs.csv"))
  
  ifelse(dir.exists(paste0("../Figures/",filepaths[i])),"",dir.create(paste0("../Figures/",filepaths[i])))
  if (nrow(cor)>50){
    width = 14
    height = 12
  }
  if (nrow(cor)<50){
    width = 8.5
    height = 7
  }
  {pdf(paste0("../Figures/",filepaths[i],"/Spearman_correlation_matrix.pdf"), width = width, height = height)
    pheatmap(cor,
             cellwidth = 8, cellheight = 8,
             border_color = NA,
             cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE, 
             annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors)
    dev.off()
    }
  {pdf(paste0("../Figures/",filepaths[i],"/Spearman_correlation_matrix_clustered.pdf"), width = width, height = height)
    pheatmap(cor,
             cellwidth = 8, cellheight = 8,
             border_color = NA,
             treeheight_row = 0, treeheight_col = 0,
             annotation_names_row = FALSE, annotation_names_col = FALSE, 
             annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors)
    dev.off()
    }
  ### PCA for visualisation ----
  mypca=PCA(expo, graph = FALSE)
  
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
  CreateScorePlot(mypca=mypca, type=as.character(covars$Family.ID),
                  mycolours=family.colours[levels(covars$Family.ID)],
                  filename=paste0("../Figures/",filepaths[i],"/PCA_score_plot.pdf"))
  if (i %in% 4:5){
    CreateScorePlot(mypca=mypca, type=as.character(covars$Batch),
                    mycolours=batch.colours[levels(covars$Batch)],
                    filename=paste0("../Figures/",filepaths[i],"/PCA_score_plot_batch.pdf"),
                    segments = FALSE, ellipse = TRUE)
  }
  if (i%in%c(2,4)){
    CreateScorePlot(mypca=mypca, type=as.character(covars$Region),
                    mycolours=region.colours[levels(covars$Region)],
                    filename=paste0("../Figures/",filepaths[i],"/PCA_score_plot_region.pdf"),
                    segments = FALSE, ellipse = TRUE)
    CreateScorePlot(mypca=mypca, type=as.character(covars$Department),
                    mycolours=depart.colours[levels(covars$Department)],
                    filename=paste0("../Figures/",filepaths[i],"/PCA_score_plot_depart.pdf"),
                    segments = FALSE, ellipse = TRUE)
  }
  
  mycor=cor(expo, mypca$ind$coord)

  {pdf(paste0("../Figures/",filepaths[i],"/PCA_correlation_circle.pdf"), width=17.5, height=5)
    par(mar = c(5,5,1,1))
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(2,2,2,1))
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
    par(mar = c(1,1,1,1))
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("left", col=annot.colours[unique(annot_sub)],
           legend=names(annot.colours[unique(annot_sub)]),
           ncol = ceiling(length(annot.colours[unique(annot_sub)])/25),
           lty = 1, pt.cex=1.2, bty = "n", cex = 1.2)
    dev.off()
  }
}

# Inspect highly correlated compounds
expo = readRDS(paste0("../Processed/",filepaths[1],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
chem = readRDS(paste0("../Processed/",filepaths[1],"/Chemical_compound_info_thresh.rds"))
high_corr = read.csv(paste0("../Exports/",filepaths[1],"/High_correlation_compound_pairs.csv"),
                     row.names = 1)
high_corr[1,]
high_corr[2,]
high_corr[3,]
chem["PCB-101",]
chem["Dimethachlor",]
cor(expo[,rownames(chem)[chem$nd_prop>0.89]], method = "spearman")

expo = readRDS(paste0("../Processed/",filepaths[2],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
chem = readRDS(paste0("../Processed/",filepaths[2],"/Chemical_compound_info_thresh.rds"))
high_corr = read.csv(paste0("../Exports/",filepaths[2],"/High_correlation_compound_pairs.csv"),
                     row.names = 1)
high_corr[1,]
high_corr[2,]
high_corr[3,]
chem["Propiconazole",]
chem["Trifloxystrobin",]
cor(expo[,rownames(chem)[chem$nd_prop+chem$NA_prop>0.84]], method = "spearman")
chem["Propiconazole","nd_prop"]+chem["Propiconazole","NA_prop"]
chem["Trifloxystrobin","nd_prop"]+chem["Trifloxystrobin","NA_prop"]

expo = readRDS(paste0("../Processed/",filepaths[3],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
chem = readRDS(paste0("../Processed/",filepaths[3],"/Chemical_compound_info_thresh.rds"))
high_corr = read.csv(paste0("../Exports/",filepaths[3],"/High_correlation_compound_pairs.csv"),
                     row.names = 1)
high_corr[1,]
chem["Pendimethalin",]
chem["Prosulfocarb",]

expo = readRDS(paste0("../Processed/",filepaths[4],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
chem = readRDS(paste0("../Processed/",filepaths[4],"/Chemical_compound_info_thresh.rds"))
high_corr = read.csv(paste0("../Exports/",filepaths[4],"/High_correlation_compound_pairs.csv"),
                     row.names = 1)
high_corr[1,]
high_corr[2,]
high_corr[3,]

expo = readRDS(paste0("../Processed/",filepaths[5],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
chem = readRDS(paste0("../Processed/",filepaths[5],"/Chemical_compound_info_thresh.rds"))
high_corr[1,]
high_corr[2,]
