## Descriptive analysis: Correlation matrix
## Rin Wada 22 June

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(grid)
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

### Correlation matrix ----
for (i in 1:length(batches)){
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  print(ncol(expo))
  
  ### Correlation matrix ----
  annot_sub = annot[colnames(expo)]
  mat_col = data.frame(Family = annot_sub)
  rownames(mat_col) = colnames(expo)
  mat_colors = list(Family = annot.colours[unique(annot_sub)])
  names(mat_colors$Family) = unique(annot_sub)
  
  cor = cor(expo, method = "spearman")
  
  ifelse(dir.exists("../Figures/Section1"),"",dir.create("../Figures/Section1"))
  if (nrow(cor)>50){
    width = 11
    height = 11
  }
  if (nrow(cor)<50){
    width = 6
    height = 6
  }
  mybreaks = seq(-1, 1, length = 100)
  {pdf(paste0("../Figures/Section1/Spearman_correlation_matrix_",suffix[i],".pdf"), width = width, height = height)
    pheatmap(cor,
             cellwidth = 8, cellheight = 8,
             border_color = NA, breaks = mybreaks,
             treeheight_row = 0, treeheight_col = 0,
             annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,
             legend = FALSE, annotation_legend = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE)
    dev.off()
  }
  
  tmp = cor(expo, method = "spearman")
  tmp[lower.tri(tmp, diag = T)] = NA
  high_corr = cbind(rownames(tmp)[which(abs(tmp) > 0.8, arr.ind = TRUE)[,1]],
                    colnames(tmp)[which(abs(tmp) > 0.8, arr.ind = TRUE)[,2]])
  ifelse(dir.exists(paste0("../Exports/",filepaths[i])),"",dir.create(paste0("../Exports/",filepaths[i])))
  write.csv(high_corr, file = paste0("../Exports/",filepaths[i],"/High_correlation_compound_pairs.csv"))

  # ifelse(dir.exists(paste0("../Figures/",filepaths[i])),"",dir.create(paste0("../Figures/",filepaths[i])))
  # if (nrow(cor)>50){
  #   width = 14
  #   height = 12
  # }
  # if (nrow(cor)<50){
  #   width = 8.5
  #   height = 7
  # }
  # {pdf(paste0("../Figures/",filepaths[i],"/Spearman_correlation_matrix.pdf"), width = width, height = height)
  #   pheatmap(cor,
  #            cellwidth = 8, cellheight = 8,
  #            border_color = NA,
  #            cluster_rows = FALSE, cluster_cols = FALSE,
  #            annotation_names_row = FALSE, annotation_names_col = FALSE, 
  #            annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors)
  #   dev.off()
  #   }
  # {pdf(paste0("../Figures/",filepaths[i],"/Spearman_correlation_matrix_clustered.pdf"), width = width, height = height)
  #   pheatmap(cor,
  #            cellwidth = 8, cellheight = 8,
  #            border_color = NA,
  #            treeheight_row = 0, treeheight_col = 0,
  #            annotation_names_row = FALSE, annotation_names_col = FALSE, 
  #            annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors)
  #   dev.off()
  #   }
}

## Legend for correlation matrices
draw_legend = function(color, breaks, legend, ...){
  color = color[!is.infinite(breaks)]
  breaks = breaks[!is.infinite(breaks)]
  
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  
  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos = height * legend_pos + (unit(1, "npc") - height)
  
  breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
  breaks = height * breaks + (unit(1, "npc") - height)
  
  h = breaks[-1] - breaks[-length(breaks)]
  
  rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
  title = textGrob("Spearman's correlation coefficient", x = 0, y = 1.02, hjust = 0.5, vjust = 0, gp = gpar(...))
  
  res = grobTree(rect, text, title,
                 vp = viewport(x = unit(0.7, "npc"), height=unit(0.8, "npc"),width=unit(4, "cm")))
  
  return(res)
}

legend = grid.pretty(range(as.vector(mybreaks)))
names(legend) = legend
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

{pdf("../Figures/Section1/Spearman_correlation_matrix_legend.pdf", height = 8.5, width = 4)
grid.newpage()
grid.draw(draw_legend(color, mybreaks, legend))
for(e in 1:length(annot.colours)) {
  grid.rect(
    2,
    e*0.55+2,
    0.5,
    0.5,
    default.units = "cm",
    gp=gpar( 
      fill=rev(annot.colours)[e],
      lwd=0
    )
  )
  grid.text(
    label = rev(names(annot.colours))[e],
    x = 2.5,
    y = e*0.55+2,
    default.units = "cm",
    hjust = 0
  )
}
dev.off()
}

# Inspect highly correlated compounds
for (i in 1:length(batches)){
  high_corr = read.csv(paste0("../Exports/",filepaths[i],"/High_correlation_compound_pairs.csv"),row.names = 1)
  for (k in 1:nrow(high_corr)){
    print(unlist(high_corr[k,]))
    for (m in 1:length(batches)){
    print(batches[m])
    expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
    if(all(unlist(high_corr[k,]) %in% colnames(expo))){print(cor(expo[,unlist(high_corr[k,])], method = "spearman")[1,2])}
    }
  }
}
