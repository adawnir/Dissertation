### Community characterisation
### 15 Aug

# Load packages
library(igraph)
library(colorspace)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

for(m in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  expo  = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  lc = readRDS(paste0("../Results/",filepaths[m],"/Graphical_network_children_community.rds"))
  # lc = membership(lc)
  lc = lc$membership
  print(all(names(lc)==rownames(covars)))
  
  
  ### Heatmap ----
  
  z = lc
  # pvals = apply(expo, 2, function(y) summary(aov(y~z))[[1]][["Pr(>F)"]][1])
  
  annot_sub = annot[colnames(expo)]
  mat_col = data.frame(Class = annot_sub)
  rownames(mat_col) = colnames(expo)
  
  # Set cluster colour
  cluster_colour = alpha(rainbow(length(unique(lc))),0.5)
  names(cluster_colour) = paste0("C",1:length(unique(lc)))
  
  if(m %in% c(2,4,5)){
    mat_row = data.frame(Communities = paste0("C",lc),
                         Regions = as.character(covars$Region))
    rownames(mat_row) = rownames(expo)
    
    mat_colors = list(Class = annot.colours[unique(annot_sub)],
                      Communities = cluster_colour,
                      Regions = region.colours[levels(covars$Region)])
    names(mat_colors$Class) = unique(annot_sub)
    names(mat_colors$Communities) = paste0("C",1:length(unique(lc)))
    names(mat_colors$Region) = levels(covars$Region)
  } else{
    mat_row = data.frame(Communities = paste0("C",lc))
    rownames(mat_row) = rownames(expo)
    
    mat_colors = list(Class = annot.colours[unique(annot_sub)],
                      Communities = cluster_colour)
    names(mat_colors$Class) = unique(annot_sub)
    names(mat_colors$Communities) = paste0("C",1:length(unique(lc)))
  }
  
  mat = expo[order(lc),]
  
  if(ncol(expo)>40){
    pdf(paste0("../Figures/Section5/Community_characterisation_heatmap_",suffix[m],".pdf"), height = 10, width = 14)
  } else {pdf(paste0("../Figures/Section5/Community_characterisation_heatmap_",suffix[m],".pdf"), height = 7, width = 10)}
  pheatmap(t(mat), color = magma(100),
           border_color = NA,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           treeheight_row = 20, treeheight_col = 0,
           annotation_row = mat_col,annotation_col = mat_row,
           annotation_colors = mat_colors,
           legend = TRUE, annotation_legend = FALSE,
           annotation_names_row = FALSE, annotation_names_col = TRUE)
  dev.off()
  
  ### Boxplots ----
  
  # z = paste0("C",lc)
  # # Group satellites
  # z[z %in% names(which(table(z) == 1))] = "Satellites"
  # order_lc = order(unique(lc[which(z!= "Satellites")]))
  # if(sum(z == "Satellites")>0){
  #   z = factor(z, levels = c(unique(z[which(z!= "Satellites")])[order_lc],"Satellites"))
  # } else {
  #   z = factor(z, levels = unique(z)[order_lc])
  # }
  # 
  # 
  # # Set cluster colour
  # cluster_colour = rainbow(length(unique(lc)))
  # names(cluster_colour) = paste0("C",1:length(unique(lc)))
  # 
  # ifelse(dir.exists(paste0("../Figures/",filepaths[m],"/Community_characterisation")),"",
  #        dir.create(paste0("../Figures/",filepaths[m],"/Community_characterisation")))
  # 
  # pvals = apply(expo, 2, function(y) summary(aov(y~z))[[1]][["Pr(>F)"]][1])
  # 
  # for (p in 1:ncol(expo)){
  #   y = 10^expo[,p]
  #   if(colnames(expo)[p] %in% c("Fipronil", "Permethrin", names(pvals)[pvals<0.05/length(pvals)])){
  #     pdf(paste0("../Figures/Supplementary/Boxplot_community_characterisation_",colnames(expo)[p],"_",suffix[m],".pdf"), height = 5, width = 10,)
  #   } else{
  #     pdf(paste0("../Figures/",filepaths[m],"/Community_characterisation/Boxplot_community_characterisation_",colnames(expo)[p],".pdf"), height = 5, width = 10)
  #   }
  #   options(scipen = 999)
  #   par(mar=c(5,5,1,1))
  #   mycol=c(cluster_colour[which(names(cluster_colour) %in% levels(z))],"grey80")
  #   names(mycol) = levels(z)
  #   boxplot(y~z,pch=19, cex=0.5, col=alpha(mycol,0.3),
  #           border=darken(mycol,0.5),
  #           las=1, xaxt="n", outline=FALSE,
  #           medcol='white', whiskcol='black', staplecol='black', outcol='black',
  #           xlab="Communities", ylab=substitute(tmp~(pg/mg),list(tmp = colnames(expo)[p])),
  #           cex.lab=1.5, log = "y")
  #   axis(side=1, at=1:length(levels(z)), labels=levels(z), cex.axis=0.8, las=2)
  #   dev.off()
  # }
  # saveRDS(pvals, paste0("../Results/",filepaths[m],"/Community_characterisation_pvals.rds"))
}

