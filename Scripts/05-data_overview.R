## Density plots
## Rin Wada 10 July

# Load packages
library(tidyverse)

### Pre log transformation ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
extract_diff = readRDS("../Results/Pooled3/Chemical_compound_info_extract_diff.rds")
suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:3){
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh.rds"))
  assign(paste0("expo_",suffix[i]),expo)
}

# Density plots
mylabels = unique(c(colnames(expo_lux),colnames(expo_fra),colnames(expo_gs)))
annot_sub = annot[mylabels]
{pdf(paste0("../Figures/Compound_dist.pdf"), width = 14, height = 8)
  par(mfrow = c(4, 7), oma = c(0.5,0.5,0,0.5), mar = c(2,2,2,1), mgp=c(1,0.05,0))
  for (k in 1:length(annot_sub)){
    var = names(annot_sub)[k]
    print(var)
    ylim = NULL
    xlim = NULL
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        ylim = c(ylim, range(density(expo[,var], na.rm = T)$y))
        xlim = c(xlim,range(density(expo[,var], na.rm = T)$x))
      }
    }
    plot(NULL, main=var, xlab="pg/mg", ylab ="Density",
         xlim=c(min(xlim),max(xlim)),ylim=c(min(ylim),max(ylim)),
         cex.main = 1,  cex.lab=0.8, cex.axis=0.8, tck=-0.01, lwd = 0.5)
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        x=density(expo[,var], na.rm = T)
        lines(x, col=batch.colours[i], lwd = 0.5)
      }
    }
    # legend("topright", lwd=1, col=batch.colours, legend = batches, cex=0.5)
  }
  dev.off()
  }

annot_sub = annot[extract_diff$Compound[extract_diff$Compound %in% mylabels]]
length(annot_sub)

ifelse(dir.exists(paste0("../Figures/",filepaths[4])), "", dir.create(paste0("../Figures/",filepaths[4])))
{pdf(paste0("../Figures/Pooled3/Compound_dist_extract_diff.pdf"))
  par(mfrow = c(3, 3), oma = c(0.5,0.5,0,0.5), mar = c(2,2,2,1), mgp=c(1,0.05,0))
  for (k in 1:length(annot_sub)){
    var = names(annot_sub)[k]
    print(var)
    ylim = NULL
    xlim = NULL
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        ylim = c(ylim, range(density(expo[,var], na.rm = T)$y))
        xlim = c(xlim,range(density(expo[,var], na.rm = T)$x))
      }
    }
    plot(NULL, main=var, xlab="pg/mg", ylab ="Density",
         xlim=c(min(xlim, na.rm = T),max(xlim, na.rm = T)),
         ylim=c(min(ylim, na.rm = T),max(ylim, na.rm = T)),
         cex.main = 1,  cex.lab=0.8, cex.axis=0.8, tck=-0.01, lwd = 0.5)
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        x=density(expo[,var], na.rm = T)
        lines(x, col=batch.colours[i], lwd = 0.5)
      }
    }
    # legend("topright", lwd=1, col=batch.colours, legend = batches, cex=0.5)
  }
  dev.off()
}

### Post log transformation ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
extract_diff = readRDS("../Results/Pooled3/Chemical_compound_info_extract_diff.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:3){
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log.rds"))
  assign(paste0("expo_",suffix[i]),expo)
}

# Density plots
mylabels = unique(c(colnames(expo_lux),colnames(expo_fra),colnames(expo_gs)))
annot_sub = annot[mylabels]
{pdf(paste0("../Figures/Compound_dist_log.pdf"), width = 14, height = 8)
  par(mfrow = c(4, 7), oma = c(0.5,0.5,0,0.5), mar = c(2,2,2,1), mgp=c(1,0.05,0))
  for (k in 1:length(annot_sub)){
    var = names(annot_sub)[k]
    print(var)
    ylim = NULL
    xlim = NULL
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        ylim = c(ylim, range(density(expo[,var], na.rm = T)$y))
        xlim = c(xlim,range(density(expo[,var], na.rm = T)$x))
      }
    }
    plot(NULL, main=var, xlab=expression(log[10](pg/mg)), ylab ="Density",
         xlim=c(min(xlim),max(xlim)),ylim=c(min(ylim),max(ylim)),
         cex.main = 1, cex.lab=0.8, cex.axis=0.8, tck=-0.01, lwd = 0.5)
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        x=density(expo[,var], na.rm = T)
        lines(x, col=batch.colours[i], lwd = 0.5)
      }
    }
    # legend("topleft", lwd=1, col=batch.colours, legend = batches, cex=0.5)
  }
  dev.off()
}
annot_sub = annot[extract_diff$Compound[extract_diff$Compound %in% mylabels]]

{pdf(paste0("../Figures/Pooled3/Compound_dist_log_extract_diff.pdf"))
  par(mfrow = c(3, 3), oma = c(0.5,0.5,0,0.5), mar = c(2,2,2,1), mgp=c(1,0.05,0))
  for (k in 1:length(annot_sub)){
    var = names(annot_sub)[k]
    print(var)
    ylim = NULL
    xlim = NULL
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        ylim = c(ylim, range(density(expo[,var], na.rm = T)$y))
        xlim = c(xlim,range(density(expo[,var], na.rm = T)$x))
      }
    }
    plot(NULL, main=var, xlab="pg/mg", ylab ="Density",
         xlim=c(min(xlim, na.rm = T),max(xlim, na.rm = T)),
         ylim=c(min(ylim, na.rm = T),max(ylim, na.rm = T)),
         cex.main = 1,  cex.lab=0.8, cex.axis=0.8, tck=-0.01, lwd = 0.5)
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        x=density(expo[,var], na.rm = T)
        lines(x, col=batch.colours[i], lwd = 0.5)
      }
    }
    # legend("topright", lwd=1, col=batch.colours, legend = batches, cex=0.5)
  }
  dev.off()
}
### Post imputation ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
extract_diff = readRDS("../Results/Pooled3/Chemical_compound_info_extract_diff.rds")
suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:3){
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  assign(paste0("expo_",suffix[i]),expo)
}

# Density plots
mylabels = unique(c(colnames(expo_lux),colnames(expo_fra),colnames(expo_gs)))
annot_sub = annot[mylabels]
{pdf(paste0("../Figures/Compound_dist_log_naimp.pdf"), width = 14, height = 8)
  par(mfrow = c(4, 7), oma = c(0.5,0.5,0,0.5), mar = c(2,2,2,1), mgp=c(1,0.05,0))
  for (k in 1:length(annot_sub)){
    var = names(annot_sub)[k]
    print(var)
    ylim = NULL
    xlim = NULL
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        ylim = c(ylim, range(density(expo[,var], na.rm = T)$y))
        xlim = c(xlim,range(density(expo[,var], na.rm = T)$x))
      }
    }
    plot(NULL, main=var, xlab=expression(log[10](pg/mg)), ylab ="Density",
         xlim=c(min(xlim),max(xlim)),ylim=c(min(ylim),max(ylim)),
         cex.main = 1, cex.lab=0.8, cex.axis=0.8, tck=-0.01, lwd = 0.5)
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        x=density(expo[,var], na.rm = T)
        lines(x, col=batch.colours[i], lwd = 0.5)
      }
    }
    # legend("topleft", lwd=1, col=batch.colours, legend = batches, cex=0.5)
  }
  dev.off()
}

annot_sub = annot[extract_diff$Compound[extract_diff$Compound %in% mylabels]]
{pdf(paste0("../Figures/Pooled3/Compound_dist_log_naimp_extract_diff.pdf"))
  par(mfrow = c(3, 3), oma = c(0.5,0.5,0,0.5), mar = c(2,2,2,1), mgp=c(1,0.05,0))
  for (k in 1:length(annot_sub)){
    var = names(annot_sub)[k]
    print(var)
    ylim = NULL
    xlim = NULL
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        ylim = c(ylim, range(density(expo[,var], na.rm = T)$y))
        xlim = c(xlim,range(density(expo[,var], na.rm = T)$x))
      }
    }
    plot(NULL, main=var, xlab="pg/mg", ylab ="Density",
         xlim=c(min(xlim, na.rm = T),max(xlim, na.rm = T)),
         ylim=c(min(ylim, na.rm = T),max(ylim, na.rm = T)),
         cex.main = 1,  cex.lab=0.8, cex.axis=0.8, tck=-0.01, lwd = 0.5)
    for (i in 1:3){
      expo = eval(parse(text = paste0("expo_", suffix[i])))
      if (var %in% colnames(expo)){
        x=density(expo[,var], na.rm = T)
        lines(x, col=batch.colours[i], lwd = 0.5)
      }
    }
    # legend("topright", lwd=1, col=batch.colours, legend = batches, cex=0.5)
  }
  dev.off()
}


