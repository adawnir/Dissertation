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
source("function.R")

# Load data sets
covar_lux = readRDS("../Processed/Luxembourg/Participant_covariate_info_subset.rds")
chem_lux = readRDS("../Processed/Luxembourg/Chemical_compound_info_subset.rds")
mat_lux = readRDS("../Processed/Luxembourg/Chemical_compound_matrix_subset_trans_imp.rds")

covar_fra = readRDS("../Processed/France/Participant_covariate_info_subset.rds")
chem_fra = readRDS("../Processed/France/Chemical_compound_info_subset.rds")
mat_fra = readRDS("../Processed/France/Chemical_compound_matrix_subset_trans_imp.rds")

covar_gs = readRDS("../Processed/GrandeSynthe/Participant_covariate_info_subset.rds")
chem_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_info_subset.rds")
mat_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_matrix_subset_trans.rds")

### Correlation matrix ----
mat_col = data.frame(Family = chem_lux$Family)
rownames(mat_col) = chem_lux$Compound
mat_colors = list(Family = family.colours)
names(mat_colors$Family) = unique(chem_lux$Family)

cor_lux = cor(mat_lux, method = "spearman")

pdf("../Figures/Luxembourg/Corr_mat.pdf", width = 11, height = 9)
pheatmap(cor_lux,
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

pdf("../Figures/Luxembourg/Corr_mat_clustered.pdf", width = 12, height = 10)
pheatmap(cor_lux,
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

mat_col = data.frame(Family = chem_fra$Family)
rownames(mat_col) = chem_fra$Compound
mat_colors = list(Family = family.colours[1:length(levels(chem_fra$Family))])

cor_fra = cor(mat_fra, method = "spearman")

names(mat_colors$Family) = unique(chem_fra$Family)
pdf("../Figures/France/Corr_mat.pdf", width = 11, height = 9)
pheatmap(cor(mat_fra, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

pdf("../Figures/France/Corr_mat_clustered.pdf", width = 12, height = 10)
pheatmap(cor(mat_fra, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

mat_col = data.frame(Family = chem_gs$Family)
rownames(mat_col) = chem_gs$Compound
mat_colors = list(Family = family.colours)
names(mat_colors$Family) = unique(chem_gs$Family)

cor_gs = cor(mat_gs, method = "spearman")

pdf("../Figures/GrandeSynthe/Corr_mat.pdf", width = 11, height = 9)
pheatmap(cor(mat_gs, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

pdf("../Figures/GrandeSynthe/Corr_mat_clustered.pdf", width = 12, height = 10)
pheatmap(cor(mat_gs, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()




