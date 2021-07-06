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
source("table_one_functions.R")

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

covar_pooled3 = readRDS("../Processed/Pooled3/Participant_covariate_info_subset.rds")
chem_pooled3 = readRDS("../Processed/Pooled3/Chemical_compound_info_subset.rds")
mat_pooled3 = readRDS("../Processed/Pooled3/Chemical_compound_matrix_subset_trans_imp.rds")

covar_pooled2 = readRDS("../Processed/Pooled2/Participant_covariate_info_subset.rds")
chem_pooled2 = readRDS("../Processed/Pooled2/Chemical_compound_info_subset.rds")
mat_pooled2 = readRDS("../Processed/Pooled2/Chemical_compound_matrix_subset_trans_imp.rds")

### (Pooled3 and Pooled 2) Table comparing means of pooled data sets (transformed and imputed) ----
mytable=NULL
for (k in 1:length(colnames(mat_pooled3))){
  print(colnames(mat_pooled3)[k])
  tmp=ContinuousTest3(x=mat_pooled3[,k], y=covar_pooled3$Batch)
  mytable=rbind(mytable, tmp)
  rownames(mytable)[nrow(mytable)]=colnames(mat_pooled3)[k]
}
colnames(mytable)=c(rep("Mean (sd)",3), "LUX vs FRA", "LUX vs GS", "FRA vs GS")
mytable[,4:6]=ReformatScientificNotation(mytable[,4:6])
ifelse(dir.exists("../Exports"), "", dir.create("../Exports"))
ifelse(dir.exists("../Exports/Pooled3"), "", dir.create("../Exports/Pooled3"))
SaveExcelWithSuperscripts(cbind(rownames(mytable),mytable), "../Exports/Pooled3/Table1.xlsx")

mytable=NULL
for (k in 1:length(colnames(mat_pooled2))){
  print(colnames(mat_pooled2)[k])
  tmp=ContinuousTest2(x=mat_pooled2[,k], y=covar_pooled2$Batch)
  mytable=rbind(mytable, tmp)
  rownames(mytable)[nrow(mytable)]=colnames(mat_pooled2)[k]
}
colnames(mytable)=c(rep("Mean (sd)",2), "LUX vs GS")
mytable[,3]=ReformatScientificNotation(mytable[,3])
ifelse(dir.exists("../Exports"), "", dir.create("../Exports"))
ifelse(dir.exists("../Exports/Pooled2"), "", dir.create("../Exports/Pooled2"))
SaveExcelWithSuperscripts(cbind(rownames(mytable),mytable), "../Exports/Pooled2/Table1.xlsx")

### Correlation matrix ----
mat_col = data.frame(Family = chem_lux$Family)
rownames(mat_col) = chem_lux$Compound
mat_colors = list(Family = family.colours[unique(chem_lux$Family)])
names(mat_colors$Family) = unique(chem_lux$Family)

cor_lux = cor(mat_lux, method = "spearman")

pdf("../Figures/Luxembourg/Corr_mat.pdf", width = 12, height = 9)
pheatmap(cor_lux,
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

pdf("../Figures/Luxembourg/Corr_mat_clustered.pdf", width = 12, height = 9)
pheatmap(cor_lux,
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         treeheight_row = 0, treeheight_col = 0,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

mat_col = data.frame(Family = chem_fra$Family)
rownames(mat_col) = chem_fra$Compound
mat_colors = list(Family = family.colours[unique(chem_fra$Family)])
names(mat_colors$Family) = unique(chem_fra$Family)

cor_fra = cor(mat_fra, method = "spearman")

pdf("../Figures/France/Corr_mat.pdf", width = 12, height = 9)
pheatmap(cor(mat_fra, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

pdf("../Figures/France/Corr_mat_clustered.pdf", width = 12, height = 9)
pheatmap(cor(mat_fra, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         treeheight_row = 0, treeheight_col = 0,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

mat_col = data.frame(Family = chem_gs$Family)
rownames(mat_col) = chem_gs$Compound
mat_colors = list(Family = family.colours[unique(chem_gs$Family)])
names(mat_colors$Family) = unique(chem_gs$Family)

cor_gs = cor(mat_gs, method = "spearman")

pdf("../Figures/GrandeSynthe/Corr_mat.pdf", width = 12, height = 9)
pheatmap(cor(mat_gs, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()

pdf("../Figures/GrandeSynthe/Corr_mat_clustered.pdf", width = 12, height = 9)
pheatmap(cor(mat_gs, method = "spearman"),
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         treeheight_row = 0, treeheight_col = 0,
         annotation_names_row = FALSE, annotation_names_col = FALSE, 
         annotation_row = mat_col,annotation_col = mat_col, annotation_colors = mat_colors,)
dev.off()





