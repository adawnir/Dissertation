## Univariate regression (including isolated children)
## Rin 12 July

# Load packages
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

### Exposure ~ Covariate ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  if (i==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else {X = covars %>% select(Age, Gender)}
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
  t0=Sys.time()
  for (k in 1:ncol(expo)){
    for (j in 1:ncol(X)){
      model1=lm(as.formula(f1))
      betas=c(betas, coefficients(model1)[2])
      pvals=c(pvals, summary(model1)$coefficients[2,4])
    }
  }
  t1=Sys.time()
  print(t1-t0)
  betas = matrix(betas, nrow = ncol(expo), ncol = ncol(X))
  pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
  rownames(pvals)=rownames(betas)=colnames(expo)
  colnames(pvals)=colnames(betas)=colnames(X)
  ifelse(dir.exists("../Results"),"",dir.create("../Results"))
  ifelse(dir.exists(paste0("../Results/",filepaths[i])),"",dir.create(paste0("../Results/",filepaths[i])))
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/univar_covar_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/univar_covar_betas.rds"))
  assign(paste0("betas_",suffix[i]),betas)
  assign(paste0("pvals_",suffix[i]),pvals)
}

## Heatmap
# Merge beta value matrices
values = merge(as.data.frame(betas_lux),as.data.frame(betas_fra), by=0, all=TRUE)
rownames(values) = values[,1]
values = merge(values[,-1], as.data.frame(betas_gs), by=0, all=TRUE)
rownames(values) = values[,1]
values = values[,-1]

# Merge p-value matrices and convert to asterisks for values smaller than Bonferroni threshold
labels = merge(as.data.frame(pvals_lux),as.data.frame(pvals_fra), by=0, all=TRUE)
rownames(labels) = labels[,1]
labels = merge(labels[,-1], as.data.frame(pvals_gs), by=0, all=TRUE)
rownames(labels) = labels[,1]
labels = labels[,-1]
labels[,1:4] = ifelse(labels[,1:4]<0.05/nrow(pvals_lux), "***", ifelse(labels[,1:4]<0.05, "*",""))
labels[,5:6] = ifelse(labels[,5:6]<0.05/nrow(pvals_fra), "***", ifelse(labels[,5:6]<0.05, "*",""))
labels[,7:8] = ifelse(labels[,7:8]<0.05/nrow(pvals_gs), "***", ifelse(labels[,7:8]<0.05, "*",""))
labels[is.na(labels)] = ""

# Sort rows
values = values[intersect(names(annot), rownames(values)),]
labels = labels[intersect(names(annot), rownames(labels)),]
annot_sub = annot[rownames(values)]

mat_row = data.frame(Batch = c(rep("Luxembourg",4), rep("France",2), rep("Grande-Synthe",2)))
rownames(mat_row) = make.unique(colnames(values))

mat_col = data.frame(Family = annot[rownames(values)])
rownames(mat_col) = rownames(values)
mat_colors = list(Batch = batch.colours[1:3], Family = annot.colours[unique(annot_sub)])
names(mat_colors$Family) = unique(annot_sub)
names(mat_colors$Batch) = batches[1:3]

{pdf("../Figures/Univariate_regression_covariates.pdf", width = 14)
  pheatmap(t(values),
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = t(labels),
           gaps_row =  c(4, 6),
           labels_row = c(colnames(betas_lux),colnames(betas_fra),colnames(betas_gs)),
           annotation_row = mat_row, annotation_col = mat_col,
           annotation_colors = mat_colors, annotation_names_col = FALSE,
           na_col = "grey80")
  dev.off()
  }

### Detection ~ Covariate ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_raw.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh.rds"))
  print(all(rownames(expo)==rownames(covars)))
  if (i==1){
    X = covars %>% select(Age, Gender, Length, Weight)
  } else {X = covars %>% select(Age, Gender)}
  
  # Convert matrix into binary
  expo = ifelse(expo=="nd",1,0)
  # Exclude chemicals with 0% non-detects
  expo = expo[,which(colSums(expo,na.rm = TRUE)!=0)]
  
  ncol(expo) * ncol(X) # number of tests
  betas = pvals = NULL
  f1='expo[,k] ~ X[,j]'
  t0=Sys.time()
  for (k in 1:ncol(expo)){
    for (j in 1:ncol(X)){
      model1=lm(as.formula(f1), family = "binomial")
      betas=c(betas, coefficients(model1)[2])
      pvals=c(pvals, summary(model1)$coefficients[2,4])
    }
  }
  t1=Sys.time()
  print(t1-t0)
  betas = matrix(betas, nrow = ncol(expo), ncol = ncol(X))
  pvals = matrix(pvals, nrow = ncol(expo), ncol = ncol(X))
  rownames(pvals)=rownames(betas)=colnames(expo)
  colnames(pvals)=colnames(betas)=colnames(X)
  saveRDS(pvals, paste0("../Results/",filepaths[i],"/univar_nd_covar_pvals.rds"))
  saveRDS(betas, paste0("../Results/",filepaths[i],"/univar_nd_covar_betas.rds"))
  assign(paste0("betas_",suffix[i]),betas)
  assign(paste0("pvals_",suffix[i]),pvals)
}

## Heatmap
# Merge beta value matrices
values = merge(as.data.frame(betas_lux),as.data.frame(betas_fra), by=0, all=TRUE)
rownames(values) = values[,1]
values = merge(values[,-1], as.data.frame(betas_gs), by=0, all=TRUE)
rownames(values) = values[,1]
values = values[,-1]

# Merge p-value matrices and convert to asterisks for values smaller than Bonferroni threshold
labels = merge(as.data.frame(pvals_lux),as.data.frame(pvals_fra), by=0, all=TRUE)
rownames(labels) = labels[,1]
labels = merge(labels[,-1], as.data.frame(pvals_gs), by=0, all=TRUE)
rownames(labels) = labels[,1]
labels = labels[,-1]
labels[,1:4] = ifelse(labels[,1:4]<0.05/nrow(pvals_lux), "***", ifelse(labels[,1:4]<0.05, "*",""))
labels[,5:6] = ifelse(labels[,5:6]<0.05/nrow(pvals_fra), "***", ifelse(labels[,5:6]<0.05, "*",""))
labels[,7:8] = ifelse(labels[,7:8]<0.05/nrow(pvals_gs), "***", ifelse(labels[,7:8]<0.05, "*",""))
labels[is.na(labels)] = ""

# Sort rows
values = values[intersect(names(annot), rownames(values)),]
labels = labels[intersect(names(annot), rownames(labels)),]
annot_sub = annot[rownames(values)]

mat_row = data.frame(Batch = c(rep("Luxembourg",4), rep("France",2), rep("Grande-Synthe",2)))
rownames(mat_row) = make.unique(colnames(values))

mat_col = data.frame(Family = annot[rownames(values)])
rownames(mat_col) = rownames(values)
mat_colors = list(Batch = batch.colours[1:3], Family = annot.colours[unique(annot_sub)])
names(mat_colors$Family) = unique(annot_sub)
names(mat_colors$Batch) = batches[1:3]

{pdf("../Figures/Univariate_nd_regression_covariates.pdf", width = 14)
  pheatmap(t(values),
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = t(labels),
           gaps_row =  c(4, 6),
           labels_row = c(colnames(betas_lux),colnames(betas_fra),colnames(betas_gs)),
           annotation_row = mat_row, annotation_col = mat_col,
           annotation_colors = mat_colors, annotation_names_col = FALSE,
           na_col = "grey80")
  dev.off()
}

### Detection frequency per participant ----
# Count non-detects for each participant
covar_lux$nd_count = apply(mat_lux, 1, function(x) sum(x=="nd", na.rm = T))
covar_fra$nd_count = apply(mat_fra, 1, function(x) sum(x=="nd", na.rm = T))
covar_gs$nd_count = apply(mat_gs, 1, function(x) sum(x=="nd", na.rm = T))

# Association with gender
mean(covar_lux$nd_count)
var(covar_lux$nd_count)
summary(covar_lux$nd_count)
glm1 <- glm(nd_count ~ Gender, poisson, covar_lux)
summary(glm1)

mean(covar_fra$nd_count)
var(covar_fra$nd_count)
summary(covar_fra$nd_count)
glm2 <- glm(nd_count ~ Gender, poisson, covar_fra)
summary(glm2)

mean(covar_gs$nd_count)
var(covar_gs$nd_count)
summary(covar_gs$nd_count)
glm3 <- glm(nd_count ~ Gender, poisson, covar_gs)
summary(glm3)

pdf("../Figures/Nd_count_Gender.pdf", width=12, height=4)
par(mfrow = c(1,3))
boxplot(nd_count ~ Gender, data = covar_lux, ylab = "Number of non-detects", main = "Luxembourg")
boxplot(nd_count ~ Gender, data = covar_fra, ylab = "Number of non-detects", main = "France")
boxplot(nd_count ~ Gender, data = covar_gs, ylab = "Number of non-detects", main = "Grande-Synthe")
dev.off()

# Association with age
glm1 <- glm(nd_count ~ Age, poisson, covar_lux)
summary(glm1)

glm2 <- glm(nd_count ~ Age, poisson, covar_fra)
summary(glm2)

glm3 <- glm(nd_count ~ Age, poisson, covar_gs)
summary(glm3)

pdf("../Figures/Nd_count_Age.pdf", width=12, height=4)
par(mfrow = c(1,3))
plot(nd_count ~ Age, data = covar_lux, pch = 19, ylab = "Number of non-detects", main = "Luxembourg")
plot(nd_count ~ Age, data = covar_fra, pch = 19, ylab = "Number of non-detects", main = "France")
plot(nd_count ~ Age, data = covar_gs, pch = 19, ylab = "Number of non-detects", main = "Grande-Synthe")
dev.off()

# Association with length of hair sample
summary(covar_lux$`Length (cm)`)
glm1 <- glm(nd_count ~ `Length (cm)`, poisson, covar_lux)
summary(glm1)

# Association with weight of hair sample
summary(covar_lux$`Weight (mg)`)
glm1 <- glm(nd_count ~ `Weight (mg)`, poisson, covar_lux)
summary(glm1)

pdf("../Figures/Nd_count_Sample.pdf", width=8, height=4)
par(mfrow = c(1,2))
plot(nd_count ~ `Length (cm)`, data = covar_lux, pch = 19, ylab = "Number of non-detects", main = "Length of hair sample (Luxembourg)")
plot(nd_count ~ `Weight (mg)`, data = covar_lux, pch = 19, ylab = "Number of non-detects", main = "Weight of hair sample (Luxembourg)")
dev.off()