## Pairwise: Characterisitng well/ill-classified families
## Rin Wada 15 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(focus)

# Initialisation
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

## Parameters
args=commandArgs(trailingOnly=TRUE)
m=as.numeric(args[1])

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

#### Fixed clustering ####
# Load data
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
fc = readRDS(paste0("../Results/",filepaths[m],"/Fixed_clusters.rds"))
expo  = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))
print(all(names(fc)==rownames(covars)))

families=unique(covars$Family.ID)

fclass=rep(NA, length(families))
for (f in 1:length(families)){
  tmp=fc[covars$Family.ID==families[f]]
  if (length(unique(tmp))==1){
    fclass[f] = 1
  } else {
    fclass[f] = 0
  }
}

X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))

### Univariate analysis ----
betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = fclass
  model1=lm(as.formula(f1))
  if(colnames(X)[k] == "gender_diff"){
    betas=c(betas, coefficients(model1)[2:length(coefficients(model1))])
    pvals=c(pvals, summary(model1)$coefficients[2:nrow(summary(model1)$coefficients),4])
  } else {
    betas=c(betas, coefficients(model1)[2])
    pvals=c(pvals, summary(model1)$coefficients[2,4])
  }
}
t1=Sys.time()
print(t1-t0)

mylabels = c(colnames(X)[1:which(colnames(X)=="gender_diff")-1],
             "All male (ref. All female)","Different gender (ref. All female)",
             colnames(X)[(which(colnames(X)=="gender_diff")+1):ncol(X)])

length(betas)==length(mylabels)

mylabels = gsub("age_mu","Mean age",mylabels)
mylabels = gsub("age_diff","Age difference",mylabels)
mylabels = gsub("length_mu","Mean sample length",mylabels)
mylabels = gsub("length_diff","Sample length difference",mylabels)
mylabels = gsub("weight_mu","Mean sample weight",mylabels)
mylabels = gsub("weight_diff","Sample weight difference",mylabels)
mylabels = gsub("region","Île-de-France (Y/N)",mylabels)


names(pvals)=names(betas)=mylabels

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Fixed_cluster_family_class_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Fixed_cluster_family_class_univar_betas.rds"))

### Multivariate analysis using stability selection sPLS regression ----

Y = fclass[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

X = scale(X)

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Fixed_cluster_family_class_multivar_output.rds"))


### Stranger metrics ----
compare_fc = sapply(fc, function(x) x==fc)
rownames(compare_fc)=colnames(compare_fc)=rownames(covars)
compare_fc = compare_fc[lower.tri(compare_fc)]

compare_family = sapply(covars$Family.ID, function(x) x==covars$Family.ID)
rownames(compare_family)=colnames(compare_family)=rownames(covars)
compare_family = compare_family[lower.tri(compare_family)]

pairs = sapply(rownames(covars), function(x) paste0(x,"--",rownames(covars)))
pairs = pairs[lower.tri(pairs)]

sclass = ifelse((!compare_fc & !compare_family), 1,
                ifelse((compare_fc & !compare_family), 0, NA))
names(sclass) = pairs
mygrep = !is.na(sclass)
sclass = sclass[mygrep]

### Univariate analysis----
X = readRDS(paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))

betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = sclass
  model1=lm(as.formula(f1))
  if(colnames(X)[k] == "gender_diff"){
    betas=c(betas, coefficients(model1)[2:length(coefficients(model1))])
    pvals=c(pvals, summary(model1)$coefficients[2:nrow(summary(model1)$coefficients),4])
  } else {
    betas=c(betas, coefficients(model1)[2])
    pvals=c(pvals, summary(model1)$coefficients[2,4])
  }
}
t1=Sys.time()
print(t1-t0)

mylabels = c(colnames(X)[1:which(colnames(X)=="gender_diff")-1],
             "All male (ref. All female)","Different gender (ref. All female)",
             colnames(X)[(which(colnames(X)=="gender_diff")+1):ncol(X)])

length(betas)==length(mylabels)

mylabels = gsub("age_mu","Mean age",mylabels)
mylabels = gsub("age_diff","Age difference",mylabels)
mylabels = gsub("length_mu","Mean sample length",mylabels)
mylabels = gsub("length_diff","Sample length difference",mylabels)
mylabels = gsub("weight_mu","Mean sample weight",mylabels)
mylabels = gsub("weight_diff","Sample weight difference",mylabels)
mylabels = gsub("region","Same region (Y/N)",mylabels)

names(pvals)=names(betas)=mylabels

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Fixed_cluster_stranger_class_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Fixed_cluster_stranger_class_univar_betas.rds"))

### Multivariate analysis using stability selection sPLS regression ----

Y = fclass[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

X = scale(X)

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Fixed_cluster_stranger_class_multivar_output.rds"))