## Multivariate clustering - Association with comembership proportion
## Rin Wada 7 Aug

# Load packages
library(focus)
library(colorspace)
library(RColorBrewer)

# Initialise
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

## Parameters
args=commandArgs(trailingOnly=TRUE)
m=as.numeric(args[1])

#### Between siblings ####
# Load data sets
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
out = readRDS(paste0("../Results/",filepaths[m],"/Consensus_clustering_output.rds"))

selprop = SelectionProportions(out)
families=levels(covars$Family.ID)

fprop=NULL
for (f in 1:length(families)){
  tmpmat=selprop[covars$Family.ID==families[f],covars$Family.ID==families[f]]
  fprop=c(fprop,mean(tmpmat[upper.tri(tmpmat)]))
}
names(fprop)=families

saveRDS(fprop, paste0("../Results/",filepaths[m],"/Comembership_prop.rds"))

### Univariate analysis ---
betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = fprop
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

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Comembership_prop_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Comembership_prop_univar_betas.rds"))

### Multivariate analysis using stability selection sPLS regression ----

Y = fprop[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

X = scale(X)

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Comembership_prop_multivar_output.rds"))

#### Between Strangers ####
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
X = readRDS(paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))
out = readRDS(paste0("../Results/",filepaths[m],"/Consensus_clustering_output.rds"))

compare_family = sapply(covars$Family.ID, function(x) x==covars$Family.ID)
rownames(compare_family)=colnames(compare_family)=rownames(covars)
compare_family = compare_family[lower.tri(compare_family)]

pairs = sapply(rownames(covars), function(x) paste0(x,"--",rownames(covars)))
pairs = pairs[lower.tri(pairs)]

mygrep = which(!compare_family)
pairs = pairs[mygrep]

selprop = SelectionProportions(out)

tmp=as.vector(combn(rownames(covars), 2, function(x) selprop[x[1],x[2]]))
sprop=tmp[mygrep]
names(sprop)=pairs

saveRDS(sprop, paste0("../Results/",filepaths[m],"/Comembership_prop_stranger.rds"))

### Univariate analysis ----
betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = sprop
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

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_univar_betas.rds"))

### Multivariate analysis using stability selection sPLS regression ----

Y = sprop[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

X = scale(X)

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Comembership_prop_stranger_multivar_output.rds"))
