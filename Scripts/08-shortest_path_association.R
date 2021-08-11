## Multivariate clustering - Association with shortest path
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

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
fsp = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))

betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = fsp
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

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Shortest_path_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Shortest_path_univar_betas.rds"))

### Multivariate analysis using stability selection sPLS regression ----
Y = fsp[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Shortest_path_multivar_output.rds"))

