## Characterisitng well/ill-classified families
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

fclass=NULL
for (f in 1:length(families)){
  mygrep = c(f, f*2)
  tmp=fc[covars$Family.ID==families[f]]
  if (length(unique(tmp))==1){
    fclass = c(fclass,rep(1,2))
  } else {
    fclass = c(fclass,rep(0,2))
  }
}

if (m==1){
  X = covars %>% select(Age, Gender, Length, Weight)
} else if(m %in% c(2,4)) {
  X = covars %>% select(Age, Gender, Region) %>%
    mutate(Region = ifelse(Region=="Île-de-France","Yes","No"))
} else {X = covars %>% select(Age, Gender)}

X = cbind(X, expo)

### Univariate analysis ----
betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = fclass
  model1=lm(as.formula(f1))
  betas=c(betas, coefficients(model1)[2])
  pvals=c(pvals, summary(model1)$coefficients[2,4])
}
t1=Sys.time()
print(t1-t0)

mylabels = colnames(X)

length(betas)==length(mylabels)

names(pvals)=names(betas)=mylabels

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Fixed_cluster_class_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Fixed_cluster_class_univar_betas.rds"))

### Multivariate analysis using stability selection sPLS regression ----

Y = fclass[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Fixed_cluster_class_multivar_output.rds"))


