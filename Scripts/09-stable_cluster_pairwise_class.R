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

#### Stable clustering ####
# Load data
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
sc = readRDS(paste0("../Results/",filepaths[m],"/Stable_clusters.rds"))
expo  = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))
print(all(names(sc)==rownames(covars)))

families=unique(covars$Family.ID)

sclass=rep(NA, length(families))
for (f in 1:length(families)){
  tmp=sc[covars$Family.ID==families[f]]
  if (length(unique(tmp))==1){
    sclass[f] = 1
  } else {
    sclass[f] = 0
  }
}

X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))

### Univariate analysis ----
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

names(pvals)=names(betas)=mylabels

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Stable_cluster_family_class_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Stable_cluster_family_class_univar_betas.rds"))

### Plotting ----
# annot_sub = annot[which(names(annot) %in% mylabels)]
# 
# mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
# 
# {pdf(paste0("../Figures/",filepaths[m],"/Stable_cluster_family_class_univariate.pdf"))
#   par(mar=c(5,5,1,1))
#   plot(betas, -log10(pvals), pch=19,
#        col=ifelse(pvals < 0.05, mycolours, "grey"),
#        cex.lab=1.5, cex = 0.7,
#        ylim = c(0, max(-log10(pvals))+0.25),
#        xlim = c(-max(abs(betas))-0.05, max(abs(betas))+0.05),
#        ylab=expression(-log[10](p)), 
#        xlab=expression(beta))
#   for (k in 1:length(mylabels)){
#     if(mylabels[k] %in% names(annot_sub)){
#       if(pvals[k] < 0.05){
#         label = substitute(Delta~tmp, list(tmp=mylabels[k]))
#         text(betas[k]+sign(betas[k])*0.01, -log10(pvals[k])+0.1,
#              labels = label, col = mycolours[k])
#       }
#     }
#     else {
#       text(betas[k]+sign(betas[k])*0.01, -log10(pvals[k])+0.1,
#            labels = mylabels[k], col = mycolours[k]) 
#     }
#   }
#   text(max(abs(betas))+0.04, -log10(0.05/ncol(X))+0.01,
#        paste0("Bonferroni threshold = ",formatC(0.05/ncol(X), digits = 2, format = "e")),
#        adj = c(1,0), col = "darkred")
#   abline(h = -log10(0.05/ncol(X)), lty = 2, col = "darkred", cex = 0.6)
#   dev.off()
# }

### Multivariate analysis using stability selection sPLS regression ----

Y = sclass[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Stable_cluster_class_multivar_output.rds"))


### Pair-wise covariates and delta exposures ----
compare_sc = sapply(sc, function(x) x==sc)
rownames(compare_sc)=colnames(compare_sc)=rownames(covars)
compare_sc = compare_sc[lower.tri(compare_sc)]

compare_family = sapply(covars$Family.ID, function(x) x==covars$Family.ID)
rownames(compare_family)=colnames(compare_family)=rownames(covars)
compare_family = compare_family[lower.tri(compare_family)]

pairs = sapply(rownames(covars), function(x) paste0(x,"--",rownames(covars)))
pairs = pairs[lower.tri(pairs)]

sclass = ifelse((!compare_sc & !compare_family), 1,
                ifelse((compare_sc & !compare_family), 0, NA))
names(sclass) = pairs
mygrep = !is.na(sclass)
sclass = sclass[mygrep]

# Load covariates and exposures
X = readRDS(paste0("../Results/",filepaths[m],"/Stranger_covariates_delta_exposures.rds"))


### Univariate analysis----
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

names(pvals)=names(betas)=mylabels

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Stable_cluster_stranger_class_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Stable_cluster_stranger_class_univar_betas.rds"))

### Plotting ----
annot_sub = annot[which(names(annot) %in% mylabels)]

mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])

{pdf(paste0("../Figures/",filepaths[m],"/Stable_cluster_class_univariate.pdf"))
  par(mar=c(5,5,1,1))
  plot(betas, -log10(pvals), pch=19,
       col=ifelse(pvals < 0.05, mycolours, "grey"),
       cex.lab=1.5, cex = 0.7,
       ylim = c(0, max(-log10(pvals))+0.25),
       xlim = c(-max(abs(betas))-0.05, max(abs(betas))+0.05),
       ylab=expression(-log[10](p)),
       xlab=expression(beta))
  for (k in 1:length(mylabels)){
    if(mylabels[k] %in% names(annot_sub)){
      if(pvals[k] < 0.05/ncol(X)){
        label = substitute(Delta~tmp, list(tmp=mylabels[k]))
        text(betas[k]+sign(betas[k])*0.01, -log10(pvals[k])+0.1,
             labels = label, col = mycolours[k])
      }
    }
    else {
      text(betas[k]+sign(betas[k])*0.01, -log10(pvals[k])+0.1,
           labels = mylabels[k], col = mycolours[k])
    }
  }
  text(max(abs(betas))+0.04, -log10(0.05/ncol(X))+0.01,
       paste0("Bonferroni threshold = ",formatC(0.05/ncol(X), digits = 2, format = "e")),
       adj = c(1,0), col = "darkred", cex = 0.6)
  abline(h = -log10(0.05/ncol(X)), lty = 2, col = "darkred")
  dev.off()
}
