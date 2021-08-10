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

#### Community clustering ####
# Load data
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
lc = readRDS(paste0("../Results/",filepaths[m],"/Graphical_network_children_community.rds"))
expo  = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))
print(all(names(lc)==rownames(covars)))

families=unique(covars$Family.ID)
lc = membership(lc)

lclass=NULL
for (f in 1:length(families)){
  mygrep = c(f, f*2)
  tmp=lc[covars$Family.ID==families[f]]
  if (length(unique(tmp))==1){
    lclass = c(lclass,rep(1,2))
  } else {
    lclass = c(lclass,rep(0,2))
  }
}

if (m==1){
  X = covars %>% select(Age, Gender, Length, Weight)
} else {X = covars %>% select(Age, Gender)}

X = cbind(X, expo)

### Univariate analysis ----
betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = lclass
  model1=lm(as.formula(f1))
  betas=c(betas, coefficients(model1)[2])
  pvals=c(pvals, summary(model1)$coefficients[2,4])
}
t1=Sys.time()
print(t1-t0)

mylabels = colnames(X)

length(betas)==length(mylabels)

names(pvals)=names(betas)=mylabels

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Community_cluster_class_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Community_cluster_class_univar_betas.rds"))

### Plotting ----
# annot_sub = annot[which(names(annot) %in% mylabels)]
# 
# mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])
# 
# {pdf(paste0("../Figures/",filepaths[m],"/Community_cluster_family_class_univariate.pdf"))
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

Y = lclass[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Comembership_prop_multivar_output.rds"))


