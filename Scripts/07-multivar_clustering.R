## Multivariate clustering
## Rin Wada 12 July

# Load packages
library(ape)
library(igraph)
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
# Load data
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))

families=levels(covars$Family.ID)

expo = scale(expo)
d=dist(expo)
h=hclust(d, method = "complete")
print(all(covars$Indiv.ID==h$labels))
h$labels=paste0(covars$Family.ID, "-",h$labels)
myphylo=as.phylo(h)

{pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all.pdf"),width=14)
  par(mar=c(0,0,0,0))
  plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
       tip.color=family.colours[as.character(covars$Family.ID)])
  dev.off()
}

{pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_graph_expo_all.pdf"))
  g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[as.character(covars$Family.ID)])
  dev.off()
}

nobs=nrow(expo)
sp=shortest.paths(g)[1:nobs,1:nobs]
rownames(sp)=colnames(sp)=myphylo$tip.label

fsp=NULL
for (f in 1:length(families)){
  tmpmat=sp[covars$Family.ID==families[f],covars$Family.ID==families[f]]
  fsp=c(fsp,mean(tmpmat[upper.tri(tmpmat)]))
}
names(fsp)=families

# Save shortest path per family
saveRDS(fsp, paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))

{pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_cont_all.pdf"),width=14)
  par(mar=c(5,5,1,1))
  plot(fsp, pch=19, col=family.colours[families], xaxt="n", las=1,
       panel.first=abline(v=1:length(families), lty=3, col="grey"),
       xlab="Family ID", cex=2,
       ylab="Shortest path length between siblings", cex.lab=1.5)
  axis(side=1, at=1:length(families), labels=families, las = 2)
  dev.off()
}

### Family-wise covariates and delta exposures ----
X=matrix(NA, nrow = length(families))
for (f in 1:length(families)){
  X[f]=mean(covars$Age[covars$Family.ID==families[f]])
}
colnames(X)[ncol(X)] = "age_mu"
rownames(X) = families

tmp = NULL
for (f in 1:length(families)){
  tmp = c(tmp, abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2]))
}
X=cbind(as.data.frame(X), tmp)
colnames(X)[ncol(X)] = "age_diff"

if (m==1){
  tmp = NULL
  for (f in 1:length(families)){
    tmp = c(tmp, mean(covars$Weight[covars$Family.ID==families[f]]))
  }
  X=cbind(X, tmp)
  colnames(X)[ncol(X)] = "weight_mu"
  
  tmp = NULL
  for (f in 1:length(families)){
    
    tmp = c(tmp, abs(covars$Weight[covars$Family.ID==families[f]][1]-covars$Weight[covars$Family.ID==families[f]][2]))
  }
  X=cbind(X, tmp)
  colnames(X)[ncol(X)] = "weight_diff"
  
  tmp = NULL
  for (f in 1:length(families)){
    tmp = c(tmp, mean(covars$Length[covars$Family.ID==families[f]]))
  }
  X=cbind(X, tmp)
  colnames(X)[ncol(X)] = "length_mu"
  
  tmp = NULL
  for (f in 1:length(families)){
    
    tmp = c(tmp, abs(covars$Length[covars$Family.ID==families[f]][1]-covars$Length[covars$Family.ID==families[f]][2]))
  }
  X=cbind(X, tmp)
  colnames(X)[ncol(X)] = "length_diff"
}

tmp = NULL
for (f in 1:length(families)){
  gender=covars$Gender[covars$Family.ID==families[f]]
  if (sum(is.na(gender))>0) {
    tmp[f] = NA
  } else if (length(unique(gender))==1){
    if (unique(gender)=="Male"){
      tmp[f]="All male"
    }
    if (unique(gender)=="Female"){
      tmp[f]="All female"
    }
  } else {
    tmp[f]="Different gender"
  }
}
X=cbind(X, tmp)
colnames(X)[ncol(X)] = "gender_diff"

if (m %in% c(2,4)){
  tmp = NULL
  for (f in 1:length(families)){
    region=as.character(covars$Region)[covars$Family.ID==families[f]]
    if (unique(region)=="Île-de-France"){
      tmp[f] = "Yes"
    } else {tmp[f] = "No"}
    }
  X=cbind(X, tmp)
  colnames(X)[ncol(X)] = "region"
}

## Delta exposure
for (k in 1:ncol(expo)){
  for (f in 1:length(families)){
    delta = expo[covars$Family.ID==families[f],k]
    tmp[f] = abs(delta[1] - delta [2])
  }
  X=cbind(X, as.numeric(tmp))
  colnames(X)[ncol(X)] = colnames(expo)[k]
}

# Save covariates and exposures
saveRDS(X, paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))

# ### Geographical ----
# {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_batch.pdf"),width=14)
#   par(mar=c(0,0,0,0))
#   plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
#        tip.color=batch.colours[as.character(covars$Batch)])
#   dev.off()
# }
# {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_graph_expo_all_batch.pdf"))
#   g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = batch.colours[as.character(covars$Batch)])
#   dev.off()
# }
# if (m == 4){
#   families=unique(covars$Region)
#   mycolours = region.colours
#   names(mycolours)=families
#   
#   {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_region.pdf"),width=14)
#     par(mar=c(0,0,0,0))
#     plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
#          tip.color=region.colours[as.character(covars$Region)])
#     dev.off()
#   }
#   
#   {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_graph_expo_all_region.pdf"))
#     g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = region.colours[as.character(covars$Region)])
#     legend("bottomright", pch=19, col=region.colours[levels(covars$Region)],
#            legend=levels(covars$Region), cex = 0.4, ncol = 1)
#     dev.off()
#   }
#   
#   {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_depart.pdf"),width=14)
#     par(mar=c(0,0,0,0))
#     plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
#          tip.color=depart.colours[as.character(covars$Department)])
#     dev.off()
#   }
#   
#   {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_graph_expo_all_depart.pdf"))
#     g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = depart.colours[as.character(covars$Department)])
#     legend("bottomright", pch=19, col=depart.colours[levels(covars$Department)],
#            legend=levels(covars$Department), cex = 0.4, ncol = 1)
#     dev.off()
#   }
# }
