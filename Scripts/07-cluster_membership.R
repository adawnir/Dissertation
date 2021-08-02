## Clustering membership
## Rin Wada 22 July

# Load packages
library(tidyverse)
library(RColorBrewer)
library(focus)
library(summarise)

### Cluster membership ascertainment ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))

  # Standardisation
  expo = scale(expo)
  
  # Dissimilarity matrix
  d=dist(expo)
  
  # Hierarchical clustering using Complete Linkage
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))
  h$labels=paste0(covars$Family.ID, "-",h$labels)
  
  ### Focus: Stability selection based HC
  # out = Clustering(expo, K = 100, tau = 0.5, seed = 290621)
  # 
  # # Save outputs
  # saveRDS(out, paste0("../Results/",filepaths[i],"/Stability_clustering_output.rds"))
  # 
  # pdf(paste0("../Figures/",filepaths[i],"/Stability_clustering_output.pdf"))
  # par(mar=c(7, 5, 7, 6))
  # CalibrationPlot(out)
  # dev.off()
  
  out = readRDS(paste0("../Results/",filepaths[i],"/Stability_clustering_output.rds"))
  covars$stab_cluster = Clusters(out)
  
  ### Summarise: HC with F-test based cluster score
  score = NULL
  for(k in seq(2,nrow(expo)-1)){
    member = cutree(h, k=k)
    score = c(score,ClusteringScore(expo, member)$score)
  }
  
  pdf(paste0("../Figures/",filepaths[i],"/F-test_clustering_score.pdf"))
  plot(y = score, x = seq(2,nrow(expo)-1),
       pch = 19,
       col = "navy", type = "b",
       ylab = "F-test based clustering score",
       xlab = "Number of clusters")
  dev.off()
  
  # set k
  lambda = seq(2,nrow(expo)-1)[which.max(score)]
  
  # cut trees
  covars$score_cluster = cutree(h, k=lambda)
  
  ### Fixed: k = # of families
  covars$fixed_cluster = cutree(h, k=length(unique(covars$Family.ID)))
  
  if (i %in% 4:5){
    covars$fixed_cluster_batch = cutree(h, k=length(unique(covars$Batch)))
  }
  if (i == 4){
    covars$fixed_cluster_country = cutree(h, k=length(unique(covars$Country)))
    covars$fixed_cluster_region = cutree(h, k=length(unique(covars$Region)))
    covars$fixed_cluster_department = cutree(h, k=length(unique(covars$Department)))
  }
  
  # Save memberships
  saveRDS(covars, paste0("../Results/",filepaths[i],"/Cluster_memberships.rds"))
}

### Cluster membership evaluation ----
# Load packages
library(fossil)
library(tidyverse)
library(RColorBrewer)

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
k_stab = k_score = k_fixed = NULL
prop_stab = prop_score = prop_fixed = NULL
rand_stab = rand_score = rand_fixed = NULL
adjrand_stab = adjrand_score = adjrand_fixed = NULL
for (i in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Results/",filepaths[i],"/Cluster_memberships.rds"))
  # Number of clusters
  k_stab = c(k_stab, length(unique(covars$stab_cluster)))
  k_score = c(k_score, length(unique(covars$score_cluster)))
  k_fixed = c(k_fixed, length(unique(covars$fixed_cluster)))
  
  # Adjusted Rand index
  adjrand_stab = c(adjrand_stab, adj.rand.index(as.numeric(covars$Family.ID),covars$stab_cluster))
  adjrand_score = c(adjrand_score, adj.rand.index(as.numeric(covars$Family.ID),covars$score_cluster))
  adjrand_fixed = c(adjrand_fixed, adj.rand.index(as.numeric(covars$Family.ID),covars$fixed_cluster))
  
  # Families reconstructed
  tmp = covars %>%
    group_by(Family.ID) %>%
    summarise(count = length(unique(stab_cluster))) %>%
    filter(count==1) %>%
    nrow(.)/length(unique(covars$Family.ID))
  prop_stab = c(prop_stab, tmp)
  
  tmp = covars %>%
    group_by(Family.ID) %>%
    summarise(count = length(unique(score_cluster))) %>%
    filter(count==1) %>%
    nrow(.)/length(unique(covars$Family.ID))
  prop_score = c(prop_score, tmp)
  
  tmp = covars %>%
    group_by(Family.ID) %>%
    summarise(count = length(unique(fixed_cluster))) %>%
    filter(count==1) %>%
    nrow(.)/length(unique(covars$Family.ID))
  prop_fixed = c(prop_fixed, tmp)
}

adjrand_stab = round(adjrand_stab,2)
prop_stab = round(prop_stab*100,1)

adjrand_score = round(adjrand_score,2)
prop_score = round(prop_score*100,1)

adjrand_fixed = round(adjrand_fixed,2)
prop_fixed = round(prop_fixed*100,1)

### PLS for visualisation ----
library(sgPLS)
library(plotrix)
library(ellipse)
for (i in 1:length(batches)){
  # Load data
  covars = readRDS(paste0("../Results/",filepaths[i],"/Cluster_memberships.rds"))
  assign(paste0("covars_",suffix[i]), covars)
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  assign(paste0("expo_",suffix[i]), expo)
  
  families=unique(covars$stab_cluster)[table(covars$stab_cluster)!=1]
  mycolours=brewer.pal(n=12,name='Paired')
  mycolours=colorRampPalette(mycolours)(length(families))
  names(mycolours)=families
  tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars$stab_cluster), families)))
  names(tmp) = setdiff(unique(covars$stab_cluster), families)
  mycolours = c(mycolours, tmp)
  mycolours = mycolours[as.character(unique(covars$stab_cluster))]
  
  myplsda = plsda(expo, as.factor(covars$stab_cluster), ncomp = 3)
  CreateScorePlot.plsda(myplsda=myplsda, type1=covars$stab_cluster, type2=covars$Family.ID,
                        legend_text = c(paste0("k=",k_stab[i]),
                                        paste0("Adj. Rand index: ",adjrand_stab[i]),
                                        paste0(prop_stab[i],"% recovered")),
                        mycolours=mycolours, filename=paste0("../Figures/",filepaths[i],"/PLSDA_score_plot_stab_cluster.pdf"))
                   
  
  families=unique(covars$score_cluster)[table(covars$score_cluster)!=1]
  mycolours=brewer.pal(n=12,name='Paired')
  mycolours=colorRampPalette(mycolours)(length(families))
  names(mycolours)=families
  tmp = rep(alpha(alpha("grey80",0.3),0.3),length(setdiff(unique(covars$score_cluster), families)))
  names(tmp) = setdiff(unique(covars$score_cluster), families)
  mycolours = c(mycolours, tmp)
  mycolours = mycolours[as.character(unique(covars$score_cluster))]
  
  myplsda = plsda(expo, as.factor(covars$score_cluster), ncomp = 3)
  CreateScorePlot.plsda(myplsda=myplsda, type1=covars$score_cluster, type2=covars$Family.ID,
                        legend_text = c(paste0("k=",k_score[i]),
                                        paste0("Adj. Rand index: ",adjrand_score[i]),
                                        paste0(prop_score[i],"% recovered")),
                        mycolours=mycolours, filename=paste0("../Figures/",filepaths[i],"/PLSDA_score_plot_score_cluster.pdf"))
  
  families=unique(covars$fixed_cluster)[table(covars$fixed_cluster)!=1]
  mycolours=brewer.pal(n=12,name='Paired')
  mycolours=colorRampPalette(mycolours)(length(families))
  names(mycolours)=families
  tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars$fixed_cluster), families)))
  names(tmp) = setdiff(unique(covars$fixed_cluster), families)
  mycolours = c(mycolours, tmp)
  mycolours = mycolours[as.character(unique(covars$fixed_cluster))]
  
  myplsda = plsda(expo, as.factor(covars$fixed_cluster), ncomp = 3)
  CreateScorePlot.plsda(myplsda=myplsda, type1=covars$fixed_cluster, type2=covars$Family.ID,
                        legend_text = c(paste0("k=",k_fixed[i]),
                                        paste0("Adj. Rand index: ",adjrand_fixed[i]),
                                        paste0(prop_fixed[i],"% recovered")),
                        mycolours=mycolours, filename=paste0("../Figures/",filepaths[i],"/PLSDA_score_plot_fixed_cluster.pdf"))
  
}

### Pooled 2: Batch vs Stability----
families=unique(covars_pooled2$stab_cluster)[table(covars_pooled2$stab_cluster)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled2$stab_cluster), families)))
names(tmp) = setdiff(unique(covars_pooled2$stab_cluster), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled2$stab_cluster))]

myplsda = plsda(expo_pooled2, as.factor(covars_pooled2$stab_cluster), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled2$stab_cluster, type2=covars_pooled2$Batch,
                      legend_text = c(paste0("k=",k_stab[5]),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled2$Batch),covars_pooled2$stab_cluster),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled2$Batch),covars_pooled2$stab_cluster),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[5],"/PLSDA_score_plot_stab_cluster_batch.pdf"))


### Pooled 3: Batch vs Stability----
families=unique(covars_pooled3$stab_cluster)[table(covars_pooled3$stab_cluster)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$stab_cluster), families)))
names(tmp) = setdiff(unique(covars_pooled3$stab_cluster), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$stab_cluster))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$stab_cluster), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$stab_cluster, type2=covars_pooled3$Batch,
                      legend_text = c(paste0("k=",k_stab[4]),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Batch),covars_pooled3$stab_cluster),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Batch),covars_pooled3$stab_cluster),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_stab_cluster_batch.pdf"))

### Pooled 3: Country vs Stability----
families=unique(covars_pooled3$stab_cluster)[table(covars_pooled3$stab_cluster)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$stab_cluster), families)))
names(tmp) = setdiff(unique(covars_pooled3$stab_cluster), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$stab_cluster))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$stab_cluster), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$stab_cluster, type2=covars_pooled3$Country,
                      legend_text = c(paste0("k=",k_stab[4]),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Country),covars_pooled3$stab_cluster),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Country),covars_pooled3$stab_cluster),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_stab_cluster_country.pdf"))

### Pooled 3: Region vs Stability----
families=unique(covars_pooled3$stab_cluster)[table(covars_pooled3$stab_cluster)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$stab_cluster), families)))
names(tmp) = setdiff(unique(covars_pooled3$stab_cluster), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$stab_cluster))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$stab_cluster), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$stab_cluster, type2=covars_pooled3$Region,
                      legend_text = c(paste0("k=",k_stab[4]),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Region),covars_pooled3$stab_cluster),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Region),covars_pooled3$stab_cluster),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_stab_cluster_region.pdf"))

### Pooled 3: Department vs Stability----
families=unique(covars_pooled3$stab_cluster)[table(covars_pooled3$stab_cluster)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$stab_cluster), families)))
names(tmp) = setdiff(unique(covars_pooled3$stab_cluster), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$stab_cluster))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$stab_cluster), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$stab_cluster, type2=covars_pooled3$Department,
                      legend_text = c(paste0("k=",k_stab[4]),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Department),covars_pooled3$stab_cluster),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Department),covars_pooled3$stab_cluster),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_stab_cluster_depart.pdf"))

### Pooled 3: Batch vs Fixed----
families=unique(covars_pooled3$fixed_cluster_batch)[table(covars_pooled3$fixed_cluster_batch)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$fixed_cluster_batch), families)))
names(tmp) = setdiff(unique(covars_pooled3$fixed_cluster_batch), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$fixed_cluster_batch))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$fixed_cluster_batch), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$fixed_cluster_batch, type2=covars_pooled3$Batch,
                      legend_text = c(paste0("k=",length(unique(covars_pooled3$Batch))),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Batch),covars_pooled3$fixed_cluster_batch),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Batch),covars_pooled3$fixed_cluster_batch),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_fixed_cluster_batch.pdf"))

### Pooled 3: Country vs Fixed----
families=unique(covars_pooled3$fixed_cluster_country)[table(covars_pooled3$fixed_cluster_country)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$fixed_cluster_country), families)))
names(tmp) = setdiff(unique(covars_pooled3$fixed_cluster_country), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$fixed_cluster_country))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$fixed_cluster_country), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$fixed_cluster_country, type2=covars_pooled3$Country,
                      legend_text = c(paste0("k=",length(unique(covars_pooled3$Country))),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Country),covars_pooled3$fixed_cluster_country),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Country),covars_pooled3$fixed_cluster_country),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_fixed_cluster_country.pdf"))

### Pooled 3: Region vs Fixed----
families=unique(covars_pooled3$fixed_cluster_region)[table(covars_pooled3$fixed_cluster_region)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$fixed_cluster_region), families)))
names(tmp) = setdiff(unique(covars_pooled3$fixed_cluster_region), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$fixed_cluster_region))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$fixed_cluster_region), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$fixed_cluster_region, type2=covars_pooled3$Region,
                      legend_text = c(paste0("k=",length(unique(covars_pooled3$Region))),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Region),covars_pooled3$fixed_cluster_region),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Region),covars_pooled3$fixed_cluster_region),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_fixed_cluster_region.pdf"))

### Pooled 3: Department vs Fixed----
families=unique(covars_pooled3$fixed_cluster_depart)[table(covars_pooled3$fixed_cluster_depart)!=1]
mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families
tmp = rep(alpha("grey80",0.3),length(setdiff(unique(covars_pooled3$fixed_cluster_depart), families)))
names(tmp) = setdiff(unique(covars_pooled3$fixed_cluster_depart), families)
mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(covars_pooled3$fixed_cluster_depart))]

myplsda = plsda(expo_pooled3, as.factor(covars_pooled3$fixed_cluster_depart), ncomp = 3)
CreateScorePlot.plsda(myplsda=myplsda, type1=covars_pooled3$fixed_cluster_depart, type2=covars_pooled3$Department,
                      legend_text = c(paste0("k=",length(unique(covars_pooled3$Department))),
                                      paste0("Adj. Rand index: ",round(adj.rand.index(as.numeric(covars_pooled3$Department),covars_pooled3$fixed_cluster_depart),2)),
                                      paste0("Rand index: ",round(rand.index(as.numeric(covars_pooled3$Department),covars_pooled3$fixed_cluster_depart),2))),
                      mycolours=mycolours, filename=paste0("../Figures/",filepaths[4],"/PLSDA_score_plot_fixed_cluster_depart.pdf"))


