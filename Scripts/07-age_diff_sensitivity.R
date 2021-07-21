## Age difference sensitivity analysis
## Rin Wada 15 July

# Load packages
library(ape)
library(igraph)
library(tidyverse)

# Initialise
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
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  families=unique(covars$Family.ID)
  
  expo = scale(expo)
  d=dist(expo)
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))
  h$labels=paste0(covars$Family.ID, "-",h$labels)
  myphylo=as.phylo(h)
  
  g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[levels(covars$Family.ID)], verbose = FALSE)
  
  nobs=nrow(expo)
  sp=shortest.paths(g)[1:nobs,1:nobs]
  rownames(sp)=colnames(sp)=myphylo$tip.label
  
  fsp=NULL
  for (f in 1:length(families)){
    tmpmat=sp[covars$Family.ID==families[f],covars$Family.ID==families[f]]
    fsp=c(fsp,mean(tmpmat[upper.tri(tmpmat)]))
  }
  names(fsp)=families
  
  # Age group
  if (i==1){
    covars$Age_group = ifelse(covars$Age<=5, "<=5",
                              ifelse(covars$Age>5 & covars$Age<=7, "6-7",
                                     ifelse(covars$Age>7 & covars$Age<=9, "8-9",
                                            ifelse(covars$Age>9 & covars$Age<=11, "10-11",
                                                   ifelse(covars$Age>11 & covars$Age<=15, "12-15", ">15")))))
    }
  if (i %in% c(2,3)){
    covars$Age_group = ifelse(covars$Age<=5, "<=5",
                              ifelse(covars$Age>5 & covars$Age<=10, "6-10",
                                     ifelse(covars$Age>10 & covars$Age<=15, "11-15", ">15")))
  }
  if (i %in% c(4,5)){
    covars$Age_group[which(covars$Country == "France")] = ifelse(covars$Age[which(covars$Country == "France")]<=5, "<=5",
                                                                 ifelse(covars$Age[which(covars$Country == "France")]>5 & covars$Age[which(covars$Country == "France")]<=10, "6-10",
                                                                        ifelse(covars$Age[which(covars$Country == "France")]>10 & covars$Age[which(covars$Country == "France")]<=15, "11-15", ">15")))
    covars$Age_group[which(covars$Country == "Luxembourg")] = ifelse(covars$Age[which(covars$Country == "Luxembourg")]<=5, "<=5",
                                                                     ifelse(covars$Age[which(covars$Country == "Luxembourg")]>5 & covars$Age[which(covars$Country == "Luxembourg")]<=7, "6-7",
                                                                            ifelse(covars$Age[which(covars$Country == "Luxembourg")]>7 & covars$Age[which(covars$Country == "Luxembourg")]<=9, "8-9",
                                                                                   ifelse(covars$Age[which(covars$Country == "Luxembourg")]>9 & covars$Age[which(covars$Country == "Luxembourg")]<=11, "10-11",
                                                                                          ifelse(covars$Age[which(covars$Country == "Luxembourg")]>11 & covars$Age[which(covars$Country == "Luxembourg")]<=15, "12-15", ">15")))))
  }
  age_group_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    tmp=covars$Age_group[covars$Family.ID==families[f]]
    if (sum(is.na(tmp))>0){
      age_group_diff[f] = NA
    } else if (length(unique(tmp))==1){
      age_group_diff[f] = 1
    } else {
      age_group_diff[f] = 2
    }
  }
  names(age_group_diff)=families
  
  gender_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    tmp=covars$Gender[covars$Family.ID==families[f]]
    if (sum(is.na(tmp))>0) {
      gender_diff[f] = NA
    } else if (length(unique(tmp))==1){
      if (unique(tmp)=="Male"){
        gender_diff[f]=1
      }
      if (unique(tmp)=="Female"){
        gender_diff[f]=2
      }
    } else {
      gender_diff[f]=3
    }
  }
  names(gender_diff)=families
  
  model1=lm(fsp~as.factor(age_group_diff))
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_group_diff.pdf"))
    par(mar=c(5,5,1,1))
    plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
         col=c("springgreen", "red")[age_group_diff],
         xlab="Family", ylab="Shortest path length between siblings")
    text(fsp, labels=ifelse(!is.na(age_group_diff), names(age_group_diff), ""), pos=3, col=darken(family.colours[names(age_group_diff)], amount=0.5))
    legend("topleft", pch=19, col=c("springgreen", "red"), ncol=1,
           legend=c("Same age group", "Different age group"))
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model1)$coefficients[2,4], format="e", digits=2)))
    dev.off()
  }
  
  if (i %in% c(2,4,5)){
    model1=lm(fsp[which(gender_diff=="1")]~as.factor(age_group_diff)[which(gender_diff=="1")])
    model2=lm(fsp[which(gender_diff=="2")]~as.factor(age_group_diff)[which(gender_diff=="2")])
    model3=lm(fsp[which(gender_diff=="3")]~as.factor(age_group_diff)[which(gender_diff=="3")])
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_group_diff_by_gender_diff.pdf"), width = 15, height = 5)
      par(mar=c(5,5,2,2), mfrow = c(1,3))
      for (k in 1:3){
        plot(fsp[which(gender_diff==k)], pch=19, cex=1, cex.lab=1.5, las=1,
             col=c("springgreen", "red")[age_group_diff[which(gender_diff==k)]],
             xlab="Family", ylab="Shortest path length between siblings",
             main = c("All male", "All female", "Different gender")[k])
        text(fsp[which(gender_diff==k)], labels=ifelse(!is.na(age_group_diff[which(gender_diff==k)]), names(age_group_diff[which(gender_diff==k)]), ""), pos=3, col=darken(family.colours[names(age_group_diff[which(gender_diff==k)])], amount=0.5))
        legend("topleft", pch=19, col=c("springgreen", "red"), ncol=1,
               legend=c("Same age group", "Different age group"))
        legend("topright", bty="n", cex=1.5,
               legend=paste0("p=",formatC(summary(eval(parse(text = paste0("model",k))))$coefficients[2,4], format="e", digits=2))) 
      }
      dev.off()
    }
    model0.1=lm(fsp[which(age_group_diff=="1")]~1)
    model0.2=lm(fsp[which(age_group_diff=="2")]~1)
    model1=lm(fsp[which(age_group_diff=="1")]~as.factor(gender_diff)[which(age_group_diff=="1")])
    model2=lm(fsp[which(age_group_diff=="2")]~as.factor(gender_diff)[which(age_group_diff=="2")])
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_gender_diff_by_age_group_diff.pdf"), width = 10, height = 5)
      par(mar=c(5,5,2,2), mfrow = c(1,2))
      for (k in 1:2){
        plot(fsp[which(age_group_diff==k)], pch=19, cex=1, cex.lab=1.5, las=1,
             col=c("skyblue", "pink", "tan")[gender_diff[which(age_group_diff==k)]],
             xlab="Family", ylab="Shortest path length between siblings",
             main = c("Same age group", "Different age group")[k])
        text(fsp[which(age_group_diff==k)], labels=ifelse(!is.na(gender_diff[which(age_group_diff==k)]), names(gender_diff[which(age_group_diff==k)]), ""), pos=3, col=darken(family.colours[names(age_group_diff[which(gender_diff==k)])], amount=0.5))
        legend("topleft", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
               legend=c("All Male", "All Female", "Different genders"))
        legend("topright", bty="n", cex=1.5,
               legend=paste0("p=",
                             formatC(anova(eval(parse(text = paste0("model",k))),
                                           eval(parse(text = paste0("model0.",k))),
                                           test = "Chisq")$`Pr(>Chi)`[2], format="e", digits=2))) 
      }
      dev.off()
    }
    }
  }
