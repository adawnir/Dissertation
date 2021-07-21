## Multivariate clustering analysis
## Rin Wada 12 July

# Load packages
library(ape)
library(igraph)
library(colorspace)

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
  
  ## Interation between age mean, age difference and gender difference
  age_mu=rep(NA, length(families))
  for (f in 1:length(families)){
    age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
  }
  names(age_mu)=families
  
  age_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    age_diff[f]=abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2])
  }
  names(age_diff)=families
  
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
  
  model0 = lm(fsp ~ age_mu + as.factor(gender_diff))
  model1 = lm(fsp ~ age_mu * as.factor(gender_diff))
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_mean_gender_diff.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1,
         col=c("skyblue", "pink", "tan")[gender_diff],
         xlab="Within-family age mean", ylab="Shortest path length between siblings")
    text(age_mu, fsp, labels=names(age_mu), pos=3, col=darken(family.colours[names(age_mu)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("topleft", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
           legend=c("All Male", "All Female", "Different genders"))
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(anova(model0,model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
    dev.off()
  }
  
  model1=lm(fsp[which(gender_diff=="1")]~age_mu[which(gender_diff=="1")])
  model2=lm(fsp[which(gender_diff=="2")]~age_mu[which(gender_diff=="2")])
  model3=lm(fsp[which(gender_diff=="3")]~age_mu[which(gender_diff=="3")])

  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_mean_by_gender_diff.pdf"), width = 15, height = 5)
    par(mar=c(5,5,2,2), mfrow = c(1,3))
    for (k in 1:3){
      plot(age_mu[which(gender_diff==k)], fsp[which(gender_diff==k)], pch=19, cex=1, cex.lab=1.5, las=1,
           col=c("skyblue", "pink", "tan")[k],
           xlab="Within-family age mean", ylab="Shortest path length between siblings",
           main = c("All male", "All female", "Different gender")[k])
      text(age_mu[which(gender_diff==k)], fsp[which(gender_diff==k)], labels=ifelse(!is.na(age_mu[which(gender_diff==k)]), names(age_mu[which(gender_diff==k)]), ""), pos=3, col=darken(family.colours[names(age_mu[which(gender_diff==k)])], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(eval(parse(text = paste0("model",k))))$coefficients[2,4], format="e", digits=2))) 
    }
    dev.off()
    }
  
  model0=lm(fsp~age_diff + as.factor(gender_diff))
  model1=lm(fsp~age_diff * as.factor(gender_diff))
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_diff_gender_diff.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_diff, fsp, pch=19, cex=1, cex.lab=1.5, las=1,
         col=c("skyblue", "pink", "tan")[gender_diff],
         xlab="Within-family age absolute difference", ylab="Shortest path length between siblings")
    text(age_diff, fsp, labels=names(age_diff), pos=3, col=darken(family.colours[names(age_diff)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("topleft", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
           legend=c("All Male", "All Female", "Different genders"))
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(anova(model0,model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
    dev.off()
  }
  
  model1=lm(fsp[which(gender_diff=="1")]~age_diff[which(gender_diff=="1")])
  model2=lm(fsp[which(gender_diff=="2")]~age_diff[which(gender_diff=="2")])
  model3=lm(fsp[which(gender_diff=="3")]~age_diff[which(gender_diff=="3")])
  
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_diff_by_gender_diff.pdf"), width = 15, height = 5)
    par(mar=c(5,5,2,2), mfrow = c(1,3))
    for (k in 1:3){
      plot(age_diff[which(gender_diff==k)], fsp[which(gender_diff==k)], pch=19, cex=1, cex.lab=1.5, las=1,
           col=c("skyblue", "pink", "tan")[k],
           xlab="Within-family age absolute difference", ylab="Shortest path length between siblings",
           main = c("All male", "All female", "Different gender")[k])
      text(age_diff[which(gender_diff==k)], fsp[which(gender_diff==k)], labels=ifelse(!is.na(age_diff[which(gender_diff==k)]), names(age_diff[which(gender_diff==k)]), ""), pos=3, col=darken(family.colours[names(age_diff[which(gender_diff==k)])], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(eval(parse(text = paste0("model",k))))$coefficients[2,4], format="e", digits=2))) 
    }
    dev.off()
  }
}

