## Univariate clustering analysis
## Rin Wada 12 July

# Load packages
library(ape)
library(igraph)
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
  
  mycolours=brewer.pal(n=12,name='Paired')
  mycolours=colorRampPalette(mycolours)(length(families))
  names(mycolours)=families
  
  expo = scale(expo)
  d=dist(expo)
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))
  h$labels=paste0(covars$Family.ID, "-",h$labels)
  myphylo=as.phylo(h)
  
  {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_expo_all.pdf"),width=14)
    par(mar=c(0,0,0,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=mycolours[covars$Family.ID])
    dev.off()
  }
  
  {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_graph_expo_all.pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo)
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
  
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all.pdf"),width=14)
    par(mar=c(5,5,1,1))
    plot(fsp, pch=19, col=mycolours[families], xaxt="n", las=1,
         panel.first=abline(v=1:length(families), lty=3, col="grey"),
         xlab="Family ID", cex=2,
         ylab="Shortest path length between siblings", cex.lab=1.5)
    axis(side=1, at=1:length(families), labels=families, las = 2)
    dev.off()}
  
  ## Relationship with age
  age_mu=rep(NA, length(families))
  for (f in 1:length(families)){
    age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
  }
  names(age_mu)=families
  
  model=lm(fsp~age_mu)
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_mean.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(age_mu)],
         xlab="Within-family age mean", ylab="Shortest path length between siblings")
    text(age_mu, fsp, labels=names(age_mu), pos=3, col=darken(mycolours[names(age_mu)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
    dev.off()
  }
  
  age_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    age_diff[f]=abs(covars$Age[covars$Family.ID==families[f]][1]-covars$Age[covars$Family.ID==families[f]][2])
  }
  names(age_diff)=families
  
  model=lm(fsp~age_diff)
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_diff.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_diff, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(age_diff)],
         xlab="Within-family age absolute difference", ylab="Shortest path length between siblings")
    text(age_diff, fsp, labels=names(age_diff), pos=3, col=darken(mycolours[names(age_diff)], amount=0.5))
    abline(h=axTicks(2), col="grey", lty=3)
    abline(v=axTicks(1), col="grey", lty=3)
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
    dev.off()
  }
  
  if (i==1){
    weight_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_mu[f]=mean(covars$Weight[covars$Family.ID==families[f]])
    }
    names(weight_mu)=families

    model=lm(fsp~weight_mu)
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_weight_mean.pdf"))
      par(mar=c(5,5,1,1))
      plot(weight_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(weight_mu)],
           xlab="Within-family weight mean", ylab="Shortest path length between siblings")
      text(weight_mu, fsp, labels=names(weight_mu), pos=3, col=darken(mycolours[names(weight_mu)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topleft", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    weight_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      weight_diff[f]=abs(covars$Weight[covars$Family.ID==families[f]][1]-covars$Weight[covars$Family.ID==families[f]][2])
    }
    names(weight_diff)=families

    model=lm(fsp~weight_diff)
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_weight_diff.pdf"))
      par(mar=c(5,5,1,1))
      plot(weight_diff, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(weight_diff)],
           xlab="Within-family weight absolute difference", ylab="Shortest path length between siblings")
      text(weight_diff, fsp, labels=names(weight_diff), pos=3, col=darken(mycolours[names(weight_diff)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
      }

    length_mu=rep(NA, length(families))
    for (f in 1:length(families)){
      length_mu[f]=mean(covars$Length[covars$Family.ID==families[f]])
    }
    names(length_mu)=families

    model=lm(fsp~length_mu)
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_length_mean.pdf"))
      par(mar=c(5,5,1,1))
      plot(length_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(length_mu)],
           xlab="Within-family length mean", ylab="Shortest path length between siblings")
      text(length_mu, fsp, labels=names(length_mu), pos=3, col=darken(mycolours[names(length_mu)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
      }

    length_diff=rep(NA, length(families))
    for (f in 1:length(families)){
      length_diff[f]=abs(covars$Length[covars$Family.ID==families[f]][1]-covars$Length[covars$Family.ID==families[f]][2])
    }
    names(length_diff)=families
    
    model=lm(fsp~length_diff)
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_length_diff.pdf"))
      par(mar=c(5,5,1,1))
      plot(length_diff, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(length_diff)],
           xlab="Within-family length absolute difference", ylab="Shortest path length between siblings")
      text(length_diff, fsp, labels=names(length_diff), pos=3, col=darken(mycolours[names(length_diff)], amount=0.5))
      abline(h=axTicks(2), col="grey", lty=3)
      abline(v=axTicks(1), col="grey", lty=3)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
      }
  }
  
  gender_diff=rep(NA, length(families))
  for (f in 1:length(families)){
    tmp=covars$Gender[covars$Family.ID==families[f]]
    if (length(unique(tmp))==1){
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
  
  model1 = lm(fsp ~ gender_diff)
  model0 = lm(fsp ~ 1)
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_by_gender.pdf"))
    par(mar=c(5,5,1,1))
    plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
         col=c("skyblue", "pink", "tan")[gender_diff],
         xlab="Family ID", ylab="Shortest path length between siblings")
    text(fsp, labels=names(gender_diff), pos=3, col=darken(mycolours[names(gender_diff)], amount=0.5))
    legend("topleft", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
           legend=c("All Male", "All Female", "Different genders"))
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
    dev.off()
  }
  if (i %in% 4:5){
    Batch=rep(NA, length(families))
    for (f in 1:length(families)){
      tmp=covars$Batch[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Batch[f] = unique(tmp)
      }
    }
    names(Batch)=families
    
    model1 = lm(fsp ~ Batch)
    model0 = lm(fsp ~ 1)
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_by_Batch.pdf"))
      par(mar=c(5,5,1,1))
      plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=batch.colours[Batch],
           xlab="Family ID", ylab="Shortest path length between siblings")
      text(fsp, labels=names(Batch), pos=3, col=darken(mycolours[names(Batch)], amount=0.5))
      legend("topleft", pch=19, col=batch.colours[unique(Batch)],
             legend=levels(covars$Batch))
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
      dev.off()
    }
  }
  if (i == 4){
    Region=rep(NA, length(families))
    Department=rep(NA, length(families))
    for (f in 1:length(families)){
      tmp=covars$Region[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Region[f] = unique(tmp)
      }
      tmp=covars$Department[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Department[f] = unique(tmp)
      }
    }
    names(Region)=families
    names(Department)=families
    
    region.colours=brewer.pal(n=12,name='Paired')
    region.colours=colorRampPalette(region.colours)(length(levels(covars$Region)))
    names(region.colours)=levels(covars$Region)
    
    depart.colours=brewer.pal(n=12,name='Paired')
    depart.colours=colorRampPalette(depart.colours)(length(levels(covars$Department)))
    names(depart.colours)=levels(covars$Department)
    
    model1 = lm(fsp ~ Region)
    model0 = lm(fsp ~ 1)
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_by_Region.pdf"))
      par(mar=c(5,5,1,1))
      plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=region.colours[Region],
           xlab="Family ID", ylab="Shortest path length between siblings")
      text(fsp, labels=names(Region), pos=3, col=darken(mycolours[names(Region)], amount=0.5))
      legend("topleft", pch=19, col=region.colours[unique(Region)],
             legend=levels(covars$Region), cex = 0.6, ncol = 2)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
      dev.off()
    }
    
    model1 = lm(fsp ~ Department)
    model0 = lm(fsp ~ 1)
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_by_Department.pdf"))
      par(mar=c(5,5,1,1))
      plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=depart.colours[Department],
           xlab="Family ID", ylab="Shortest path length between siblings")
      text(fsp, labels=names(Department), pos=3, col=darken(mycolours[names(Department)], amount=0.5))
      legend("topleft", pch=19, col=depart.colours[unique(Department)],
             legend=levels(covars$Department), cex = 0.6, ncol = 2)
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
      dev.off()
    }
  }
}

