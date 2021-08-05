## Multivariate clustering analysis
## Rin Wada 12 July

# Load packages
library(ape)
library(igraph)
library(focus)
library(colorspace)
library(RColorBrewer)

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
  
  families=levels(covars$Family.ID)
  
  expo = scale(expo)
  d=dist(expo)
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))
  h$labels=paste0(covars$Family.ID, "-",h$labels)
  myphylo=as.phylo(h)
  
  {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_expo_all.pdf"),width=14)
    par(mar=c(0,0,0,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=family.colours[as.character(covars$Family.ID)])
    dev.off()
  }
  
  {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_graph_expo_all.pdf"))
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
  
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all.pdf"),width=14)
    par(mar=c(5,5,1,1))
    plot(fsp, pch=19, col=family.colours[families], xaxt="n", las=1,
         panel.first=abline(v=1:length(families), lty=3, col="grey"),
         xlab="Family ID", cex=2,
         ylab="Shortest path length between siblings", cex.lab=1.5)
    axis(side=1, at=1:length(families), labels=families, las = 2)
    dev.off()
    }
  
  ## Relationship with age
  age_mu=rep(NA, length(families))
  for (f in 1:length(families)){
    age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
  }
  names(age_mu)=families
  
  model=lm(fsp~age_mu)
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_vs_age_mean.pdf"))
    par(mar=c(5,5,1,1))
    plot(age_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(age_mu)],
         xlab="Within-family age mean", ylab="Shortest path length between siblings")
    text(age_mu, fsp, labels=names(age_mu), pos=3, col=darken(family.colours[names(age_mu)], amount=0.5))
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
    plot(age_diff, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(age_diff)],
         xlab="Within-family age absolute difference", ylab="Shortest path length between siblings")
    text(age_diff, fsp, labels=names(age_diff), pos=3, col=darken(family.colours[names(age_diff)], amount=0.5))
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
      plot(weight_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(weight_mu)],
           xlab="Within-family weight mean", ylab="Shortest path length between siblings")
      text(weight_mu, fsp, labels=names(weight_mu), pos=3, col=darken(family.colours[names(weight_mu)], amount=0.5))
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
      plot(weight_diff, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(weight_diff)],
           xlab="Within-family weight absolute difference", ylab="Shortest path length between siblings")
      text(weight_diff, fsp, labels=names(weight_diff), pos=3, col=darken(family.colours[names(weight_diff)], amount=0.5))
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
      plot(length_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(length_mu)],
           xlab="Within-family length mean", ylab="Shortest path length between siblings")
      text(length_mu, fsp, labels=names(length_mu), pos=3, col=darken(family.colours[names(length_mu)], amount=0.5))
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
      plot(length_diff, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=family.colours[names(length_diff)],
           xlab="Within-family length absolute difference", ylab="Shortest path length between siblings")
      text(length_diff, fsp, labels=names(length_diff), pos=3, col=darken(family.colours[names(length_diff)], amount=0.5))
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
    if (sum(is.na(tmp))>0) {
      gender_diff[f] = NA
    } else if (length(unique(tmp))==1){
      if (unique(tmp)=="Male"){
        gender_diff[f]="All male"
      }
      if (unique(tmp)=="Female"){
        gender_diff[f]="All female"
      }
    } else {
      gender_diff[f]="Different genders"
    }
  }
  names(gender_diff)=families
  
  gender_col = c("skyblue", "pink", "tan")
  names(gender_col) = c("All male", "All female", "Different genders")
  
  model1 = lm(fsp[!is.na(gender_diff)] ~ gender_diff[!is.na(gender_diff)])
  print(summary(model1))
  model0 = lm(fsp[!is.na(gender_diff)] ~ 1)
  {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_by_gender.pdf"))
    par(mar=c(5,5,1,1))
    plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
         col=gender_col[gender_diff],
         xlab="Family", ylab="Shortest path length between siblings")
    text(fsp, labels=ifelse(!is.na(gender_diff),names(gender_diff),""), pos=3,
         col=darken(family.colours[names(gender_diff)], amount=0.5))
    legend("topleft", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
           legend=c("All Male", "All Female", "Different genders"))
    legend("topright", bty="n", cex=1.5,
           legend=paste0("p=",formatC(anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2], format="e", digits=2)))
    dev.off()
  }
  if (i %in% c(2,4)){
    Region=rep(NA, length(families))
    Department=rep(NA, length(families))
    for (f in 1:length(families)){
      tmp=as.character(covars$Region)[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Region[f] = unique(tmp)
      }
      tmp=as.character(covars$Department)[covars$Family.ID==families[f]]
      if (length(unique(tmp))==1){
        Department[f] = unique(tmp)
      }
    }
    names(Region)=families
    names(Department)=families
    
    Department = ifelse(Department=="Paris",1,0)
    Region = ifelse(Region=="Île-de-France",1,0)
    
    model = lm(fsp ~ as.factor(Region))
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_by_Region.pdf"))
      par(mar=c(5,5,1,1))
      plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=ifelse(Region==1,region.colours["Île-de-France"],"grey"),
           xlab="Family", ylab="Shortest path length between siblings")
      text(fsp, labels=names(Region), pos=3, col=darken(family.colours[names(Region)], amount=0.5))
      legend("top", pch=19, col=c(region.colours["Île-de-France"],"grey"),
             legend=c("Île-de-France","Other region"))
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
    
    model = lm(fsp ~ as.factor(Department))
    {pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_cont_all_by_Department.pdf"))
      par(mar=c(5,5,1,1))
      plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
           col=ifelse(Department==1,depart.colours["Paris"],"grey"),
           xlab="Family", ylab="Shortest path length between siblings")
      text(fsp, labels=names(Department), pos=3, col=darken(family.colours[names(Department)], amount=0.5))
      legend("top", pch=19, col=c(depart.colours["Paris"],"grey"),
             legend=c("Paris","Other department"))
      legend("topright", bty="n", cex=1.5,
             legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
      dev.off()
    }
  }
  
  ## Delta exposure
  delta_mat=matrix(NA, nrow=length(families), ncol=ncol(expo))
  for (f in 1:length(families)){
    for (k in 1:ncol(expo)){
      tmp = expo[covars$Family.ID==families[f],k]
      delta_mat[f,k] = abs(tmp[1] - tmp [2])
    }
  }
  rownames(delta_mat)=families
  colnames(delta_mat)=colnames(expo)
  
  ### Multivariate analysis using stability selection LASSO regression ----
  if (i == 1){
    X = cbind(age_mu, age_diff, gender_diff,
              length_mu, length_diff, weight_mu, weight_diff,
              as.data.frame(delta_mat))
  }
  if (i %in% c(2,4)){
    X = cbind(age_mu, age_diff, gender_diff,
              Region,
              as.data.frame(delta_mat))
  }
  if (i %in% c(3,5)){
    X = cbind(age_mu, age_diff, gender_diff,
              as.data.frame(delta_mat))
  }
  
  Y = fsp[complete.cases(X)]
  X = model.matrix(~., X)[,-1]
  
  stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)
  
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_output.pdf"))
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_beta_consistency.pdf"))
  par(mar=c(5,5,1,1))
  plot(a/(a+b), selprop, las=1, pch=19, col="navy", cex.lab=1.5,
       xlab="Proportion of positive beta among non-zero betas", 
       ylab="Selection Proportion") # Ideally no point around 0.5 with high selection proportion
  dev.off()
  
  # Extract beta
  beta = stab$Beta[ArgmaxId(stab)[1],,]
  beta_mu = rowMeans(beta[-nrow(beta),] %*% diag(beta[nrow(beta),]))
  names(beta_mu) = names(selprop)
  
  selprop_ranked <- sort(selprop, decreasing = TRUE)
  beta_mu_ranked <- beta_mu[order(-selprop)]
  print(all(names(selprop_ranked)==names(beta_mu_ranked)))
  
  mylabels = gsub("`","",names(selprop_ranked))
  mylabels = gsub("age_mu","Mean age",mylabels)
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("gender_diff","",mylabels)
  mylabels = gsub("All male","All male (ref. All female)",mylabels)
  mylabels = gsub("Different genders","Different genders (ref. All female)",mylabels)
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("Region","Île-de-France (Y/N)",mylabels)
  
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_selprop.pdf"), width = 14)
  par(mar=c(15, 5, 1, 1))
  plot(selprop_ranked,
       type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
       col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[i], no = "grey"),
       xaxt = "n", xlab = "", ylab = "Selection proportions"
  )
  abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
  axis(side = 1, at = 1:length(selprop_ranked), labels = NA)
  for (k in 1:length(selprop_ranked)){
    if (beta_mu_ranked[k] > 0){
      col = "red"
    } else if (beta_mu_ranked[k] < 0){
      col = "blue"
    } else {
      col = "black"
    }
    if (mylabels[k] %in% colnames(delta_mat)){
      label = substitute(Delta~tmp, list(tmp = mylabels[k]))
    } else {
      label = mylabels[k]
    }
    axis(side = 1, at = k, las = 2, labels = label,
         col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
  }
  dev.off()
  
  ## Only covariates ----
  # This line throws error in sgPLS::sPLS internal function, unless scale = FALSE
  stab = VariableSelection(xdata = X[,1:(ncol(X)-ncol(delta_mat))], ydata = Y,
                           implementation = SparsePLS, Lambda = 1:(ncol(X) - 1))
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_covariates_only_output.pdf"))
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(stab)
  dev.off()

  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_covariates_only_beta_consistency.pdf"))
  par(mar=c(5,5,1,1))
  plot(a/(a+b), selprop, las=1, pch=19, col="navy", cex.lab=1.5,
       xlab="Proportion of positive beta among non-zero betas",
       ylab="Selection Proportion") # Ideally no point around 0.5 with high selection proportion
  dev.off()

  # Extract beta
  beta = stab$Beta[ArgmaxId(stab)[1],,]
  beta_mu = rowMeans(beta[-nrow(beta),] %*% diag(beta[nrow(beta),]))
  names(beta_mu) = names(selprop)

  selprop_ranked <- sort(selprop, decreasing = TRUE)
  beta_mu_ranked <- beta_mu[order(-selprop)]
  print(all(names(selprop_ranked)==names(beta_mu_ranked)))

  mylabels = gsub("age_mu","Mean age",names(selprop_ranked))
  mylabels = gsub("age_diff","Age difference",mylabels)
  mylabels = gsub("gender_diff","",mylabels)
  mylabels = gsub("All male","All male (ref. All female)",mylabels)
  mylabels = gsub("Different genders","Different genders (ref. All female)",mylabels)
  mylabels = gsub("length_mu","Mean sample length",mylabels)
  mylabels = gsub("length_diff","Sample length difference",mylabels)
  mylabels = gsub("weight_mu","Mean sample weight",mylabels)
  mylabels = gsub("weight_diff","Sample weight difference",mylabels)
  mylabels = gsub("Region","Île-de-France (Y/N)",mylabels)

  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_covariates_only_selprop.pdf"), width = 14)
  par(mar=c(15, 5, 1, 1))
  plot(selprop_ranked,
       type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
       col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[i], no = "grey"),
       xaxt = "n", xlab = "", ylab = "Selection proportions"
  )
  abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
  axis(side = 1, at = 1:length(selprop_ranked), labels = NA)
  for (k in 1:length(selprop_ranked)){
    if (beta_mu_ranked[k] > 0){
      col = "red"
    } else if (beta_mu_ranked[k] < 0){
      col = "blue"
    } else {
      col = "black"
    }
    axis(side = 1, at = k, las = 2, labels = mylabels[k],
         col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
  }
  dev.off()

  ## Only exposures ----
  stab = VariableSelection(xdata = X[,(ncol(X)-ncol(delta_mat)+1):ncol(X)], ydata = Y, implementation = SparsePLS)
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_exposures_only_output.pdf"))
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(stab)
  dev.off()
  
  # Checking consistency in sign of the beta coefficients for the variables with high selprop
  # Ideally no point around 0.5 with high selection proportion
  selprop=SelectionProportions(stab)
  a=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x>0)})
  a = a[-length(a)]
  b=apply(stab$Beta[ArgmaxId(stab)[1],,],1,FUN=function(x){sum(x<0)})
  b = b[-length(b)]
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_exposures_only_beta_consistency.pdf"))
  par(mar=c(5,5,1,1))
  plot(a/(a+b), selprop, las=1, pch=19, col="navy", cex.lab=1.5,
       xlab="Proportion of positive beta among non-zero betas", 
       ylab="Selection Proportion") # Ideally no point around 0.5 with high selection proportion
  dev.off()
  
  # Extract beta
  beta = stab$Beta[ArgmaxId(stab)[1],,]
  beta_mu = rowMeans(beta[-nrow(beta),] %*% diag(beta[nrow(beta),]))
  names(beta_mu) = names(selprop)
  
  selprop_ranked <- sort(selprop, decreasing = TRUE)
  beta_mu_ranked <- beta_mu[order(-selprop)]
  print(all(names(selprop_ranked)==names(beta_mu_ranked)))
  
  mylabels = gsub("`","",names(selprop_ranked))
  
  pdf(paste0("../Figures/",filepaths[i],"/Shortest_path_multivariate_exposures_only_selprop.pdf"), width = 14)
  par(mar=c(15, 5, 1, 1))
  plot(selprop_ranked,
       type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
       col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[i], no = "grey"),
       xaxt = "n", xlab = "", ylab = "Selection proportions"
  )
  abline(h = Argmax(stab)[2], lty = 2, col = "darkred")
  axis(side = 1, at = 1:length(selprop_ranked), labels = NA)
  for (k in 1:length(selprop_ranked)){
    if (beta_mu_ranked[k] > 0){
      col = "red"
    } else if (beta_mu_ranked[k] < 0){
      col = "blue"
    } else {
      col = "black"
    }
    axis(side = 1, at = k, las = 2, labels = substitute(Delta~tmp, list(tmp=mylabels[k])),
         col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
  }
  dev.off()
}

### Geographical ----
for (i in 4:5){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[i],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[i],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  families=unique(covars$Batch)
  
  # mycolours=brewer.pal(n=12,name='Paired')
  # mycolours=colorRampPalette(mycolours)(length(families))
  # names(mycolours)=families

  expo = scale(expo)
  d=dist(expo)
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))
  h$labels=paste0(covars$Family.ID, "-",h$labels)
  myphylo=as.phylo(h)
  
  {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_expo_all_batch.pdf"),width=14)
    par(mar=c(0,0,0,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=batch.colours[as.character(covars$Batch)])
    dev.off()
  }
  {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_graph_expo_all_batch.pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = batch.colours[as.character(covars$Batch)])
    dev.off()
  }
  
  if (i == 4){
    families=unique(covars$Region)
    mycolours = region.colours
    names(mycolours)=families
    
    {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_expo_all_region.pdf"),width=14)
      par(mar=c(0,0,0,0))
      plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
           tip.color=region.colours[as.character(covars$Region)])
      dev.off()
    }
    
    {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_graph_expo_all_region.pdf"))
      g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = region.colours[as.character(covars$Region)])
      legend("bottomright", pch=19, col=region.colours[levels(covars$Region)],
             legend=levels(covars$Region), cex = 0.4, ncol = 1)
      dev.off()
    }
    
    {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_expo_all_depart.pdf"),width=14)
      par(mar=c(0,0,0,0))
      plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
           tip.color=depart.colours[as.character(covars$Department)])
      dev.off()
    }
    
    {pdf(paste0("../Figures/",filepaths[i],"/Hierarchical_cont_graph_expo_all_depart.pdf"))
      g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = depart.colours[as.character(covars$Department)])
      legend("bottomright", pch=19, col=depart.colours[levels(covars$Department)],
             legend=levels(covars$Department), cex = 0.4, ncol = 1)
      dev.off()
    }
  }
}

