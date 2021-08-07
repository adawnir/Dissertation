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

betas = pvals = NULL
f1='fsp ~ delta_mat[,k]'
t0=Sys.time()
for (k in 1:ncol(delta_mat)){
  model1=lm(as.formula(f1))
  betas=c(betas, coefficients(model1)[2])
  pvals=c(pvals, summary(model1)$coefficients[2,4])
}
t1=Sys.time()
print(t1-t0)
names(pvals)=names(betas)=colnames(delta_mat)

{pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_vs_delta_expo.pdf"))
  par(mar=c(5,5,1,1))
  plot(betas, -log10(pvals), pch=19,
       col=ifelse(pvals < 0.05/length(betas), annot.colours[annot_sub], "grey"),
       cex.lab=1, cex = 0.7,
       ylim = c(0, max(-log10(pvals))+0.25),
       xlim = c(-max(abs(betas))-0.4, max(abs(betas))+0.4),
       ylab=expression(-log[10](p)), 
       xlab=expression(beta))
  text(betas+sign(betas)*0.1, -log10(pvals)+0.25,
       labels = ifelse(pvals < 0.05/length(betas), names(betas), ""),
       col = annot.colours[annot_sub])
  abline(h = -log10(0.05/length(betas)), lty = 2)
  dev.off()
}

### Multivariate analysis using stability selection LASSO regression ----
if (m == 1){
  X = cbind(age_mu, age_diff, gender_diff,
            length_mu, length_diff, weight_mu, weight_diff,
            as.data.frame(delta_mat))
}
if (m %in% c(2,4)){
  X = cbind(age_mu, age_diff, gender_diff,
            Region,
            as.data.frame(delta_mat))
}
if (m %in% c(3,5)){
  X = cbind(age_mu, age_diff, gender_diff,
            as.data.frame(delta_mat))
}

Y = fsp[complete.cases(X)]
X = model.matrix(~., X)[,-1]

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_output.pdf"))
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
pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_beta_consistency.pdf"))
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

pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_selprop.pdf"), width = 14)
par(mar=c(15, 5, 1, 1))
plot(selprop_ranked,
     type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
     col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[m], no = "grey"),
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
pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_covariates_only_output.pdf"))
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
pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_covariates_only_beta_consistency.pdf"))
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

pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_covariates_only_selprop.pdf"), width = 14)
par(mar=c(15, 5, 1, 1))
plot(selprop_ranked,
     type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
     col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[m], no = "grey"),
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
pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_exposures_only_output.pdf"))
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
pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_exposures_only_beta_consistency.pdf"))
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

pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_multivariate_exposures_only_selprop.pdf"), width = 14)
par(mar=c(15, 5, 1, 1))
plot(selprop_ranked,
     type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
     col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[m], no = "grey"),
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

### Geographical ----
# Load data
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
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

{pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_batch.pdf"),width=14)
  par(mar=c(0,0,0,0))
  plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
       tip.color=batch.colours[as.character(covars$Batch)])
  dev.off()
}
{pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_graph_expo_all_batch.pdf"))
  g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = batch.colours[as.character(covars$Batch)])
  dev.off()
}
if (m == 4){
  families=unique(covars$Region)
  mycolours = region.colours
  names(mycolours)=families
  
  {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_region.pdf"),width=14)
    par(mar=c(0,0,0,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=region.colours[as.character(covars$Region)])
    dev.off()
  }
  
  {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_graph_expo_all_region.pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = region.colours[as.character(covars$Region)])
    legend("bottomright", pch=19, col=region.colours[levels(covars$Region)],
           legend=levels(covars$Region), cex = 0.4, ncol = 1)
    dev.off()
  }
  
  {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_depart.pdf"),width=14)
    par(mar=c(0,0,0,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=depart.colours[as.character(covars$Department)])
    dev.off()
  }
  
  {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_graph_expo_all_depart.pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = depart.colours[as.character(covars$Department)])
    legend("bottomright", pch=19, col=depart.colours[levels(covars$Department)],
           legend=levels(covars$Department), cex = 0.4, ncol = 1)
    dev.off()
  }
}