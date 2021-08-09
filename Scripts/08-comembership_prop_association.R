## Multivariate clustering - Association with comembership proportion
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

covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
X = readRDS(paste0("../Results/",filepaths[m],"/Family_covariates_delta_exposures.rds"))
out = readRDS(paste0("../Results/",filepaths[m],"/Stability_clustering_output.rds"))

selprop = SelectionProportions(out)
families=levels(covars$Family.ID)

fprop=NULL
for (f in 1:length(families)){
  tmpmat=prop[covars$Family.ID==families[f],covars$Family.ID==families[f]]
  fprop=c(fprop,mean(tmpmat[upper.tri(tmpmat)]))
}
names(fprop)=families

betas = pvals = NULL
f1='Y ~ X[,k]'
t0=Sys.time()
for (k in 1:ncol(X)){
  Y = fprop
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
mylabels = gsub("Region","Île-de-France (Y/N)",mylabels)

names(pvals)=names(betas)=mylabels

saveRDS(pvals, paste0("../Results/",filepaths[m],"/Comembership_prop_univar_pvals.rds"))
saveRDS(betas, paste0("../Results/",filepaths[m],"/Comembership_prop_univar_betas.rds"))

annot_sub = annot[which(names(annot) %in% mylabels)]

mycolours = c(rep(batch.colours[m],length(mylabels)-length(annot_sub)),annot.colours[annot_sub])

{pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_univariate.pdf"))
  par(mar=c(5,5,1,1))
  plot(betas, -log10(pvals), pch=19,
       col=ifelse(pvals < 0.05, mycolours, "grey"),
       cex.lab=1, cex = 0.7,
       ylim = c(0, max(-log10(pvals))+0.25),
       xlim = c(-max(abs(betas))-5, max(abs(betas))+5),
       ylab=expression(-log[10](p)), 
       xlab=expression(beta))
  for (k in 1:length(mylabels)){
    if(mylabels[k] %in% names(annot_sub)){
      if(pvals[k] < 0.05){
        label = substitute(Delta~tmp, list(tmp=mylabels[k]))
        text(betas[k]+sign(betas[k])*2, -log10(pvals[k])+0.1,
             labels = label, col = mycolours[k])
      } else {
        text(betas[k]+sign(betas[k])*2, -log10(pvals[k])+0.1,
             labels = mylabels[k], col = mycolours[k]) 
      }
    }
  }
  text(max(abs(betas))+4.9, -log10(0.05/ncol(X))+0.01, paste0("Bonferroni threshold = ",formatC(-log10(0.05/ncol(X)), format = "e")), adj = c(1,0), col = "darkred")
  abline(h = -log10(0.05/ncol(X)), lty = 2, col = "darkred")
  dev.off()
}

### Multivariate analysis using stability selection sPLS regression ----

Y = fprop[complete.cases(X)]
X = X[complete.cases(X),]
X = model.matrix(~., X)[,-1]

stab = VariableSelection(xdata = X, ydata = Y, implementation = SparsePLS)

saveRDS(stab, paste0("../Results/",filepaths[m],"/Comembership_prop_univar_stab.rds"))

pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_multivariate_output.pdf"))
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
pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_multivariate_beta_consistency.pdf"))
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

pdf(paste0("../Figures/",filepaths[m],"/Comembership_prop_multivariate_selprop.pdf"), width = 14)
par(mar=c(15, 5, 1, 1))
plot(selprop_ranked,
     type = "h", lwd = 3, las = 1, cex.lab = 1.3, bty = "n", ylim = c(0, 1),
     col = ifelse(selprop_ranked >= Argmax(stab)[2], yes = batch.colours[m], no = "grey"),
     xaxt = "n", xlab = "", ylab = "Selection proportions"
)
text(length(selprop_ranked)-0.01, Argmax(stab)[2]+0.01, paste0("Selection proportion threshold = ",Argmax(stab)[2]), adj = c(1,0), col = "darkred")
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
  if (mylabels[k] %in% names(annot_sub)){
    label = substitute(Delta~tmp, list(tmp = mylabels[k]))
  } else {
    label = mylabels[k]
  }
  axis(side = 1, at = k, las = 2, labels = label,
       col.axis = ifelse(selprop_ranked[k] >= Argmax(stab)[2], col, "grey"))
}
dev.off()
