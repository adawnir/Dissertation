## Pooling data sets
## Rin Wada 21 June 2021

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(imputeLCMD)
source("function.R")

# Load data sets
covar_lux = readRDS("../Processed/Luxembourg/Participant_covariate_info_subset.rds")
chem_lux = readRDS("../Processed/Luxembourg/Chemical_compound_info.rds")
mat_lux = readRDS("../Processed/Luxembourg/Chemical_compound_matrix_raw.rds")

covar_fra = readRDS("../Processed/France/Participant_covariate_info_subset.rds")
chem_fra = readRDS("../Processed/France/Chemical_compound_info.rds")
mat_fra = readRDS("../Processed/France/Chemical_compound_matrix_raw.rds")

covar_gs = readRDS("../Processed/GrandeSynthe/Participant_covariate_info_subset.rds")
chem_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_info.rds")
mat_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_matrix_raw.rds")

### Covariate information ----
# Merge covariate information
covar = bind_rows(covar_lux, covar_fra, covar_gs) %>%
  select(Batch, Indiv.ID, Age, Gender, Weight, Length, Family.ID,
        Area, Department, Region, Country)
rownames(covar) = covar$Indiv.ID
str(covar)

ifelse(dir.exists("../Figures/Pooled3"),"",dir.create("../Figures/Pooled3"))
# Age
pdf("../Figures/Pooled3/Age_dist.pdf", width = 5, height = 5)
hist(covar$Age, main = "Pooled (LUX/FRA/GS)", xlab = "Age (years)")
dev.off()

# Table 1
# N
nrow(covar)
# Families & number of isolated
table(table(covar$Family.ID))
# Age
cat(round(mean(covar$Age, na.rm = T),2),
    " (", round(sd(covar$Age, na.rm = T),2), ")\n", sep = "")
cat(sum(is.na(covar$Age)),
    " (", round(sum(is.na(covar$Age))/length(covar$Age)*100,2), "%)\n", sep = "")
# Gender
table(covar$Gender, useNA = "ifany")
round(prop.table(table(covar$Gender, useNA = "ifany"))*100,2)
cat(sum(is.na(covar$Gender)),
    " (", round(sum(is.na(covar$Gender))/length(covar$Gender)*100,2), "%)\n", sep = "")

### Chemical compound information ----
# Remove Fipronil and Fipronil-sulfonate extracted using SPME
chem_lux = chem_lux[-which(chem_lux$Compound %in% c("Fipronil.SPME", "Fipronil-sulfone.SPME")),]
mat_lux = mat_lux[,-which(colnames(mat_lux) %in% c("Fipronil.SPME", "Fipronil-sulfone.SPME"))]

# Check extraction consistency
extract_diff = full_join(chem_lux, chem_fra, by = "Compound", suffix = c(".lux", ".fra")) %>%
  full_join(chem_gs, by = "Compound") %>%
  mutate(Family = coalesce(Family,Family)) %>%
  select(Compound, Family, starts_with("Extraction")) %>%
  .[!is.equal(list(.$Extraction.lux, .$Extraction.fra, .$Extraction)),]

# Merge chemical compound information
mat = bind_rows(mat_lux, mat_fra, mat_gs)
mat = mat[covar$Indiv.ID,]

chem = full_join(chem_lux, chem_fra, by = "Compound", suffix = c(".lux", ".fra")) %>%
  full_join(chem_gs, by = "Compound") %>%
  mutate(Family = coalesce(Family,Family)) %>%
  select(Compound, Family, starts_with("LOD"), starts_with("LOQ"),starts_with("Extraction"))
colnames(chem) = ifelse(grepl("\\.", colnames(chem))|colnames(chem)%in%c("Compound", "Family"),
                        colnames(chem), paste0(colnames(chem),".gs"))

# Proportion of NAs per chemical compounds
all(colnames(mat)==chem$Compound)
chem$NA_prop = apply(mat, 2, function(x) sum(is.na(x))/nrow(mat))

# Proportion observed (not missing)
tmp = chem$NA_prop
names(tmp) = chem$Compound
prop = 1-tmp
prop[which(is.na(prop))] = 0 # Replace NA with 0
prop = prop[order(chem$Family,-prop)]
mylabels = names(prop)

background = TRUE
family = rep(levels(chem$Family), as.numeric(table(chem$Family)))
mycolours = c("navy")
background_colour = rep(family.colours,
                        times = as.numeric(table(chem$Family)))
myspacing = 1
xseq = seq(1:length(prop))

{pdf("../Figures/Pooled3/Observed_prop.pdf", width=14, height=7)
  par(mar=c(7,4,9,1))
  plot(prop,
       col=mycolours,
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1))
  xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
  if (background){
    for (k in seq(1,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(background_colour[k],0.95), border=NA)
    }
    for (k in seq(2,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(background_colour[k],0.99), border=NA)
    }
    box()
  }
  abline(v=xseqgreysep,lty=1,lwd=0.1,col="grey")
  par(new = TRUE)
  plot(prop,
       col= mycolours,
       xaxt="n", ylab="Proportion observed", xlab = "",
       type="h", lwd=2, ylim = c(0,1))
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.5)}
  xseqblack=c(xseq[!duplicated(family)]-myspacing/2, max(xseq)+myspacing/2)
  abline(v=xseqblack,lty=3,col="black")
  for (k in 1:(length(xseqblack)-1)){
    axis(side=3, at=xseqblack[c(k,k+1)]+c(2,-2), line=0.5, labels=NA)
  }
  for (k in 1:(length(xseqblack)-1)){
    axis(side=3, at=mean(xseqblack[c(k,k+1)]), line=0.2, tick=FALSE,
         labels=unique(family)[k], cex.axis = 0.7, las = 2,
         col.axis = darken(family.colours[k],0.5))
  }
  dev.off()
}

# Proportion of detected per chemical compounds
all(colnames(mat)==chem$Compound)
chem$nd_prop = apply(mat, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(mat))

tmp = chem$nd_prop + chem$NA_prop
names(tmp) = chem$Compound
prop = 1-tmp
prop = prop[order(chem$Family,-prop)]
mylabels = names(prop)

background = TRUE
family = rep(levels(chem$Family), as.numeric(table(chem$Family)))
mycolours = c("navy")
background_colour = rep(family.colours,
                        times = as.numeric(table(chem$Family)))
myspacing = 1
xseq = seq(1:length(prop))

{pdf("../Figures/Pooled3/Detection_prop.pdf", width=14, height=7)
  par(mar=c(7,4,9,1))
  plot(prop,
       col=mycolours,
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1))
  xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
  if (background){
    for (k in seq(1,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(background_colour[k],0.95), border=NA)
    }
    for (k in seq(2,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(background_colour[k],0.99), border=NA)
    }
    box()
  }
  abline(v=xseqgreysep,lty=1,lwd=0.1,col="grey")
  par(new = TRUE)
  plot(prop,
       col= mycolours,
       xaxt="n", ylab="Proportion observed", xlab = "",
       type="h", lwd=2, ylim = c(0,1))
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.5)}
  xseqblack=c(xseq[!duplicated(family)]-myspacing/2, max(xseq)+myspacing/2)
  abline(v=xseqblack,lty=3,col="black")
  for (k in 1:(length(xseqblack)-1)){
    axis(side=3, at=xseqblack[c(k,k+1)]+c(2,-2), line=0.5, labels=NA)
  }
  for (k in 1:(length(xseqblack)-1)){
    axis(side=3, at=mean(xseqblack[c(k,k+1)]), line=0.2, tick=FALSE,
         labels=unique(family)[k], cex.axis = 0.7, las = 2,
         col.axis = darken(family.colours[k],0.5))
  }
  dev.off()
}

# Setting detection threshold
sorted_prop = sapply(prop, sort)
thrseq=seq(0,1,by=0.05)
N = NULL
for (k in thrseq){
  N=c(N,sum(sorted_prop>k))
}

{pdf("../Figures/Pooled3/Thresholds_proportion_detection.pdf", width=10, height=7)
  par(mar=c(5,5,1,1))
  plot(thrseq, N, type="b", pch=19, col=mycolours[1], las=1, cex.lab=1.5,
       xlab="Threshold in proportion of detection",
       ylab="Number of kept variables", ylim = c(0,152))
  abline(v = 0.1, lty = 2)
  dev.off()
}

{pdf("../Figures/Pooled3/Detection_prop_threshold.pdf", width=14, height=7)
  par(mar=c(7,4,9,1))
  plot(prop,
       col=ifelse(prop>0.1, mycolours, lighten(mycolours, 0.5)),
       xaxt="n", ylab="", xlab = "",
       type="n", lwd=2, ylim = c(0,1))
  xseqgreysep=c(min(xseq)-myspacing/2,apply(rbind(xseq[-1],xseq[-length(xseq)]),2,mean),max(xseq)+myspacing/2)
  if (background){
    for (k in seq(1,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(background_colour[k],0.95), border=NA)
    }
    for (k in seq(2,length(xseqgreysep),by=2)){
      polygon(x=c(xseqgreysep[k],xseqgreysep[k+1],xseqgreysep[k+1],xseqgreysep[k]),
              y=c(-10,-10,10,10), col=lighten(background_colour[k],0.99), border=NA)
    }
    box()
  }
  abline(v=xseqgreysep,lty=1,lwd=0.1,col="grey")
  par(new = TRUE)
  plot(prop,
       col= ifelse(prop>0.1, mycolours, lighten(mycolours, 0.5)),
       xaxt="n", ylab="Proportion observed", xlab = "",
       type="h", lwd=2, ylim = c(0,1))
  abline(h=0.1,lty=2)
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = mylabels[i], las=2, cex.axis = 0.5)}
  xseqblack=c(xseq[!duplicated(family)]-myspacing/2, max(xseq)+myspacing/2)
  abline(v=xseqblack,lty=3,col="black")
  for (k in 1:(length(xseqblack)-1)){
    axis(side=3, at=xseqblack[c(k,k+1)]+c(2,-2), line=0.5, labels=NA)
  }
  for (k in 1:(length(xseqblack)-1)){
    axis(side=3, at=mean(xseqblack[c(k,k+1)]), line=0.2, tick=FALSE,
         labels=unique(family)[k], cex.axis = 0.7, las = 2,
         col.axis = darken(family.colours[k],0.5))
  }
  dev.off()
}

### Inclusion/exclusion criteria (Chemical compound) ----
# Number of chemical compounds with 90% or more nd or NA
sum(chem$nd_prop + chem$NA_prop>=0.9)

# Filter out >=90% nd or NA
mat = mat[,which(chem$nd_prop + chem$NA_prop<0.9)]
chem = chem[which(chem$nd_prop + chem$NA_prop<0.9),]
ncol(mat)

max(rowSums(is.na(mat))/nrow(mat))

# Save data sets
saveRDS(covar, "../Processed/Pooled3/Participant_covariate_info_subset.rds")
saveRDS(chem, "../Processed/Pooled3/Chemical_compound_info_subset.rds")
saveRDS(mat, "../Processed/Pooled3/Chemical_compound_matrix_raw_subset.rds")

### Recoding nd ----
# Replace nd with random values from 0 to minimum detection (Gaussian)
min.lod = chem[,grepl("LOD",colnames(chem))] %>%
  na_if(0) %>%
  apply(., 1, function(x) min(x, na.rm = T))

# Set LOD for each data set
# If LOD is 0 or missing set it as minimum of three data sets
tmp = mat[covar_lux$Indiv.ID,]
lod = ifelse(chem$LOD.lux==0|is.na(chem$LOD.lux), min.lod, chem$LOD.lux)
names(lod) = chem$Compound

tmp2 = mat[covar_fra$Indiv.ID,]
lod2 = ifelse(chem$LOD.fra==0|is.na(chem$LOD.fra), min.lod, chem$LOD.fra)
names(lod2) = chem$Compound

tmp3 = mat[covar_gs$Indiv.ID,]
lod3 = ifelse(chem$LOD.gs==0|is.na(chem$LOD.gs), min.lod, chem$LOD.gs)
names(lod3) = chem$Compound
for(k in 1:ncol(mat)){
  set.seed(150621)
  tmp[tmp[,k]=="nd"&!is.na(tmp[,k]),k] = rtruncnorm(sum(tmp[,k]=="nd"&!is.na(tmp[,k])),
                                                    min=0,
                                                    max=lod[colnames(mat)[k]]) 
  tmp2[tmp2[,k]=="nd"&!is.na(tmp2[,k]),k] = rtruncnorm(sum(tmp2[,k]=="nd"&!is.na(tmp2[,k])),
                                                    min=0,
                                                    max=lod2[colnames(mat)[k]]) 
  tmp3[tmp3[,k]=="nd"&!is.na(tmp3[,k]),k] = rtruncnorm(sum(tmp3[,k]=="nd"&!is.na(tmp3[,k])),
                                                    min=0,
                                                    max=lod3[colnames(mat)[k]]) 
}

# Merge
mat_pooled = rbind(tmp, tmp2, tmp3)

# Convert to numeric (Replace any characters with NaN and add rownames)
mat_pooled = apply(mat_pooled, 2, as.numeric)
rownames(mat_pooled) = covar$Indiv.ID

summary(mat_pooled)
sum(is.infinite(mat_pooled))
# Save data sets
saveRDS(mat_pooled, "../Processed/Pooled3/Chemical_compound_matrix_subset.rds")

### Transformation ----
# Density plots
ncol(mat_pooled)
ifelse(dir.exists("../Figures/Pooled3"),"",dir.create("../Figures/Pooled3"))
pdf("../Figures/Pooled3/Compound_dist.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_pooled)){
  plot(density(mat_pooled[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_pooled)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

# Log10 Transformation
mat_pooled = log10(mat_pooled)

sum(is.infinite(mat_pooled))

# Density plots
ncol(mat_pooled)
pdf("../Figures/Pooled3/Compound_dist_trans.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_pooled)){
  plot(density(mat_pooled[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_pooled)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

# Save data sets
saveRDS(mat_pooled, "../Processed/Pooled3/Chemical_compound_matrix_subset_trans.rds")

### Imputation ----
# impute.QRILC
set.seed(7)
mat_pooled_imp=impute.QRILC(t(mat_pooled))
mat_pooled_imp=t(mat_pooled_imp[[1]])
rownames(mat_pooled_imp) = rownames(mat_pooled)
colnames(mat_pooled_imp) = colnames(mat_pooled)

# Save data sets
saveRDS(mat_pooled_imp, "../Processed/Pooled3/Chemical_compound_matrix_subset_trans_imp.rds")

# Compounds with different methods
nrow(extract_diff)
tmp = mat_pooled_imp[covar_lux$Indiv.ID,]
tmp2 = mat_pooled_imp[covar_fra$Indiv.ID,]
tmp3 = mat_pooled_imp[covar_gs$Indiv.ID,]

diff = intersect(colnames(mat_pooled_imp), extract_diff$Compound)
pdf("../Figures/Pooled3/Extract_diff_Compound_dist_trans_imp.pdf", height = 8, width = 11)
par(mfcol = c(3, 6), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:length(diff)){
  plot(density(tmp[,which(colnames(mat_pooled_imp)==diff[k])], na.rm = T),
       col="navy", xlab="", ylab="", main=paste0(diff[k], " (LUX)"),
       cex.main = 0.5, cex.axis=0.5, tck=-0.05, xlim = c(-10,10))
  plot(density(tmp2[,which(colnames(mat_pooled_imp)==diff[k])], na.rm = T),
       col="navy", xlab="", ylab="", main=paste0(diff[k], " (FRA)"),
       cex.main = 0.5, cex.axis=0.5, tck=-0.05, xlim = c(-10,10))
  plot(density(tmp3[,which(colnames(mat_pooled_imp)==diff[k])], na.rm = T),
       col="navy", xlab="", ylab="", main=paste0(diff[k], " (GS)"),
       cex.main = 0.5, cex.axis=0.5, tck=-0.05, xlim = c(-10,10))
}
dev.off()
