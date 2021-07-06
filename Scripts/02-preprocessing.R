## Data exploration and pre-processing
## Rin Wada 15 June

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(colorspace)
library(imputeLCMD)
source("function.R")

# Load data sets
covar_lux = readRDS("../Processed/Luxembourg/Participant_covariate_info.rds")
chem_lux = readRDS("../Processed/Luxembourg/Chemical_compound_info.rds")
mat_lux = readRDS("../Processed/Luxembourg/Chemical_compound_matrix_raw.rds")

covar_fra = readRDS("../Processed/France/Participant_covariate_info.rds")
chem_fra = readRDS("../Processed/France/Chemical_compound_info.rds")
mat_fra = readRDS("../Processed/France/Chemical_compound_matrix_raw.rds")

covar_gs = readRDS("../Processed/GrandeSynthe/Participant_covariate_info.rds")
chem_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_info.rds")
mat_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_matrix_raw.rds")

### Inclusion/exclusion criteria (Participant) ----
min(covar_lux$Age, na.rm = T)
max(covar_lux$Age, na.rm = T)

min(covar_fra$Age, na.rm = T)
max(covar_fra$Age, na.rm = T)

min(covar_gs$Age, na.rm = T)
max(covar_gs$Age, na.rm = T)

# Exclude adults aged 18 or over
covar_gs = covar_gs[-which(covar_gs$Age >= 18),]
max(covar_gs$Age, na.rm = T)

pdf("../Figures/Age_dist_subset.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
hist(covar_lux$Age, main = "Luxembourg", xlab = "Age (years)")
hist(covar_fra$Age, main = "France", xlab = "Age (years)")
hist(covar_gs$Age, main = "Grande-Synthe", xlab = "Age (years)")
dev.off()

# Exclude hair samples with <35mg
covar_lux = covar_lux[-which(covar_lux$Weight < 35),]
min(covar_lux$Weight, na.rm = T)

# Adjust chemical matrices
mat_lux = mat_lux[covar_lux$Indiv.ID,]
mat_gs = mat_gs[covar_gs$Indiv.ID,]

pdf("../Figures/Sample_length_weight_subset.pdf", width=8, height=4)
par(mfrow = c(1,2))
hist(covar_lux$Length, xlab = "Length (cm)", main = "Length of hair sample (Luxembourg)")
hist(covar_lux$Weight, xlab = "Weight (mg)", main = "Weight of hair sample (Luxembourg)")
dev.off()

# Randomly sampling 2 children out of 3-4 child groups
set.seed(280621)
covar_lux = covar_lux[which(group_sample(covar_lux$Family.ID, covar_lux$Family.ID, 2)),]
covar_fra = covar_fra[which(group_sample(covar_fra$Siblings.Groups, covar_fra$Siblings.Groups, 2, "Isolated")),]
covar_gs = covar_gs[which(group_sample(covar_gs$Siblings.Groups, covar_gs$Siblings.Groups, 2, "Isolated")),]
nrow(covar_lux)
nrow(covar_fra)
nrow(covar_gs)

# Adjust chemical matrices
mat_lux = mat_lux[covar_lux$Indiv.ID,]
mat_fra = mat_fra[covar_fra$Indiv.ID,]
mat_gs = mat_gs[covar_gs$Indiv.ID,]

### Covariate infromation ---- 
# Recoding Family Sibling ID
covar_lux$Family.ID = ifelse(covar_lux$Sibling.ID=="Only child","Isolated",covar_lux$Family.ID)
covar_lux$Sibling.ID = ifelse(covar_lux$Sibling.ID=="Only child","Isolated",covar_lux$Sibling.ID)
table(covar_lux$Family.ID)

# Adjust family id
tmp = covar_lux$Family.ID[which(covar_lux$Family.ID!="Isolated")]
t = table(covar_lux$Family.ID)[unique(tmp)]
family.number = rep(paste0("L",seq(1, length.out = sum(!duplicated(tmp)))),t)
covar_lux$Family.ID[which(covar_lux$Family.ID!="Isolated")] = family.number
covar_lux$Sibling.ID[which(covar_lux$Family.ID!="Isolated")] = ave(tmp, tmp, FUN=function(a) 1:length(a))
table(covar_lux$Family.ID)
table(covar_lux$Sibling.ID)
covar_lux$Area = NA
covar_lux$Department = NA
covar_lux$Region = "Luxembourg"
covar_lux$Country = "Luxembourg"
covar_lux$Batch = "LUX"
str(covar_lux)
covar_lux = covar_lux %>% select(Indiv.ID, Age, Gender, Family.ID, Sibling.ID,
                                 Weight, Length,
                                 Area, Department, Region, Country, Batch)

# Recoding Family Sibling ID
table(covar_fra$Siblings.Groups)
# Adjust family id
tmp = covar_fra$Siblings.Groups[which(covar_fra$Siblings.Groups!="Isolated")]
t = table(covar_fra$Siblings.Groups)[unique(tmp)]
family.number = rep(paste0("F",seq(1, length.out = sum(!duplicated(tmp)))),t)
covar_fra$Family.ID = covar_fra$Siblings.Groups
covar_fra$Sibling.ID = covar_fra$Siblings.Groups
covar_fra$Family.ID[which(covar_fra$Family.ID!="Isolated")] = family.number
covar_fra$Sibling.ID[which(covar_fra$Family.ID!="Isolated")] = ave(tmp, tmp, FUN=function(a) 1:length(a))
table(covar_fra$Family.ID)
table(covar_fra$Sibling.ID)
covar_fra$Country = "France"
covar_fra$Batch = "FRA"
covar_fra$Siblings.Groups = NULL
str(covar_fra)
covar_fra = covar_fra %>% select(Indiv.ID, Age, Gender, Family.ID, Sibling.ID,
                                 Area, Department, Region, Country, Batch)

# Recoding Family Sibling ID
table(covar_gs$Siblings.Groups)
covar_gs$Siblings.Groups = ifelse(covar_gs$Siblings.Groups%in%c("Family 2", "Family 3"),
                                  "Isolated",covar_gs$Siblings.Groups)
# Adjust family id
tmp = covar_gs$Siblings.Groups[which(covar_gs$Siblings.Groups!="Isolated")]
t = table(covar_gs$Siblings.Groups)[unique(tmp)]
family.number = rep(paste0("G",seq(1, length.out = sum(!duplicated(tmp)))),t)
covar_gs$Family.ID = covar_gs$Siblings.Groups
covar_gs$Sibling.ID = covar_gs$Siblings.Groups
covar_gs$Family.ID[which(covar_gs$Family.ID!="Isolated")] = family.number
covar_gs$Sibling.ID[which(covar_gs$Family.ID!="Isolated")] = ave(tmp, tmp, FUN=function(a) 1:length(a))
table(covar_gs$Family.ID)
table(covar_gs$Sibling.ID)
covar_gs$Area= "Grande-Synthe"
covar_gs$Department = "Nord"
covar_gs$Region = "Hauts-de-France"
covar_gs$Country = "France"
covar_gs$Batch = "GS"
covar_gs$Siblings.Groups = NULL
str(covar_gs)

# Total number of families
sum(covar_lux$Family.ID=="Isolated") + sum(table(table(covar_lux$Family.ID[which(covar_lux$Family.ID!="Isolated")])))
sum(covar_fra$Family.ID=="Isolated") + sum(table(table(covar_fra$Family.ID[which(covar_fra$Family.ID!="Isolated")])))
sum(covar_gs$Family.ID=="Isolated") + sum(table(table(covar_gs$Family.ID[which(covar_gs$Family.ID!="Isolated")])))

# Number of each type of family
table(table(covar_lux$Family.ID))
table(table(covar_fra$Family.ID))
table(table(covar_gs$Family.ID))

# Age (Mean, SD)
cat(round(mean(covar_lux$Age, na.rm = T),2),
    " (", round(sd(covar_lux$Age, na.rm = T),2), ")\n", sep = "")
cat(round(mean(covar_fra$Age, na.rm = T),2),
    " (", round(sd(covar_fra$Age, na.rm = T),2), ")\n", sep = "")
cat(round(mean(covar_gs$Age, na.rm = T),2),
    " (", round(sd(covar_gs$Age, na.rm = T),2), ")\n", sep = "")

# Age (NA %)
cat(sum(is.na(covar_lux$Age)),
    " (", round(sum(is.na(covar_lux$Age))/length(covar_lux$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_fra$Age)),
    " (", round(sum(is.na(covar_fra$Age))/length(covar_fra$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_gs$Age)),
    " (", round(sum(is.na(covar_gs$Age))/length(covar_gs$Age)*100,2), "%)\n", sep = "")

# Gender (N, %)
table(covar_lux$Gender, useNA = "ifany")
round(prop.table(table(covar_lux$Gender, useNA = "ifany"))*100,2)

table(covar_fra$Gender, useNA = "ifany")
round(prop.table(table(covar_fra$Gender, useNA = "ifany"))*100,2)

table(covar_gs$Gender, useNA = "ifany")
round(prop.table(table(covar_gs$Gender, useNA = "ifany"))*100,2)

# Gender (NA %)
cat(sum(is.na(covar_lux$Gender)),
    " (", round(sum(is.na(covar_lux$Gender))/length(covar_lux$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_fra$Gender)),
    " (", round(sum(is.na(covar_fra$Gender))/length(covar_fra$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_gs$Gender)),
    " (", round(sum(is.na(covar_gs$Gender))/length(covar_gs$Gender)*100,2), "%)\n", sep = "")

### Detection frequency per chemical----
# Number of chemical compounds measured
nrow(chem_lux) # Fipronil and Fipronil-sulfone are in duplicate
nrow(chem_fra)
nrow(chem_gs)

# Remove Fipronil and Fipronil-sulfonate extracted using SPME
chem_lux = chem_lux[-which(chem_lux$Compound %in% c("Fipronil.SPME", "Fipronil-sulfone.SPME")),]
mat_lux = mat_lux[,-which(colnames(mat_lux) %in% c("Fipronil.SPME", "Fipronil-sulfone.SPME"))]

# Proportion of NAs per chemical compounds
all(colnames(mat_lux)==chem_lux$Compound)
chem_lux$NA_prop = apply(mat_lux, 2, function(x) sum(is.na(x))/nrow(mat_lux))

all(colnames(mat_fra)==chem_fra$Compound)
chem_fra$NA_prop = apply(mat_fra, 2, function(x) sum(is.na(x))/nrow(mat_fra))

all(colnames(mat_gs)==chem_gs$Compound)
chem_gs$NA_prop = apply(mat_gs, 2, function(x) sum(is.na(x))/nrow(mat_gs))

tmp = chem_lux$NA_prop
names(tmp) = chem_lux$Compound
tmp2 = chem_fra$NA_prop
names(tmp2) = chem_fra$Compound
tmp3 = chem_gs$NA_prop
names(tmp3) = chem_gs$Compound

# Proportion observed (not missing)
prop = 1-t(bind_rows(tmp, tmp2, tmp3))
prop[which(is.na(prop))] = 0 # Replace NA with 0
prop = prop[order(chem_lux$Family,-prop[,1],-prop[,2],-prop[,3]),]
prop = cbind(prop,rep(NA,length(tmp)),rep(NA,length(tmp)))
mylabels = rownames(prop)
prop = as.vector(t(prop))
prop=prop[c(length(prop),1:length(prop)-1)]

background = TRUE
family = rep(levels(chem_lux$Family), as.numeric(table(chem_lux$Family)))
mycolours = c("tomato","royalblue","forestgreen")
background_colour = rep(family.colours,
                        times = as.numeric(table(chem_lux$Family)))
myspacing = 5
xseq = seq(myspacing/2,length(tmp)*myspacing, by=myspacing)

{pdf("../Figures/Observed_prop.pdf", width=14, height=7)
  par(mar=c(7,4,9,1))
  plot(prop,
       col=c(NA, mycolours, NA),
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
       col=c(NA, mycolours, NA),
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
         col.axis = darken(family.colours, 0.5)[k])
  }
  legend("topright", lty=1, lwd=2, col=mycolours,
         legend = c("Luxembourg","France","Grande-Synthe"), cex=0.7, bg="white")
  dev.off()
}

# Proportion of detected per chemical compounds
all(colnames(mat_lux)==chem_lux$Compound)
chem_lux$nd_prop = apply(mat_lux, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(mat_lux))

all(colnames(mat_fra)==chem_fra$Compound)
chem_fra$nd_prop = apply(mat_fra, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(mat_fra))

all(colnames(mat_gs)==chem_gs$Compound)
chem_gs$nd_prop = apply(mat_gs, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(mat_gs))

tmp = chem_lux$nd_prop + chem_lux$NA_prop
names(tmp) = chem_lux$Compound
tmp2 = chem_fra$nd_prop + chem_fra$NA_prop
names(tmp2) = chem_fra$Compound
tmp3 = chem_gs$nd_prop + chem_gs$NA_prop
names(tmp3) = chem_gs$Compound

# Setting detection threshold
prop = list(tmp, tmp2, tmp3)
sorted_prop = lapply(prop, sort)
N=list(NULL, NULL, NULL)
thrseq=seq(0,1,by=0.05)
for (i in 1:length(sorted_prop)){
  for (k in thrseq){
    N[[i]]=c(N[[i]],sum(sorted_prop[[i]]>k))
  }
}

{pdf("../Figures/Thresholds_proportion_detection.pdf", width=10, height=7)
  par(mar=c(5,5,1,1))
  plot(thrseq, N[[1]], type="b", pch=19, col=mycolours[1], las=1, cex.lab=1.5,
       xlab="Threshold in proportion of detection",
       ylab="Number of kept variables", ylim = c(0,152))
  points(thrseq, N[[2]], type="b", pch=19, col=mycolours[2])
  points(thrseq, N[[3]], type="b", pch=19, col=mycolours[3])
  abline(v = 0.1, lty = 2)
  legend("topright", lty=1, lwd=2, col=mycolours,
         legend = c("Luxembourg","France","Grande-Synthe"), cex=0.7, bg="white")
  dev.off()
  }

prop = 1-t(bind_rows(tmp, tmp2, tmp3))
prop = prop[order(chem_lux$Family,-prop[,1],-prop[,2],-prop[,3]),]
prop = cbind(prop,rep(NA,length(tmp)),rep(NA,length(tmp)))
mylabels = rownames(prop)
prop = as.vector(t(prop))
prop=prop[c(length(prop),1:length(prop)-1)]

background = TRUE
family = rep(levels(chem_lux$Family), as.numeric(table(chem_lux$Family)))
mycolours = c("tomato","royalblue","forestgreen")
background_colour = rep(family.colours,
                        times = as.numeric(table(chem_lux$Family)))

myspacing = 5
xseq = seq(myspacing/2,length(tmp)*myspacing, by=myspacing)

{pdf("../Figures/Detection_prop.pdf", width=14, height=7)
  par(mar=c(7,4,9,1))
  plot(prop,
       col=c(NA, mycolours, NA),
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
       col=c(NA, mycolours, NA),
       xaxt="n", ylab="Proportion detected", xlab = "",
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
         col.axis = darken(family.colours, 0.5)[k])
  }
  legend("topright", lty=1, lwd=2, col=mycolours,
         legend = c("Luxembourg","France","Grande-Synthe"), cex=0.7, bg="white")
  dev.off()
}

{pdf("../Figures/Detection_prop_threshold.pdf", width=14, height=7)
  par(mar=c(7,4,9,1))
  plot(prop,
       col=c(NA, mycolours, NA),
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
       col=ifelse(prop<0.1, lighten(c(NA, mycolours, NA),0.5), c(NA, mycolours, NA)),
       xaxt="n", ylab="Proportion detected", xlab = "",
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
         col.axis = darken(family.colours, 0.5)[k])
  }
  abline(h = 0.1, lty = 2)
  legend("topright", lty=1, lwd=2, col=mycolours,
         legend = c("Luxembourg","France","Grande-Synthe"), cex=0.7, bg="white")
  dev.off()
}

### Inclusion/exclusion criteria (Chemical compound) ----
# Number of chemical compounds with 90% or more nd or NA
sum(chem_lux$nd_prop + chem_lux$NA_prop>=0.9)
sum(chem_fra$nd_prop + chem_fra$NA_prop>=0.9)
sum(chem_gs$nd_prop + chem_gs$NA_prop>=0.9)

# Filter out >=90% nd or NA
mat_lux = mat_lux[,which(chem_lux$nd_prop + chem_lux$NA_prop<0.9)]
mat_fra = mat_fra[,which(chem_fra$nd_prop + chem_fra$NA_prop<0.9)]
mat_gs = mat_gs[,which(chem_gs$nd_prop + chem_gs$NA_prop<0.9)]

chem_lux = chem_lux[which(chem_lux$nd_prop + chem_lux$NA_prop<0.9),]
chem_fra = chem_fra[which(chem_fra$nd_prop + chem_fra$NA_prop<0.9),]
chem_gs = chem_gs[which(chem_gs$nd_prop + chem_gs$NA_prop<0.9),]

ncol(mat_lux)
ncol(mat_fra)
ncol(mat_gs)

# Save data sets
saveRDS(covar_lux, "../Processed/Luxembourg/Participant_covariate_info_subset.rds")
saveRDS(covar_fra, "../Processed/France/Participant_covariate_info_subset.rds")
saveRDS(covar_gs, "../Processed/GrandeSynthe/Participant_covariate_info_subset.rds")

saveRDS(chem_lux, "../Processed/Luxembourg/Chemical_compound_info_subset.rds")
saveRDS(chem_fra, "../Processed/France/Chemical_compound_info_subset.rds")
saveRDS(chem_gs, "../Processed/GrandeSynthe/Chemical_compound_info_subset.rds")

saveRDS(mat_lux, "../Processed/Luxembourg/Chemical_compound_matrix_raw_subset.rds")
saveRDS(mat_fra, "../Processed/France/Chemical_compound_matrix_raw_subset.rds")
saveRDS(mat_gs, "../Processed/GrandeSynthe/Chemical_compound_matrix_raw_subset.rds")

### Recoding nd ----
# Replace nd with random values from 0 to minimum detection (Gaussian)
for(k in 1:ncol(mat_lux)){
  set.seed(150621)
  mat_lux[mat_lux[,k]=="nd"&!is.na(mat_lux[,k]),k] = rtruncnorm(sum(mat_lux[,k]=="nd"&!is.na(mat_lux[,k])), min=0, max=chem_lux$LOD[k]) 
}
for(k in 1:ncol(mat_fra)){
  set.seed(150621)
  mat_fra[mat_fra[,k]=="nd"&!is.na(mat_fra[,k]),k] = rtruncnorm(sum(mat_fra[,k]=="nd"&!is.na(mat_fra[,k])), min=0, max=chem_fra$LOD[k]) 
}
for(k in 1:ncol(mat_gs)){
  set.seed(150621)
  mat_gs[mat_gs[,k]=="nd"&!is.na(mat_gs[,k]),k] = rtruncnorm(sum(mat_gs[,k]=="nd"&!is.na(mat_gs[,k])), min=0, max=chem_gs$LOD[k]) 
}

# Convert to numeric (Replace any characters with NaN and add rownames)
mat_lux = apply(mat_lux, 2, as.numeric)
rownames(mat_lux) = covar_lux$Indiv.ID

mat_fra = apply(mat_fra, 2, as.numeric)
rownames(mat_fra) = covar_fra$Indiv.ID

mat_gs = apply(mat_gs, 2, as.numeric)
rownames(mat_gs) = covar_gs$Indiv.ID

# Save data sets
saveRDS(mat_lux, "../Processed/Luxembourg/Chemical_compound_matrix_subset.rds")
saveRDS(mat_fra, "../Processed/France/Chemical_compound_matrix_subset.rds")
saveRDS(mat_gs, "../Processed/GrandeSynthe/Chemical_compound_matrix_subset.rds")

### Transformation ----
# Density plots
ncol(mat_lux)
ifelse(dir.exists("../Figures/Luxembourg"),"",dir.create("../Figures/Luxembourg"))
pdf("../Figures/Luxembourg/Compound_dist.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_lux)){
  plot(density(mat_lux[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_lux)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

ncol(mat_fra)
ifelse(dir.exists("../Figures/France"),"",dir.create("../Figures/France"))
pdf("../Figures/France/Compound_dist.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_fra)){
  plot(density(mat_fra[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_fra)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

ncol(mat_gs)
ifelse(dir.exists("../Figures/GrandeSynthe"),"",dir.create("../Figures/GrandeSynthe"))
pdf("../Figures/GrandeSynthe/Compound_dist.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_gs)){
  plot(density(mat_gs[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_gs)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

# Log10 Transformation
mat_lux = log10(mat_lux)
mat_fra = log10(mat_fra)
mat_gs = log10(mat_gs)

# Density plots
ncol(mat_lux)
pdf("../Figures/Luxembourg/Compound_dist_trans.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_lux)){
  plot(density(mat_lux[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_lux)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

ncol(mat_fra)
pdf("../Figures/France/Compound_dist_trans.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_fra)){
  plot(density(mat_fra[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_fra)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

ncol(mat_gs)
pdf("../Figures/GrandeSynthe/Compound_dist_trans.pdf", height = 11, width = 8)
par(mfrow = c(7, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_gs)){
  plot(density(mat_gs[,k], na.rm = T),
       col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_gs)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

# Save data sets
saveRDS(mat_lux, "../Processed/Luxembourg/Chemical_compound_matrix_subset_trans.rds")
saveRDS(mat_fra, "../Processed/France/Chemical_compound_matrix_subset_trans.rds")
saveRDS(mat_gs, "../Processed/GrandeSynthe/Chemical_compound_matrix_subset_trans.rds")

### Imputation ----
# Participants with high missing rates
NA_prop = list(apply(mat_lux, 1, function(x) sum(is.na(x))/ncol(mat_lux)),
               apply(mat_fra, 1, function(x) sum(is.na(x))/ncol(mat_fra)),
               apply(mat_gs, 1, function(x) sum(is.na(x))/ncol(mat_gs)))
pdf("../Figures/Participant_missing_prop.pdf", width = 9, height = 5)
par(mfrow = c(1,3))
boxplot(NA_prop[1], ylab = "Participant missing rate", main = "Luxembourg", ylim = c(0,1))
boxplot(NA_prop[2], ylab = "Participant missing rate", main = "France", ylim = c(0,1))
boxplot(NA_prop[3], ylab = "Participant missing rate", main = "Grande-Synthe", ylim = c(0,1))
dev.off()

# impute.QRILC
set.seed(7)
mat_lux_imp=impute.QRILC(t(mat_lux))
mat_lux_imp=t(mat_lux_imp[[1]])
rownames(mat_lux_imp) = rownames(mat_lux)
colnames(mat_lux_imp) = colnames(mat_lux)

set.seed(7)
mat_fra_imp=impute.QRILC(t(mat_fra))
mat_fra_imp=t(mat_fra_imp[[1]])
rownames(mat_fra_imp) = rownames(mat_fra)
colnames(mat_fra_imp) = colnames(mat_fra)

# Save data sets
saveRDS(mat_lux_imp, "../Processed/Luxembourg/Chemical_compound_matrix_subset_trans_imp.rds")
saveRDS(mat_fra_imp, "../Processed/France/Chemical_compound_matrix_subset_trans_imp.rds")
