## Data exploration
## Rin Wada 24 May

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(openxlsx)
library(tidyverse)

# Load data
lux = read.xlsx("../Data/Data_sample.xlsx")

# Covariate information
covar = dat[1:39,1:7]
rownames(covar) = covar$ID
colnames(covar)[6:7] = c("Weight (mg)","Length (cm)")

# tabe 1
covar$Gender[which(covar$Gender=="Jules")] = "Male"
table(covar$Gender)
mean(covar$Age)
sd(covar$Age)
aggregate(Age ~ Gender, covar, mean)
aggregate(Age ~ Gender, covar, sd)
aggregate(Gender ~ FAMILY.ID, covar, function(x) sum(x=="Female")/length(x))

# Chemical component matrix
mat = dat[1:39,8:ncol(dat)]
rownames(mat) = covar$ID

# Chemical group, extraction method and LOQ
comp = data.frame(Component = colnames(mat),
                  LOQ = as.numeric(gsub("\\*","",dat[40,8:ncol(dat)])),
                  Group = unlist(dat[42,8:ncol(dat)]),
                  Extraction = unlist(dat[41,8:ncol(dat)]))
comp = comp %>% fill(Group)
mylist = c("DMP", "DMDTP", "DEP", "DEDTP", "ClCF3CA", "Propargite",
           "Spinosyn.A", "Bisphenol.A", "Bisphenol.S**")
comp$Quantitative = ifelse(comp$Component %in% mylist,
                           "No","Yes")
table(comp$Quantitative)

# Fine print
dat[43:nrow(dat),1]

# Quality control - Replace nd and No result with NA
mat = na_if(mat, "nd")
mat = na_if(mat, "No result")
mat = apply(mat, 2, as.numeric)
summary(mat)
NA_prop = apply(mat, 2, function(x) sum(is.na(x))/nrow(mat))
comp$NA_prop = NA_prop

# Quality control - Replace values below LOQ with NA
LOQ_prop = NULL
for(k in 1:ncol(mat)){
  LOQ_prop = c(LOQ_prop, sum(mat[,k]<comp$LOQ[k], na.rm = T)/nrow(mat))
  mat[,k]=ifelse(mat[,k]<comp$LOQ[k],NA,mat[,k])
}
comp$LOQ_prop = LOQ_prop
prop = apply(mat, 2, function(x) sum(is.na(x))/nrow(mat))
comp$Total_prop = prop
prop = sort(prop, decreasing = TRUE)

# Plot
ifelse(dir.exists("../Figures"),"",dir.create("../Figures"))
pdf(paste0("../Figures/Missing_Prop.pdf"), height=5, width=14)
par(mar=c(9,5,1,1))
plot(prop,
     col=ifelse(prop<0.5,"navy","red"),
     xaxt="n", ylab="Proportion missing", xlab = "",
     type="h", lwd=3)
abline(h=0.5,lty=2,col="black")
for(i in 1:length(prop)){
  axis(1, at=i, labels = names(prop)[i], las=2, cex.axis = 0.5,
       col.axis = ifelse(names(prop)[i] %in% mylist, "red", "black")) 
}
dev.off()

# Filter out low quality components
comp_sub = comp[which(comp$Total_prop<0.5),]
comp_sub = comp_sub[which(comp_sub$Quantitative=="Yes"),]

mat_sub = mat[,which(colnames(mat) %in% comp_sub$Component)]

# Density plots
pdf("../Figures/Component_dist.pdf")
par(mfrow = c(3, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_sub)){
  xfull=density(mat_sub[,k], na.rm = T)
  plot(xfull, col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_sub)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

# Transformation
mat_sub_trans = log10(mat_sub)
pdf("../Figures/Component_dist_transformed.pdf")
par(mfrow = c(3, 5), oma = c(0.5,0.5,0,0.5), mar = c(1,1,1,1), mgp=c(0,0.2,0))
for (k in 1:ncol(mat_sub_trans)){
  xfull=density(mat_sub_trans[,k], na.rm = T)
  plot(xfull, col="navy", xlab="", ylab="", xaxt="n", main=colnames(mat_sub_trans)[k],
       cex.main = 0.5, cex.axis=0.5, tck=-0.05)
}
dev.off()

# Participants with high missing rates
NA_prop = apply(mat_sub, 1, function(x) sum(is.na(x))/ncol(mat_sub))
pdf("../Figures/Participant_missing_prop_dist.pdf")
boxplot(NA_prop, ylab = "Proportion missing", main = "Participant missing rate")
dev.off()

# Save processed data sets
rownames(mat) = rownames(mat_sub) = rownames(mat_sub_trans) = rownames(covar)

ifelse(dir.exists("../Processed"),"",dir.create("../Processed"))
saveRDS(comp, "../Processed/Chemical_component_info.rds")
saveRDS(comp_sub, "../Processed/Chemical_component_info_subset.rds")
saveRDS(covar, "../Processed/Participant_covariate_info.rds")
saveRDS(mat, "../Processed/Chemical_component_matrix.rds")
saveRDS(mat_sub, "../Processed/Chemical_component_matrix_subset.rds")
saveRDS(mat_sub_trans, "../Processed/Chemical_component_matrix_subset_transformed.rds")


