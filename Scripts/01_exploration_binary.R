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
dat = read.xlsx("../Data/Data_sample.xlsx")

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

# BINARY FEATURES: Quality control - Replace No result with NA and nd with 0
mat = na_if(mat, "No result")
mat[mat=="nd"] = 0
mat = apply(mat, 2, as.numeric)
summary(mat)
NA_prop = apply(mat, 2, function(x) sum(is.na(x))/nrow(mat))
comp$NA_prop = NA_prop

# Quality control - Replace values below LOQ with 0 and above LOQ with 1 
for(k in 1:ncol(mat)){
  mat[,k]=ifelse(mat[,k]<comp$LOQ[k],0,1)
}
LOQ_prop = apply(mat, 2, function(x) sum(x == 0)/nrow(mat))
comp$LOQ_prop = LOQ_prop

# Plot
prop = sort(LOQ_prop, decreasing = TRUE)
pdf(paste0("../Figures/Low_prop_binary.pdf"), height=5, width=14)
par(mar=c(9,5,1,1))
plot(prop,
     col=ifelse(prop==1,"red","navy"),
     xaxt="n", ylab="Proportion not detected/<LOQ", xlab = "",
     type="h", lwd=3)
for(i in 1:length(prop)){
  axis(1, at=i, labels = names(prop)[i], las=2, cex.axis = 0.5,
       col.axis = ifelse(names(prop)[i] %in% mylist, "red", "black")) 
}
dev.off()

# Plot
prop = sort(NA_prop, decreasing = TRUE)
pdf(paste0("../Figures/Missing_Prop_binary.pdf"), height=5, width=14)
par(mar=c(9,5,1,1))
plot(prop,
     col=ifelse(prop<0.5,"navy","red"),
     xaxt="n", ylab="Proportion no results", xlab = "",
     type="h", lwd=3)
for(i in 1:length(prop)){
  axis(1, at=i, labels = names(prop)[i], las=2, cex.axis = 0.5,
       col.axis = ifelse(names(prop)[i] %in% mylist, "red", "black")) 
}
dev.off()

# Filter out low quality components
comp_sub = comp[-which(comp$LOQ_prop==1),]
comp_sub = comp_sub[-which(comp_sub$NA_prop>0.5),]
comp_sub = comp_sub[which(comp_sub$Quantitative=="Yes"),]

mat_sub = mat[,which(colnames(mat) %in% comp_sub$Component)]

# Save processed data sets
rownames(mat) = rownames(mat_sub)

ifelse(dir.exists("../Processed"),"",dir.create("../Processed"))
saveRDS(comp_sub, "../Processed/Chemical_component_info_subset_binary.rds")
saveRDS(mat_sub, "../Processed/Chemical_component_matrix_subset_binary.rds")

