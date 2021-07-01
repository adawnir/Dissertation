## Dendrogram visualisation (Pooled LUX/FRA/GS) 
## Rin Wada 30 June

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)
library(RColorBrewer)
library(colorspace)
source("function.R")

### Including isolated----
# Load data sets
mat = readRDS("../Processed/Pooled3/Chemical_compound_matrix_subset_trans_imp.rds")
covar = readRDS("../Results/Pooled3/Clustering_membership.rds")

# Build dataset
sample = covar$Indiv.ID
family = covar$Family.ID
cluster = covar$hc_fs
data <- cbind(sample,family,cluster, mat)
rownames(data) = data[,1]

# Compute Euclidean distance between samples
dist=dist(data[ ,c(4:ncol(data))], method = "euclidean", diag=TRUE)

# Perfor clustering with hclust
hc <- hclust(dist, method = "complete")
dhc <- as.dendrogram(hc)

#So if I Want to color each leaf of the Tree, I have to change the attribute of each leaf. This can be done using the dendrapply function. So I create a function that # # add 3 attributes to the leaf : one for the color (“lab.col”) ,one for the font “lab.font” and one for the size (“lab.cex”).
i=0
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    line=match(attributes(n)$label,data[,1])
    mycolours = random.colour.pal(length(unique(data[,2])), seed = 100)
    mycolours2 = random.colour.pal(length(unique(data[,3])), seed = 300621)
    col_family= mycolours[match(data[line,2], levels(factor(data[,2])))];
    col_cluster=mycolours2[match(data[line,3], levels(factor(data[,3])))];
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,pch=20,col=col_family))
    attr(n,"edgePar")<-c(a$edgePar,list(col=darken(col_cluster, 0.2), lwd = 1)) }
  return(n)
}

# Finally I just have to apply this to my dendrogram
dL <- dendrapply(dhc, colLab)

pdf("../Figures/Pooled3/HC_complete_dendrogram.pdf")
plot(dL, leaflab = "none", main="Pooled3 (N=36)")
dev.off()

