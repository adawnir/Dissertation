## Dendrogram visualisation (excluding isolated children) 
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

# Load data sets
mat_lux = readRDS("../Results/Luxembourg/Chemical_compound_matrix_no_isolated.rds")
mat_fra = readRDS("../Results/France/Chemical_compound_matrix_no_isolated.rds")
mat_gs = readRDS("../Results/GrandeSynthe/Chemical_compound_matrix_no_isolated.rds")

covar_lux = readRDS("../Results/Luxembourg/Clustering_membership_no_isolated.rds")
covar_fra = readRDS("../Results/France/Clustering_membership_no_isolated.rds")
covar_gs = readRDS("../Results/GrandeSynthe/Clustering_membership_no_isolated.rds")

### Luxembourg----
# Build dataset
sample = covar_lux$Indiv.ID
family = covar_lux$Family.ID
cluster = covar_lux$hc_fs
data <- cbind(sample,family,cluster, mat_lux)
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

pdf("../Figures/Luxembourg/HC_complete_dendrogram.pdf")
plot(dL, leaflab = "none", main="Luxembourg (N=36)")
dev.off()

### France----
# Build dataset
sample = covar_fra$Indiv.ID
family = covar_fra$Family.ID
cluster = covar_fra$hc_fs
data <- cbind(sample,family,cluster, mat_fra)
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

pdf("../Figures/France/HC_complete_dendrogram.pdf")
plot(dL, leaflab = "none", main="France (N=70)")
dev.off()

### GrandeSynthe----
# Build dataset
sample = covar_gs$Indiv.ID
family = covar_gs$Family.ID
cluster = covar_gs$hc_fs
data <- cbind(sample,family,cluster, mat_gs)
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

pdf("../Figures/GrandeSynthe/HC_complete_dendrogram.pdf")
plot(dL, leaflab = "none", main="GrandeSynthe (N=26)")
dev.off()
