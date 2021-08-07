# Community detection: PFER
# Rin 4 Aug

# Load packages
library(tidyverse)
library(igraph)
library(colorspace)
library(focus)
library(parallel)

# Initialise
rm(list=ls())
path="~/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

## Parameters
args=commandArgs(trailingOnly=TRUE)
m=as.numeric(args[2])
nchunks=as.numeric(args[1])

### PFER calibration ----
suffix = c("lux","fra","gs","pooled3","pooled2")
# Load data
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))

X = t(expo)
start.type = ifelse(m %in% c(4,5),"cold","warm")
n = c(Inf, seq(1000,200,-100), seq(100,10,-10))

t0=Sys.time()

no_cores=min(detectCores(), nchunks)
print("Number of cores:")
print(no_cores)
cl <- makeCluster(no_cores, type="FORK")
zero=parSapply(cl=cl, 1:nchunks, FUN=function(k){
  out = GraphicalModel(xdata = X, PFER_thr = n[k], start = start.type)
  return(sum(rowSums(Adjacency(out))==0))
})
stopCluster(cl)

t1=Sys.time()
print(t1-t0)

zero=unlist(zero)

pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_PFER_thres.pdf"))
par(mar = c(5, 5, 1, 1))
plot(1:length(zero), zero,
     ylab = "Number of nodes with zero edges",
     xlab = "PFER threshold",
     xaxt = "n",
     cex.lab = 1.5,
     pch = 19, type = "b", col = batch.colours[m])
axis(side = 1, at = 1:length(n), labels = n)
dev.off()

if(sum(zero != min(zero))!=0){
  start = c(Inf, seq(1000,200,-100), seq(100,10,-10))[which(zero == min(zero))[length(which(zero == min(zero)))]]
  if(start != 10){
    stop = c(Inf, seq(1000,200,-100), seq(100,10,-10))[which(zero == min(zero))[length(which(zero == min(zero)))]+1]
  }

  if(is.infinite(start)){
    n = seq(20000,1000,-1000)
    
    t0=Sys.time()
    
    no_cores=min(detectCores(), nchunks)
    print("Number of cores:")
    print(no_cores)
    cl <- makeCluster(no_cores, type="FORK")
    zero=parSapply(cl=cl, 1:nchunks, FUN=function(k){
      out = GraphicalModel(xdata = X, PFER_thr = n[k], start = start.type)
      return(sum(rowSums(Adjacency(out))==0))
    })
    stopCluster(cl)
    
    t1=Sys.time()
    print(t1-t0)
    
    zero=unlist(zero)
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_PFER_thres_broad.pdf"))
    par(mar = c(5, 5, 1, 1))
    plot(1:length(zero), zero,
         ylab = "Number of nodes with zero edges",
         xlab = "PFER threshold",
         xaxt = "n",
         cex.lab = 1.5,
         pch = 19, type = "b", col = batch.colours[m])
    axis(side = 1, at = 1:length(n), labels = n)
    dev.off()
  } else if (start <=100){
    print(paste0("Set threshold to ", stop))
    } else{
      n = seq(start,stop,-10)
      
      t0=Sys.time()
      
      no_cores=min(detectCores(), nchunks)
      zero=sapply(1:length(n),FUN=function(k){
        out = GraphicalModel(xdata = X, PFER_thr = n[k], start = start.type)
        return(sum(rowSums(Adjacency(out))==0))
      })
      t1=Sys.time()
      print(t1-t0)
      
      zero=unlist(zero)
      pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_PFER_thres_fine.pdf"))
      par(mar = c(5, 5, 1, 1))
      plot(1:length(zero), zero,
           ylab = "Number of nodes with zero edges",
           xlab = "PFER threshold",
           xaxt = "n",
           cex.lab = 1.5,
           pch = 19, type = "b", col = batch.colours[m])
      axis(side = 1, at = 1:length(n), labels = n)
      dev.off()
    }
} else {
  print("Set threshold to 10")
}
