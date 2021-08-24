## Univariate clustering analysis
## Rin Wada 12 July

# Load packages
library(ape)
library(igraph)
library(colorspace)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

# ## Parameters
# args=commandArgs(trailingOnly=TRUE)
# m=as.numeric(args[1])

# Load data sets
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

suffix = c("lux","fra","gs","pooled3","pooled2")

m = 4
# Load data
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
print(all(rownames(expo)==rownames(covars)))

families=unique(covars$Family.ID)

withinf_sp=matrix(NA,nrow=length(families),ncol=ncol(expo))
null_sp=matrix(NA,nrow=(nrow(expo)*(nrow(expo)-1))/2,ncol=ncol(expo))

ifelse(dir.exists(paste0("../Figures/",filepaths[m],"/Univar_clustering")),"",
       dir.create(paste0("../Figures/",filepaths[m],"/Univar_clustering")))
for (p in 1:ncol(expo)){
  d=dist(expo[,p])
  h=hclust(d, method = "complete")
  # print(all(covars$Indiv.ID==h$labels))
  h$labels=paste0(covars$Family.ID, "-",h$labels)
  myphylo=as.phylo(h)
  
  {pdf(paste0("../Figures/",filepaths[m],"/Univar_clustering/Hierarchical_cont_expo_",p,".pdf"),width=14)
    par(mar=c(0,0,3,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=family.colours[as.character(covars$Family.ID)], main=colnames(expo)[p])
    dev.off()}

  {pdf(paste0("../Figures/",filepaths[m],"/Univar_clustering/Hierarchical_cont_graph_expo_",p,".pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[levels(covars$Family.ID)])
    dev.off()}

  g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[levels(covars$Family.ID)],
                      verbose = FALSE)
  
  nobs=nrow(expo)
  sp=shortest.paths(g)[1:nobs,1:nobs]
  rownames(sp)=colnames(sp)=myphylo$tip.label
  
  fsp=NULL
  for (f in 1:length(families)){
    tmpmat=sp[covars$Family.ID==families[f],covars$Family.ID==families[f]]
    fsp=c(fsp,mean(tmpmat[upper.tri(tmpmat)]))
  }
  withinf_sp[,p]=fsp
  
  null_sp[,p]=sp[upper.tri(sp)]
}
colnames(withinf_sp)=colnames(expo)
rownames(withinf_sp)=families

for (p in 1:ncol(withinf_sp)){
  {pdf(paste0("../Figures/",filepaths[m],"/Univar_clustering/Shortest_path_cont_",p,".pdf"),width=14)
    par(mar=c(5,5,1,1))
    plot(withinf_sp[,p], pch=19, col=family.colours[rownames(withinf_sp)], xaxt="n", las=2,
         panel.first=abline(v=1:nrow(withinf_sp), lty=3, col="grey"),
         xlab="Family ID", cex=2,
         ylab="Shortest path length between siblings", cex.lab=1.5)
    axis(side=1, at=1:nrow(withinf_sp), labels=families, las = 2)
    dev.off()}
}

x=as.vector(row(withinf_sp))
y=as.vector(withinf_sp)
ybis=as.vector(null_sp)
z=as.vector(col(withinf_sp))*2
zbis=as.vector(col(null_sp))*2-1
y=c(y,ybis)
z=c(z,zbis)

annot_sub = annot[colnames(expo)]
ifelse(dir.exists("../Figures/Section3"),"",dir.create("../Figures/Section3"))
{pdf(paste0("../Figures/Section3/Boxplot_univariate_cont_shortest_path_",suffix[m],".pdf"), width=12, height=8)
  par(mar=c(20,5,1,1))
  mycol=rep(annot.colours[annot_sub], each=2)
  mycol[seq(1,ncol(expo)*2,by=2)]="grey80"
  boxplot(y~z,pch=19, cex=0.5, col=mycol, border=mycol,
          las=1, xaxt="n", outline=FALSE,
          medcol='white', whiskcol='black', staplecol='black', outcol='black',
          xlab="", ylab="Shortest path", cex.lab=1.5)
  for (k in 1:ncol(withinf_sp)){
    axis(side=1, at=2*k-0.5, labels=colnames(withinf_sp)[k], cex.axis=0.8, las=2)
  }
  xgroup=2*c(which(!duplicated(annot_sub))-0.75, length(unique(z))/2+0.25)
  abline(v=xgroup, lty=3)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
         col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis=0.9)
  }
  dev.off()
}

# 90% interval
q = apply(withinf_sp, 2, function(x) quantile(x,c(0.05,0.95)))
null_q = apply(null_sp, 2, function(x) quantile(x,c(0.05,0.95)))
if(all((q[1,] <= null_q[2]) & (q[2,] >= null_q[1]))){
  print("All compounds overlap with null distribution")
}

## Overall distribution of the shortest paths by family ID

withinf_sp=t(withinf_sp)
# ids=sort.list(apply(withinf_sp,2,median))
# withinf_sp=withinf_sp[,ids]

x=as.vector(row(withinf_sp))
y=as.vector(withinf_sp)
z=as.vector(col(withinf_sp))
ybis=as.vector(null_sp)
zbis=rep(0,length(ybis))
y=c(y,ybis)
z=c(z,zbis)

{pdf(paste0("../Figures/Section3/Boxplot_univariate_cont_shortest_path_by_family_",suffix[m],".pdf"), width=14)
  par(mar=c(5,5,1,1))
  boxplot(y~z,pch=19, cex=0.5, col=c("grey80", family.colours[colnames(withinf_sp)]),
          border=c("grey80", family.colours[colnames(withinf_sp)]), #outline=FALSE,
          las=1, xaxt="n", outline=FALSE,
          medcol='white', whiskcol='black', staplecol='black', outcol='black',
          xlab="Family ID", ylab="Shortest path", cex.lab=1.5)
  abline(h=median(ybis), lty=3)
  axis(side=1, at=1, labels="Ref", cex.axis=0.8, las = 2)
  for (k in 1:ncol(withinf_sp)){
    axis(side=1, at=k+1, labels=colnames(withinf_sp)[k], cex.axis=0.8, las = 2)
  }
  dev.off()
}

# 90% interval
q = apply(withinf_sp, 2, function(x) quantile(x,c(0.05,0.95)))
null_q = quantile(null_sp,c(0.05,0.95))
if(all((q[1,] <= null_q[2]) & (q[2,] >= null_q[1]))){
  print("All families overlap with null distribution")
}

for (m in c(1:3,5)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  families=unique(covars$Family.ID)
  
  withinf_sp=matrix(NA,nrow=length(families),ncol=ncol(expo))
  null_sp=matrix(NA,nrow=(nrow(expo)*(nrow(expo)-1))/2,ncol=ncol(expo))
  
  ifelse(dir.exists(paste0("../Figures/",filepaths[m],"/Univar_clustering")),"",
         dir.create(paste0("../Figures/",filepaths[m],"/Univar_clustering")))
  for (p in 1:ncol(expo)){
    d=dist(expo[,p])
    h=hclust(d, method = "complete")
    # print(all(covars$Indiv.ID==h$labels))
    h$labels=paste0(covars$Family.ID, "-",h$labels)
    myphylo=as.phylo(h)
    
    {pdf(paste0("../Figures/",filepaths[m],"/Univar_clustering/Hierarchical_cont_expo_",p,".pdf"),width=14)
      par(mar=c(0,0,3,0))
      plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
           tip.color=family.colours[as.character(covars$Family.ID)], main=colnames(expo)[p])
      dev.off()}
    
    {pdf(paste0("../Figures/",filepaths[m],"/Univar_clustering/Hierarchical_cont_graph_expo_",p,".pdf"))
      g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[levels(covars$Family.ID)])
      dev.off()}
    
    g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[levels(covars$Family.ID)],
                        verbose = FALSE)
    
    nobs=nrow(expo)
    sp=shortest.paths(g)[1:nobs,1:nobs]
    rownames(sp)=colnames(sp)=myphylo$tip.label
    
    fsp=NULL
    for (f in 1:length(families)){
      tmpmat=sp[covars$Family.ID==families[f],covars$Family.ID==families[f]]
      fsp=c(fsp,mean(tmpmat[upper.tri(tmpmat)]))
    }
    withinf_sp[,p]=fsp
    
    null_sp[,p]=sp[upper.tri(sp)]
  }
  colnames(withinf_sp)=colnames(expo)
  rownames(withinf_sp)=families
  
  for (p in 1:ncol(withinf_sp)){
    {pdf(paste0("../Figures/",filepaths[m],"/Univar_clustering/Shortest_path_cont_",p,".pdf"),width=14)
      par(mar=c(5,5,1,1))
      plot(withinf_sp[,p], pch=19, col=family.colours[rownames(withinf_sp)], xaxt="n", las=2,
           panel.first=abline(v=1:nrow(withinf_sp), lty=3, col="grey"),
           xlab="Family ID", cex=2,
           ylab="Shortest path length between siblings", cex.lab=1.5)
      axis(side=1, at=1:nrow(withinf_sp), labels=families, las = 2)
      dev.off()}
  }
  
  x=as.vector(row(withinf_sp))
  y=as.vector(withinf_sp)
  ybis=as.vector(null_sp)
  z=as.vector(col(withinf_sp))*2
  zbis=as.vector(col(null_sp))*2-1
  y=c(y,ybis)
  z=c(z,zbis)
  
  annot_sub = annot[colnames(expo)]
  ifelse(dir.exists("../Figures/Supplementary"),"",dir.create("../Figures/Supplementary"))
  {pdf(paste0("../Figures/Supplementary/Boxplot_univariate_cont_shortest_path_",suffix[m],".pdf"), width=12, height=8)
    par(mar=c(20,5,1,1))
    mycol=rep(annot.colours[annot_sub], each=2)
    mycol[seq(1,ncol(expo)*2,by=2)]="grey80"
    boxplot(y~z,pch=19, cex=0.5, col=mycol, border=mycol,
            las=1, xaxt="n", outline=FALSE,
            medcol='white', whiskcol='black', staplecol='black', outcol='black',
            xlab="", ylab="Shortest path", cex.lab=1.5)
    for (k in 1:ncol(withinf_sp)){
      axis(side=1, at=2*k-0.5, labels=colnames(withinf_sp)[k], cex.axis=0.8, las=2)
    }
    xgroup=2*c(which(!duplicated(annot_sub))-0.75, length(unique(z))/2+0.25)
    abline(v=xgroup, lty=3)
    axis(side=1, line=8, at=xgroup, labels=NA)
    tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
    for (k in 1:length(unique(annot_sub))){
      axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2,
           col.axis=darken(annot.colours[unique(annot_sub)[k]], amount=0.5),
           cex.axis=0.9)
    }
    dev.off()
  }
  
  # 90% interval
  q = apply(withinf_sp, 2, function(x) quantile(x,c(0.05,0.95)))
  null_q = apply(null_sp, 2, function(x) quantile(x,c(0.05,0.95)))
  if(all((q[1,] <= null_q[2]) & (q[2,] >= null_q[1]))){
    print("All compounds overlap with null distribution")
  }
  
  ## Overall distribution of the shortest paths by family ID
  
  withinf_sp=t(withinf_sp)
  # ids=sort.list(apply(withinf_sp,2,median))
  # withinf_sp=withinf_sp[,ids]
  
  x=as.vector(row(withinf_sp))
  y=as.vector(withinf_sp)
  z=as.vector(col(withinf_sp))
  ybis=as.vector(null_sp)
  zbis=rep(0,length(ybis))
  y=c(y,ybis)
  z=c(z,zbis)
  
  {pdf(paste0("../Figures/Supplementary/Boxplot_univariate_cont_shortest_path_by_family_",suffix[m],".pdf"), width=14)
    par(mar=c(5,5,1,1))
    boxplot(y~z,pch=19, cex=0.5, col=c("grey80", family.colours[colnames(withinf_sp)]),
            border=c("grey80", family.colours[colnames(withinf_sp)]), #outline=FALSE,
            las=1, xaxt="n", outline=FALSE,
            medcol='white', whiskcol='black', staplecol='black', outcol='black',
            xlab="Family ID", ylab="Shortest path", cex.lab=1.5)
    abline(h=median(ybis), lty=3)
    axis(side=1, at=1, labels="Ref", cex.axis=0.8, las = 2)
    for (k in 1:ncol(withinf_sp)){
      axis(side=1, at=k+1, labels=colnames(withinf_sp)[k], cex.axis=0.8, las = 2)
    }
    dev.off()
  }
  
  # 90% interval
  q = apply(withinf_sp, 2, function(x) quantile(x,c(0.05,0.95)))
  null_q = quantile(null_sp,c(0.05,0.95))
  if(all((q[1,] <= null_q[2]) & (q[2,] >= null_q[1]))){
    print("All families overlap with null distribution")
  }
}

# print(which(colnames(expo)=="Fipronil"))
