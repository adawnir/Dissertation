## Univariate clustering analysis
## Rin Wada 12 July

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ape)
library(igraph)

### Univariate clustering ----

# p=1
# delta_mat=matrix(NA, nrow=nrow(expo), ncol=nrow(expo))
# for (i in 1:nrow(expo)){
#   for (j in 1:nrow(expo)){
#     delta=expo[i,p]-expo[j,p]
#     delta_mat[i,j]=delta
#   }
# }

families=unique(covars$Family.ID)

mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families

withinf_sp=matrix(NA,nrow=length(families),ncol=ncol(expo))
null_sp=matrix(NA,nrow=(nrow(expo)*(nrow(expo)-1))/2,ncol=ncol(expo))

for (p in 1:ncol(expo)){
  d=dist(expo[,p])
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))
  h$labels=paste0(covars$Family.ID, "-",h$labels)
  myphylo=as.phylo(h)
  
  {pdf(paste0("../Figures/",filepath,"/Hierarchical_cont_expo_",p,".pdf"),width=14)
    par(mar=c(0,0,3,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=mycolours[covars$Family.ID], main=colnames(expo)[p])
    dev.off()}
  
  {pdf(paste0("../Figures/",filepath,"/Hierarchical_cont_graph_expo_",p,".pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo)
    dev.off()}
  
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
  {pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_",p,".pdf"),width=14)
    par(mar=c(5,5,1,1))
    plot(withinf_sp[,p], pch=19, col=mycolours[rownames(withinf_sp)], xaxt="n", las=1,
         panel.first=abline(v=1:ncol(withinf_sp), lty=3, col="grey"), 
         xlab="Family ID", cex=2,
         ylab="Shortest path length between siblings", cex.lab=1.5)
    axis(side=1, at=1:nrow(withinf_sp), labels=families, las = 2)
    dev.off()}
}

which(colnames(expo)=="Fipronil")

## Overall distribution of the shortest paths by pollutant

# ids=sort.list(apply(withinf_sp,2,median))
# withinf_sp=withinf_sp[,ids]

# x=as.vector(row(withinf_sp))
# y=as.vector(withinf_sp)
# z=as.vector(col(withinf_sp))

x=as.vector(row(withinf_sp))
y=as.vector(withinf_sp)
ybis=as.vector(null_sp)
z=as.vector(col(withinf_sp))*2
zbis=as.vector(col(null_sp))*2-1
y=c(y,ybis)
z=c(z,zbis)

mycolours = family.colours[unique(annot)]
annot_sub = annot[colnames(expo)]

{pdf(paste0("../Figures/",filepath,"/Boxplot_univariate_cont_shortest_path.pdf"), width=12, height=8)
  par(mar=c(20,5,1,1))
  mycol=rep(mycolours[as.character(annot_sub)], each=2)
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
         col.axis=darken(mycolours[as.character(unique(annot_sub)[k])], amount=0.5),
         cex.axis=0.9)
  }
  dev.off()
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

mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families

{pdf(paste0("../Figures/",filepath,"/Boxplot_univariate_cont_shortest_path_by_family.pdf"), width=14)
  par(mar=c(5,5,1,1))
  boxplot(y~z,pch=19, cex=0.5, col=c("grey80", mycolours[colnames(withinf_sp)]), 
          border=c("grey80", mycolours[colnames(withinf_sp)]), #outline=FALSE, 
          las=1, xaxt="n", outline=FALSE,
          medcol='white', whiskcol='black', staplecol='black', outcol='black',
          xlab="Family ID", ylab="Shortest path", cex.lab=1.5)
  abline(h=median(ybis), lty=3)
  axis(side=1, at=1, labels="Ref", cex.axis=0.8)
  for (k in 1:ncol(withinf_sp)){
    axis(side=1, at=k+1, labels=colnames(withinf_sp)[k], cex.axis=0.8)
  }
  dev.off()}

