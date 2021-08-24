## Multivariate clustering
## Rin Wada 12 July

# Load packages
library(ape)
library(igraph)
library(focus)
library(colorspace)
library(RColorBrewer)

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

families=levels(covars$Family.ID)

expo = scale(expo)
d=dist(expo)
h=hclust(d, method = "complete")
print(all(covars$Indiv.ID==h$labels))
h$labels=as.character(covars$Family.ID)
myphylo=as.phylo(h)

{pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all.pdf"),width=14)
  par(mar=c(0,0,0,0))
  plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
       tip.color=family.colours[as.character(covars$Family.ID)])
  dev.off()
}

{pdf(paste0("../Figures/Section3/Hierarchical_cont_graph_expo_all_",suffix[m],".pdf"))
  g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[as.character(covars$Family.ID)])
  dev.off()
}

nobs=nrow(expo)
sp=shortest.paths(g)[1:nobs,1:nobs]
rownames(sp)=colnames(sp)=myphylo$tip.label

fsp=NULL
for (f in 1:length(families)){
  tmpmat=sp[covars$Family.ID==families[f],covars$Family.ID==families[f]]
  fsp=c(fsp,mean(tmpmat[upper.tri(tmpmat)]))
}
names(fsp)=families

# Save shortest path per family
saveRDS(fsp, paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))

# {pdf(paste0("../Figures/Section3/Shortest_path_cont_all_",suffix[m],".pdf"),width=14,height=6)
#   par(mar=c(5,5,1,1))
#   plot(fsp, pch=19, col=family.colours[families], xaxt="n", las=1,
#        panel.first=abline(v=1:length(families), lty=3, col="grey"),
#        xlab="Family ID", cex=2,
#        ylab="Shortest path length between siblings", cex.lab=1.5)
#   axis(side=1, at=1:length(families), labels=families, las = 2)
#   dev.off()
# }

families=unique(covars$Region)
  mycolours = region.colours
  names(mycolours)=families

  {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_region.pdf"),width=14)
    par(mar=c(0,0,0,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=region.colours[as.character(covars$Region)])
    dev.off()
  }

  {pdf(paste0("../Figures/Section3/Hierarchical_cont_graph_expo_all_region_",suffix[m],".pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = region.colours[as.character(covars$Region)])
    legend("topright", pch=19, col=region.colours[levels(covars$Region)],
           legend=levels(covars$Region), pt.cex = 2, cex = 0.8, ncol = 2, bty = "n")
    dev.off()
  }

for (m in c(1:3,5)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  families=levels(covars$Family.ID)
  
  expo = scale(expo)
  d=dist(expo)
  h=hclust(d, method = "complete")
  print(all(covars$Indiv.ID==h$labels))
  h$labels=as.character(covars$Family.ID)
  myphylo=as.phylo(h)
  
  {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all.pdf"),width=14)
    par(mar=c(0,0,0,0))
    plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
         tip.color=family.colours[as.character(covars$Family.ID)])
    dev.off()
  }
  
  {pdf(paste0("../Figures/Supplementary/Hierarchical_cont_graph_expo_all_",suffix[m],".pdf"))
    g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = family.colours[as.character(covars$Family.ID)])
    dev.off()
  }
  
  nobs=nrow(expo)
  sp=shortest.paths(g)[1:nobs,1:nobs]
  rownames(sp)=colnames(sp)=myphylo$tip.label
  
  fsp=NULL
  for (f in 1:length(families)){
    tmpmat=sp[covars$Family.ID==families[f],covars$Family.ID==families[f]]
    fsp=c(fsp,mean(tmpmat[upper.tri(tmpmat)]))
  }
  names(fsp)=families
  
  # Save shortest path per family
  saveRDS(fsp, paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))
  
  # if(m==5){
  #   pdf(paste0("../Figures/Supplementary/Shortest_path_cont_all_",suffix[m],".pdf"),width=14,height=6)
  # }
  # pdf(paste0("../Figures/",filepaths[m],"/Shortest_path_cont_all.pdf"),width=14,height=6)
  #   par(mar=c(5,5,1,1))
  #   plot(fsp, pch=19, col=family.colours[families], xaxt="n", las=1,
  #        panel.first=abline(v=1:length(families), lty=3, col="grey"),
  #        xlab="Family ID", cex=2,
  #        ylab="Shortest path length between siblings", cex.lab=1.5)
  #   axis(side=1, at=1:length(families), labels=families, las = 2)
  #   dev.off()
  
  if(m==2){
    families=unique(covars$Region)
    mycolours = region.colours
    names(mycolours)=families
    
    {pdf(paste0("../Figures/",filepaths[m],"/Hierarchical_cont_expo_all_region.pdf"),width=14)
      par(mar=c(0,0,0,0))
      plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
           tip.color=region.colours[as.character(covars$Region)])
      dev.off()
    }
    
    {pdf(paste0("../Figures/Supplementary/Hierarchical_cont_graph_expo_all_region_",suffix[m],".pdf"))
      g=ClusteringToGraph(covars=covars, myphylo=myphylo, mycol = region.colours[as.character(covars$Region)])
      legend("topright", pch=19, col=region.colours[levels(covars$Region)],
             legend=levels(covars$Region), cex = 0.8, ncol = 1, bty = "n")
      dev.off()
    }
  }
}

# ### Combined ----
# 
# for(m in 1:3){
#   fsp = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))
#   assign(paste0("fsp_",suffix[m]),fsp)
# }
# 
# 
# {pdf(paste0("../Figures/Section3/Shortest_path_cont_all_combined.pdf"),width=10,height=4.5)
#   par(oma=c(7,5,1,0), mar=c(0,0,0,1), xpd = TRUE)
#   layout(matrix(c(1,2,3),1,3, byrow = TRUE), width = c(18,35,13))
#   plot(fsp_lux, pch=19, col=family.colours[names(fsp_lux)], xaxt="n", las=1,
#        panel.first=abline(v=1:length(fsp_lux), lty=3, col="grey"),
#        xlab="", cex=2,
#        ylab="", cex.lab=1.5)
#   axis(side=1, at=1:length(fsp_lux), labels=names(fsp_lux), las = 2)
#   mtext("Family ID", side = 1, line = 4.5, outer = TRUE)
#   mtext("Shortest path length between siblings", side = 2, line = 3, outer = TRUE, las = 3)
#   mtext(batches[1:3], side = 1, line = 3, at = c(0.12,0.5,0.9), adj = c(0.5,0),outer = TRUE)
#   plot(fsp_fra, pch=19, col=family.colours[names(fsp_fra)], xaxt="n", las=1,
#        panel.first=abline(v=1:length(fsp_fra), lty=3, col="grey"),
#        yaxt = "n", xlab="", cex=2, ylab="", cex.lab=1.5)
#   axis(side=1, at=1:length(fsp_fra), labels=names(fsp_fra), las = 2)
#   plot(fsp_gs, pch=19, col=family.colours[names(fsp_gs)], xaxt="n", las=1,
#        panel.first=abline(v=1:length(fsp_gs), lty=3, col="grey"),
#        yaxt = "n",
#        xlab="", cex=2,
#        ylab="", cex.lab=1.5)
#   axis(side=1, at=1:length(fsp_gs), labels=names(fsp_gs), las = 2)
#   dev.off()
# }

### Compare ----
covars = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh_no_isolated.rds"))
  
for(m in 1:4){
    fsp = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))
    assign(paste0("fsp_",suffix[m]),fsp)
    }

  {pdf(paste0("../Figures/Section3/Shortest_path_cont_all_compare.pdf"),width=7,height=7)
    par(mar=c(2,5,3,6), mgp=c(4,1,0), pty="s")
    plot(fsp_pooled3, c(fsp_lux, fsp_fra, fsp_gs), pch=19,
         yaxt = "n", xaxt = "n",
         ylim = c(2,max(c(fsp_lux, fsp_fra, fsp_gs, fsp_pooled3))),
         xlim = c(2,max(c(fsp_lux, fsp_fra, fsp_gs, fsp_pooled3))),
         col=family.colours[names(fsp_pooled3)], las=1,
         panel.first=list(abline(v=seq(2,20,4), lty=3, col="grey"),
                          abline(h=seq(2,20,4), lty=3, col="grey"),
                          abline(0,1, lty=2, col="grey")),
         xlab="", ylab = "",
         cex=3, cex.lab=1.5)
    axis(4, at = seq(2,20,4))
    mtext("Shortest path length between siblings \n(By cohort)", 4, 4, cex = 1.5)
    axis(3, at = seq(2,20,4))
    mtext("Shortest path length between siblings \n(Pooled)", 3, 2, cex = 1.5)
    par(xpd = TRUE)
    jitter = seq(-1,25,length.out = sum(fsp_pooled3<5 & c(fsp_lux, fsp_fra, fsp_gs) < 5))
    names(jitter) = names(fsp_pooled3)[fsp_pooled3<5 & c(fsp_lux, fsp_fra, fsp_gs) < 5][order(fsp_pooled3[fsp_pooled3<5 & c(fsp_lux, fsp_fra, fsp_gs) < 5])]
    
    for(k in 1:length(fsp_pooled3)){
      if(fsp_pooled3[k]<5 & c(fsp_lux, fsp_fra, fsp_gs)[k]<5){
        y = jitter[names(fsp_pooled3)[k]]
        x = -1
        segments(fsp_pooled3[k], c(fsp_lux, fsp_fra, fsp_gs)[k], x,  y, col = family.colours[names(fsp_pooled3)][k])
        points(x, y,
               col = family.colours[names(fsp_pooled3)][k],
               pch = 19,
               cex = 3)
        text(x, y, labels = names(fsp_pooled3)[k], cex = 0.8)
      } else{
        text(fsp_pooled3[k], c(fsp_lux, fsp_fra, fsp_gs)[k], labels = names(fsp_pooled3)[k], cex=0.8)
      }
    }

    dev.off()
  }
  
  mycolours = region.colours[as.character(covars$Region[!duplicated(covars$Family.ID)])]
  {pdf(paste0("../Figures/Section3/Shortest_path_cont_all_compare_region.pdf"),width=9,height=7)
    par(mar=c(2,5,3,16), mgp=c(4,1,0), pty="s")
    plot(fsp_pooled3, c(fsp_lux, fsp_fra, fsp_gs), pch=19,
         yaxt = "n", xaxt = "n",
         ylim = c(2,max(c(fsp_lux, fsp_fra, fsp_gs, fsp_pooled3))),
         xlim = c(2,max(c(fsp_lux, fsp_fra, fsp_gs, fsp_pooled3))),
         col=mycolours, las=1,
         panel.first=list(abline(v=seq(2,20,4), lty=3, col="grey"),
                          abline(h=seq(2,20,4), lty=3, col="grey"),
                          abline(0,1, lty=2, col="grey")),
         xlab="", ylab = "",
         cex=3, cex.lab=1.5)
    axis(4, at = seq(2,20,4))
    mtext("Shortest path length between siblings \n(By cohort)", 4, 4, cex = 1.5)
    axis(3, at = seq(2,20,4))
    mtext("Shortest path length between siblings \n(Pooled)", 3, 2, cex = 1.5)
    par(xpd = TRUE)
    jitter = seq(-1,25,length.out = sum(fsp_pooled3<5 & c(fsp_lux, fsp_fra, fsp_gs) < 5))
    names(jitter) = names(fsp_pooled3)[fsp_pooled3<5 & c(fsp_lux, fsp_fra, fsp_gs) < 5][order(fsp_pooled3[fsp_pooled3<5 & c(fsp_lux, fsp_fra, fsp_gs) < 5])]
    
    for(k in 1:length(fsp_pooled3)){
      if(fsp_pooled3[k]<5 & c(fsp_lux, fsp_fra, fsp_gs)[k]<5){
        y = jitter[names(fsp_pooled3)[k]]
        x = -1
        segments(fsp_pooled3[k], c(fsp_lux, fsp_fra, fsp_gs)[k], x,  y, col = mycolours[k])
        points(x, y,
               col = mycolours[k],
               pch = 19,
               cex = 3)
        text(x, y, labels = names(fsp_pooled3)[k], cex = 0.8)
      } else{
        text(fsp_pooled3[k], c(fsp_lux, fsp_fra, fsp_gs)[k], labels = names(fsp_pooled3)[k], cex=0.8)
      }
    }
    coord = par("usr")
    legend(x = coord[2]+4, y = coord[4],
           legend = levels(covars$Region),
           pch = 19,
           pt.cex = 2,
           col = region.colours[levels(covars$Region)],
           bty = "n")
    dev.off()
  }
  
  ### Compare (Sensitivity)----
  for(m in c(1,3,5)){
    fsp = readRDS(paste0("../Results/",filepaths[m],"/Shortest_path_cont_all.rds"))
    assign(paste0("fsp_",suffix[m]),fsp)
  }
  
  {pdf(paste0("../Figures/Supplementary/Shortest_path_cont_all_compare_pooled2.pdf"),width=7,height=6.5)
    par(mar=c(1,5,3,6), mgp=c(4,1,0), pty="s")
    plot(fsp_pooled2, c(fsp_lux, fsp_gs), pch=19,
         yaxt = "n", xaxt = "n",
         ylim = c(2,max(c(fsp_lux, fsp_gs, fsp_pooled2))),
         xlim = c(2,max(c(fsp_lux, fsp_gs, fsp_pooled2))),
         col=family.colours[names(fsp_pooled2)], las=1,
         panel.first=list(abline(v=seq(2,15,4), lty=3, col="grey"),
                          abline(h=seq(2,15,4), lty=3, col="grey"),
                          abline(0,1, lty=2, col="grey")),
         xlab="", ylab = "",
         cex=3, cex.lab=1.5)
    axis(4, at = seq(2,15,4))
    mtext("Shortest path length between siblings \n(By cohort)", 4, 4, cex = 1.5)
    axis(3, at = seq(2,15,4))
    mtext("Shortest path length between siblings \n(Pooled)", 3, 2, cex = 1.5)
    par(xpd = TRUE)
    jitter = seq(2,16,length.out = sum(fsp_pooled2<5 & c(fsp_lux, fsp_gs) < 5))
    names(jitter) = names(fsp_pooled2)[fsp_pooled2<5 & c(fsp_lux, fsp_gs) < 5][order(fsp_pooled2[fsp_pooled2<5 & c(fsp_lux, fsp_gs) < 5])]
    
    for(k in 1:length(fsp_pooled2)){
      if(fsp_pooled2[k]<5 & c(fsp_lux, fsp_gs)[k]<5){
        y = jitter[names(fsp_pooled2)[k]]
        x = 0
        segments(fsp_pooled2[k], c(fsp_lux, fsp_gs)[k], x,  y, col = family.colours[names(fsp_pooled2)][k])
        points(x, y,
               col = family.colours[names(fsp_pooled2)][k],
               pch = 19,
               cex = 3)
        text(x, y, labels = names(fsp_pooled2)[k])
      } else{
        text(fsp_pooled2[k], c(fsp_lux, fsp_gs)[k], labels = names(fsp_pooled2)[k], cex=0.8)
      }
    }
    
    dev.off()
  }  
  
  