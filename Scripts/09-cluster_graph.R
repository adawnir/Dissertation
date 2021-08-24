## Clustering visualisations
## 11 Aug

### Consensus clustering ----
# Load packages
# library(sgPLS)
# library(plotrix)
# library(ellipse)
library(igraph)
library(focus)
library(colorspace)

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")
for (m in 1:length(batches)){
  out = readRDS(paste0("../Results/",filepaths[m],"/Consensus_clustering_output.rds"))
  pdf(paste0("../Figures/",filepaths[m],"/Consensus_clustering_output.pdf")) 
  par(mar = c(7, 5, 7, 6))
  CalibrationPlot(out)
  dev.off()
}

perf = NULL
ifelse(dir.exists("../Figures/Section4"),"",dir.create("../Figures/Section4"))
for (m in 1:length(batches)){
  # Load data
sc = readRDS(paste0("../Results/",filepaths[m],"/Stable_clusters.rds"))
stab = readRDS(paste0("../Results/",filepaths[m],"/Consensus_clustering_output.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))

k = max(sc)

perf = rbind(perf, cbind(k,ClusteringPerformance(sc, covars$Family.ID)))

# families=unique(sc)[table(sc)!=1]
# mycolours=brewer.pal(n=12,name='Paired')
# mycolours=colorRampPalette(mycolours)(length(families))
# names(mycolours)=families
# 
# tmp = rep(alpha("grey80",0.3),length(setdiff(unique(sc), families)))
# names(tmp) = setdiff(unique(sc), families)
# 
# mycolours = c(mycolours, tmp)
# mycolours = mycolours[as.character(unique(sc))]
#   
# myplsda = plsda(expo, as.factor(sc), ncomp = 3)
# 
# CreateScorePlot.plsda(myplsda=myplsda, type1=sc, type2=covars$Family.ID,
#                       mycolours=mycolours,
#                       filename=paste0("../Figures/",filepaths[m],"/Consensus_clustering_PLSDA_score_plot.pdf"))

myadjacency = Adjacency(stab)
mynode_colour = family.colours[as.character(covars$Family.ID)]
names(mynode_colour) = rownames(covars)

communities = lapply(unique(sc), function(i) rownames(covars)[sc==i])
names(communities) = unique(sc)

cluster_colour = rainbow(length(unique(sc)), alpha = 0.1)[sc]
names(cluster_colour) = sc

mygraph = Graph(myadjacency, node_colour = mynode_colour,
                node_label = as.character(covars$Family.ID),
                satellites = TRUE)

E(mygraph)$color = darken(cluster_colour[as.character(sc[ends(mygraph, E(mygraph))[,1]])],0.5)  

E(mygraph)$width = 1
V(mygraph)$size = 10
V(mygraph)$label.cex = 0.8

# Plot clusters in graph
if (m ==4){
  pdf(paste0("../Figures/Section4/Consensus_clustering_graph_",suffix[m],".pdf"))
} else {pdf(paste0("../Figures/Supplementary/Consensus_clustering_graph_",suffix[m],".pdf"))}
par(mar=c(0,0,0,0))
set.seed(1)
plot(mygraph,
     vertex.frame.color = darken(cluster_colour,0.7),
     mark.groups = communities,
     mark.col = cluster_colour[as.character(1:length(communities))],
     mark.border = darken(cluster_colour[as.character(1:length(communities))], 0.5))
dev.off()

if (m %in% c(2,4)){
  mynode_colour = region.colours[as.character(covars$Region)]
  names(mynode_colour) = rownames(covars)
  if (m==4){
    pdf(paste0("../Figures/Section4/Consensus_clustering_graph_region_",suffix[m],".pdf"))
  } else {
    pdf(paste0("../Figures/Supplementary/Consensus_clustering_graph_region_",suffix[m],".pdf"))
    
  }
  # Plot clusters in graph
  par(mar=c(0,0,5,0), xpd = TRUE)
  set.seed(1)
  plot(mygraph,
       vertex.color = mynode_colour,
       vertex.frame.color = darken(cluster_colour,0.7),
       mark.groups = communities,
       mark.col = cluster_colour[as.character(1:length(communities))],
       mark.border = darken(cluster_colour[as.character(1:length(communities))], 0.5))
  legend("top",pch=19, inset = c(0,-0.15),
         col = region.colours[levels(covars$Region)],
         legend = levels(covars$Region),
         ncol = ceiling(length(levels(covars$Region))/3), bty = "n")
  dev.off()
}
}
saveRDS(perf, "../Results/Consensus_clustering_performance.rds")

### Community detection ----
suffix = c("lux","fra","gs","pooled3","pooled2")

perf = NULL
ifelse(dir.exists("../Figures/Section5"),"",dir.create("../Figures/Section5"))
for (m in 1:length(batches)){
  A = readRDS(paste0("../Results/",filepaths[m],"/Graphical_network_children_adjacency_matrix.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  
  mynode_colour = family.colours[as.character(covars$Family.ID)]
  names(mynode_colour) = rownames(covars)
  
  # Getting igraph object
  mygraph=Graph(adjacency=A,
                node_label = as.character(covars$Family.ID),
                node_colour=mynode_colour,
                satellites = TRUE)
  
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children.pdf"))
  # par(mar=c(0,0,5,0), xpd = TRUE)
  par(mar=c(0,0,0,0))
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       vertex.label.cex = 0.7,
       edge.width = 1,
       vertex.size = 9)
  # legend("top",pch=19, inset = c(0,-0.15),
  #        col = family.colours[levels(covars$Family.ID)],
  #        legend = levels(covars$Family.ID), cex = 0.6,
  #        ncol = ceiling(length(levels(covars$Family.ID))/5), bty = "n")
  dev.off()
  
  ### Community detection
  # Louvain
  lc = cluster_louvain(mygraph)
  
  k = max(lc$membership)
  
  perf = rbind(perf, cbind(k,ClusteringPerformance(lc$membership, covars$Family.ID)))
  
  # Set cluster colour
  cluster_colour = rainbow(length(unique(lc$membership)))[lc$membership]
  names(cluster_colour) = lc$membership
  
  E(mygraph)$color = ifelse(igraph::crossing(lc, mygraph),
                            alpha('grey',0.6),
                            darken(cluster_colour[as.character(membership(lc)[ends(mygraph, E(mygraph))[,1]])],0.5))  
  
  # # Plot clusters in graph
  # pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community.pdf"))
  # par(mar=c(0,0,0,5), xpd = TRUE)
  # set.seed(1)
  # plot(mygraph, layout=layout_with_fr(mygraph),
  #      edge.width = 1,
  #      vertex.size = 9,
  #      vertex.label.cex = 0.7,
  #      vertex.frame.color = darken(cluster_colour,0.5),
  #      mark.groups = communities(lc),
  #      mark.col = cluster_colour[as.character(1:length(communities(lc)))],
  #      mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
  #      vertex.color = mynode_colour)
  # legend("right",pch=19, inset = c(-0.15,0),
  #        col = unique(cluster_colour[order(names(cluster_colour))]),
  #        legend = 1:length(unique(cluster_colour)), cex = 0.6,
  #        ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
  # dev.off()
  
  mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
  # Plot clusters in graph separately
  if (m ==4){
    pdf(paste0("../Figures/Section5/Graphical_network_children_community_",suffix[m],".pdf"))
  } else {pdf(paste0("../Figures/Supplementary/Graphical_network_children_community_",suffix[m],".pdf"))}
  par(mar=c(0,0,0,0), xpd = TRUE)
  set.seed(1)
  plot(mygraph_split, layout=layout_with_fr(mygraph_split),
       edge.width = 1,
       vertex.size = 10,
       vertex.label.cex = 0.8,
       vertex.frame.color = darken(cluster_colour,0.5),
       mark.groups = communities(lc),
       mark.col = alpha(cluster_colour[as.character(1:length(communities(lc)))],0.2),
       mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
       vertex.color = mynode_colour)
  # legend("right",pch=21, inset = c(-0.15,0),
  #        pt.cex = 2,
  #        pt.bg = alpha(unique(cluster_colour[order(as.numeric(names(cluster_colour)))]),0.2),
  #        col = darken(unique(cluster_colour[order(as.numeric(names(cluster_colour)))]),0.5),
  #        legend = paste0("C",1:length(unique(cluster_colour))),
  #        ncol = ceiling(length(unique(cluster_colour))/30), bty = "n")
  dev.off()
  
  if(m %in% c(2,4)){
    mynode_colour = region.colours[as.character(covars$Region)]
    names(mynode_colour) = rownames(covars)
    mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
    
    if (m==4){
      pdf(paste0("../Figures/Section5/Graphical_network_children_community_region_",suffix[m],".pdf"))
    } else {
      pdf(paste0("../Figures/Supplementary/Graphical_network_children_region_",suffix[m],".pdf"))
    }
    par(mar=c(0,0,5,0), xpd = TRUE)
    set.seed(1)
    plot(mygraph_split, layout=layout_with_fr(mygraph_split),
         edge.width = 1,
         vertex.size = 10,
         vertex.label.cex = 0.8,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = alpha(cluster_colour[as.character(1:length(communities(lc)))],0.2),
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("top",pch=19, inset = c(0,-0.15),
           col = region.colours[levels(covars$Region)],
           legend = levels(covars$Region),
           ncol = ceiling(length(levels(covars$Region))/3), bty = "n")
    # legend("right",pch=21, inset = c(-0.15,0),
    #        pt.cex = 2,
    #        pt.bg = alpha(unique(cluster_colour[order(names(cluster_colour))]),0.2),
    #        col = darken(unique(cluster_colour[order(names(cluster_colour))]),0.5),
    #        legend = paste0("C",1:length(unique(cluster_colour))),
    #        ncol = ceiling(length(unique(cluster_colour))/30), bty = "n")
    dev.off()
  }
}
saveRDS(perf, "../Results/Community_detection_performance.rds")


