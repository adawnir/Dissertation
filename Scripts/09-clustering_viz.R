## Clustering visualisations
## Rin Wada 11 Aug

### Consensus clustering (PLS-DA score plot) ----
# Load packages
library(sgPLS)
library(plotrix)
library(ellipse)
library(igraph)


# Initialisation
rm(list=ls())
path="~/HDAML/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

perf = NULL
for (m in 1:length(batches)){# Load data
sc = readRDS(paste0("../Results/",filepaths[m],"/Stable_clusters.rds"))
stab = readRDS(paste0("../Results/",filepaths[m],"/Consensus_clustering_output.rds"))
covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))

families=unique(sc)[table(sc)!=1]

k = max(sc)

perf = rbind(perf, cbind(k,ClusteringPerformance(sc, covars$Family.ID)))

mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families

tmp = rep(alpha("grey80",0.3),length(setdiff(unique(sc), families)))
names(tmp) = setdiff(unique(sc), families)

mycolours = c(mycolours, tmp)
mycolours = mycolours[as.character(unique(sc))]
  
myplsda = plsda(expo, as.factor(sc), ncomp = 3)

CreateScorePlot.plsda(myplsda=myplsda, type1=sc, type2=covars$Family.ID,
                      mycolours=mycolours,
                      filename=paste0("../Figures/",filepaths[m],"/Consensus_clustering_PLSDA_score_plot.pdf"))

myadjacency = Adjacency(stab)
mynode_colour = family.colours[as.character(covars$Family.ID)]
names(mynode_colour) = rownames(covars)

communities = lapply(unique(sc), function(i) rownames(covars)[sc==i])
names(communities) = unique(sc)

cluster_colour = rainbow(length(unique(sc)), alpha = 0.6)[sc]
names(cluster_colour) = sc

mygraph = Graph(myadjacency, node_colour = mynode_colour,
                node_label = as.character(covars$Family.ID),
                satellites = TRUE)

E(mygraph)$color = darken(cluster_colour[as.character(sc[ends(mygraph, E(mygraph))[,1]])],0.5)  

E(mygraph)$width = 1
V(mygraph)$size = 9
V(mygraph)$label.cex = 0.7



# Plot clusters in graph
pdf(paste0("../Figures/",filepaths[m],"/Consensus_clustering_graph.pdf"))
par(mar=c(0,0,0,5), xpd = TRUE)
set.seed(1)
plot(mygraph,
     vertex.frame.color = darken(cluster_colour,0.5),
     mark.groups = communities,
     mark.col = cluster_colour[as.character(1:length(communities))],
     mark.border = darken(cluster_colour[as.character(1:length(communities))], 0.5))
legend("right",pch=19, inset = c(-0.15,0),
       col = cluster_colour[as.character(unique(sc))],
       legend = unique(sc), cex = 0.6,
       ncol = ceiling(length(unique(cluster_colour))/25), bty = "n")

dev.off()
}

mynode_colour = region.colours[as.character(covars$Region)]
names(mynode_colour) = rownames(covars)

# Plot clusters in graph
pdf(paste0("../Figures/",filepaths[m],"/Consensus_clustering_graph_region.pdf"))
par(mar=c(0,0,5,5), xpd = TRUE)
set.seed(1)
plot(mygraph,
     vertex.frame.color = darken(cluster_colour,0.5),
     mark.groups = communities,
     mark.col = cluster_colour[as.character(1:length(communities))],
     mark.border = darken(cluster_colour[as.character(1:length(communities))], 0.5))
legend("right",pch=19, inset = c(-0.15,0),
       col = cluster_colour[as.character(unique(sc))],
       legend = unique(sc), cex = 0.6,
       ncol = ceiling(length(unique(cluster_colour))/25), bty = "n")
legend("top",pch=19, inset = c(0,-0.15),
       col = region.colours[levels(covars$Region)],
       legend = levels(covars$Region), cex = 0.6,
       ncol = ceiling(length(levels(covars$Region))/3), bty = "n")
dev.off()

saveRDS(perf, "../Results/Consensus_clustering_performance.rds")

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

  # Set cluster colour
  cluster_colour = rainbow(length(unique(lc$membership)), alpha = 0.6)[membership(lc)]
  names(cluster_colour) = membership(lc)
  
  E(mygraph)$color = ifelse(igraph::crossing(lc, mygraph),
                            alpha('grey',0.6),
                            darken(cluster_colour[as.character(membership(lc)[ends(mygraph, E(mygraph))[,1]])],0.5))  
  
  # Plot clusters in graph
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community.pdf"))
  par(mar=c(0,0,0,5), xpd = TRUE)
  set.seed(1)
  plot(mygraph, layout=layout_with_fr(mygraph),
       edge.width = 1,
       vertex.size = 9,
       vertex.label.cex = 0.7,
       vertex.frame.color = darken(cluster_colour,0.5),
       mark.groups = communities(lc),
       mark.col = cluster_colour[as.character(1:length(communities(lc)))],
       mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
       vertex.color = mynode_colour)
  legend("right",pch=19, inset = c(-0.15,0),
         col = unique(cluster_colour[order(names(cluster_colour))]),
         legend = 1:length(unique(cluster_colour)), cex = 0.6,
         ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
  dev.off()
  
  mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
  # Plot clusters in graph separately
  pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_split.pdf"))
  par(mar=c(0,0,0,5), xpd = TRUE)
  set.seed(1)
  plot(mygraph_split, layout=layout_with_fr(mygraph_split),
       edge.width = 1,
       vertex.size = 9,
       vertex.label.cex = 0.7,
       vertex.frame.color = darken(cluster_colour,0.5),
       mark.groups = communities(lc),
       mark.col = cluster_colour[as.character(1:length(communities(lc)))],
       mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
       vertex.color = mynode_colour)
  legend("right",pch=19, inset = c(-0.15,0),
         col = unique(cluster_colour[order(names(cluster_colour))]),
         legend = 1:length(unique(cluster_colour)), cex = 0.6,
         ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
  dev.off()
  
  if(m %in% c(4,5)){
    mynode_colour = batch.colours[as.character(covars$Batch)]
    names(mynode_colour) = rownames(covars)
    
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_batch_community.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph, layout=layout_with_fr(mygraph),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = batch.colours[levels(covars$Batch)],
           legend = levels(covars$Batch), cex = 0.6,
           horiz = TRUE, bty = "n")
    dev.off()
    
    mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
    # Plot clusters in graph separately
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_batch_split.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph_split, layout=layout_with_fr(mygraph_split),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = batch.colours[covars$Batch])
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = batch.colours[levels(covars$Batch)],
           legend = levels(covars$Batch), cex = 0.6,
           horiz = TRUE, bty = "n")
    dev.off()  
  }
  
  if(m %in% c(2,4)){
    mynode_colour = region.colours[as.character(covars$Region)]
    names(mynode_colour) = rownames(covars)
    
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_region_community.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph, layout=layout_with_fr(mygraph),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = region.colours[levels(covars$Region)],
           legend = levels(covars$Region), cex = 0.6,
           ncol = ceiling(length(levels(covars$Region))/3), bty = "n")
    dev.off()
    
    mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
    
    # Plot clusters in graph separately
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_region_split.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph_split, layout=layout_with_fr(mygraph_split),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = region.colours[levels(covars$Region)],
           legend = levels(covars$Region), cex = 0.6,
           ncol = ceiling(length(levels(covars$Region))/3), bty = "n")
    dev.off()  
    
    # Department
    mynode_colour = depart.colours[as.character(covars$Department)]
    names(mynode_colour) = rownames(covars)
    
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_depart_community.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph, layout=layout_with_fr(mygraph),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = depart.colours[levels(covars$Department)],
           legend = levels(covars$Department), cex = 0.6,
           ncol = ceiling(length(levels(covars$Department))/3), bty = "n")
    dev.off()
    
    mygraph_split = delete_edges(mygraph, E(mygraph)[igraph::crossing(lc, mygraph)])
    
    # Plot clusters in graph separately
    pdf(paste0("../Figures/",filepaths[m],"/Graphical_network_children_community_depart_split.pdf"))
    par(mar=c(0,0,5,5), xpd = TRUE)
    set.seed(1)
    plot(mygraph_split, layout=layout_with_fr(mygraph_split),
         edge.width = 1,
         vertex.size = 9,
         vertex.label.cex = 0.7,
         vertex.frame.color = darken(cluster_colour,0.5),
         mark.groups = communities(lc),
         mark.col = cluster_colour[as.character(1:length(communities(lc)))],
         mark.border = darken(cluster_colour[as.character(1:length(communities(lc)))], 0.5),
         vertex.color = mynode_colour)
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour[order(names(cluster_colour))]),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
    legend("top",pch=19, inset = c(0,-0.15),
           col = depart.colours[levels(covars$Department)],
           legend = levels(covars$Department), cex = 0.6,
           ncol = ceiling(length(levels(covars$Department))/3), bty = "n")
    dev.off()  
  }
}
