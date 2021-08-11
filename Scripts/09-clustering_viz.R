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
       col = unique(cluster_colour[order(names(cluster_colour))]),
       legend = 1:length(unique(cluster_colour)), cex = 0.6,
       ncol = ceiling(length(unique(cluster_colour))/10), bty = "n")
dev.off()
}
saveRDS(perf, "../Results/Consensus_clustering_performance.rds")

# sclc = cluster_louvain(mygraph)
# set.seed(1)
# plot(sclc, mygraph)


lc = readRDS(paste0("../Results/",filepaths[m],"/Graphical_network_children_community.rds"))
all(membership(lc) == membership(sclc))
ClusteringPerformance(sc, covars$Family.ID)
ClusteringPerformance(membership(lc), covars$Family.ID)
ClusteringPerformance(membership(sclc), covars$Family.ID)
