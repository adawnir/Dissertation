# Load packages
library(ggplot2)
# library(grid)
# library(rworldmap)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

### Plotting ----
# Initialisation
rm(list=ls())
path="~/HDAML/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

suffix = c("lux","fra","gs","pooled3","pooled2")

for (m in 1:length(batches)){
  # Load data
  expo = readRDS(paste0("../Processed/",filepaths[m],"/Exposure_matrix_ndimp_thresh_log_naimp_no_isolated.rds"))
  covars = readRDS(paste0("../Processed/",filepaths[m],"/Participant_covariate_info_thresh_no_isolated.rds"))
  print(all(rownames(expo)==rownames(covars)))
  
  out = readRDS(paste0("../Results/",filepaths[m],"/Stability_clustering_output.rds"))
  sc = readRDS(paste0("../Results/",filepaths[m],"/Cluster_memberships.rds"))
  sc = sc$stab_cluster
  names(sc) = rownames(covars)

  if (m == 4){
    lux = map_data("world", "Luxembourg")
    france = map_data("france")
    
    mycoords = rbind(lux, france)
    
    # area = read.csv("../Dictionaries/French_area_codes.csv")
    # 
    # covars$Department = area$Department[match(covars$Area, area$Code.Commune)]
    # covars$Department[is.na(covars$Department)] = "Luxembourg"
    
    setdiff(covars$Department, unique(mycoords$region))
    
    
    # # Add some data for each member
    # value = aggregate(expo[,"Propiconazole"] ~ Country, covars, mean)
    # mytable <- data.frame(country = levels(covars$Country), value = value[,2])
    # mycoords$value <- mytable$value[match(mycoords$region,mytable$country)]
    
    world <- ne_countries(scale = "medium", returnclass = "sf")
    fr_states <- ne_states(country = "france", returnclass = "sf")
    lux <- ne_states(country = "luxembourg", returnclass = "sf")
    
    setdiff(covars$Department, fr_states$name)
    
    fr_states$name[grep("Marne",fr_states$name)]
    
    fr_states$name[which(fr_states$name=="Seien-et-Marne")] = "Seine-et-Marne"
    
    n = as.vector(table(covars$Department))
    names(n) = names(table(covars$Department))
    n = n[as.character(fr_states$name)]
    n[is.na(n)] = 0
    names(n) = fr_states$name
    
    centroids = as.tibble(st_coordinates(st_centroid(fr_states))) %>%
      rbind(st_coordinates(st_centroid(lux))) %>%
      mutate(region = c(fr_states$name,lux$name)) %>%
      filter(region %in% unique(covars$Department))
    
    samples = NULL
    for (i in 1:length(centroids$region)){
      samples = rbind(samples, centroids[rep(i, sum(centroids$region[i]==covars$Department)), ])  
    }
    samples = samples[order(factor.order(covars$Department)),]
    samples$cluster = sc[order(factor.order(covars$Department))]

    regbydepart = region.colours[area$Region]
    names(regbydepart) = area$Department
    regbydepart = regbydepart[!duplicated(names(regbydepart))]
    
    ggplot() +
      geom_sf(data = world) +
      geom_sf(data = lux, fill = region.colours["Luxembourg"])) + 
      geom_sf(data = fr_states, fill = regbydepart[fr_states$name]) +
      coord_sf(xlim = range(mycoords$long),  ylim = range(mycoords$lat)) +
      theme_void() +
      theme(legend.position = "none")
    
    ggplot() +
      geom_sf(data = world) +
      geom_sf(data = lux, fill = region.colours["Luxembourg"]) + 
      geom_sf(data = fr_states, fill = depart[fr_states$name]) +
      coord_sf(xlim = range(mycoords$long),  ylim = range(mycoords$lat)) +
      theme_void() +
      theme(legend.position = "none")
    
  
  
  mynode_colour = family.colours[as.character(covars$Family.ID)]
  names(mynode_colour) = rownames(covars)
  
  cluster_colour = rainbow(length(unique(sc)), alpha = 0.6)[sc]
  names(cluster_colour) = sc
  
  communities = lapply(unique(sc), function(i) rownames(covars)[sc==i])
  names(communities) = unique(sc)
  
  {pdf(paste0("../Figures/",filepaths[m],"/Stab_clustering_graph.pdf"))
    par(mar=c(0,5,0,5))
    set.seed(1)
    plot(mygraph,
         mark.col = cluster_colour[as.character(1:length(communities))],
         mark.border = darken(cluster_colour[as.character(1:length(communities))], 0.5))
    legend("right",pch=19, inset = c(-0.15,0),
           col = unique(cluster_colour),
           legend = 1:length(unique(cluster_colour)), cex = 0.6,
           ncol = ceiling(length(unique(cluster_colour))/20), bty = "n")
    # legend("top",pch=19, inset = c(0,-0.15),
    #        col = depart.colours[levels(covars$Department)],
    #        legend = levels(covars$Department), cex = 0.6,
    #        ncol = ceiling(length(levels(covars$Department))/3), bty = "n")
    dev.off()
  }
}
