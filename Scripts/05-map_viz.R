# Load packages
# library(grid)
# library(rworldmap)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(tidyverse)

### Plotting ----
# Initialisation
rm(list=ls())
path="~/HDAML/Dissertation/Scripts"
setwd(path)

# Load custom
source("functions.R")
source("graph_param.R")

### Region ----
covars = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))

lux = map_data("world", "Luxembourg")
france = map_data("france")

mycoords = rbind(lux, france)

world <- ne_countries(scale = "large", returnclass = "sf")
fr_states <- ne_states(country = "france", returnclass = "sf")

lux <- ne_states(country = "luxembourg", returnclass = "sf")

setdiff(covars$Department, fr_states$name)

fr_states$name[grep("Marne",fr_states$name)]

fr_states$name[which(fr_states$name=="Seien-et-Marne")] = "Seine-et-Marne"

centroids = as.tibble(st_coordinates(st_centroid(fr_states))) %>%
  rbind(st_coordinates(st_centroid(lux))) %>%
  mutate(region = c(fr_states$name,lux$name)) %>%
  filter(region %in% unique(covars$Department))
centroids$X[centroids$region=="Nord"] = 2.3026
centroids$Y[centroids$region=="Nord"] = 51.0135
centroids = rbind(centroids, c(5.963604356638484,44.41791982331576,"Hautes-Alpes/Alpes-de-Haute-Provence",0))
centroids = centroids %>% mutate(X = as.numeric(X), Y = as.numeric(Y))

region = area$Region
names(region) = area$Department
region = c(region[!duplicated(names(region))], "Luxembourg")
names(region)[length(region)] = "Luxembourg"
centroids$group = region[centroids$region]

new_centroids = centroids %>%
  group_by(group) %>%
  summarise(X = mean(X),
            Y = mean(Y))
new_centroids$pop = as.vector(table(covars$Region)[new_centroids$group])
new_centroids$label = new_centroids$group
new_centroids$label[new_centroids$label=="Hauts-de-France"] = "Grande-Synthe"
pdf("../Figures/Section1/Map_of_France_Luxembourg_region.pdf")
  ggplot() +
    geom_sf(data = world, alpha = 0.5, lwd = 0.2, color = "grey") +
    # geom_sf(data = france, alpha = 0.5) +
    geom_point(data=centroids, aes(x=X, y=Y), pch = 17,
               color = region.colours[centroids$group])+
    geom_point(data=new_centroids, aes(x=X, y=Y, size=pop),
               color = region.colours[new_centroids$group], alpha = 0.7)+
    geom_text(data=new_centroids, aes(x=X, y=Y, label=label), nudge_y = -0.4) +
    scale_size_continuous(name="Population", range=c(2,16), breaks=c(1,10,20,40,60)) +
    coord_sf(xlim = range(mycoords$long),  ylim = range(mycoords$lat)+1) +
    theme_void() +
    guides(colour = guide_legend())+
    theme(legend.position = c(0.9,0.8))
dev.off()


### Department ----
covars = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))

lux = map_data("world", "Luxembourg")
france = map_data("france")

mycoords = rbind(lux, france)

world <- ne_countries(scale = "large", returnclass = "sf")
fr_states <- ne_states(country = "france", returnclass = "sf")

lux <- ne_states(country = "luxembourg", returnclass = "sf")

setdiff(covars$Department, fr_states$name)

fr_states$name[grep("Marne",fr_states$name)]

fr_states$name[which(fr_states$name=="Seien-et-Marne")] = "Seine-et-Marne"

centroids = as.tibble(st_coordinates(st_centroid(fr_states))) %>%
  rbind(st_coordinates(st_centroid(lux))) %>%
  mutate(region = c(fr_states$name,lux$name)) %>%
  filter(region %in% unique(covars$Department))
centroids$X[centroids$region=="Nord"] = 2.3026
centroids$Y[centroids$region=="Nord"] = 51.0135
centroids = rbind(centroids, c(5.963604356638484,44.41791982331576,"Hautes-Alpes/Alpes-de-Haute-Provence",0))
centroids = centroids %>% mutate(X = as.numeric(X), Y = as.numeric(Y))
centroids$pop = as.vector(table(covars$Department)[centroids$region])

depart = area$Department
names(depart) = area$Region
depart = depart[!duplicated(depart)]
mycolours = region.colours[c(names(depart),"Luxembourg")]
names(mycolours) = c(depart,"Luxembourg")

pdf("../Figures/Map_of_France_Luxembourg.pdf")
ggplot() +
  geom_sf(data = world, alpha = 0.5) +
  # geom_sf(data = france, alpha = 0.5) +
  geom_label(aes(x=centroids$X[centroids$region=="Luxembourg"],
                 y=centroids$Y[centroids$region=="Luxembourg"]-.5, label="Luxembourg"),
             fill = "white",label.size = NA) +
  geom_label(aes(x=centroids$X[centroids$region=="Nord"],
                 y=centroids$Y[centroids$region=="Nord"]-.5, label="Grande-Synthe"),
             fill = "white",label.size = NA) +
  geom_point(data=centroids, aes(x=X, y=Y, size=pop),
             color = depart.colours[centroids$region], alpha = 0.9) +
  scale_size_continuous(name="Population", range=c(2,16), breaks=c(1,5,10,25,50)) +
  scale_alpha_continuous(trans="log") +
  coord_sf(xlim = range(mycoords$long),  ylim = range(mycoords$lat)+1) +
  theme_void() +
  guides(colour = guide_legend())+
  theme(legend.position = c(0.9,0.8))
dev.off()
