## Univariate clustering analysis
## Rin Wada 12 July

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ape)
library(igraph)