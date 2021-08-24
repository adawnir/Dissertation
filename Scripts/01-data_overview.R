## Raw data overview
## 8 July

# Load packages
library(openxlsx)
library(tidyverse)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

covars_lux = readRDS(paste0("../Data/",filepaths[1],"/Participant_covariate_info.rds"))
chem_lux = readRDS(paste0("../Data/",filepaths[1],"/Chemical_compound_info.rds"))
expo_lux = readRDS(paste0("../Data/",filepaths[1],"/Exposure_matrix_raw.rds"))

covars_fra = readRDS(paste0("../Data/",filepaths[2],"/Participant_covariate_info.rds"))
chem_fra = readRDS(paste0("../Data/",filepaths[2],"/Chemical_compound_info.rds"))
expo_fra = readRDS(paste0("../Data/",filepaths[2],"/Exposure_matrix_raw.rds"))

covars_gs = readRDS(paste0("../Data/",filepaths[3],"/Participant_covariate_info.rds"))
chem_gs = readRDS(paste0("../Data/",filepaths[3],"/Chemical_compound_info.rds"))
expo_gs = readRDS(paste0("../Data/",filepaths[3],"/Exposure_matrix_raw.rds"))

## Family sibling ID
# Descriptive
table(covars_lux$annot.ID)
table(table(covars_lux$annot.ID))

table(covars_fra$Siblings.Groups)
table(table(covars_fra$Siblings.Groups))

table(covars_gs$Siblings.Groups)
table(table(covars_gs$Siblings.Groups))

# Missing gender information
cat(sum(is.na(covars_lux$Gender)),
    " (", round(sum(is.na(covars_lux$Gender))/length(covars_lux$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_fra$Gender)),
    " (", round(sum(is.na(covars_fra$Gender))/length(covars_fra$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_gs$Gender)),
    " (", round(sum(is.na(covars_gs$Gender))/length(covars_gs$Gender)*100,2), "%)\n", sep = "")

# Age distribution
ifelse(dir.exists("../Figures/"),"",dir.create("../Figures/"))
{pdf("../Figures/Age_dist.pdf", width = 12, height = 4)
  par(mfrow=c(1,3), mar = c(5,5,2,1))
  hist(covars_lux$Age, main = batches[1], xlab = "Age (years)", col = batch.colours[1])
  hist(covars_fra$Age, main = batches[2], xlab = "Age (years)", col = batch.colours[2])
  hist(covars_gs$Age, main = batches[3], xlab = "Age (years)", col = batch.colours[3])
  dev.off()
  }

# Missing age information
cat(sum(is.na(covars_lux$Age)),
    " (", round(sum(is.na(covars_lux$Age))/length(covars_lux$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_fra$Age)),
    " (", round(sum(is.na(covars_fra$Age))/length(covars_fra$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_gs$Age)),
    " (", round(sum(is.na(covars_gs$Age))/length(covars_gs$Age)*100,2), "%)\n", sep = "")

## Hair sample length and weight (Luxembourg)
ifelse(dir.exists(paste0("../Figures/",filepaths[1])),"",dir.create(paste0("../Figures/",filepaths[1])))
{pdf(paste0("../Figures/",filepaths[1],"/Sample_length_weight_dist.pdf"), width=8, height=4)
  par(mfrow = c(1,2), mar = c(5,5,1,1))
  hist(covars_lux$Length,
       xlab = "Length of hair sample (cm)", main = NULL, col = batch.colours[1])
  hist(covars_lux$Weight,
       xlab = "Weight of hair sample (mg)", main = NULL, col = batch.colours[1])
  dev.off()
}

# Check consistency of colnames
setdiff(colnames(expo_lux),colnames(expo_gs))
setdiff(colnames(expo_fra),colnames(expo_lux))
setdiff(colnames(expo_gs),colnames(expo_lux))

# Check chemical family consistency
all(chem_lux$Family[which(chem_lux$Compound%in%chem_gs$Compound)]==as.character(chem_gs$Family))
all(chem_lux$Family[which(chem_lux$Compound%in%chem_fra$Compound)]==as.character(chem_fra$Family))

