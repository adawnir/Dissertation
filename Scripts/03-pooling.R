## Pooling data sets
## Rin Wada 21 June 2021

# Load packages
library(tidyverse)

### Pooled3 (LUX/FRA/GS) ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

covars_lux = readRDS(paste0("../Processed/",filepaths[1],"/Participant_covariate_info_thresh.rds"))
chem_lux = readRDS(paste0("../Processed/",filepaths[1],"/Chemical_compound_info.rds"))
expo_lux = readRDS(paste0("../Processed/",filepaths[1],"/Exposure_matrix_raw.rds"))

covars_fra = readRDS(paste0("../Processed/",filepaths[2],"/Participant_covariate_info_thresh.rds"))
chem_fra = readRDS(paste0("../Processed/",filepaths[2],"/Chemical_compound_info.rds"))
expo_fra = readRDS(paste0("../Processed/",filepaths[2],"/Exposure_matrix_raw.rds"))

covars_gs = readRDS(paste0("../Processed/",filepaths[3],"/Participant_covariate_info_thresh.rds"))
chem_gs = readRDS(paste0("../Processed/",filepaths[3],"/Chemical_compound_info.rds"))
expo_gs = readRDS(paste0("../Processed/",filepaths[3],"/Exposure_matrix_raw.rds"))

# Merge covariate information
covars = bind_rows(covars_lux, covars_fra, covars_gs) %>%
  select(Batch, Indiv.ID, Age, Gender, Weight, Length, Family.ID,
        Area, Department, Region, Country)
rownames(covars) = covars$Indiv.ID
str(covars)

# Check extraction consistency
extract_diff = full_join(chem_lux, chem_fra, by = "Compound", suffix = c(".lux", ".fra")) %>%
  full_join(chem_gs, by = "Compound") %>%
  mutate(Family = coalesce(Family.lux,Family)) %>%
  select(Compound, Family, starts_with("Extraction")) %>%
  .[!is.equal(list(.$Extraction.lux, .$Extraction.fra, .$Extraction)),]
dim(extract_diff)

ifelse(dir.exists("../Results/"),"",dir.create("../Results/"))
ifelse(dir.exists(paste0("../Results/",filepaths[4])),"",dir.create(paste0("../Results/",filepaths[4])))
saveRDS(extract_diff, paste0("../Results/",filepaths[4],"/Chemical_compound_info_extract_diff.rds"))

# Merge chemical compound information
expo = bind_rows(expo_lux, expo_fra, expo_gs)
expo = expo[rownames(covars),]
all(rownames(expo)==rownames(covars))

chem = full_join(chem_lux, chem_fra, by = "Compound", suffix = c(".lux", ".fra")) %>%
  full_join(chem_gs, by = "Compound") %>%
  mutate(Family = coalesce(Family.lux,Family)) %>%
  select(Compound, Family, starts_with("LOD"), starts_with("LOQ"),starts_with("Extraction"))
colnames(chem) = ifelse(grepl("\\.", colnames(chem))|colnames(chem)%in%c("Compound", "Family"),
                        colnames(chem), paste0(colnames(chem),".gs"))

# Proportion of NAs per chemical compounds
all(colnames(expo)==chem$Compound)
chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))

# Proportion of detected per chemical compounds
all(colnames(expo)==chem$Compound)
chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))

# Number of chemical compounds with 90% or more nd or NA
sum(chem$nd_prop + chem$NA_prop>=0.9)

# Save data sets
ifelse(dir.exists(paste0("../Processed/",filepaths[4])),"",dir.create(paste0("../Processed/",filepaths[4])))
saveRDS(covars, paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))
saveRDS(chem, paste0("../Processed/",filepaths[4],"/Chemical_compound_info.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[4],"/Exposure_matrix_raw.rds"))

# Filter out >=90% nd or NA
expo = expo[,which(chem$nd_prop + chem$NA_prop<0.9)]
chem = chem[which(chem$nd_prop + chem$NA_prop<0.9),]
ncol(expo)

max(rowSums(is.na(expo))/nrow(expo))

# Save data sets
saveRDS(chem, paste0("../Processed/",filepaths[4],"/Chemical_compound_info_thresh.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[4],"/Exposure_matrix_raw_thresh.rds"))

## Recoding nd
# Replace nd with random values from 0 to minimum detection (Gaussian)
min.lod = apply(chem[,grepl("LOD",colnames(chem))], 1, function(x) min(x, na.rm = T))

# Set LOD for each data set
# If LOD is 0 or missing set it as minimum of three data sets
tmp = expo[covars_lux$Indiv.ID,]
lod = ifelse(is.na(chem$LOD.lux), min.lod, chem$LOD.lux)
names(lod) = chem$Compound

tmp2 = expo[covars_fra$Indiv.ID,]
lod2 = ifelse(is.na(chem$LOD.fra), min.lod, chem$LOD.fra)
names(lod2) = chem$Compound

tmp3 = expo[covars_gs$Indiv.ID,]
lod3 = ifelse(is.na(chem$LOD.gs), min.lod, chem$LOD.gs)
names(lod3) = chem$Compound

for(k in 1:ncol(expo)){
  set.seed(150621)
  tmp[tmp[,k]=="nd"&!is.na(tmp[,k]),k] = rtruncnorm(sum(tmp[,k]=="nd"&!is.na(tmp[,k])),
                                                    min=0,
                                                    max=lod[colnames(expo)[k]]) 
  tmp2[tmp2[,k]=="nd"&!is.na(tmp2[,k]),k] = rtruncnorm(sum(tmp2[,k]=="nd"&!is.na(tmp2[,k])),
                                                    min=0,
                                                    max=lod2[colnames(expo)[k]]) 
  tmp3[tmp3[,k]=="nd"&!is.na(tmp3[,k]),k] = rtruncnorm(sum(tmp3[,k]=="nd"&!is.na(tmp3[,k])),
                                                    min=0,
                                                    max=lod3[colnames(expo)[k]]) 
}

# Merge
expo = rbind(tmp, tmp2, tmp3)

# Convert to numeric (Replace any characters with NaN and add rownames)
expo = apply(expo, 2, as.numeric)
rownames(expo) = rownames(covars)

summary(expo)
# Save data sets
saveRDS(expo, paste0("../Processed/",filepaths[4],"/Exposure_matrix_ndimp_thresh.rds"))


### Pooled2 (LUX/GS) ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

covars_lux = readRDS(paste0("../Processed/",filepaths[1],"/Participant_covariate_info_thresh.rds"))
chem_lux = readRDS(paste0("../Data/",filepaths[1],"/Chemical_compound_info.rds"))
expo_lux = readRDS(paste0("../Data/",filepaths[1],"/Exposure_matrix_raw.rds"))

covars_gs = readRDS(paste0("../Processed/",filepaths[3],"/Participant_covariate_info_thresh.rds"))
chem_gs = readRDS(paste0("../Data/",filepaths[3],"/Chemical_compound_info.rds"))
expo_gs = readRDS(paste0("../Data/",filepaths[3],"/Exposure_matrix_raw.rds"))

# Merge covariate information
covars = bind_rows(covars_lux, covars_gs) %>%
  select(Batch, Indiv.ID, Age, Gender, Weight, Length, Family.ID,
         Area, Department, Region, Country)
rownames(covars) = covars$Indiv.ID
str(covars)

# Check extraction consistency
extract_diff = full_join(chem_lux, chem_gs, by = "Compound", suffix = c(".lux", ".gs")) %>%
  mutate(Family = coalesce(Family.lux,Family.gs)) %>%
  select(Compound, Family, starts_with("Extraction")) %>%
  .[!is.equal(list(.$Extraction.lux, .$Extraction.gs)),]
dim(extract_diff)

# Merge chemical compound information
expo = bind_rows(expo_lux, expo_gs)
expo = expo[rownames(covars),]
all(rownames(expo)==rownames(covars))

chem = full_join(chem_lux, chem_gs, by = "Compound", suffix = c(".lux", ".gs")) %>%
  mutate(Family = coalesce(Family.lux,Family.gs)) %>%
  select(Compound, Family, starts_with("LOD"), starts_with("LOQ"),starts_with("Extraction"))

# Proportion of NAs per chemical compounds
all(colnames(expo)==chem$Compound)
chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))

# Proportion of detected per chemical compounds
all(colnames(expo)==chem$Compound)
chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))

# Number of chemical compounds with 90% or more nd or NA
sum(chem$nd_prop + chem$NA_prop>=0.9)

# Save data sets
ifelse(dir.exists(paste0("../Processed/",filepaths[5])),"",dir.create(paste0("../Processed/",filepaths[5])))
saveRDS(covars, paste0("../Processed/",filepaths[5],"/Participant_covariate_info_thresh.rds"))
saveRDS(chem, paste0("../Processed/",filepaths[5],"/Chemical_compound_info.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[5],"/Exposure_matrix_raw.rds"))

# Filter out >=90% nd or NA
expo = expo[,which(chem$nd_prop + chem$NA_prop<0.9)]
chem = chem[which(chem$nd_prop + chem$NA_prop<0.9),]
ncol(expo)

max(rowSums(is.na(expo))/nrow(expo))

# Save data sets
saveRDS(chem, paste0("../Processed/",filepaths[5],"/Chemical_compound_info_thresh.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[5],"/Exposure_matrix_raw_thresh.rds"))

## Recoding nd
# Replace nd with random values from 0 to minimum detection (Gaussian)
min.lod = apply(chem[,grepl("LOD",colnames(chem))], 1, function(x) min(x, na.rm = T))

# Set LOD for each data set
# If LOD is missing set it as minimum of three data sets
tmp = expo[covars_lux$Indiv.ID,]
lod = ifelse(is.na(chem$LOD.lux), min.lod, chem$LOD.lux)
names(lod) = chem$Compound

tmp3 = expo[covars_gs$Indiv.ID,]
lod3 = ifelse(is.na(chem$LOD.gs), min.lod, chem$LOD.gs)
names(lod3) = chem$Compound

for(k in 1:ncol(expo)){
  set.seed(150621)
  tmp[tmp[,k]=="nd"&!is.na(tmp[,k]),k] = rtruncnorm(sum(tmp[,k]=="nd"&!is.na(tmp[,k])),
                                                    min=0,
                                                    max=lod[colnames(expo)[k]])
  tmp3[tmp3[,k]=="nd"&!is.na(tmp3[,k]),k] = rtruncnorm(sum(tmp3[,k]=="nd"&!is.na(tmp3[,k])),
                                                       min=0,
                                                       max=lod3[colnames(expo)[k]]) 
}

# Merge
expo = rbind(tmp, tmp3)

# Convert to numeric (Replace any characters with NaN and add rownames)
expo = apply(expo, 2, as.numeric)
rownames(expo) = rownames(covars)

summary(expo)
# Save data sets
saveRDS(expo, paste0("../Processed/",filepaths[5],"/Exposure_matrix_ndimp_thresh.rds"))
