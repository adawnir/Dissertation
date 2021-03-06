## Pooling data sets
## 21 June 2021

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
covars$Family.ID = last.level(covars$Family.ID, "Isolated")
levels(covars$Family.ID)
str(covars)

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
rownames(chem) = chem$Compound

# Reorder columns
all(colnames(expo)==rownames(chem))
expo = expo[,match(rownames(chem),colnames(expo))]
all(colnames(expo)==rownames(chem))

# Proportion of NAs per chemical compounds
chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))

# Proportion of non-detects per chemical compounds
all(colnames(expo)==rownames(chem))
chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))

# Detection rate
chem$detect_rate = 1-(chem$nd_prop + chem$NA_prop)

# Number of chemical compounds with more than 10% detected
sum(chem$detect_rate>0.1)

# Number of chemical compounds with 10% or less detected
sum(chem$detect_rate<=0.1)

# Save data sets
ifelse(dir.exists(paste0("../Processed/",filepaths[4])),"",dir.create(paste0("../Processed/",filepaths[4])))
saveRDS(covars, paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))
saveRDS(chem, paste0("../Processed/",filepaths[4],"/Chemical_compound_info.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[4],"/Exposure_matrix_raw.rds"))

# Filter out compounds with 10% or less detected
all(colnames(expo)==rownames(chem))
expo = expo[,which(chem$detect_rate>0.1)]
chem = chem[which(chem$detect_rate>0.1),]
ncol(expo)

# Number of chemical compounds not detected in one data set
grep = intersect(rownames(chem_lux)[which(chem_lux$detect_rate!=0)],
                 rownames(chem_fra)[which(chem_fra$detect_rate!=0)]) %>%
  intersect(rownames(chem_gs)[which(chem_gs$detect_rate!=0)])

sum(rownames(chem) %in% grep)

expo = expo[,which(rownames(chem) %in% grep)]
chem = chem[which(rownames(chem) %in% grep),]
ncol(expo)

table(chem$Family)
chem$Family = droplevels(chem$Family)
table(chem$Family)

# Save data sets
saveRDS(chem, paste0("../Processed/",filepaths[4],"/Chemical_compound_info_thresh.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[4],"/Exposure_matrix_raw_thresh.rds"))

# Check extraction consistency
extract_diff = chem %>%
  select(Compound, Family, starts_with("Extraction")) %>%
  .[!is.equal(list(.$Extraction.lux, .$Extraction.fra, .$Extraction.gs)),]
dim(extract_diff)

ifelse(dir.exists("../Results/"),"",dir.create("../Results/"))
ifelse(dir.exists(paste0("../Results/",filepaths[4])),"",dir.create(paste0("../Results/",filepaths[4])))
saveRDS(extract_diff, paste0("../Results/",filepaths[4],"/Chemical_compound_info_extract_diff.rds"))

## Recoding nd
# Replace nd with random values from 0 to minimum detection (Gaussian)
# Set LOD for each data set
tmp1 = expo[covars$Batch=="LUX",]
lod1 = chem$LOD.lux
names(lod1) = rownames(chem)

tmp2 = expo[covars$Batch=="FRA",]
lod2 = chem$LOD.fra
names(lod2) = rownames(chem)

tmp3 = expo[covars$Batch=="GS",]
lod3 = chem$LOD.gs
names(lod3) = rownames(chem)

for(k in 1:ncol(expo)){
  tmp1[which(tmp1[,k]=="nd"),k] = rtruncnorm(sum(tmp1[,k]=="nd", na.rm = T),min=0,max=lod1[colnames(expo)[k]]) 
  tmp2[which(tmp2[,k]=="nd"),k] = rtruncnorm(sum(tmp2[,k]=="nd", na.rm = T),min=0,max=lod2[colnames(expo)[k]]) 
  tmp3[which(tmp3[,k]=="nd"),k] = rtruncnorm(sum(tmp3[,k]=="nd", na.rm = T),min=0,max=lod3[colnames(expo)[k]]) 
}

# Merge
expo = rbind(tmp1, tmp2, tmp3)

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
chem_lux = readRDS(paste0("../Processed/",filepaths[1],"/Chemical_compound_info.rds"))
expo_lux = readRDS(paste0("../Processed/",filepaths[1],"/Exposure_matrix_raw.rds"))

covars_gs = readRDS(paste0("../Processed/",filepaths[3],"/Participant_covariate_info_thresh.rds"))
chem_gs = readRDS(paste0("../Processed/",filepaths[3],"/Chemical_compound_info.rds"))
expo_gs = readRDS(paste0("../Processed/",filepaths[3],"/Exposure_matrix_raw.rds"))

# Merge covariate information
covars = bind_rows(covars_lux, covars_gs) %>%
  select(Batch, Indiv.ID, Age, Gender, Weight, Length, Family.ID,
         Area, Department, Region, Country)
rownames(covars) = covars$Indiv.ID
covars$Family.ID = last.level(covars$Family.ID, "Isolated")
levels(covars$Family.ID)
str(covars)

# Check extraction consistency
extract_diff = full_join(chem_lux, chem_gs, by = "Compound", suffix = c(".lux", ".gs")) %>%
  mutate(Family = coalesce(Family.lux,Family.gs)) %>%
  select(Compound, Family, starts_with("Extraction")) %>%
  .[!is.equal(list(.$Extraction.lux, .$Extraction.gs)),]
dim(extract_diff)

if(nrow(extract_diff)>0){
  ifelse(dir.exists(paste0("../Results/",filepaths[5])),"",dir.create(paste0("../Results/",filepaths[5])))
  saveRDS(extract_diff, paste0("../Results/",filepaths[5],"/Chemical_compound_info_extract_diff.rds"))
}

# Merge chemical compound information
expo = bind_rows(expo_lux, expo_gs)
expo = expo[rownames(covars),]
all(rownames(expo)==rownames(covars))

chem = full_join(chem_lux, chem_gs, by = "Compound", suffix = c(".lux", ".gs")) %>%
  mutate(Family = coalesce(Family.lux,Family.gs)) %>%
  select(Compound, Family, starts_with("LOD"), starts_with("LOQ"),starts_with("Extraction"))

rownames(chem) = chem$Compound

# Reorder columns
all(colnames(expo)==rownames(chem))
expo = expo[,match(rownames(chem),colnames(expo))]
all(colnames(expo)==rownames(chem))

# Proportion of NAs per chemical compounds
chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))

# Proportion of detected per chemical compounds
all(colnames(expo)==rownames(chem))
chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))

# Detection rate
chem$detect_rate = 1-(chem$nd_prop + chem$NA_prop)

# Number of chemical compounds with more than 10% detected
sum(chem$detect_rate>0.1)

# Number of chemical compounds with 10% or less detected
sum(chem$detect_rate<=0.1)

# Save data sets
ifelse(dir.exists(paste0("../Processed/",filepaths[5])),"",dir.create(paste0("../Processed/",filepaths[5])))
saveRDS(covars, paste0("../Processed/",filepaths[5],"/Participant_covariate_info_thresh.rds"))
saveRDS(chem, paste0("../Processed/",filepaths[5],"/Chemical_compound_info.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[5],"/Exposure_matrix_raw.rds"))

# Filter out compounds with 10% or less detected
expo = expo[,which(chem$detect_rate>0.1)]
chem = chem[which(chem$detect_rate>0.1),]
ncol(expo)

# Number of chemical compounds not detected in one data set
grep = intersect(chem_lux$Compound[which(chem_lux$detect_rate!=0)], chem_gs$Compound[which(chem_gs$detect_rate!=0)])

sum(rownames(chem) %in% grep)

expo = expo[,which(rownames(chem) %in% grep)]
chem = chem[which(rownames(chem) %in% grep),]
ncol(expo)

table(chem$Family)
chem$Family = droplevels(chem$Family)
table(chem$Family)

# Save data sets
saveRDS(chem, paste0("../Processed/",filepaths[5],"/Chemical_compound_info_thresh.rds"))
saveRDS(expo, paste0("../Processed/",filepaths[5],"/Exposure_matrix_raw_thresh.rds"))

## Recoding nd
# Replace nd with random values from 0 to minimum detection (Gaussian)
# Set LOD for each data set
tmp1 = expo[covars$Batch=="LUX",]
lod1 = chem$LOD.lux
names(lod1) = rownames(chem)

tmp3 = expo[covars$Batch=="GS",]
lod3 = chem$LOD.gs
names(lod3) = rownames(chem)

for(k in 1:ncol(expo)){
  set.seed(150621)
  tmp1[which(tmp1[,k]=="nd"),k] = rtruncnorm(sum(tmp1[,k]=="nd", na.rm = T),min=0,max=lod1[colnames(expo)[k]]) 
  set.seed(150621)
  tmp3[which(tmp3[,k]=="nd"),k] = rtruncnorm(sum(tmp3[,k]=="nd", na.rm = T),min=0,max=lod3[colnames(expo)[k]]) 
}


# Merge
expo = rbind(tmp1, tmp3)

# Convert to numeric (Replace any characters with NaN and add rownames)
expo = apply(expo, 2, as.numeric)
rownames(expo) = rownames(covars)

summary(expo)
# Save data sets
saveRDS(expo, paste0("../Processed/",filepaths[5],"/Exposure_matrix_ndimp_thresh.rds"))
