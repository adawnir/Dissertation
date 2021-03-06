## Data exploration and pre-processing
## 15 June

# Load packages
library(tidyverse)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

covars_lux = readRDS(paste0("../Data/",filepaths[1],"/Participant_covariate_info.rds"))
chem_lux = readRDS(paste0("../Data/",filepaths[1],"/Chemical_compound_info.rds"))
expo_lux = readRDS(paste0("../Data/",filepaths[1],"/Exposure_matrix_raw.rds"))

covars_fra = readRDS(paste0("../Data/",filepaths[2],"/Participant_covariate_info.rds"))
chem_fra = readRDS(paste0("../Data/",filepaths[2],"/Chemical_compound_info.rds"))
expo_fra = readRDS(paste0("../Data/",filepaths[2],"/Exposure_matrix_raw.rds"))

covars_gs = readRDS(paste0("../Data/",filepaths[3],"/Participant_covariate_info.rds"))
chem_gs = readRDS(paste0("../Data/",filepaths[3],"/Chemical_compound_info.rds"))
expo_gs = readRDS(paste0("../Data/",filepaths[3],"/Exposure_matrix_raw.rds"))

### Inclusion/exclusion criteria (Participant) ----
min(covars_lux$Age, na.rm = T)
max(covars_lux$Age, na.rm = T)

min(covars_fra$Age, na.rm = T)
max(covars_fra$Age, na.rm = T)

min(covars_gs$Age, na.rm = T)
max(covars_gs$Age, na.rm = T)

## Exclude adults aged 18 or over
covars_gs = covars_gs[-which(covars_gs$Age >= 18),]
max(covars_gs$Age, na.rm = T)

# Exclude hair samples with <35mg
min(covars_lux$Weight, na.rm = T)
covars_lux = covars_lux[-which(covars_lux$Weight < 35),]
min(covars_lux$Weight, na.rm = T)

# Randomly sampling 2 children out of 3-4 child groups
set.seed(280621)
covars_lux = covars_lux[which(group_sample(covars_lux$Family.ID, covars_lux$Family.ID, 2)),]
covars_fra = covars_fra[which(group_sample(covars_fra$Siblings.Groups, covars_fra$Siblings.Groups, 2, "Isolated")),]
covars_gs = covars_gs[which(group_sample(covars_gs$Siblings.Groups, covars_gs$Siblings.Groups, 2, "Isolated")),]
nrow(covars_lux)
nrow(covars_fra)
nrow(covars_gs)

# Recoding Family Sibling ID
covars_lux$Family.ID = ifelse(covars_lux$Sibling.ID=="Only child","Isolated",covars_lux$Family.ID)
covars_lux$Sibling.ID = ifelse(covars_lux$Sibling.ID=="Only child","Isolated",covars_lux$Sibling.ID)
table(covars_lux$Family.ID)

# Adjust family id
tmp = covars_lux$Family.ID[which(covars_lux$Family.ID!="Isolated")]
t = table(covars_lux$Family.ID)[unique(tmp)]
family.number = rep(paste0("L",seq(1, length.out = sum(!duplicated(tmp)))),t)
covars_lux$Family.ID[which(covars_lux$Family.ID!="Isolated")] = family.number
covars_lux$Sibling.ID[which(covars_lux$Family.ID!="Isolated")] = ave(tmp, tmp, FUN=function(a) 1:length(a))
table(covars_lux$Family.ID)
table(covars_lux$Sibling.ID)
str(covars_lux)

covars_lux = covars_lux %>%
  select(Indiv.ID, Age, Gender, Family.ID, Sibling.ID, Weight, Length,
         Area, Department, Region, Country, Batch) %>%
  mutate_if(is.character, factor.order)
str(covars_lux)

# Recoding Family Sibling ID
table(covars_fra$Siblings.Groups)
# Adjust family id
tmp = covars_fra$Siblings.Groups[which(covars_fra$Siblings.Groups!="Isolated")]
t = table(covars_fra$Siblings.Groups)[unique(tmp)]
family.number = rep(paste0("F",seq(1, length.out = sum(!duplicated(tmp)))),t)
covars_fra$Family.ID = covars_fra$Siblings.Groups
covars_fra$Sibling.ID = covars_fra$Siblings.Groups
covars_fra$Family.ID[which(covars_fra$Family.ID!="Isolated")] = family.number
covars_fra$Sibling.ID[which(covars_fra$Family.ID!="Isolated")] = ave(tmp, tmp, FUN=function(a) 1:length(a))
table(covars_fra$Family.ID)
table(covars_fra$Sibling.ID)
str(covars_fra)
covars_fra = covars_fra = covars_fra %>%
  select(Indiv.ID, Age, Gender, Family.ID, Sibling.ID,
         Area, Department, Region, Country, Batch) %>%
  mutate_if(is.character, factor.order)
covars_fra$Family.ID = last.level(covars_fra$Family.ID, "Isolated")
levels(covars_fra$Family.ID)
str(covars_fra)

# Recoding Family Sibling ID
table(covars_gs$Siblings.Groups)
covars_gs$Siblings.Groups = ifelse(covars_gs$Siblings.Groups%in%c("Family 2", "Family 3"),
                                  "Isolated",covars_gs$Siblings.Groups)
# Adjust family id
tmp = covars_gs$Siblings.Groups[which(covars_gs$Siblings.Groups!="Isolated")]
t = table(covars_gs$Siblings.Groups)[unique(tmp)]
family.number = rep(paste0("G",seq(1, length.out = sum(!duplicated(tmp)))),t)
covars_gs$Family.ID = covars_gs$Siblings.Groups
covars_gs$Sibling.ID = covars_gs$Siblings.Groups
covars_gs$Family.ID[which(covars_gs$Family.ID!="Isolated")] = family.number
covars_gs$Sibling.ID[which(covars_gs$Family.ID!="Isolated")] = ave(tmp, tmp, FUN=function(a) 1:length(a))
table(covars_gs$Family.ID)
table(covars_gs$Sibling.ID)
str(covars_gs)
covars_gs = covars_gs = covars_gs %>%
  select(Indiv.ID, Age, Gender, Family.ID, Sibling.ID,
         Area, Department, Region, Country, Batch) %>%
  mutate_if(is.character, factor.order)
str(covars_gs)
covars_gs$Family.ID = last.level(covars_gs$Family.ID, "Isolated")
levels(covars_gs$Family.ID)
# Relevel gender
covars_gs$Gender = relevel(covars_gs$Gender, "Female")

# Add rownames
rownames(covars_lux) = covars_lux$Indiv.ID
rownames(covars_fra) = covars_fra$Indiv.ID
rownames(covars_gs) = covars_gs$Indiv.ID

# Adjust chemical matrices
expo_lux = expo_lux[rownames(covars_lux),]
expo_fra = expo_fra[rownames(covars_fra),]
expo_gs = expo_gs[rownames(covars_gs),]

all(rownames(covars_lux)==rownames(rownames(expo_lux)))
all(rownames(covars_fra)==rownames(rownames(expo_fra)))
all(rownames(covars_gs)==rownames(rownames(expo_gs)))


suffix = c("lux","fra","gs")
for(i in 1:3){
  expo = eval(parse(text = paste0("expo_",suffix[i])))
  chem = eval(parse(text = paste0("chem_",suffix[i])))
  all(colnames(expo)==rownames(chem))
  
  # Proportion of NAs per chemical compounds
  chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))
  # Proportion of detected per chemical compounds
  chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))
  
  # Detection rate
  if(sum(is.na(chem$nd_prop))+sum(is.na(chem$NA_prop))==0){
    chem$detect_rate = 1-(chem$nd_prop + chem$NA_prop)
  } else {warning("NA in NA_prop or nd_prop")}
  
  assign(paste0("chem_",suffix[i]), chem)
  assign(paste0("expo_",suffix[i]), expo)
}

ifelse(dir.exists("../Processed/"),"",dir.create("../Processed/"))
ifelse(dir.exists(paste0("../Processed/",filepaths[1])),"",dir.create(paste0("../Processed/",filepaths[1])))
ifelse(dir.exists(paste0("../Processed/",filepaths[2])),"",dir.create(paste0("../Processed/",filepaths[2])))
ifelse(dir.exists(paste0("../Processed/",filepaths[3])),"",dir.create(paste0("../Processed/",filepaths[3])))

saveRDS(covars_lux, paste0("../Processed/",filepaths[1],"/Participant_covariate_info_thresh.rds"))
saveRDS(covars_fra,  paste0("../Processed/",filepaths[2],"/Participant_covariate_info_thresh.rds"))
saveRDS(covars_gs,  paste0("../Processed/",filepaths[3],"/Participant_covariate_info_thresh.rds"))

saveRDS(chem_lux, paste0("../Processed/",filepaths[1],"/Chemical_compound_info.rds"))
saveRDS(chem_fra,  paste0("../Processed/",filepaths[2],"/Chemical_compound_info.rds"))
saveRDS(chem_gs,  paste0("../Processed/",filepaths[3],"/Chemical_compound_info.rds"))

saveRDS(expo_lux, paste0("../Processed/",filepaths[1],"/Exposure_matrix_raw.rds"))
saveRDS(expo_fra,  paste0("../Processed/",filepaths[2],"/Exposure_matrix_raw.rds"))
saveRDS(expo_gs,  paste0("../Processed/",filepaths[3],"/Exposure_matrix_raw.rds"))

### Inclusion/exclusion criteria (Chemical compound) ----
# Number of chemical compounds with more than 10% detected
sum(chem_lux$detect_rate>0.1)
sum(chem_fra$detect_rate>0.1)
sum(chem_gs$detect_rate>0.1)

# Number of chemical compounds with 10% or less detected
sum(chem_lux$detect_rate<=0.1)
sum(chem_fra$detect_rate<=0.1)
sum(chem_gs$detect_rate<=0.1)

# Filter out compounds with 10% or less detected
expo_lux = expo_lux[,which(chem_lux$detect_rate>0.1)]
expo_fra = expo_fra[,which(chem_fra$detect_rate>0.1)]
expo_gs = expo_gs[,which(chem_gs$detect_rate>0.1)]

chem_lux = chem_lux[which(chem_lux$detect_rate>0.1),]
chem_fra = chem_fra[which(chem_fra$detect_rate>0.1),]
chem_gs = chem_gs[which(chem_gs$detect_rate>0.1),]

ncol(expo_lux)
ncol(expo_fra)
ncol(expo_gs)

table(chem_lux$Family)
table(chem_fra$Family)
table(chem_gs$Family)
chem_gs$Family = droplevels(chem_gs$Family)


saveRDS(chem_lux, paste0("../Processed/",filepaths[1],"/Chemical_compound_info_thresh.rds"))
saveRDS(chem_fra,  paste0("../Processed/",filepaths[2],"/Chemical_compound_info_thresh.rds"))
saveRDS(chem_gs,  paste0("../Processed/",filepaths[3],"/Chemical_compound_info_thresh.rds"))

saveRDS(expo_lux, paste0("../Processed/",filepaths[1],"/Exposure_matrix_raw_thresh.rds"))
saveRDS(expo_fra,  paste0("../Processed/",filepaths[2],"/Exposure_matrix_raw_thresh.rds"))
saveRDS(expo_gs,  paste0("../Processed/",filepaths[3],"/Exposure_matrix_raw_thresh.rds"))

### Recoding nd ----
# Replace nd with random values from 0 to minimum detection (Gaussian)
for(k in 1:ncol(expo_lux)){
  expo_lux[which(expo_lux[,k]=="nd"),k] = rtruncnorm(sum(expo_lux[,k]=="nd", na.rm = T),min=0,max=chem_lux$LOD[k]) 
}
for(k in 1:ncol(expo_fra)){
  expo_fra[which(expo_fra[,k]=="nd"),k] = rtruncnorm(sum(expo_fra[,k]=="nd", na.rm = T),min=0,max=chem_fra$LOD[k]) 
}
for(k in 1:ncol(expo_gs)){
  expo_gs[which(expo_gs[,k]=="nd"),k] = rtruncnorm(sum(expo_gs[,k]=="nd", na.rm = T),min=0,max=chem_gs$LOD[k]) 
}

# Convert to numeric (Replace any characters with NaN and add rownames)
expo_lux = apply(expo_lux, 2, as.numeric)
rownames(expo_lux) = rownames(covars_lux)

expo_fra = apply(expo_fra, 2, as.numeric)
rownames(expo_fra) = rownames(covars_fra)

expo_gs = apply(expo_gs, 2, as.numeric)
rownames(expo_gs) = rownames(covars_gs)

# Check for NA prop
summary(chem_lux$NA_prop)
summary(chem_fra$NA_prop)
summary(chem_gs$NA_prop)

# Save data sets
saveRDS(expo_lux, paste0("../Processed/",filepaths[1],"/Exposure_matrix_ndimp_thresh.rds"))
saveRDS(expo_fra,  paste0("../Processed/",filepaths[2],"/Exposure_matrix_ndimp_thresh.rds"))
saveRDS(expo_gs,  paste0("../Processed/",filepaths[3],"/Exposure_matrix_ndimp_thresh.rds"))

