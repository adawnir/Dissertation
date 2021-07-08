## Data exploration and aata cleaning
## Rin Wada 24 May

# Load packages
library(openxlsx)
library(tidyverse)

### Prepare data (LUX) ----
# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
mydata = read.xlsx("../Data/Data_sample.xlsx")
filepath = filepaths[1]

## Covariate information
covars = mydata[1:39,1:7]
colnames(covars) = c("Family.ID", "Sibling.ID", "Indiv.ID", "Age", "Gender","Weight","Length")

## Gender
table(covars$Gender, useNA = "ifany")
# Recoding gender
covars$Gender[which(covars$Gender=="Jules")] = "Male"
table(covars$Gender, useNA = "ifany")

## Area
covar_lux$Area = "Luxembourg"
covar_lux$Department = "Luxembourg"
covar_lux$Region = "Luxembourg"
covar_lux$Country = "Luxembourg"
covar_lux$Batch = "LUX"

## Exposure matrix
expo = mydata[1:39,8:ncol(mydata)]
rownames(expo) = expo$Indiv.ID
expo = na_if(expo, "No result")

# Translate chemical names to English
colnames(expo) = readLines(paste0("../Dictionaries/Chemical_names_",filepath,".txt"))

## Chemical compound information
chem = data.frame(Compound = colnames(expo),
                  Family = unlist(mydata[42,8:ncol(mydata)]),
                  LOQ = as.numeric(gsub("\\*","",mydata[40,8:ncol(mydata)])),
                  Extraction = unlist(mydata[41,8:ncol(mydata)]))
chem = chem %>% fill(Family)

# Translate chemical family names from French to English
table(chem$Family, useNA = "ifany")
chem$Family[which(chem$Family=="PHENYLPYRAZOLES LC")] = "PHENYLPYRAZOLES"
names(chem_family) = unique(chem$Family)
chem$Family = factor(chem_family[chem$Family], levels = chem_family)
table(chem$Family, useNA = "ifany")

# Translate extraction method from French to English
table(chem$Extraction, useNA = "ifany")
chem$Extraction[which(chem$Extraction=="Inj. Liquide")] = "Liquid Injection"
table(chem$Extraction, useNA = "ifany")

# Semi-quantitative chemicals
mylist = c("DMP", "DMDTP", "DEP", "DEDTP", "ClCF3CA", "Propargite", "Spinosyn-A", "Bisphenol-A", "Bisphenol-S")
chem$Quantitative = ifelse(chem$Compound %in% mylist,"No","Yes")
table(chem$Quantitative, useNA = "ifany")

# Set rownames as compound name
rownames(chem) = chem$Compound
# Sort by family then compound name
chem = chem[order(chem$Family,chem$Compound),]
# Sort exposure matrix in same order
expo = expo[,chem$Compound]
# Check consistency
all(colnames(expo)==chem$Compound)

# Number of chemical compounds measured
nrow(chem) # Fipronil and Fipronil-sulfone are in duplicate

# Remove Fipronil and Fipronil-sulfonate extracted using SPME
chem = chem[-which(chem$Compound %in% c("Fipronil.SPME", "Fipronil-sulfone.SPME")),]
expo = expo[,-which(colnames(expo) %in% c("Fipronil.SPME", "Fipronil-sulfone.SPME"))]

# Proportion of NAs per chemical compounds
all(colnames(expo)==chem$Compound)
chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))

# Proportion of detected per chemical compounds
all(colnames(expo)==chem$Compound)
chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))

# Minimum detection per chemical compound
all(colnames(expo)==chem$Compound)
chem$LOD = apply(expo, 2, function(x) min(as.numeric(unlist(x)), na.rm = T))
# Check for infinite values
summary(chem$nd_prop[is.infinite(chem$LOD)]) # never-detected compounds have infinite values
# Replace LOD of never-detected compounds with NA
chem$LOD[which(is.infinite(chem$LOD))] = NA

## Save chemical family annotation as named list
annot = chem$Family
names(annot) = chem$Compound
ifelse(dir.exists("../Data"),"",dir.create("../Data"))
saveRDS(annot, "../Data/Chemical_compound_family_annotation.rds")

ifelse(dir.exists(paste0("../Data/",filepath)),"",dir.create(paste0("../Data/",filepath)))
saveRDS(covars, paste0("../Data/",filepath,"/Participant_covariate_info.rds"))
saveRDS(chem, paste0("../Data/",filepath,"/Chemical_compound_info.rds"))
saveRDS(expo, paste0("../Data/",filepath,"/Exposure_matrix_raw.rds"))

### Prepare data (FRA) ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
mydata = read.xlsx("../Data/Raw_Data_French_Children.xlsx")
filepath = filepaths[2]

## Covariate information
covars = mydata[1:142,1:5]
colnames(covars) = c("Siblings.Groups", "Indiv.ID", "Age", "Gender","Area")
covars$Age = as.numeric(covars$Age)

## Gender
table(covars$Gender, useNA = "ifany")
# Recoding gender
covars$Gender = ifelse(covars$Gender=="M", "Male", "Female")
table(covars$Gender, useNA = "ifany")

## Area
table(covars$Area, useNA = "ifany")

area = read.csv("../Dictionaries/French_area_codes.csv")
covars$Region = area$Region[which(covars$Area %in% area$Code.Commune)]
covars$Department = area$Department[which(covars$Area %in% area$Code.Commune)]
covars$Country = "France"
covars$Batch = "FRA"

## Exposure matrix
expo = mydata[1:142,6:ncol(mydata)]
rownames(expo) = covars$Indiv.ID
expo = na_if(expo, "No result") %>% na_if("No result (interf.)")

# Translate chemical names to English
colnames(expo) = readLines(paste0("../Dictionaries/Chemical_names_",filepath,".txt"))

## Chemical compound information
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
chem = data.frame(Compound = colnames(expo),
                  Family = annot[colnames(expo)],
                  LOQ = as.numeric(gsub("\\*","",mydata[143,6:ncol(mydata)])),
                  Extraction = unlist(mydata[144,6:ncol(mydata)]))

# Translate chemical family names from French to English
table(chem$Family, useNA = "ifany")
chem$Family = factor(chem_family[chem$Family], levels = chem_family)
table(chem$Family, useNA = "ifany")

# Translate extraction method from French to English
table(chem$Extraction, useNA = "ifany")
chem$Extraction[which(chem$Extraction=="Inj. Liquide")] = "Liquid Injection"
table(chem$Extraction, useNA = "ifany")

# Set rownames as compound name
rownames(chem) = chem$Compound
# Sort by family then compound name
chem = chem[order(chem$Family,chem$Compound),]
# Sort exposure matrix in same order
expo = expo[,chem$Compound]
# Check consistency
all(colnames(expo)==chem$Compound)

# Number of chemical compounds measured
nrow(chem)

# Proportion of NAs per chemical compounds
all(colnames(expo)==chem$Compound)
chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))

# Proportion of detected per chemical compounds
all(colnames(expo)==chem$Compound)
chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))

# Minimum detection per chemical compound
all(colnames(expo)==chem$Compound)
chem$LOD = apply(expo, 2, function(x) min(as.numeric(unlist(x)), na.rm = T))
# Check for infinite values
summary(chem$nd_prop[is.infinite(chem$LOD)]) # never-detected compounds have infinite values
# Replace LOD of never-detected compounds with NA
chem$LOD[which(is.infinite(chem$LOD))] = NA

ifelse(dir.exists(paste0("../Data/",filepath)),"",dir.create(paste0("../Data/",filepath)))
saveRDS(covars, paste0("../Data/",filepath,"/Participant_covariate_info.rds"))
saveRDS(chem, paste0("../Data/",filepath,"/Chemical_compound_info.rds"))
saveRDS(expo, paste0("../Data/",filepath,"/Exposure_matrix_raw.rds"))

### Prepare data (GS) ----
# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("graph_param.R")

# Load data
mydata = read.xlsx("../Data/Raw_Data_Grande-Synthe_Children.xlsx")
filepath = filepaths[3]

## Covariate information
covars = mydata[1:44,1:4]
colnames(covars) = c("Siblings.Groups", "Indiv.ID", "Age", "Gender")

# Recoding family ID
table(covars$Siblings.Groups, useNA = "ifany")
covars$Siblings.Groups = gsub("Famiy", "Family",covars$Siblings.Groups)
covars$Siblings.Groups = gsub("Isolate$", "Isolated",covars$Siblings.Groups)

## Gender
table(covars$Gender, useNA = "ifany")
# Recoding gender
covars$Gender = ifelse(covars$Gender=="M", "Male", ifelse(covars$Gender=="H", "Male", "Female"))
table(covars$Gender, useNA = "ifany")

## Area
covars$Area = "Grande-Synthe"
area = read.csv("../Dictionaries/French_area_codes.csv")
covars$Region = area$Region[which(covars$Area %in% area$Code.Commune)]
covars$Department = area$Department[which(covars$Area %in% area$Code.Commune)]
covars$Country = "France"
covars$Batch = "GS"

## Exposure matrix
expo = mydata[1:44,5:ncol(mydata)]
rownames(expo) = covars$Indiv.ID
expo = na_if(expo, "No result") %>% na_if("No result (interf.)")

# Translate chemical names to English
colnames(expo) = readLines(paste0("../Dictionaries/Chemical_names_",filepath,".txt"))

## Chemical compound information
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")
chem = data.frame(Compound = colnames(expo),
                  Family = unlist(mydata[47,5:ncol(mydata)]),
                  LOQ = as.numeric(gsub("\\*","",mydata[45,5:ncol(mydata)])),
                  Extraction = unlist(mydata[46,5:ncol(mydata)]))
chem = chem %>% fill(Family)

# Translate chemical family names from French to English
table(chem$Family, useNA = "ifany")
chem$Family[which(chem$Family=="PHENYLPYRAZOLES LC")] = "PHENYLPYRAZOLES"
names(chem_family) = unique(chem$Family)
chem$Family = factor(chem_family[chem$Family], levels = chem_family)
table(chem$Family, useNA = "ifany")

# Translate extraction method from French to English
table(chem$Extraction, useNA = "ifany")
chem$Extraction[which(chem$Extraction=="Inj. Liquide")] = "Liquid Injection"
table(chem$Extraction, useNA = "ifany")

# Set rownames as compound name
rownames(chem) = chem$Compound
# Sort by family then compound name
chem = chem[order(chem$Family,chem$Compound),]
# Sort exposure matrix in same order
expo = expo[,chem$Compound]
# Check consistency
all(colnames(expo)==chem$Compound)

# Number of chemical compounds measured
nrow(chem)

# Proportion of NAs per chemical compounds
all(colnames(expo)==chem$Compound)
chem$NA_prop = apply(expo, 2, function(x) sum(is.na(x))/nrow(expo))

# Proportion of detected per chemical compounds
all(colnames(expo)==chem$Compound)
chem$nd_prop = apply(expo, 2, function(x) sum(x=="nd", na.rm = TRUE)/nrow(expo))

# Minimum detection per chemical compound
all(colnames(expo)==chem$Compound)
chem$LOD = apply(expo, 2, function(x) min(as.numeric(unlist(x)), na.rm = T))
# Check for infinite values
summary(chem$nd_prop[is.infinite(chem$LOD)]) # never-detected compounds have infinite values
# Replace LOD of never-detected compounds with NA
chem$LOD[which(is.infinite(chem$LOD))] = NA

ifelse(dir.exists(paste0("../Data/",filepath)),"",dir.create(paste0("../Data/",filepath)))
saveRDS(covars, paste0("../Data/",filepath,"/Participant_covariate_info.rds"))
saveRDS(chem, paste0("../Data/",filepath,"/Chemical_compound_info.rds"))
saveRDS(expo, paste0("../Data/",filepath,"/Exposure_matrix_raw.rds"))
