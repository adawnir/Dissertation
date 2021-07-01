## Data exploration and aata cleaning
## Rin Wada 24 May

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

source("function.R")

# Load packages
library(openxlsx)
library(tidyverse)

# Load data
lux = read.xlsx("../Data/Data_sample.xlsx")
fra = read.xlsx("../Data/Raw_Data_French_Children.xlsx")
gs = read.xlsx("../Data/Raw_Data_Grande-Synthe_Children.xlsm")

### Covariate information----
# Luxembourg children
covar_lux = lux[1:39,1:7]
colnames(covar_lux) = c("Family.ID", "Sibling.ID", "Indiv.ID", "Age", "Gender",
                  "Weight (mg)","Length (cm)")

# French children
covar_fra = fra[1:142,1:5]
colnames(covar_fra) = c("Siblings.Groups", "Indiv.ID", "Age", "Gender","Area")
covar_fra$Age = as.numeric(covar_fra$Age)

covar_gs = gs[1:44,1:4]
colnames(covar_gs) = c("Siblings.Groups", "Indiv.ID", "Age", "Gender")
covar_gs$Siblings.Groups = gsub("Famiy", "Family",covar_gs$Siblings.Groups)
covar_gs$Siblings.Groups = gsub("Isolate$", "Isolated",covar_gs$Siblings.Groups)

### Hair sample length and weight ----
pdf("../Figures/Sample_length_weight.pdf", width=8, height=4)
par(mfrow = c(1,2))
hist(covar_lux$`Length (cm)`, xlab = "Length (cm)", main = "Length of hair sample (Luxembourg)")
hist(covar_lux$`Weight (mg)`, xlab = "Weight (mg)", main = "Weight of hair sample (Luxembourg)")
dev.off()

### Family sibling ID----
# Descriptive
table(covar_lux$Family.ID)
table(table(covar_lux$Family.ID))

table(covar_fra$Siblings.Groups)
table(table(covar_fra$Siblings.Groups))

table(covar_gs$Siblings.Groups)
table(table(covar_gs$Siblings.Groups))

### Gender ----
table(covar_lux$Gender, useNA = "ifany")
table(covar_fra$Gender, useNA = "ifany")
table(covar_gs$Gender, useNA = "ifany")

# Recoding gender
covar_lux$Gender[which(covar_lux$Gender=="Jules")] = "Male"

covar_fra$Gender[which(covar_fra$Gender=="M")] = "Male"
covar_fra$Gender[which(covar_fra$Gender=="F")] = "Female"

covar_gs$Gender[which(covar_gs$Gender=="M")] = "Male"
covar_gs$Gender[which(covar_gs$Gender=="H")] = "Male"
covar_gs$Gender[which(covar_gs$Gender=="F")] = "Female"

# Missing gender information
cat(sum(is.na(covar_lux$Gender)),
    " (", round(sum(is.na(covar_lux$Gender))/length(covar_lux$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_fra$Gender)),
    " (", round(sum(is.na(covar_fra$Gender))/length(covar_fra$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_gs$Gender)),
    " (", round(sum(is.na(covar_gs$Gender))/length(covar_gs$Gender)*100,2), "%)\n", sep = "")

### Age ----
# Age distribution
pdf("../Figures/Age_dist.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
hist(covar_lux$Age, main = "Luxembourg", xlab = "Age (years)")
hist(covar_fra$Age, main = "France", xlab = "Age (years)")
hist(covar_gs$Age, main = "Grande-Synthe", xlab = "Age (years)")
dev.off()

# Missing age information
cat(sum(is.na(covar_lux$Age)),
    " (", round(sum(is.na(covar_lux$Age))/length(covar_lux$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_fra$Age)),
    " (", round(sum(is.na(covar_fra$Age))/length(covar_fra$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covar_gs$Age)),
    " (", round(sum(is.na(covar_gs$Age))/length(covar_gs$Age)*100,2), "%)\n", sep = "")

### Chemical compound matrix ----
mat_lux = lux[1:39,8:ncol(lux)]
rownames(mat_lux) = covar_lux$Indiv.ID

mat_fra = fra[1:142,6:ncol(fra)]
rownames(mat_fra) = covar_fra$Indiv.ID

mat_gs = gs[1:44,5:ncol(gs)]
rownames(mat_gs) = covar_gs$Indiv.ID

# Replace no results (due to technical issues) with NA
mat_lux = na_if(mat_lux, "No result")
mat_fra = na_if(mat_fra, "No result")
mat_fra = na_if(mat_fra, "No result (interf.)")
mat_gs = na_if(mat_gs, "No result")

# Translate chemical names to English
colnames(mat_lux) = readLines("../Dictionaries/Chemical_names_Luxembourg.txt")
colnames(mat_fra) = readLines("../Dictionaries/Chemical_names_France.txt")
colnames(mat_gs) = readLines("../Dictionaries/Chemical_names_GrandeSynthe.txt")

setdiff(colnames(mat_lux),colnames(mat_gs))
setdiff(colnames(mat_fra),colnames(mat_lux))
setdiff(colnames(mat_gs),colnames(mat_lux))

### Chemical group, extraction method and LOQ----
## Luxembourg
chem_lux = data.frame(Compound = colnames(mat_lux),
                      Family = str_to_title(unlist(lux[42,8:ncol(lux)])),
                      LOQ = as.numeric(gsub("\\*","",lux[40,8:ncol(lux)])),
                      Extraction = unlist(lux[41,8:ncol(lux)]))
chem_lux = chem_lux %>% fill(Family)

# Translate chemical family names from French to English
unique(str_to_title(chem_lux$Family))
chem_lux$Family[which(chem_lux$Family=="Organochlores")] = "Organochlorines"
chem_lux$Family[which(chem_lux$Family=="Organophosphores Et Metabolites")] = "Organophosphate metabolites"
chem_lux$Family[which(chem_lux$Family=="Pyrethrinoides")] = "Pyrethroids"
chem_lux$Family[which(chem_lux$Family=="Metabolites De Pyrethrinoides")] = "Pyrethroid metabolites"
chem_lux$Family[which(chem_lux$Family=="Pcbs")] = "PCBs"
chem_lux$Family[which(chem_lux$Family=="Pbdes")] = "PBDEs"
chem_lux$Family[which(chem_lux$Family=="Herbicides Acides")] = "Acidic herbicides"
chem_lux$Family[which(chem_lux$Family=="Anilino-Pyrimidines")] = "Anilinopyrimidines"
chem_lux$Family[which(chem_lux$Family=="Neonicotinoides")] = "Neonicotinoids"
chem_lux$Family[which(chem_lux$Family=="Phenylpyrazoles Lc")] = "Phenylpyrazoles"
chem_lux$Family[which(chem_lux$Family=="Strobilurines")] = "Strobilurins"
chem_lux$Family[which(chem_lux$Family=="Triazines / Triazinones / Diazines")] = "Triazines/Triazinones/Diazines"
chem_lux$Family[which(chem_lux$Family=="Urees Substituees")] = "Subtituted ureas"
chem_lux$Family[which(chem_lux$Family=="Divers")] = "Other"
unique(chem_lux$Family)

# Translate extraction method from French to English
chem_lux$Extraction[which(chem_lux$Extraction=="Inj. Liquide")] = "Liquid Injection"

# Semi-quantitative chemicals
mylist = c("DMP", "DMDTP", "DEP", "DEDTP", "ClCF3CA", "Propargite", "Spinosyn-A", "Bisphenol-A", "Bisphenol-S")
chem_lux$Quantitative = ifelse(chem_lux$Compound %in% mylist,"No","Yes")
table(chem_lux$Quantitative)

# Set rownames as compound name
rownames(chem_lux) = chem_lux$Compound

## France
chem_fra = data.frame(Compound = colnames(mat_fra),
                      Family = chem_lux[colnames(mat_fra),"Family"],
                      LOQ = as.numeric(gsub("\\*","",fra[143,6:ncol(fra)])),
                      Extraction = unlist(fra[144,6:ncol(fra)]))

# Translate extraction method from French to English
chem_fra$Extraction[which(chem_fra$Extraction=="Inj. Liquide")] = "Liquid Injection"

# Set rownames as compound name
rownames(chem_fra) = chem_fra$Compound

## Grand-Synthe
chem_gs = data.frame(Compound = colnames(mat_gs),
                     Family = str_to_title(unlist(gs[47,5:ncol(gs)])),
                     LOQ = as.numeric(gsub("\\*","",gs[45,5:ncol(gs)])),
                     Extraction = unlist(gs[46,5:ncol(gs)]))
chem_gs = chem_gs %>% fill(Family)

# Translate chemical family names from French to English
unique(str_to_title(chem_gs$Family))
chem_gs$Family[which(chem_gs$Family=="Organochlores")] = "Organochlorines"
chem_gs$Family[which(chem_gs$Family=="Organophosphores Et Metabolites")] = "Organophosphate metabolites"
chem_gs$Family[which(chem_gs$Family=="Pyrethrinoides")] = "Pyrethroids"
chem_gs$Family[which(chem_gs$Family=="Metabolites De Pyrethrinoides")] = "Pyrethroid metabolites"
chem_gs$Family[which(chem_gs$Family=="Pcbs")] = "PCBs"
chem_gs$Family[which(chem_gs$Family=="Pbdes")] = "PBDEs"
chem_gs$Family[which(chem_gs$Family=="Herbicides Acides")] = "Acidic herbicides"
chem_gs$Family[which(chem_gs$Family=="Anilino-Pyrimidines")] = "Anilinopyrimidines"
chem_gs$Family[which(chem_gs$Family=="Neonicotinoides")] = "Neonicotinoids"
chem_gs$Family[which(chem_gs$Family=="Phenylpyrazoles Lc")] = "Phenylpyrazoles"
chem_gs$Family[which(chem_gs$Family=="Strobilurines")] = "Strobilurins"
chem_gs$Family[which(chem_gs$Family=="Triazines / Triazinones / Diazines")] = "Triazines/Triazinones/Diazines"
chem_gs$Family[which(chem_gs$Family=="Urees Substituees")] = "Subtituted ureas"
chem_gs$Family[which(chem_gs$Family=="Divers")] = "Other"
unique(chem_gs$Family)

# Translate extraction method from French to English
chem_gs$Extraction[which(chem_gs$Extraction=="Inj. Liquide")] = "Liquid Injection"

# Set rownames as compound name
rownames(chem_gs) = chem_gs$Compound

# Check chemical family consistency
all(chem_lux[rownames(chem_gs),"Family"]==chem_gs$Family)

# Set chemical family as factor
chem_lux$Family = factor(chem_lux$Family, levels = unique(chem_lux$Family))
levels(chem_lux$Family)
chem_fra$Family = factor(chem_fra$Family,
                             levels = intersect(unique(chem_lux$Family),unique(chem_fra$Family)))
levels(chem_fra$Family)
chem_gs$Family = factor(chem_gs$Family, levels = unique(chem_lux$Family))
levels(chem_gs$Family)

# Sort by family then compound name
chem_lux = chem_lux[order(chem_lux$Family,chem_lux$Compound),]
chem_fra = chem_fra[order(chem_fra$Family,chem_fra$Compound),]
chem_gs = chem_gs[order(chem_gs$Family,chem_gs$Compound),]

# Sort matrix in same order
mat_lux = mat_lux[,chem_lux$Compound]
mat_fra = mat_fra[,chem_fra$Compound]
mat_gs = mat_gs[,chem_gs$Compound]

# Minimum detection per chemical compound
all(colnames(mat_lux)==chem_lux$Compound)
all(colnames(mat_fra)==chem_fra$Compound)
all(colnames(mat_gs)==chem_gs$Compound)

chem_lux$LOD = apply(mat_lux, 2, function(x) min(as.numeric(unlist(x)), na.rm = T))
chem_fra$LOD = apply(mat_fra, 2, function(x) min(as.numeric(unlist(x)), na.rm = T))
chem_gs$LOD = apply(mat_gs, 2, function(x) min(as.numeric(unlist(x)), na.rm = T))

chem_lux$LOD = ifelse(is.infinite(chem_lux$LOD),0,chem_lux$LOD)
chem_fra$LOD = ifelse(is.infinite(chem_fra$LOD),0,chem_fra$LOD)
chem_gs$LOD = ifelse(is.infinite(chem_gs$LOD),0,chem_gs$LOD)

ifelse(dir.exists("../Processed"),"",dir.create("../Processed"))
ifelse(dir.exists("../Processed/Luxembourg"),"",dir.create("../Processed/Luxembourg"))
saveRDS(covar_lux, "../Processed/Luxembourg/Participant_covariate_info.rds")
saveRDS(chem_lux, "../Processed/Luxembourg/Chemical_compound_info.rds")
saveRDS(mat_lux, "../Processed/Luxembourg//Chemical_compound_matrix_raw.rds")

ifelse(dir.exists("../Processed/France"),"",dir.create("../Processed/France"))
saveRDS(covar_fra, "../Processed/France/Participant_covariate_info.rds")
saveRDS(chem_fra, "../Processed/France/Chemical_compound_info.rds")
saveRDS(mat_fra, "../Processed/France//Chemical_compound_matrix_raw.rds")

ifelse(dir.exists("../Processed/GrandeSynthe"),"",dir.create("../Processed/GrandeSynthe"))
saveRDS(covar_gs, "../Processed/GrandeSynthe/Participant_covariate_info.rds")
saveRDS(chem_gs, "../Processed/GrandeSynthe/Chemical_compound_info.rds")
saveRDS(mat_gs, "../Processed/GrandeSynthe//Chemical_compound_matrix_raw.rds")
