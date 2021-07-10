## Table 1s
## Rin Wada 8 July

# Load packages
library(tidyverse)

# Initialise
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load custom functions and parameters
source("functions.R")
source("table_functions.R")
source("graph_param.R")

# Load data
annot = readRDS("../Data/Chemical_compound_family_annotation.rds")

covars_lux = readRDS(paste0("../Processed/",filepaths[1],"/Participant_covariate_info_thresh.rds"))
chem_lux = readRDS(paste0("../Processed/",filepaths[1],"/Chemical_compound_info_thresh.rds"))
expo_lux = readRDS(paste0("../Processed/",filepaths[1],"/Exposure_matrix_ndimp_thresh.rds"))

covars_fra = readRDS(paste0("../Processed/",filepaths[2],"/Participant_covariate_info_thresh.rds"))
chem_fra = readRDS(paste0("../Processed/",filepaths[2],"/Chemical_compound_info_thresh.rds"))
expo_fra = readRDS(paste0("../Processed/",filepaths[2],"/Exposure_matrix_ndimp_thresh.rds"))

covars_gs = readRDS(paste0("../Processed/",filepaths[3],"/Participant_covariate_info_thresh.rds"))
chem_gs = readRDS(paste0("../Processed/",filepaths[3],"/Chemical_compound_info_thresh.rds"))
expo_gs = readRDS(paste0("../Processed/",filepaths[3],"/Exposure_matrix_ndimp_thresh.rds"))

covars_pooled3 = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))
chem_pooled3 = readRDS(paste0("../Processed/",filepaths[4],"/Chemical_compound_info_thresh.rds"))
expo_pooled3 = readRDS(paste0("../Processed/",filepaths[4],"/Exposure_matrix_ndimp_thresh.rds"))

covars_pooled2 = readRDS(paste0("../Processed/",filepaths[5],"/Participant_covariate_info_thresh.rds"))
chem_pooled2 = readRDS(paste0("../Processed/",filepaths[5],"/Chemical_compound_info_thresh.rds"))
expo_pooled2 = readRDS(paste0("../Processed/",filepaths[5],"/Exposure_matrix_ndimp_thresh.rds"))

### LOD, Detection rate and concentration range ----
ifelse(dir.exists("../Exports"), "", dir.create("../Exports"))
suffix = c("lux","fra","gs","pooled3","pooled2")
for (i in 1:length(batches)){
  expo = eval(parse(text = paste0("expo_",suffix[i])))
  chem = eval(parse(text = paste0("chem_",suffix[i])))
  annot_sub = annot[colnames(expo)]
  all(colnames(expo)[k] == rownames(chem)[k])
  mytable=NULL
  for (k in 1:length(colnames(expo))){
    print(colnames(expo)[k])
    lod = formatC(chem$LOD[k], format="e", digits=2)
    detect_rate = formatC(1-(chem$nd_prop[k]+chem$NA_prop[k]), format="f", digits=2)
    q = quantile(expo[,k], na.rm = T)
    q = ifelse(sapply(q, function(x) round(x,2)) <= 0.1,
               formatC(q, format="e", digits=2),
               formatC(q, format="f", digits=2))
    if(i %in% c(4,5)){
      missing_rate = formatC(chem$NA_prop[k], format="f", digits=2)
      tmp=c(lod,detect_rate,missing_rate,q)
    } else {tmp=c(lod,detect_rate,q)}
    if(!duplicated(annot_sub)[k]){
      mytable=rbind(mytable, rep(NA,length(tmp)))
      rownames(mytable)[nrow(mytable)]=as.character(annot_sub[k])
    }
    mytable=rbind(mytable, tmp)
    rownames(mytable)[nrow(mytable)]=colnames(expo)[k]
  }
  if(i %in% c(4,5)){
    colnames(mytable)=c("LOD (pg/mg)","Detection rate","Missing rate","Min", "Q1","Q2","Q3","Max")
  } else {colnames(mytable)=c("LOD (pg/mg)","Detection rate","Min", "Q1","Q2","Q3","Max")}
  mytable[,-2]=ReformatScientificNotation(mytable[,-2])
  ifelse(dir.exists(paste0("../Exports/",filepaths[i])), "", dir.create(paste0("../Exports/",filepaths[i])))
  SaveExcelWithSuperscripts(cbind(rownames(mytable),mytable), paste0("../Exports/",filepaths[i],"/Table1_chemical_compound.xlsx"))
}

### Covariate----
# Total number of families
sum(covars_lux$Family.ID=="Isolated") + sum(table(table(covars_lux$Family.ID[which(covars_lux$Family.ID!="Isolated")])))
sum(covars_fra$Family.ID=="Isolated") + sum(table(table(covars_fra$Family.ID[which(covars_fra$Family.ID!="Isolated")])))
sum(covars_gs$Family.ID=="Isolated") + sum(table(table(covars_gs$Family.ID[which(covars_gs$Family.ID!="Isolated")])))

# Number of each type of family
table(table(covars_lux$Family.ID))
table(table(covars_fra$Family.ID))
table(table(covars_gs$Family.ID))

# Age (Mean, SD)
cat(round(mean(covars_lux$Age, na.rm = T),2),
    " (", round(sd(covars_lux$Age, na.rm = T),2), ")\n", sep = "")
cat(round(mean(covars_fra$Age, na.rm = T),2),
    " (", round(sd(covars_fra$Age, na.rm = T),2), ")\n", sep = "")
cat(round(mean(covars_gs$Age, na.rm = T),2),
    " (", round(sd(covars_gs$Age, na.rm = T),2), ")\n", sep = "")

# Age (NA %)
cat(sum(is.na(covars_lux$Age)),
    " (", round(sum(is.na(covars_lux$Age))/length(covars_lux$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_fra$Age)),
    " (", round(sum(is.na(covars_fra$Age))/length(covars_fra$Age)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_gs$Age)),
    " (", round(sum(is.na(covars_gs$Age))/length(covars_gs$Age)*100,2), "%)\n", sep = "")

# Gender (N, %)
table(covars_lux$Gender, useNA = "ifany")
round(prop.table(table(covars_lux$Gender, useNA = "ifany"))*100,2)

table(covars_fra$Gender, useNA = "ifany")
round(prop.table(table(covars_fra$Gender, useNA = "ifany"))*100,2)

table(covars_gs$Gender, useNA = "ifany")
round(prop.table(table(covars_gs$Gender, useNA = "ifany"))*100,2)

# Gender (NA %)
cat(sum(is.na(covars_lux$Gender)),
    " (", round(sum(is.na(covars_lux$Gender))/length(covars_lux$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_fra$Gender)),
    " (", round(sum(is.na(covars_fra$Gender))/length(covars_fra$Gender)*100,2), "%)\n", sep = "")
cat(sum(is.na(covars_gs$Gender)),
    " (", round(sum(is.na(covars_gs$Gender))/length(covars_gs$Gender)*100,2), "%)\n", sep = "")

# Table 1
# N
nrow(covars)
# Families & number of isolated
table(table(covars$Family.ID))
# Age
cat(round(mean(covars$Age, na.rm = T),2),
    " (", round(sd(covars$Age, na.rm = T),2), ")\n", sep = "")
cat(sum(is.na(covars$Age)),
    " (", round(sum(is.na(covars$Age))/length(covars$Age)*100,2), "%)\n", sep = "")
# Gender
table(covars$Gender, useNA = "ifany")
round(prop.table(table(covars$Gender, useNA = "ifany"))*100,2)
cat(sum(is.na(covars$Gender)),
    " (", round(sum(is.na(covars$Gender))/length(covars$Gender)*100,2), "%)\n", sep = "")

mytable=NULL
for (k in 1:length(colnames(mat_pooled3))){
  print(colnames(mat_pooled3)[k])
  tmp=ContinuousTest3(x=mat_pooled3[,k], y=covar_pooled3$Batch)
  mytable=rbind(mytable, tmp)
  rownames(mytable)[nrow(mytable)]=colnames(mat_pooled3)[k]
}
colnames(mytable)=c(rep("Mean (sd)",3), "LUX vs FRA", "LUX vs GS", "FRA vs GS")
mytable[,4:6]=ReformatScientificNotation(mytable[,4:6])
ifelse(dir.exists("../Exports"), "", dir.create("../Exports"))
ifelse(dir.exists("../Exports/Pooled3"), "", dir.create("../Exports/Pooled3"))
SaveExcelWithSuperscripts(cbind(rownames(mytable),mytable), "../Exports/Pooled3/Table1.xlsx")

