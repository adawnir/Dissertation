## Table 2s
## Rin Wada 12 July

# Load packages
library(tidyverse)

### Table 2: Pooled3----
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
expo = readRDS(paste0("../Processed/",filepaths[4],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
chem = readRDS(paste0("../Processed/",filepaths[4],"/Chemical_compound_info_thresh.rds"))
covars = readRDS(paste0("../Processed/",filepaths[4],"/Participant_covariate_info_thresh.rds"))
annot_sub = annot[colnames(expo)]
mytable=NULL
for (k in 1:length(colnames(expo))){
  print(colnames(expo)[k])
  tmp=ContinuousTest3(x=expo[,k], y=covars$Batch)
  missing = formatC(chem$NA_prop[k], format="f", digits=2)
  tmp = c(tmp,missing)
  if(!duplicated(annot_sub)[k]){
    mytable=rbind(mytable, rep(NA,length(tmp)))
    rownames(mytable)[nrow(mytable)]=as.character(annot_sub[k])
  }
  mytable=rbind(mytable, tmp)
  rownames(mytable)[nrow(mytable)]=colnames(expo)[k]
}
colnames(mytable)=c(rep("Mean (sd)",3), "P-value", "% missing")
# mytable=ReformatScientificNotation(mytable)
SaveExcelWithSuperscripts(cbind(rownames(mytable),mytable), paste0("../Exports/",filepaths[4],"/Table2.xlsx"))

### Table 2: Pooled2----
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
expo = readRDS(paste0("../Processed/",filepaths[5],"/Exposure_matrix_ndimp_thresh_log_naimp.rds"))
chem = readRDS(paste0("../Processed/",filepaths[5],"/Chemical_compound_info_thresh.rds"))
covars = readRDS(paste0("../Processed/",filepaths[5],"/Participant_covariate_info_thresh.rds"))
annot_sub = annot[colnames(expo)]

mytable=NULL
for (k in 1:length(colnames(expo))){
  print(colnames(expo)[k])
  tmp=ContinuousTest2(x=expo[,k], y=covars$Batch)
  missing = formatC(chem$NA_prop[k], format="f", digits=2)
  tmp = c(tmp,missing)
  if(!duplicated(annot_sub)[k]){
    mytable=rbind(mytable, rep(NA,length(tmp)))
    rownames(mytable)[nrow(mytable)]=as.character(annot_sub[k])
  }
  mytable=rbind(mytable, tmp)
  rownames(mytable)[nrow(mytable)]=colnames(expo)[k]
}
colnames(mytable)=c(rep("Mean (sd)",2), "P-value", "% missing")
# mytable=ReformatScientificNotation(mytable)
SaveExcelWithSuperscripts(cbind(rownames(mytable),mytable), paste0("../Exports/",filepaths[5],"/Table2.xlsx"))
