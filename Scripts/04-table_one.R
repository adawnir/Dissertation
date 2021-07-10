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
covars = covars_pooled3
mytable=NULL
out = sum(unique(covars$Family.ID)!="Isolated")
for (k in c("LUX","FRA","GS")){
  myn = sum(unique(covars$Family.ID[covars$Batch==k])!="Isolated")
  out=c(out, myn)
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Families with siblings, N"

myn = sum(covars$Family.ID=="Isolated")
myp = formatC(myn/length(covars$Family.ID)*100, format = "f", digits = 1)
out = paste0(myn, " (", myp, "%)")
for (k in c("LUX","FRA","GS")){
  myn = sum(covars$Family.ID[covars$Batch==k]=="Isolated")
  myp = formatC(myn/sum(covars$Batch==k), format = "f", digits = 1)
  out=c(out, paste0(myn, " (", myp, "%)"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Only child, N (%)"

mymean=formatC(mean(covars$Age, na.rm=TRUE), format="f", digits=2)
mysd=formatC(sd(covars$Age, na.rm=TRUE), format="f", digits=2)
out=paste0(mymean, " (", mysd, ")")
for (k in c("LUX","FRA","GS")){
  mymean=formatC(mean(covars$Age[covars$Batch==k], na.rm=TRUE), format="f", digits=2)
  mysd=formatC(sd(covars$Age[covars$Batch==k], na.rm=TRUE), format="f", digits=2)
  out=c(out, paste0(mymean, " (", mysd, ")"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Age (years), Mean (SD)"

myn = sum(is.na(covars$Age))
myp = formatC(myn/length(covars$Age)*100, format = "f", digits = 1)
out = paste0(myn, " (", myp, "%)")
for (k in c("LUX","FRA","GS")){
  myn = sum(is.na(covars$Age[covars$Batch==k]))
  myp = formatC(myn/sum(covars$Batch==k)*100, format = "f", digits = 1)
  out=c(out, paste0(myn, " (", myp, "%)"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Missing, N (%)"

mytable=rbind(mytable, rep(NA, ncol(mytable)))
rownames(mytable)[nrow(mytable)]="Gender, N (%)"

myn = sum(covars$Gender=="Female", na.rm = T)
myp = formatC(myn/length(covars$Gender)*100, format = "f", digits = 1)
out = paste0(myn, " (", myp, "%)")
for (k in c("LUX","FRA","GS")){
  myn = sum(covars$Gender[covars$Batch==k]=="Female", na.rm = T)
  myp = formatC(myn/sum(covars$Batch==k)*100, format = "f", digits = 1)
  out=c(out, paste0(myn, " (", myp, "%)"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Female"

myn = sum(covars$Gender=="Male", na.rm = T)
myp = formatC(myn/length(covars$Gender)*100, format = "f", digits = 1)
out = paste0(myn, " (", myp, "%)")
for (k in c("LUX","FRA","GS")){
  myn = sum(covars$Gender[covars$Batch==k]=="Male", na.rm = T)
  myp = formatC(myn/sum(covars$Batch==k)*100, format = "f", digits = 1)
  out=c(out, paste0(myn, " (", myp, "%)"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Male"

myn = sum(is.na(covars$Gender))
myp = formatC(myn/length(covars$Gender)*100, format = "f", digits = 1)
out = paste0(myn, " (", myp, "%)")
for (k in c("LUX","FRA","GS")){
  myn = sum(is.na(covars$Gender[covars$Batch==k]))
  myp = formatC(myn/sum(covars$Batch==k)*100, format = "f", digits = 1)
  out=c(out, paste0(myn, " (", myp, "%)"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Missing, N (%)"

mymean=formatC(mean(covars$Weight, na.rm=TRUE), format="f", digits=2)
mysd=formatC(sd(covars$Weight, na.rm=TRUE), format="f", digits=2)
out=paste0(mymean, " (", mysd, ")")
for (k in c("LUX","FRA","GS")){
  if (k = "LUX"){
    mymean=formatC(mean(covars$Age[covars$Batch==k], na.rm=TRUE), format="f", digits=2)
    mysd=formatC(sd(covars$Age[covars$Batch==k], na.rm=TRUE), format="f", digits=2)
    out=c(out, paste0(mymean, " (", mysd, ")"))
  } else {out = c(out, "-")}
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Sample weight (mg), Mean (SD)"

myn = sum(is.na(covars$Weight))
myp = formatC(myn/length(covars$Weight)*100, format = "f", digits = 1)
out = paste0(myn, " (", myp, "%)")
for (k in c("LUX","FRA","GS")){
  myn = sum(is.na(covars$Weight[covars$Batch==k]))
  myp =  formatC(myn/sum(covars$Batch==k)*100, format = "f", digits = 1)
  out=c(out, paste0(myn, " (", myp, "%)"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Missing, N (%)"

mymean=formatC(mean(covars$Length, na.rm=TRUE), format="f", digits=2)
mysd=formatC(sd(covars$Length, na.rm=TRUE), format="f", digits=2)
out=paste0(mymean, " (", mysd, ")")
for (k in c("LUX","FRA","GS")){
  if (k == "LUX"){
    mymean=formatC(mean(covars$Age[covars$Batch==k], na.rm=TRUE), format="f", digits=2)
    mysd=formatC(sd(covars$Age[covars$Batch==k], na.rm=TRUE), format="f", digits=2)
    out=c(out, paste0(mymean, " (", mysd, ")"))
  } else {out = c(out, "-")}
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Sample length (cm), Mean (SD)"

myn = sum(is.na(covars$Length))
myp = formatC(myn/length(covars$Length)*100, format = "f", digits = 1)
out = paste0(myn, " (", myp, "%)")
for (k in c("LUX","FRA","GS")){
  myn = sum(is.na(covars$Length[covars$Batch==k]))
  myp = formatC(myn/sum(covars$Batch==k)*100, format = "f", digits = 1)
  out=c(out, paste0(myn, " (", myp, "%)"))
}
mytable=rbind(mytable, out)
rownames(mytable)[nrow(mytable)]="Missing, N (%)"
out = nrow(covars)
for (k in c("LUX","FRA","GS")){
  out=c(out, sum(covars$Batch==k))
}
colnames(mytable)=paste0(c("Pooled","Luxembourg","France","Grande-Synthe"),
                         " (N=",out,")")
SaveExcelWithSuperscripts(cbind(rownames(mytable),mytable), "../Exports/Table1_covariates.xlsx")
