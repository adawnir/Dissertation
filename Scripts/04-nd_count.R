## Non-detects count analysis
## Rin Wada 15 June

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(tidyverse)

# Load data sets
# Save data sets
covar_lux = readRDS("../Processed/Luxembourg/Participant_covariate_info_subset.rds")
covar_fra = readRDS("../Processed/France/Participant_covariate_info_subset.rds")
covar_gs = readRDS("../Processed/GrandeSynthe/Participant_covariate_info_subset.rds")

chem_lux = readRDS("../Processed/Luxembourg/Chemical_compound_info_subset.rds")
chem_fra = readRDS("../Processed/France/Chemical_compound_info_subset.rds")
chem_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_info_subset.rds")

mat_lux = readRDS("../Processed/Luxembourg/Chemical_compound_matrix_raw_subset.rds")
mat_fra = readRDS("../Processed/France/Chemical_compound_matrix_raw_subset.rds")
mat_gs = readRDS("../Processed/GrandeSynthe/Chemical_compound_matrix_raw_subset.rds")


### Detection frequency per participant ----
# Count non-detects for each participant
covar_lux$nd_count = apply(mat_lux, 1, function(x) sum(x=="nd", na.rm = T))
covar_fra$nd_count = apply(mat_fra, 1, function(x) sum(x=="nd", na.rm = T))
covar_gs$nd_count = apply(mat_gs, 1, function(x) sum(x=="nd", na.rm = T))

# Association with gender
mean(covar_lux$nd_count)
var(covar_lux$nd_count)
summary(covar_lux$nd_count)
glm1 <- glm(nd_count ~ Gender, poisson, covar_lux)
summary(glm1)

mean(covar_fra$nd_count)
var(covar_fra$nd_count)
summary(covar_fra$nd_count)
glm2 <- glm(nd_count ~ Gender, poisson, covar_fra)
summary(glm2)

mean(covar_gs$nd_count)
var(covar_gs$nd_count)
summary(covar_gs$nd_count)
glm3 <- glm(nd_count ~ Gender, poisson, covar_gs)
summary(glm3)

pdf("../Figures/Nd_count_Gender.pdf", width=12, height=4)
par(mfrow = c(1,3))
boxplot(nd_count ~ Gender, data = covar_lux, ylab = "Number of non-detects", main = "Luxembourg")
boxplot(nd_count ~ Gender, data = covar_fra, ylab = "Number of non-detects", main = "France")
boxplot(nd_count ~ Gender, data = covar_gs, ylab = "Number of non-detects", main = "Grande-Synthe")
dev.off()

# Association with age
glm1 <- glm(nd_count ~ Age, poisson, covar_lux)
summary(glm1)

glm2 <- glm(nd_count ~ Age, poisson, covar_fra)
summary(glm2)

glm3 <- glm(nd_count ~ Age, poisson, covar_gs)
summary(glm3)

pdf("../Figures/Nd_count_Age.pdf", width=12, height=4)
par(mfrow = c(1,3))
plot(nd_count ~ Age, data = covar_lux, pch = 19, ylab = "Number of non-detects", main = "Luxembourg")
plot(nd_count ~ Age, data = covar_fra, pch = 19, ylab = "Number of non-detects", main = "France")
plot(nd_count ~ Age, data = covar_gs, pch = 19, ylab = "Number of non-detects", main = "Grande-Synthe")
dev.off()

# Association with length of hair sample
summary(covar_lux$`Length (cm)`)
glm1 <- glm(nd_count ~ `Length (cm)`, poisson, covar_lux)
summary(glm1)

# Association with weight of hair sample
summary(covar_lux$`Weight (mg)`)
glm1 <- glm(nd_count ~ `Weight (mg)`, poisson, covar_lux)
summary(glm1)

pdf("../Figures/Nd_count_Sample.pdf", width=8, height=4)
par(mfrow = c(1,2))
plot(nd_count ~ `Length (cm)`, data = covar_lux, pch = 19, ylab = "Number of non-detects", main = "Length of hair sample (Luxembourg)")
plot(nd_count ~ `Weight (mg)`, data = covar_lux, pch = 19, ylab = "Number of non-detects", main = "Weight of hair sample (Luxembourg)")
dev.off()
