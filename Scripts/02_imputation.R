## Data imputation
## Rin Wada 24 May

# Initialisation
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Load packages
library(imputeLCMD)

# Load data
mat=readRDS("../Processed/Chemical_component_matrix_subset_transformed.rds")

# impute.QRILC
print("Imputing")
set.seed(7)
mat_imp=impute.QRILC(t(mat))
mat_imp=data.frame(t(mat_imp[[1]]))
print("Imputation completed!")
rownames(mat_imp) = rownames(mat)
colnames(mat_imp) = colnames(mat)

## Save
saveRDS(mat_imp, "../Processed/Chemical_component_matrix_subset_transformed_imputed.rds")







