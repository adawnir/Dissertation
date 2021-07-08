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

### Detection of chemical compound vs family id ----
# Subset matrix families with siblings only
nd_lux_sub = nd_lux[which(covar_lux$Family.ID!="Isolated"),]
nd_fra_sub = nd_fra[which(covar_fra$Family.ID!="Isolated"),]
nd_gs_sub = nd_gs[which(covar_gs$Family.ID!="Isolated"),]

# Luxembourg
pval_lux_sub = NULL
f1='nd_lux[,i] ~ covar_lux$Family.ID'
f0='nd_lux[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(nd_lux_sub)){
  model1=glm(as.formula(f1), family = "binomial")
  model0=glm(as.formula(f0), family = "binomial")
  pval_lux_sub=c(pval_lux_sub, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pval_lux_sub) = colnames(nd_lux)

# France
pval_fra_sub = NULL
f1='nd_fra[,i] ~ covar_fra$Family.ID'
f0='nd_fra[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(nd_fra_sub)){
  model1=glm(as.formula(f1), family = "binomial")
  model0=glm(as.formula(f0), family = "binomial")
  pval_fra_sub=c(pval_fra_sub, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pval_fra_sub) = colnames(nd_fra)

# Grande-Synthe
pval_gs_sub = NULL
f1='nd_gs[,i] ~ covar_gs$Family.ID'
f0='nd_gs[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(nd_gs_sub)){
  model1=glm(as.formula(f1), family = "binomial")
  model0=glm(as.formula(f0), family = "binomial")
  pval_gs_sub=c(pval_gs_sub, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pval_gs_sub) = colnames(nd_gs)

saveRDS(pval_gs_sub, "../Results/nd_pval_family.rds")

# Manhattan plot
mycolours = c("tomato","royalblue","forestgreen")
values = cbind(pval_lux_sub, pval_fra_sub[match(names(pval_lux_sub), names(pval_fra_sub))]) %>%
  cbind(., pval_gs_sub[match(rownames(.), names(pval_gs_sub))])
values = -log10(values)

bonf = c(0.05/length(pval_lux_sub), 0.05/length(pval_fra_sub), 0.05/length(pval_gs_sub))
bonf = -log10(bonf)

xseq = seq(1, nrow(values))
{pdf("../Figures/Nd_Family.pdf", width=14, height=7)
  par(mar=c(10,3,3,3))
  plot(values[,1],
       col=mycolours[1], cex = 0.7,
       xaxt="n", yaxt="n", ylab="", xlab = "",
       type="n", ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)))
  axis(4, axTicks(2))
  abline(h = -log10(0.05), lty = 2, col = "grey")
  abline(h = bonf[1], lty = 2, col = darken(mycolours[1], 0.5))
  abline(h = bonf[2], lty = 2, col = darken(mycolours[2], 0.5))
  abline(h = bonf[3], lty = 2, col = darken(mycolours[3], 0.5))
  abline(v = xseq, lty = 3, col = "grey", lwd = 0.7)
  mtext(side = 4, text = expression(-log[10](italic(p))), line = 2)
  points(values[,1], pch = 17, col = mycolours[1], cex = 0.7)
  points(values[,2], pch = 19, col = mycolours[2], cex = 0.7)
  points(values[,3], pch = 15, col = mycolours[3], cex = 0.7)
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.7)}
  legend("topright", pch=c(17,19,15, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
         col=c(mycolours,darken(mycolours, 0.5), "grey"),
         legend = c("Luxembourg","France","Grande-Synthe",
                    "Bonferroni threshold (LUX)", "Bonferroni threshold (FRA)", "Bonferroni threshold (GS)",
                    "Nominal threshold"),
         cex=0.7, bg="white")
  dev.off()
}




