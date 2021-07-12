### Remove isolated children (no siblings) ----
covar_lux = covar_lux[which(covar_lux$Family.ID != "Isolated"),]
mat_lux = mat_lux[covar_lux$Indiv.ID,]
nd_lux = nd_lux[covar_lux$Indiv.ID,]

covar_fra = covar_fra[which(covar_fra$Family.ID != "Isolated"),]
mat_fra = mat_fra[covar_fra$Indiv.ID,]
nd_fra = nd_fra[covar_fra$Indiv.ID,]

covar_gs = covar_gs[which(covar_gs$Family.ID != "Isolated"),]
mat_gs = mat_gs[covar_gs$Indiv.ID,]
nd_gs = nd_gs[covar_gs$Indiv.ID,]

saveRDS(covar_lux, "../Processed/Luxembourg/Participant_covariate_info_subset_no_isolated.rds")
saveRDS(covar_fra, "../Processed/France/Participant_covariate_info_subset_no_isolated.rds")
saveRDS(covar_gs, "../Processed/GrandeSynthe/Participant_covariate_info_subset_no_isolated.rds")

saveRDS(mat_lux, "../Processed/Luxembourg/Chemical_compound_matrix_subset_trans_imp_no_isolated.rds")
saveRDS(mat_fra, "../Processed/France/Chemical_compound_matrix_subset_trans_imp_no_isolated.rds")
saveRDS(mat_gs, "../Processed/GrandeSynthe/Chemical_compound_matrix_subset_trans_no_isolated.rds")

saveRDS(nd_lux, "../Processed/Luxembourg/Chemical_compound_nd_matrix_subset_no_isolated.rds")
saveRDS(nd_fra, "../Processed/France/Chemical_compound_nd_matrix_subset_no_isolated.rds")
saveRDS(nd_gs, "../Processed/GrandeSynthe/Chemical_compound_nd_matdrix_subset_no_isolated.rds")

### Univariate regression (family id) ----
## Exposure vs family ID by compound
# fraembourg
pvals_lux = NULL
f1='mat_lux[,i] ~ covar_lux$Family.ID'
f0='mat_lux[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(mat_lux)){
  model1=lm(as.formula(f1))
  model0=lm(as.formula(f0))
  pvals_lux=c(pvals_lux, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pvals_lux) = colnames(mat_lux)

saveRDS(pvals_lux, "../Results/Luxembourg/univar_family_pvals.rds")

# France
pvals_fra = NULL
f1='mat_fra[,i] ~ covar_fra$Family.ID'
f0='mat_fra[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(mat_fra)){
  model1=lm(as.formula(f1))
  model0=lm(as.formula(f0))
  pvals_fra=c(pvals_fra, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pvals_fra) = colnames(mat_fra)

saveRDS(pvals_fra, "../Results/France/univar_family_pvals.rds")

# Grande-Synthe
pvals_gs = NULL
f1='mat_gs[,i] ~ covar_gs$Family.ID'
f0='mat_gs[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(mat_gs)){
  model1=lm(as.formula(f1))
  model0=lm(as.formula(f0))
  pvals_gs=c(pvals_gs, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pvals_gs) = colnames(mat_gs)

saveRDS(pvals_gs, "../Results/GrandeSynthe/univar_family_pvals.rds")

# Manhattan plot
mycolours = c("tomato","royalblue","forestgreen")
values = merge(as.data.frame(pvals_lux),as.data.frame(pvals_fra), by=0, all=TRUE)
rownames(values) = values[,1]
values = merge(values[,-1], as.data.frame(pvals_gs), by=0, all=TRUE)
rownames(values) = values[,1]
values = values[,-1]
values = -log10(values)

# Sort rows
values = values[intersect(names(annot), rownames(values)),]

bonf = c(0.05/length(pvals_lux), 0.05/length(pvals_fra), 0.05/length(pvals_gs))
bonf = -log10(bonf)

xseq = seq(1, nrow(values))
annot_sub=annot[rownames(values)]
{pdf("../Figures/Univariate_regression_family.pdf", width=14, height=8)
  par(mar=c(20,3,3,3))
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
  points(values[,3], pch = 15, col = mycolours[3], cex = 0.8)
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.7)}
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(family.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("topright", pch=c(17,19,15, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
         col=c(mycolours,darken(mycolours, 0.5), "grey"),
         legend = c("Luxembourg","France","Grande-Synthe",
                    "Bonferroni threshold (LUX)", "Bonferroni threshold (FRA)", "Bonferroni threshold (GS)",
                    "Nominal threshold"),
         cex=0.7, bg="white")
  dev.off()
}

## Detection status vs family ID by compound
# Luxembourg
pvals_lux = NULL
f1='nd_lux[,i] ~ covar_lux$Family.ID'
f0='nd_lux[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(nd_lux)){
  model1=glm(as.formula(f1), family = "binomial")
  model0=glm(as.formula(f0), family = "binomial")
  pvals_lux=c(pvals_lux, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pvals_lux) = colnames(nd_lux)

saveRDS(pvals_lux, "../Results/Luxembourg/univar_nd_family_pvals.rds")

# France
pvals_fra = NULL
f1='nd_fra[,i] ~ covar_fra$Family.ID'
f0='nd_fra[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(nd_fra)){
  model1=glm(as.formula(f1), family = "binomial")
  model0=glm(as.formula(f0), family = "binomial")
  pvals_fra=c(pvals_fra, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pvals_fra) = colnames(nd_fra)

saveRDS(pvals_fra, "../Results/France/univar_nd_family_pvals.rds")

# Grande-Synthe
pvals_gs = NULL
f1='nd_gs[,i] ~ covar_gs$Family.ID'
f0='nd_gs[,i] ~ 1'
t0=Sys.time()
for (i in 1:ncol(nd_gs)){
  model1=glm(as.formula(f1), family = "binomial")
  model0=glm(as.formula(f0), family = "binomial")
  pvals_gs=c(pvals_gs, anova(model0, model1, test = 'Chisq')$`Pr(>Chi)`[2])
}
t1=Sys.time()
print(t1-t0)
names(pvals_gs) = colnames(nd_gs)

saveRDS(pvals_gs, "../Results/GrandeSynthe/univar_nd_family_pvals.rds")

# Manhattan plot
mycolours = c("tomato","royalblue","forestgreen")
values = merge(as.data.frame(pvals_lux),as.data.frame(pvals_fra), by=0, all=TRUE)
rownames(values) = values[,1]
values = merge(values[,-1], as.data.frame(pvals_gs), by=0, all=TRUE)
rownames(values) = values[,1]
values = values[,-1]
values = -log10(values)

# Sort rows
values = values[intersect(names(annot), rownames(values)),]

bonf = c(0.05/length(pvals_lux), 0.05/length(pvals_fra), 0.05/length(pvals_gs))
bonf = -log10(bonf)

xseq = seq(1, nrow(values))
annot_sub=annot[rownames(values)]
{pdf("../Figures/Univariate_regression_nd_family.pdf", width=14, height=8)
  par(mar=c(20,3,3,3))
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
  points(values[,3], pch = 15, col = mycolours[3], cex = 0.8)
  for(i in 1:length(xseq)){
    axis(1, at=xseq[i], labels = rownames(values)[i], las=2, cex.axis = 0.7)}
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(family.colours[unique(annot_sub)[k]], amount=0.5), cex.axis = 0.9)
  }
  legend("topright", pch=c(17,19,15, rep(NA,4)), lty = c(rep(NA,3),rep(2,4)),
         col=c(mycolours,darken(mycolours, 0.5), "grey"),
         legend = c("Luxembourg","France","Grande-Synthe",
                    "Bonferroni threshold (LUX)", "Bonferroni threshold (FRA)", "Bonferroni threshold (GS)",
                    "Nominal threshold"),
         cex=0.7, bg="white")
  dev.off()
}

### Within-family variation ----
families=unique(covars$Family.ID)
family_sd=matrix(NA, nrow=length(families), ncol=ncol(expo))
for (p in 1:ncol(expo)){
  for (f in 1:length(families)){
    family_sd[f,p]=sd(expo[covars$Family.ID==families[f],p])
  }
}
colnames(family_sd)=colnames(expo)
rownames(family_sd)=families
overall_sd=apply(expo,2,sd)

mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families

x=as.vector(row(family_sd))
y=as.vector(family_sd)
z=as.vector(col(family_sd))
mycolours = family.colours[unique(annot)]
annot_sub=annot[colnames(family_sd)]

{pdf(paste0("../Figures/",filepath,"/Univariate_sd_family_overall_cont.pdf"), width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(z,y, pch=19, cex=0.5, col=mycolours[x], las=1, xaxt="n",
       xlab="", ylab="Standard deviation", cex.lab=1.5,
       panel.first=abline(v=unique(z),lty=3,col="grey"),
       ylim=range(c(family_sd,overall_sd)))
  points(overall_sd, pch=15)
  for (k in 1:ncol(family_sd)){
    axis(side=1, at=k, labels=colnames(family_sd)[k], cex.axis=0.8, las=2)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(colnames(family_sd))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(mycolours[as.character(unique(annot_sub)[k])], amount=0.5), cex.axis = 0.9)
  }
  legend("top", col=c(mycolours, "black"), 
         pch=c(rep(19,length(families)),15), 
         pt.cex=c(rep(0.5,length(families)),1), 
         legend=c(families, "Overall"), ncol=10, bg="white")
  dev.off()
}

### Intra-class correlation ----
library(lme4)

expo = mat_lux
covars = covar_lux

### Univariate linear mixed models
icc_lux=NULL
for (p in 1:ncol(expo)){
  x=expo[,p]
  model=lmer(x~(1|covars$Family.ID))
  vcov = as.data.frame(VarCorr(model))$vcov
  icc_lux=c(icc_lux,vcov[1]/sum(vcov))
}
names(icc_lux) = colnames(expo)

expo = mat_fra
covars = covar_fra

### Univariate linear mixed models
icc_fra=NULL
for (p in 1:ncol(expo)){
  x=expo[,p]
  model=lmer(x~(1|covars$Family.ID))
  vcov = as.data.frame(VarCorr(model))$vcov
  icc_fra=c(icc_fra,vcov[1]/sum(vcov))
}
names(icc_fra) = colnames(expo)

expo = mat_gs
covars = covar_gs

### Univariate linear mixed models
icc_gs=NULL
for (p in 1:ncol(expo)){
  x=expo[,p]
  model=lmer(x~(1|covars$Family.ID))
  vcov = as.data.frame(VarCorr(model))$vcov
  icc_gs=c(icc_gs,vcov[1]/sum(vcov))
}
names(icc_gs) = colnames(expo)

# Manhattan plot
mycolours = c("tomato","royalblue","forestgreen")
values = merge(as.data.frame(icc_lux),as.data.frame(icc_fra), by=0, all=TRUE)
rownames(values) = values[,1]
values = merge(values[,-1], as.data.frame(icc_gs), by=0, all=TRUE)
rownames(values) = values[,1]
values = values[,-1]

# Sort rows
values = values[intersect(names(annot), rownames(values)),]

mycolours=c("tomato","royalblue","forestgreen")
annot_sub=annot[rownames(values)]
xseq = seq(1, nrow(values))

{pdf(paste0("../Figures/Intra_class_correlation_univariate_expo_cont.pdf"), width=14, height=8)
  par(mar=c(20,5,1,1))
  plot(values[,1], pch=17, cex=0.7, las=1, xaxt="n", type = "n",
       ylim = c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)),
       xlab="", ylab="Intra-Class Correlation", cex.lab=1.5,
       panel.first=abline(v=xseq,lty=3,col="grey"),
       col=mycolours[1])
  points(values[,1], pch = 17, col = mycolours[1], cex = 0.7)
  points(values[,2], pch = 19, col = mycolours[2], cex = 0.7)
  points(values[,3], pch = 15, col = mycolours[3], cex = 0.8)
  for (k in 1:length(xseq)){
    axis(side=1, at=xseq[k], labels=rownames(values)[k], cex.axis=0.8, las=2)
  }
  xgroup=c(which(!duplicated(annot_sub))-0.5, length(rownames(values))+0.5)
  axis(side=1, line=8, at=xgroup, labels=NA)
  tmp=apply(rbind(xgroup[-length(xgroup)],xgroup[-1]),2,mean)
  for (k in 1:length(unique(annot_sub))){
    axis(side=1, line=8, at=tmp[k], labels=unique(annot_sub)[k], tick=FALSE, las=2, 
         col.axis=darken(family.colours[unique(annot_sub)[k]], amount=0.5),
         cex.axis=0.9)
  }
  legend("topright", pch=c(17,19,15),
         col=mycolours,
         legend = c("Luxembourg","France","Grande-Synthe"),
         cex=0.7, bg="white")
  dev.off()
}
