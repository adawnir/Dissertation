### Multivariate clustering----

families=unique(covars$Family.ID)

mycolours=brewer.pal(n=12,name='Paired')
mycolours=colorRampPalette(mycolours)(length(families))
names(mycolours)=families

d=dist(expo)
h=hclust(d, method = "complete")
print(all(covars$Indiv.ID==h$labels))
h$labels=paste0(covars$Family.ID, "-",h$labels)
myphylo=as.phylo(h)

{pdf(paste0("../Figures/",filepath,"/Hierarchical_cont_expo_all.pdf"),width=14)
  par(mar=c(0,0,0,0))
  plot(myphylo, direction="downwards", cex=0.5, #srt=180, adj=1,
       tip.color=mycolours[covars$Family.ID])
  dev.off()
}

{pdf(paste0("../Figures/",filepath,"/Hierarchical_cont_graph_expo_all.pdf"))
  g=ClusteringToGraph(covars=covars, myphylo=myphylo)
  dev.off()
}

nobs=nrow(expo)
sp=shortest.paths(g)[1:nobs,1:nobs]
rownames(sp)=colnames(sp)=myphylo$tip.label

fsp=NULL
for (f in 1:length(families)){
  tmpmat=sp[covars$Family.ID==families[f],covars$Family.ID==families[f]]
  fsp=c(fsp,mean(tmpmat[upper.tri(tmpmat)]))
}
names(fsp)=families

{pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all.pdf"),width=14)
  par(mar=c(5,5,1,1))
  plot(fsp, pch=19, col=mycolours[families], xaxt="n", las=1,
       panel.first=abline(v=1:length(families), lty=3, col="grey"),
       xlab="Family ID", cex=2,
       ylab="Shortest path length between siblings", cex.lab=1.5)
  axis(side=1, at=1:length(families), labels=families, las = 2)
  dev.off()}

## Relationship with age

age_mu=rep(NA, length(families))
for (f in 1:length(families)){
  age_mu[f]=mean(covars$Age[covars$Family.ID==families[f]])
}
names(age_mu)=families

model=lm(fsp~age_mu)
{pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all_vs_age_mean.pdf"))
  par(mar=c(5,5,1,1))
  plot(age_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(age_mu)],
       xlab="Within-family age mean", ylab="Shortest path length between siblings")
  text(age_mu, fsp, labels=names(age_mu), pos=3, col=darken(mycolours[names(age_mu)], amount=0.5))
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  legend("topright", bty="n", cex=1.5,
         legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
  dev.off()
}

age_sd=rep(NA, length(families))
for (f in 1:length(families)){
  age_sd[f]=sd(covars$Age[covars$Family.ID==families[f]])
}
names(age_sd)=families

model=lm(fsp~age_sd)
{pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all_vs_age_sd.pdf"))
  par(mar=c(5,5,1,1))
  plot(age_sd, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(age_sd)],
       xlab="Within-family age standard deviation", ylab="Shortest path length between siblings")
  text(age_sd, fsp, labels=names(age_sd), pos=3, col=darken(mycolours[names(age_sd)], amount=0.5))
  abline(h=axTicks(2), col="grey", lty=3)
  abline(v=axTicks(1), col="grey", lty=3)
  legend("topright", bty="n", cex=1.5,
         legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
  dev.off()
}

# weight_mu=rep(NA, length(families))
# for (f in 1:length(families)){
#   weight_mu[f]=mean(covars$Weight[covars$Family.ID==families[f]])
# }
# names(weight_mu)=families
# 
# model=lm(fsp~weight_mu)
# {pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all_vs_weight_mean.pdf"))
#   par(mar=c(5,5,1,1))
#   plot(weight_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(weight_mu)],
#        xlab="Within-family weight mean", ylab="Shortest path length between siblings")
#   text(weight_mu, fsp, labels=names(weight_mu), pos=3, col=darken(mycolours[names(weight_mu)], amount=0.5))
#   abline(h=axTicks(2), col="grey", lty=3)
#   abline(v=axTicks(1), col="grey", lty=3)
#   legend("topleft", bty="n", cex=1.5,
#          legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
#   dev.off()
#   }
# 
# 
# weight_sd=rep(NA, length(families))
# for (f in 1:length(families)){
#   weight_sd[f]=sd(covars$Weight[covars$Family.ID==families[f]])
# }
# names(weight_sd)=families
# 
# model=lm(fsp~weight_sd)
# {pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all_vs_weight_sd.pdf"))
#   par(mar=c(5,5,1,1))
#   plot(weight_sd, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(weight_sd)],
#        xlab="Within-family weight standard deviation", ylab="Shortest path length between siblings")
#   text(weight_sd, fsp, labels=names(weight_sd), pos=3, col=darken(mycolours[names(weight_sd)], amount=0.5))
#   abline(h=axTicks(2), col="grey", lty=3)
#   abline(v=axTicks(1), col="grey", lty=3)
#   legend("topright", bty="n", cex=1.5,
#          legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
#   dev.off()
#   }
# 
# length_mu=rep(NA, length(families))
# for (f in 1:length(families)){
#   length_mu[f]=mean(covars$Length[covars$Family.ID==families[f]])
# }
# names(length_mu)=families
# 
# model=lm(fsp~length_mu)
# {pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all_vs_length_mean.pdf"))
#   par(mar=c(5,5,1,1))
#   plot(length_mu, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(length_mu)],
#        xlab="Within-family length mean", ylab="Shortest path length between siblings")
#   text(length_mu, fsp, labels=names(length_mu), pos=3, col=darken(mycolours[names(length_mu)], amount=0.5))
#   abline(h=axTicks(2), col="grey", lty=3)
#   abline(v=axTicks(1), col="grey", lty=3)
#   legend("topright", bty="n", cex=1.5,
#          legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
#   dev.off()
#   }
# 
# length_sd=rep(NA, length(families))
# for (f in 1:length(families)){
#   length_sd[f]=sd(covars$Length[covars$Family.ID==families[f]])
# }
# names(length_sd)=families
# 
# model=lm(fsp~length_sd)
# {pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all_vs_length_sd.pdf"))
#   par(mar=c(5,5,1,1))
#   plot(length_sd, fsp, pch=19, cex=1, cex.lab=1.5, las=1, col=mycolours[names(length_sd)],
#        xlab="Within-family length mean", ylab="Shortest path length between siblings")
#   text(length_sd, fsp, labels=names(length_sd), pos=3, col=darken(mycolours[names(length_sd)], amount=0.5))
#   abline(h=axTicks(2), col="grey", lty=3)
#   abline(v=axTicks(1), col="grey", lty=3)
#   legend("topright", bty="n", cex=1.5,
#          legend=paste0("p=",formatC(summary(model)$coefficients[2,4], format="e", digits=2)))
#   dev.off()
#   }


gender_diff=rep(NA, length(families))
for (f in 1:length(families)){
  tmp=covars$Gender[covars$Family.ID==families[f]]
  if (length(unique(tmp))==1){
    if (unique(tmp)=="Male"){
      gender_diff[f]=1
    }
    if (unique(tmp)=="Female"){
      gender_diff[f]=2
    }
  } else {
    gender_diff[f]=3
  } 
}
names(gender_diff)=families

model=lm(fsp~gender_diff)
{pdf(paste0("../Figures/",filepath,"/Shortest_path_cont_all_by_gender.pdf"))
  par(mar=c(5,5,1,1))
  plot(fsp, pch=19, cex=1, cex.lab=1.5, las=1, 
       col=c("skyblue", "pink", "tan")[gender_diff],
       xlab="Family ID", ylab="Shortest path length between siblings")
  text(fsp, labels=names(gender_diff), pos=3, col=darken(mycolours[names(gender_diff)], amount=0.5))
  legend("topleft", pch=19, col=c("skyblue", "pink", "tan"), ncol=1,
         legend=c("All Male", "All Female", "Different genders"))
  dev.off()
}
