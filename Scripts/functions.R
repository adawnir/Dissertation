# Group sampler
group_sample = function(x,group,n, exclude = NULL){
  if(is.null(exclude)){
    tmp = ave(x, group,
              FUN = function(a)
                if(length(a)>n)
                {sample(c(rep(TRUE, n), rep(FALSE, length(a)-n)),length(a))} else{TRUE})
  }
  else{
    tmp = ave(x, group,
              FUN = function(a)
                if(length(a)>n & all(a!=exclude))
                {sample(c(rep(TRUE, n), rep(FALSE, length(a)-n)),length(a))} else{TRUE})
  }
  return(as.logical(tmp))
}

# Make factor in the order of original data
factor.order = function(x) {factor(x, levels = unique(x))}

# Reorder last level
last.level = function(f, last) {
  if (!is.factor(f)) stop("f must be a factor")
  orig.levels = levels(f)
  if (! last %in% orig.levels) stop("last must be a level of f")
  new.levels = c(setdiff(orig.levels, last), last)
  factor(f, levels = new.levels)
}

# Compare multiple vectors and ignore NA
is.equal = function(mylist) {
  check.eq = sapply(mylist[-1], function(x) {x == mylist[[1]]})
  check.eq[is.na(check.eq)] = TRUE
  as.logical(apply(check.eq, 1, prod))
}

# Randomly sample from truncated normal distribution
rtruncnorm = function(n, mean = 0, sd = 1, min = -Inf, max = Inf){
  if (min > max) stop('Error: Truncation range is empty')
  x = runif(n, pnorm(min, mean, sd), pnorm(max, mean, sd))
  qnorm(x, mean, sd)
}

# Flexible formatting
flex_format = function(f, digits = 2, thresh = 0.1){
  ifelse(sapply(f, function(x) round(x,digits)) <= thresh,
         formatC(f, format="e", digits=2),
         formatC(f, format="f", digits=2))
}

# Create PCA score plot
CreateScorePlot=function(mypca, filename=NULL, type, mycolours, comp=NULL, pch=19, segments = TRUE, legend = TRUE){
  if (is.null(comp)){
    comp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
  }
  if (isTRUE(segments)){
    compare = sapply(type, function(x) x==type)
    compare[which(type=="Isolated"),which(type=="Isolated")] = FALSE
    compare[lower.tri(compare, diag = TRUE)] = NA
    start = which(compare, arr.ind=TRUE)[,1]
    end = which(compare, arr.ind=TRUE)[,2]
  }
  if (isTRUE(legend)){
    extra = 2.5
  } else {
    extra = 0
  }
  if (!is.null(filename)){
    pdf(paste0(filename), width=15+extra, height=5) 
  }
  if (isTRUE(legend)){
    par(mar = c(5,5,1,1))
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(2,2,2,1))
  } else {
    par(mar = c(5,5,1,1), mfrow = c(1,3))
  }
  for (k in 1:nrow(comp)){
    xcomp=comp[k,1]
    ycomp=comp[k,2]
    S = mypca$ind$coord[,c(xcomp,ycomp)]
    plot(S, pch=pch, cex=1.2, las=1, type = "n",
         col=mycolours[type],
         xlab=paste0("Comp ",xcomp," (", round(ev[xcomp], digits=2), "% e.v.)"),
         ylab=paste0("Comp ",ycomp," (", round(ev[ycomp], digits=2), "% e.v.)"),
         cex.lab = 1.5)
    if (isTRUE(segments)){
      segments(S[start,1],S[start,2],S[end,1],S[end,2], col = alpha("grey",0.5))
    }
    points(S, pch=pch, cex=1.2, col=mycolours[type])
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
  }
  if (isTRUE(legend)){
    par(mar = c(1,1,1,1))
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("left", col=mycolours, ncol = ceiling(length(mycolours)/25),
           pch=pch, pt.cex=1.2, legend=names(mycolours), bty = "n", cex = 1.2)
  }
  if (!is.null(filename)){
    print("Saved to filename")
    dev.off()
  }
}

# Creat PLS-DA score plot
CreateScorePlot.plsda=function(myplsda, filename=NULL, type1, type2, mycolours, comp=NULL, cex = 1.5, pch=19,
                               legend = TRUE, legend_text = NULL, legend_type1 = TRUE){
  if (is.null(comp)){
    comp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
  }
  compare = sapply(type2, function(x) x==type2)
  compare[which(type2=="Isolated"),which(type2=="Isolated")] = FALSE
  compare[lower.tri(compare, diag = TRUE)] = NA
  start = which(compare, arr.ind=TRUE)[,1]
  end = which(compare, arr.ind=TRUE)[,2]
  
  compare2 = sapply(type1, function(x) x==type1)
  compare2[lower.tri(compare2, diag = TRUE)] = NA
  start2 = which(compare2, arr.ind=TRUE)[,1]
  end2 = which(compare2, arr.ind=TRUE)[,2]
  if (isTRUE(legend)){
    extra = 2.5
  } else {
    extra = 0
  }
  if (!is.null(filename)){
    pdf(paste0(filename), width=15+extra, height=5) 
  }
  if (isTRUE(legend)){
    par(mar = c(5,5,1,1))
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths=c(2,2,2,1))
  } else {
    par(mar = c(5,5,1,1), mfrow = c(1,3))
  }
  for (k in 1:nrow(comp)){
    xcomp=comp[k,1]
    ycomp=comp[k,2]
    ev = myplsda$explained_variance$X*100
    S = myplsda$variates$X[,c(xcomp, ycomp)]
    plot(S, pch=pch, cex=cex, las=1, type = "n",
         col=mycolours[type1],
         xlab=paste0("Comp ",xcomp," (", round(ev[xcomp], digits=2), "% e.v.)"),
         ylab=paste0("Comp ",ycomp," (", round(ev[ycomp], digits=2), "% e.v.)"),
         cex.lab = 1.5)
    segments(S[start,1],S[start,2],S[end,1],S[end,2], lwd = 2,
             col = ifelse(paste(start,end) %in% paste(start2,end2),mycolours[type1[start]],alpha("grey80",0.3)))
    points(S, pch=pch, cex=cex,col=mycolours[type1])
    families=unique(type1)
    for (f in 1:length(families)){
      tmpmat=S[type1==families[f],]
      if (!is.null(nrow(tmpmat))){
        if (nrow(tmpmat)>2){
          lines(ellipse(cov(tmpmat), centre = c(mean(tmpmat[,1]),mean(tmpmat[,2])), level = 0.9),
                lty = 3, col = mycolours[f])
        }
      }
      # else {
      #   draw.ellipse(tmpmat[1], tmpmat[2], 0.5, 0.5,lty = 3, border = mycolours[f])
      # }
    }
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
  }
  if (isTRUE(legend)){
    par(mar = c(1,1,1,1))
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    if (isTRUE(legend_type1)){
    legend("left", col=mycolours, ncol = ceiling(length(mycolours)/25),
           pch=pch, pt.cex=cex, legend=names(mycolours), bty = "n", cex = 1.2)
      }
    legend("topleft", legend=legend_text, bty = "n", cex = 1.2)
  }
  if (!is.null(filename)){
    print("Saved to filename")
    dev.off()
  }
}

# Create clustering dendrogram graph
ClusteringToGraph=function(covars,myphylo,mycol=NULL,verbose = TRUE){
  myphylo=as.phylo(h)
  graph_edges = myphylo$edge
  graph_net=graph.edgelist(graph_edges, directed=FALSE)
  set.seed(1)
  graph_layout = layout_with_kk(graph_net)
  nobs = nrow(expo)
  
  if (isTRUE(verbose)){
    par(mar=c(0,1,0,1))
    plot(graph_layout[,1], graph_layout[,2], type = "n", axes = FALSE,
         xlab = "", ylab = "")
    # draw tree branches
    segments(
      x0 = graph_layout[graph_edges[,1],1], 
      y0 = graph_layout[graph_edges[,1],2],
      x1 = graph_layout[graph_edges[,2],1],
      y1 = graph_layout[graph_edges[,2],2],
      col = "darkgrey", lwd = 1)
    # add labels
    points(graph_layout[1:nobs,1], graph_layout[1:nobs,2], pch=19,
           col=mycol)
    text(graph_layout[1:nobs,1], graph_layout[1:nobs,2],
         col=mycol,
         myphylo$tip.label, cex = 0.5, xpd = TRUE, font = 1)
  }
  return(graph_net)
}
