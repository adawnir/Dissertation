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

factor.order = function(x) {factor(x, levels = unique(x))}

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
CreateScorePlot=function(mypca, filename=NULL, type, mycolours, comp=NULL,pch=19){
  if (is.null(comp)){
    comp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
  }
  compare = sapply(type, function(x) x==type)
  start = which(compare, arr.ind=TRUE)[,1]
  end = which(compare, arr.ind=TRUE)[,2]
  if (!is.null(filename)){
    pdf(paste0(filename), width=14, height=5) 
  }
  par(mfrow=c(1,3))
  for (k in 1:nrow(comp)){
    xcomp=comp[k,1]
    ycomp=comp[k,2]
    S = mypca$ind$coord[,c(xcomp,ycomp)]
    plot(S, pch=pch, cex=0.7, las=1, 
         col=ifelse(is.na(type), "grey", mycolours[type]),
         xlab=paste0("Comp ",xcomp," (", round(ev[xcomp], digits=2), "% e.v.)"),
         ylab=paste0("Comp ",ycomp," (", round(ev[ycomp], digits=2), "% e.v.)"))
    segments(S[start,1],S[start,2],S[end,1],S[end,2], col = alpha("grey",0.5))
    abline(v=axTicks(1), lty=3, col="grey")
    abline(h=axTicks(2), lty=3, col="grey")
    # if (k==1){
    #   legend("top", col=c(mycolours), 
    #          pch=c(rep(19,length(families))), 
    #          # cex=0.7,
    #          pt.cex=c(rep(0.5,length(families)),1), 
    #          legend=c(families), ncol=10, bg="white")
    # }
  }
  if (!is.null(filename)){
    print("Saved to filename")
    dev.off()
  }
}

# Create clustering dendrogram graph
ClusteringToGraph=function(covars,myphylo){
  myphylo=as.phylo(h)
  graph_edges = myphylo$edge
  graph_net=graph.edgelist(graph_edges, directed=FALSE)
  set.seed(1)
  graph_layout = layout_with_kk(graph_net)
  nobs = nrow(expo)
  
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
         col=mycolours[covars$Family.ID])
  text(graph_layout[1:nobs,1], graph_layout[1:nobs,2],
       col=mycolours[covars$Family.ID],
       myphylo$tip.label, cex = 0.5, xpd = TRUE, font = 1)
  
  return(graph_net)
}
