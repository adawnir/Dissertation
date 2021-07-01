# Chemical family colours
family = c("Organochlorines", "Organophosphate metabolites", "Pyrethroids",
           "Pyrethroid metabolites", "PCBs", "PBDEs", "Acidic herbicides", "Anilinopyrimidines",
           "Azoles", "Benzamides", "Carbamates", "Carboxamides", "Neonicotinoids", "Oxadiazines",
           "Phenylpyrazoles", "Strobilurins", "Triazines/Triazinones/Diazines", "Subtituted ureas",
           "Dinitroanilines", "Other")
family.colours = c("slateblue","darkblue","deepskyblue4","deepskyblue3", "darkturquoise",
                   "darkcyan", "springgreen", "darkgreen","yellowgreen", "darkkhaki",
                   "goldenrod", "peru", "sienna", "darkorange", "orangered",
                   "salmon", "hotpink","violetred","darkorchid","grey50")

# Subset string from right (x = string value, n = length of subset)
substrRight = function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

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

# Hierarchical clustering: average Silhouette Score
avg_sil_hc = function(k, hc, dist) {
  res = cutree(hc, k=k)
  ss = silhouette(res, dist)
  mean(ss[,3])
}

# Kmeans clustering: average Silhouette Score
avg_sil_km = function(k, mat, nstart = 10) {
  km.res = kmeans(mat, centers = k, nstart = nstart)
  ss = silhouette(km.res$cluster, dist(mat))
  mean(ss[, 3])
}

# Hierarchical clustering: total within-cluster sum of square
wss_hc = function(i, hc, x){
  cl = cutree(hc, i)
  spl = split(x, cl)
  wss = sum(sapply(spl, function(d) sum(scale(d, scale = FALSE)^2)))
}

# Kmeans clustering: total within-cluster sum of square 
wss_km = function(k, mat, nstart = 10){
  kmeans(mat, k, nstart = nstart)$tot.withinss
}

# Get stable k
GetArgmax=function(calib_object){
  argmax=matrix(NA, nrow=ncol(calib_object$Lambda), ncol=2)
  for (block_id in 1:ncol(calib_object$Lambda)){
    if (ncol(calib_object$Lambda)==1){
      myS=calib_object$S
    } else {
      myS=calib_object$S_blocks[,block_id,drop=FALSE]
    }
    myS[is.na(myS)]=0
    myid=which.max(myS[,1])
    argmax[block_id,]=c(calib_object$Lambda[myid,block_id], calib_object$P[myid,block_id])
  }
  colnames(argmax)=c("lambda","pi")
  return(argmax)
}

# Random colour palette generator
random.colour.pal = function(n, seed = 1){
  set.seed(seed)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  mycolours = sample(col_vector, n)
  return(mycolours)
}
