


### functions from 'Assessing the Calibration of High-Dimensional Ensemble Forecasts Using Rank Histograms'


## Minimum spanning tree ranks 
mst.rank <- function (x) {
  l.mst <- NULL
  for(f in 1:(dim(x)[2])) {
    euc.dist <- rdist(t(x[,-f]))
    l.mst <- c(l.mst,sum(spantree(euc.dist)$dist))
  }
  x.rank <- rank(l.mst,ties="random")
  return(x.rank)
}

## Multivariate ranks 
mv.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- numeric(d[2])
  for(i in 1:d[2]) {
    x.prerank[i] <- sum(apply(x<=x[,i],2,all))
  }
  x.rank <- rank(x.prerank,ties="random")
  return(x.rank)
}

## Average ranks
avg.rank <- function(x)  {
  x.ranks <- apply(x,1,rank)
  x.preranks <- apply(x.ranks,1,mean)
  x.rank <- rank(x.preranks,ties="random")
  return(x.rank)
}

## Band depth ranks
bd.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- array(NA,dim=d)
  for(i in 1:d[1]) {
    tmp.ranks <- rank(x[i,])
    x.prerank[i,] <- (d[2] - tmp.ranks) * (tmp.ranks - 1)
  }
  x.rank <- apply(x.prerank,2,mean) + d[2] - 1
  x.rank <- rank(x.rank,ties="random")
  return(x.rank)
} 


######### rank histograms for data tables ###########

# The data table should have the key variable (most commonly YM) as first column and the ranks of the observations as second

rhist.dt <- function(B, ens.size = 9, breaks = ens.size + 1, hist_xlab="", hist_ylab="", hist_ylim=NULL)
{
  hist(as.vector(B[[2]]),breaks=seq(0, breaks, by=1), main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray80",border="gray60",ylim=hist_ylim)
  abline(a=length(B[[1]])/breaks, b=0, lty=2, col="gray30")
}




#' truncation function for freezing level of sea water
#' 
#' @author Claudio Heinrich
#' 
#' @example trc(c(-2,3))
#' 
#' @export

trc = function (x){ 
  truncation.value = -1.769995
  x = truncation.value * (x < truncation.value) + x * (x >= truncation.value)
  return(x)}
