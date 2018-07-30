
#' crps score for normal distribution, resistant to missing values and 0 standard deviation  
#'
#' @param y vector of observations.
#' @param mean,sd mean and sd of the forecasting distribution
#'                   
#' @return vector of the same length as y containing crps scores.
#'
#' @author Claudio Heinrich
#' @examples crps.na.rm(c(NA,rnorm(10)), 1,1)
#' 
#' @importFrom scoringRules crps
#' 
#' @export


crps.na.rm = function(y, mean, sd){
  
  na_loc = which( is.na(y) | is.na(mean) | is.na(sd) | sd == 0)
  
  if(identical(na_loc,integer(0))){
    x = scoringRules::crps(y, family = "normal", mean = mean, sd = sd)
  }else{
    x = rep(0,length(y))
    x[-na_loc] = scoringRules::crps(y[-na_loc], family = "normal", mean = mean[-na_loc], sd = sd[-na_loc])
    x[na_loc] = NA
  }
  return(x)
}


GneitingWeightFct = function(x,L)
  {
    
    ret_value = rep(0,length(x))
    
    ret_value[x == 0] = 1
    
    t = x[x != 0]/L
    ret_value[x != 0] = (1-t)*sin(2*pi*t)/(2*pi*t)+(1-cos(2*pi*t))/(2*pi^2*t)   
    
    return(ret_value)
  }







### functions from 'Assessing the Calibration of High-Dimensional Ensemble Forecasts Using Rank Histograms'


## Minimum spanning tree ranks 
mst.rank <- function (x) {
  l.mst <- NULL
  for(f in 1:(dim(x)[2])) {
    euc.dist <- rdist::rdist(t(x[,-f]))
    l.mst <- c(l.mst,sum(vegan::spantree(euc.dist)$dist))
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

rhist_dt <- function(B, ens_size = 9, breaks = ens_size + 1, hist_xlab="", hist_ylab="", hist_ylim=NULL)
{
  hist(as.vector(B[[2]]),breaks=seq(0, max(as.vector(B[[2]])), length.out = breaks), main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray80",border="gray60",ylim=hist_ylim)
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


