
#' crps score for normal distribution, resistant to missing values and 0 standard deviation  
#'
#' @param y vector of observations.
#' @param mean,sd mean and sd of the forecasting distribution
#'                   
#' @return vector of the same length as y containing crps scores.
#'
#' @author Claudio Heinrich
#' @examples crps_na_rm(c(NA,rnorm(10)), 1,1)
#' 
#' @importFrom scoringRules crps
#' 
#' @export 


crps_na_rm = function(y, mean, sd)
  {
    na_loc = which( is.na(y) | is.na(mean) | is.na(sd) | sd == 0)
    
    if(identical(na_loc,integer(0)))
    {
      x = scoringRules::crps(y, family = "normal", mean = mean, sd = sd)
    }else
    {
      x = rep(0,length(y))
      x[-na_loc] = scoringRules::crps(y[-na_loc], family = "normal", mean = mean[-na_loc], sd = sd[-na_loc])
      x[na_loc] = NA
    }
    return(x)
  }


#'@export

GneitingWeightFct = function(x,L)
  {
    
    ret_value = rep(0,length(x))
    
    ret_value[x == 0] = 1
    
    t = x[x != 0]/L
    ret_value[x != 0] = (1-t)*sin(2*pi*t)/(2*pi*t)+(1-cos(2*pi*t))/(2*pi^2*t)   
    
    return(ret_value)
  }





#' truncation function for freezing level of sea water
#' 
#' @author Claudio Heinrich
#' 
#' @example trc(c(-2,3))
#' 
#' @export

trc = function (x,truncation.value = -1.79){ 
  x = truncation.value * (x < truncation.value) + x * (x >= truncation.value)
  return(x)}


