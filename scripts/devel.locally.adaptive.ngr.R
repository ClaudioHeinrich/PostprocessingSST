rm(list = ls())

##--- Set up ----
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
##---------------

##-- Load Data --
dt = load_combined_wide()
setkey(dt,"grid_id","YM")
##---------------

##-- Single Example ----
YM_test = 24000
month_test = 12
dt_train = dt[(YM < YM_test)]## & (month == month_test)]
dt_train[,SST_climatology := mean(SST_bar),grid_id]
dt_train[,Ens_climatology := mean(Ens_bar),grid_id]

ff = function(Ens_cl, Ens, SST_cl, SST)
{
  if(any(is.na(c(Ens_cl,SST_cl))))
  {
    return(-Inf)
  }else{
    D = data.frame(Y = SST - SST_cl,X = Ens - Ens_cl)
    mod = lm("Y ~ X - 1", data = D)
    return(mod$coef[1])
  }
}

dt_reg = dt_train[,ff(Ens_climatology,Ens_bar,SST_climatology,SST_bar),grid_id]
##----------------------
