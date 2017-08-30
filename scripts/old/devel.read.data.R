rm(list = ls())

##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)
print_figs = TRUE
##------------------------

##------ Set up -------
dt_combined_all = load_combined()
dt = rbindlist(dt_combined_all)
dt[, YM := year * 12 + month]
setkey(dt, "YM", "Lon", "Lat")
##---------------------

##------- Loop -------------
YM_all = dt[,unique(YM)]
for(j in 1:length(YM_all))
{
  print(j)
  plot_system(dt, YM_all[j])
}
##--------------------------


##------ Let's look at overall bias --
dt[,bias := SST_bar - SST_hat_grid]
ym_index = dt[, .("Year" = min(year),"Month" = min(month)),keyby=YM]
bias_ym = dt[,.("Bias"= mean(bias,na.rm=TRUE)),keyby = YM]
yy_all = dt[,range(year)]
##if(print_figs){pdf("
plot(bias_ym, xlab="Year",ylab="Global Mean Bias (C)",axes = FALSE, type = "l")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])

##--------------------------------------


      
