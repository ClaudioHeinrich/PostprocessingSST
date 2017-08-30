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

