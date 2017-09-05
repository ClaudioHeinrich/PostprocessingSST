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
dt[,residual:=SST_bar - SST_hat_grid]
##---------------------

lons = c(-80,0)
lats = c(20,80)


##------- Loop -------------
YM_all = dt[,unique(YM)]
for(j in 1:length(YM_all))
{
  print(j)
  plot_system(dt, YM_all[j])
  plot_system(dt, YM_all[j],
              file_out = paste0("./figures/system_local_",YM_all[j]),
              lons = lons,
              lats = lats)
}
##--------------------------


