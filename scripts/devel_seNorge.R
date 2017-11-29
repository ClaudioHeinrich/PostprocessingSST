rm(list = ls())

##------ Setup -------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
##---------------------

load("~/PostClimDataNoBackup/SFE/Derived/senorge2_gcfs1_combined.RData")

setkey(dt_senorge,senorge_grid_id, month, year)

dt_senorge[,temp_climatology:=(cumsum(temp) - temp) / (year - min(year)),.(senorge_grid_id, month)]

A = dt_senorge[,.(Mean_temp = mean(temp, na.rm = TRUE), Mean_Ens = mean(Ens_bar, na.rm = TRUE)),.(month,year)]
A[,"Running_Climatology" := (cumsum(Mean_temp) - Mean_temp) / (year - min(year)),month]
A[,"Running_Ens" := (cumsum(Mean_Ens) - Mean_Ens) / (year - min(year)),month]

A[,"Diff_Climate" := Mean_temp - Running_Climatology]
A[,"Diff_Ens" := Mean_Ens - Running_Ens]

id = 1578720

dt_id = dt_senorge[senorge_grid_id == id][month == 4]

rr = range(dt_id[,.(temp, Ens_bar)])
plot(dt_id[,.(year,temp)], type = "l", ylim = rr)
lines(dt_id[,.(year, Ens_bar)], col="blue")
abline(h = mean(dt_id[,temp]), col="red")
