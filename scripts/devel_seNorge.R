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


rr = range(na.omit(A[,.(Diff_Climate,Diff_Ens)]))

for(y in 1983:2015){

plot(A[year == y,Diff_Climate], main = y,type = 'l', ylim = rr)
lines(A[year == y,Diff_Ens],col = "blue")
}
