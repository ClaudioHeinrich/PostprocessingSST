rm(list = ls())

##----- Setup -------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
print_figs = FALSE
##-------------------

##---- Load -----
load("~/PostClimDataNoBackup/SFE/Derived/senorge2_gfsc1_map.RData")
load("~/PostClimDataNoBackup/SFE/Derived/senorge2_gfsc1_combined.RData")
##---------------

##--- Monthly bias correction ----
dt_senorge[,residual:=temp - Ens_bar]
setkey(dt_senorge, senorge_grid_id, month, year)
dt_senorge[,residual_bar:= (cumsum(residual) - residual) / (year - min(year)),.(senorge_grid_id, month)]
dt_senorge[year == min(year),residual_bar:=0.0, .(senorge_grid_id, month)]
dt_senorge[,paste0("Ens_bias_",1:15) := .SD + residual_bar,.SDcols = paste0("Ens_",1:15)]
dt_senorge[,Ens_bias_bar:= Ens_bar + residual_bar]
dt_senorge[,Residual_bias := temp - Ens_bias_bar]
setkey(dt_senorge, year, month, senorge_grid_id)
##--------------------------------




