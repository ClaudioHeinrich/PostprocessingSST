rm(list = ls())

##----- Set up ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
library(rgdal)
library(parallel)
##---------------------

load("~/PostClimDataNoBackup/SFE/Derived/dt_senorge_1957_1990.RData")
load("~/PostClimDataNoBackup/SFE/Derived/dt_senorge_1991_2015.RData")

yy = 1985
mm = 1

dt_ens = load_ensemble(yy,mm)

YM_i = yy * 12 + mm


dt_ens_grid = dt_ens[Ens == 1,.(Lon,Lat)]
dt_ens_grid[,ens_grid_id := 1:.N]

dt_senorge_grid = dt_senorge_early[,.(lon = head(lon,1),lat = head(lat,1)),senorge_grid_id]
M = matrix(NA,dim(dt_senorge_grid)[1],4)

helper = function(j)
{
  if(j %% 1e2 == 0)print(paste0(j," of ",dim(dt_senorge_grid)[1]))
  p = dt_senorge_early[j,.(lon,lat)]
  d = geosphere::distHaversine(p,dt_ens_grid[,.(Lon,Lat)])
  return(order(d)[1:4])
}

l = mclapply(1:dim(dt_senorge_grid)[1], "helper", mc.cores = 30, mc.silent = FALSE)
