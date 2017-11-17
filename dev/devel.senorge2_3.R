## Now try to align ensemble
rm(list = ls())

##----- Set up ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
library(rgdal)
library(parallel)
##---------------------

##----- Get the SeNorge Data -----
load("~/PostClimDataNoBackup/SFE/Derived/dt_senorge_1957_1990.RData")
load("~/PostClimDataNoBackup/SFE/Derived/dt_senorge_1991_2015.RData")
load("~/PostClimDataNoBackup/SFE/Derived/GCFS1_wide.RData")
##--------------------------------

##------ Assemble grid objects --------
dt_ens_grid = dt_ens[YM == min(YM),.(lon,lat,GCFS1_id)]
dt_senorge_grid = dt_senorge_early[,.(lon = head(lon,1),lat = head(lat,1)),senorge_grid_id]
##--------------------------------------

##---- Loop through seNorge ----
helper = function(j)
{
  if(j %% 1e2 == 0)print(paste0(j," of ",dim(dt_senorge_grid)[1]))
  p = dt_senorge_early[j,.(lon,lat)]
  d = geosphere::distHaversine(p,dt_ens_grid[,.(lon,lat)])
  return(c(order(d)[1:4], d[ order(d)[1:4]]))
}
l = mclapply(1:dim(dt_senorge_grid)[1], "helper", mc.cores = 30, mc.silent = FALSE)
M = matrix(unlist(l),ncol= 8, byrow = TRUE)
colnames(M) = c(paste0("ens_id_",1:4),paste0("dist_ens_",1:4))
##--------------------------------

##---- Add info to grid ----------
dt_senorge_grid = cbind(dt_senorge_grid, M)
##--------------------------------

##------- Save -------
save(dt_ens_grid,dt_senorge_grid, file = "~/PostClimDataNoBackup/SFE/Derived/senorge2_gfsc1_map.RData")
##--------------------
