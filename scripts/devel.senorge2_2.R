rm(list = ls())

##----- Set up ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
library(rgdal)
library(parallel)
##---------------------

##----- Grid Restrictions -----
lon_bound = c(4.5,7.8)
lat_bound = c(57.95, 61.04)
##-----------------------------

##------ Load Observations -------------------
obsdir = "~/PostClimDataNoBackup/seNorge/"
ff = paste0(obsdir,"seNorge2_TEMP1d_dataset_1957_1990.nc")
ncobs = nc_open(ff)
##--------------------------------------------
  
##------- Get high-level info -------
##senorge_obs = ncvar_get(ncobs, "mean_temperature")
x_coord = ncvar_get(ncobs,"X")
y_coord = ncvar_get(ncobs,"Y")
days = ncvar_get(ncobs,"time")
DD = as.Date(days,origin = "1900-01-01")
year = as.numeric(format(DD,"%Y"))
month = as.numeric(format(DD,"%m"))
YM = year * 12 + month
YM_all = unique(YM)
##------------------------------------

##------ Organize grid ---------
N_grid = length(x_coord) * length(y_coord)
Coord_long = cbind(rep(x_coord,length = N_grid),rep(y_coord,each = length(x_coord)))
SP <- SpatialPoints(Coord_long, proj4string=CRS("+proj=utm +zone=33") )
SP_lon_lat = spTransform(SP, CRS("+proj=longlat"))
##-------------------------------

##------ Now construct sub grid-level YM averages---
dt_list = list()
for(i in 1:length(YM_all))
{
  print(paste0("On ",i," of ",length(YM_all)))
  w_in = which(YM == YM_all[i])
  dt_list[[i]] = data.table(year = year[w_in[1]],
                            month = month[w_in[1]],
                            YM= YM_all[i],
                            lon = SP_lon_lat$coords.x1,
                            lat = SP_lon_lat$coords.x2,
                            senorge_grid_id = 1:length(SP_lon_lat$coords.x1),
                            temp=0.0)
  n_in = length(w_in)
  for(j in 1:length(w_in))
  {
    temp_j = as.vector(ncvar_get(ncobs, "mean_temperature",start = c(1,1,w_in[j]), count = c(-1,-1,1)))
    dt_list[[i]][,temp := temp + temp_j/n_in]
  }

  dt_list[[i]] = dt_list[[i]][ (lon > lon_bound[1]) &
                               (lon < lon_bound[2]) &
                               (lat > lat_bound[1]) &
                               (lat < lat_bound[2]) ]
}
##-------------------------------------------------------

dt_senorge_early = rbindlist(dt_list)

save(dt_senorge_early,file = "~/PostClimDataNoBackup/SFE/Derived/dt_senorge_1957_1990.RData")

quit(save = "no")




