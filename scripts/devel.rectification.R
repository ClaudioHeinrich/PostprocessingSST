rm(list = ls())

##------ Set up ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
##------------------------

##----- Load Constituents ---
y = 2000
m = 9
dt_ens = load_ensemble(y,m)
dt_obs = load_observations(y,m)
##--------------------------


##-------- Now Organize the Ensemble ----
dt_ens_coord = unique(dt_ens[,.(Lon,Lat)])
setkey(dt_ens_coord, "Lat","Lon")
dt_ens_coord[,Lon_ID := f_lon(Lon)]
dt_ens_coord[,Lat_ID := f_lat(Lat)]
dt_ens_coord[, grid:=Lon_ID + Lat_ID]

A = dt_ens_coord[(Lon > 48) & (Lon < 53) & (Lat > 48) & (Lat < 53)]
B = dt_obs_coord[(Lon > 48) & (Lon < 53) & (Lat > 48) & (Lat < 53)]

plot(A, pch = 19)
points(B, pch = 19, col="blue")
