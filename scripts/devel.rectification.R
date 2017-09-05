rm(list = ls())

##------ Set up ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
print_figs = FALSE
##------------------------

##----- Load Constituents ---
y = 2000
m = 9
dt_ens = load_ensemble(y,m)
dt_obs = load_observations(y,m)
##--------------------------


##--- Now Organize the Ensemble ----
dt_ens_coord = unique(dt_ens[,.(Lon,Lat)])
setkey(dt_ens_coord, "Lat","Lon")

dt_obs_coord = unique(dt_obs[,.(Lon,Lat)])
setkey(dt_obs_coord, "Lat","Lon")
##-----------------------------------

##----- Maps of grids ---------
if(print_figs){pdf("./figures/grid_ensemble.pdf")}else{X11()}
plot(dt_ens_coord[,.(Lon,Lat)], pch = ".", col="grey50", xlab = "Longitude", ylab = "Latitude", main = "Ensemble")
map("world", add = TRUE)
if(print_figs)dev.off()

if(print_figs){pdf("./figures/grid_obs.pdf")}else{X11()}
plot(dt_obs_coord[,.(Lon,Lat)], pch = ".", col="grey50", xlab = "Longitude", ylab = "Latitude", main = "Observations")
map("world", add = TRUE)
if(print_figs)dev.off()
##------------------------------


##------ Let's restrict a bit ---
lons = c(-10,0)
lats = c(40,50)
ss = restrict_coord(lons,lats)
dt_obs_sub= dt_obs_coord[ eval(ss) ]
dt_ens_sub =dt_ens_coord[ eval(ss) ]
##--------------------------------

##------ Get all boxes ---------
lon_all = sort(dt_obs_sub[,unique(Lon)])
cutoff_lon = head(lon_all,-1) + diff(lon_all)/2
n_lon = length(cutoff_lon)

lat_all = sort(dt_obs_sub[,unique(Lat)])
cutoff_lat = head(lat_all,-1) + diff(lat_all)/2
n_lat = length(cutoff_lat)
##------------------------------

##----- Maps of grids ---------
if(print_figs){pdf("./figures/grid_sub_ensemble.pdf")}else{X11()}
plot(dt_ens_sub, pch = 19, col="black", xlab = "Longitude", ylab = "Latitude")
points(dt_obs_sub,pch=19, col="blue")
map("world", add = TRUE)
for(i in 1:length(cutoff_lon))abline(v = cutoff_lon[i], lty =2 , col="blue")
for(i in 1:length(cutoff_lon))abline(h = cutoff_lat[i], lty =2 , col="blue")
if(print_figs)dev.off()
##------------------------------
