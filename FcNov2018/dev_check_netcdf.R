rm(list = ls())

library(PostProcessing)
library(data.table)
library(ncdf4)


setwd("~/PostClimDataNoBackup/SFE/FcNov2018/")

nc = nc_open("./fcnov2018.nc")

Lon = ncvar_get(nc, "lon")
Lat = ncvar_get(nc, "lat")
months = ncvar_get(nc,"month")

plot_x = function(X,Lon,Lat,mn, rr = NULL){
    X11()
    dt = data.table("Lon" = rep(Lon,times = length(Lat)),
                    "Lat" = rep(Lat, each = length(Lon)),
                    "x" = as.vector(X))
    plot_smooth(dt, exclude_ocean = TRUE, rr = rr, mn = mn)
}

ts_q_below = ncvar_get(nc,"ts_q_below")
for(j in 1:7){
    plot_x(ts_q_below[,,1,j], Lon, Lat, mn = paste0("Q Below ",j), rr = c(0,1))
}

ts_mean = ncvar_get(nc,"ts_mean")
plot_x(ts_mean[,,1], Lon, Lat, mn = "Mean")
ts_climatology = ncvar_get(nc,"ts_mean")
plot_x(ts_climatology[,,1], Lon, Lat, mn = "Climatology")
ts_anamoly = ncvar_get(nc,"ts_anamoly")
plot_x(ts_anamoly[,,1], Lon, Lat, mn = "Anamoly")

prect_q_below = ncvar_get(nc,"prect_q_below")
for(j in 1:7){
    plot_x(prect_q_below[,,1,j], Lon, Lat, mn = paste0("Q Below ",j), rr = c(0,1))
}

prect_mean = ncvar_get(nc,"prect_mean")
plot_x(prect_mean[,,1], Lon, Lat, mn = "Mean")
prect_climatology = ncvar_get(nc,"prect_mean")
plot_x(prect_climatology[,,1], Lon, Lat, mn = "Climatology")
prect_anamoly = ncvar_get(nc,"prect_anamoly")
plot_x(prect_anamoly[,,1], Lon, Lat, mn = "Anamoly")


nc_close(nc)


