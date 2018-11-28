rm(list = ls())

library(PostProcessing)
library(data.table)
library(ncdf4)


setwd("~/PostClimDataNoBackup/SFE/FcNov2018/")

nc = nc_open("./fcnov2018.nc")

Lon = ncvar_get(nc, "lon")
Lat = ncvar_get(nc, "lat")
months = ncvar_get(nc,"month")

ts_mean = ncvar_get(nc,"ts_anamoly")

ts_dec = ts_mean[,,2]

dt_ts = data.table("Lon" = rep(Lon,times = length(Lat)),
                   "Lat" = rep(Lat, each = length(Lon)),
                   "ts" = as.vector(ts_dec))
rr = max(abs(dt_ts[,ts]),na.rm=TRUE)
plot_smooth(dt_ts, exclude_ocean = TRUE, rr = c(-rr,rr))

nc_close(nc)


