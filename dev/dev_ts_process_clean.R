rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")
load("./FcNov2018/ts_hindcast_cleaned.RData")
print_figs = TRUE

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

Lon_bounds = c(-10,40)
DT_clean = DT_clean[between(lon,Lon_bounds[1],Lon_bounds[2])]
DT_clean[,grid_id:=.GRP,.(lon,lat)]
DT_clean[,"Lon":= lon]
DT_clean[,"Lat":= lat]
DT_clean[,climatology := ff_mean(obs_erai_ts),.(grid_id,month)]
DT_clean[,ecmwf_climatology := ff_mean(ecmwf_ts_bar),.(grid_id,month)]
DT_clean[,ukmo_climatology := ff_mean(ukmo_ts_bar),.(grid_id,month)]
DT_clean[,obs_anamoly := obs_erai_ts - climatology]
DT_clean[,ecmwf_anamoly := ecmwf_ts_bar - ecmwf_climatology]
DT_clean[,ukmo_anamoly := ukmo_ts_bar - ukmo_climatology]

mod_trend_obs = mgcv::gam(ukmo_anamoly ~ te(Lon,Lat), data = DT_clean)
mod_trend_ = mgcv::gam(obs_anamoly ~ te(Lon,Lat), data = DT_clean)
mod3 = mgcv::gam(ecmwf_anamoly ~ te(Lon,Lat), data = DT_clean)

anamoly_map = DT_clean[year > 2006 & month == 12, mean(obs_anamoly), .(Lon,Lat)]
