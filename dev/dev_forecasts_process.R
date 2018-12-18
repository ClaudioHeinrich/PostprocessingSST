rm(list = ls())

library(PostProcessing)
library(ncdf4)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/Fc_201812/")

ens_names = c("cmcc", "dwd", "ecmwf", "meteo_france", "ukmo")

dt_model = list()
for(j in 1:length(ens_names)){
    nm = paste0("./forecasts/", ens_names[j], "_seasonal-monthly-single-levels_vintage-2018-12.nc")
    nc = nc_open(nm)
    Times = as.Date(ncvar_get(nc, "time") / 24, origin = "1900-01-01")
    if(ens_names[j] == "ecmwf"){
        Lons = ncvar_get(nc, "lon")
        Lats = ncvar_get(nc, "lat")
    }else{
        Lons = ncvar_get(nc, "longitude")
        Lats = ncvar_get(nc, "latitude")
    }
    A = ncvar_get(nc, "t2m")
    dt_m = list()
    for(t in 1:length(Times)){
        A_mean = apply(A[,,,t] - 272.15, c(1, 2), "mean")
        A_sd = apply(A[,,,t] - 272.15, c(1, 2), "sd")
        y_t = year(Times[t])
        m_t = month(Times[t])
        dt_m[[t]] = data.table(year = y_t,
                               month = m_t,
                               Lon = rep(Lons, times = length(Lats)),
                               Lat = rep(Lats, each = length(Lons)),
                               as.vector(A_mean),
                               as.vector(A_sd))
        names(dt_m[[t]])[5:6] = c(paste0(ens_names[j], "_bar"), paste0(ens_names[j], "_sd"))
    }
    dt_model[[j]] = rbindlist(dt_m)
    nc_close(nc)
}

DT_forecast = merge(dt_model[[1]], dt_model[[2]], by = c("year","month","Lon","Lat"))
for(j in 3:length(ens_names)){
    DT_forecast = merge(DT_forecast, dt_model[[j]], by = c("year","month","Lon","Lat"))
}

DT_forecast[Lon > 180,Lon := Lon - 360]

save(DT_forecast, file = "./Forecast_ts_processed.RData")

