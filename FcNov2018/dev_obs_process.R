rm(list = ls())

library(PostProcessing)
library(data.table)
library(ncdf4)
setwd("~/PostClimDataNoBackup/SFE/")

yy = 1993:2018
dt_y = list()

for(j in 1:length(yy)){
    y = yy[j]
    
    nc = nc_open(paste0("./Fc_201812/EraI_ts/EraI_ts_",y,"_p5.nc"))
    
    Times = as.Date(ncvar_get(nc, "time") / 24, origin = "1900-01-01")
    Lons = ncvar_get(nc, "lon")
    Lats = ncvar_get(nc, "lat")
    T_all = ncvar_get(nc, "t2m")
    
    dt_m = list()
    for(t in 1:length(Times)){
        y_t = year(Times[t])
        m_t = month(Times[t])
        A = T_all[,,t] - 272.15
        dt_m[[t]] = data.table(year = y_t,
                               month = m_t,
                               Lon = rep(Lons, times = length(Lats)),
                               Lat = rep(Lats, each = length(Lons)),
                               obs_t2m = as.vector(A))
    }
    dt_y[[j]] = rbindlist(dt_m)
}

DT_obs_t2m = rbindlist(dt_y)

save(DT_obs_t2m, file = "./Fc_201812/Obs_2tm.RData")
