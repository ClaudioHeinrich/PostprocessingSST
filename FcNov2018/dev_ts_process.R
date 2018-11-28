rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

DT_mf = readr::read_rds("./FcNov2018/mf_ts_hindcast.rds")
DT = readr::read_rds("./FcNov2018/wfc_ts_hindcast_vintage-11.rds")

DT_forecast = readr::read_rds("./FcNov2018/wfc_ts_forecast.rds")

setkey(DT,"lon","lat","month","year")
DT = DT[!is.na(lon)]
DT[,grid_id:=.GRP,.(lon,lat)]
DT[,"Lon" := lon]
DT[,"Lat" := lat]
Lon_bounds = c(-10,40)
DT_train = DT[is.finite(obs_erai_ts) & between(Lon, Lon_bounds[1], Lon_bounds[2])]
setkey(DT_train,"lon","lat","month","year")

nms_adjust = c("ecmwf_ts_bar","ukmo_ts_bar", "norcpm_ts_bar")

q_report = c(.1,.25,.33,.5,.67,.75,.9)
for(i in 1:length(q_report)){
    q = q_report[i]
    q_lab = as.integer(q * 100)
    DT_train[, eval(paste0("q_", q_lab)) := quantile(obs_erai_ts, q), .(grid_id,month)]
}

DT_merge = merge(DT_train[,.("ecmwf_climatology" = mean(ecmwf_ts_bar, na.rm = TRUE),
                       "ukmo_climatology" = mean(ukmo_ts_bar, na.rm = TRUE),
                       "norcpm_climatology" = mean(norcpm_ts_bar, na.rm = TRUE)),
                       .(Lon,Lat,month)],



DT_forecast[,"Lon" := lon]
DT_forecast[,"Lat" := lat]
A = unique(DT_train[,.(grid_id,Lon,Lat)])
DT_forecast = merge(DT_forecast, A, by = c("Lon","Lat"), all = FALSE)
setkey(DT_forecast, "Lon","Lat","year","month")
DT_forecast = merge(DT_forecast,
                    DT_train[,.("q_low" = tail(q_low, 1),
                                "q_high" = tail(q_high, 1),
                                "ecmwf_climatology" = tail(ecmwf_climatology,1),
                                "ukmo_climatology" = tail(ukmo_climatology, 1),
                                "climatology" = tail(climatology, 1)),
                             .(grid_id, month)],
                    by = c("grid_id", "month"),
                    all.x = TRUE,
                    all.y = FALSE)
DT_forecast[, ecmwf_anamoly := ecmwf_prect_bar - ecmwf_climatology]
DT_forecast[, ukmo_anamoly := ukmo_prect_bar - ukmo_climatology]

save(DT_train,DT_forecast, file = "./FcNov2018/prect_processed.RData")
