rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

DT_mf = readr::read_rds("./FcNov2018/mf_prect_hindcast.rds")
DT = readr::read_rds("./FcNov2018/wfc_prect_hindcast_vintage-11.rds")
DT_forecast = readr::read_rds("./FcNov2018/wfc_prect_forecast.rds")

setkey(DT,"lon","lat","month","year")
DT = DT[!is.na(lon)]
DT[,grid_id:=.GRP,.(lon,lat)]
DT[,"Lon" := lon]
DT[,"Lat" := lat]
Lon_bounds = c(-10,40)
DT_train = DT[is.finite(obs_erai_prect) & between(Lon, Lon_bounds[1], Lon_bounds[2])]
setkey(DT_train,"lon","lat","month","year")

nms_adjust = c("ecmwf_prect_bar","ukmo_prect_bar")
ff = function(x){return(60 * 60 * 24 * 30 * 1e3 * x)}
for(j in 1:length(nms_adjust)){
    DT_train[,eval(nms_adjust[j]) := ff(get(nms_adjust[j]))]
    DT_forecast[,eval(nms_adjust[j]) := ff(get(nms_adjust[j]))]
}

##setkey(DT_forecast,"lon","lat","month","year")
DT_train[,climatology := ff_mean(obs_erai_prect),.(grid_id,month)]
DT_train[,ecmwf_climatology := ff_mean(ecmwf_prect_bar),.(grid_id,month)]
##DT_final[,mf_climatology := ff_mean(mf_ts_bar),.(grid_id,month)]
DT_train[,ukmo_climatology := ff_mean(ukmo_prect_bar),.(grid_id,month)]

DT_train[, obs_anamoly := obs_erai_prect - climatology]
DT_train[, ecmwf_anamoly := ecmwf_prect_bar - ecmwf_climatology]
##DT_final[,mf_anamoly := mf_ts_bar - mf_climatology]
DT_train[, ukmo_anamoly := ukmo_prect_bar - ukmo_climatology]
DT_train[, obs_anamoly_1 := shift(obs_anamoly, 1, NA, "lag"), .(grid_id, month)]

DT_train[,q_low := quantile(obs_erai_prect, .25), .(grid_id,month)]
DT_train[,q_med := quantile(obs_erai_prect, .5), .(grid_id,month)]
DT_train[,q_high := quantile(obs_erai_prect, .75), .(grid_id,month)]

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
