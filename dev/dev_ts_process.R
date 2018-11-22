rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

DT = readr::read_rds("./FcNov2018/wfc_ts_hindcast.rds")
setkey(DT,"lon","lat","month","year")

DT[,grid_id:=.GRP,.(lon,lat)]
DT[,N:=.N,grid_id]
DT_info = DT[is.finite(lon)]
Lon_bounds = c(-10,40)
##Lat_bounds = c(31,80)
DT[,"Lon" := lon]
DT[, "Lat" := lat]
DT_small = DT[between(Lon,Lon_bounds[1],Lon_bounds[2]),.(Lon,Lat,
                                                        grid_id,
                                                        ecmwf_ts_bar,
                                                        ##         mf_ts_bar,
                                                        ukmo_ts_bar,
                                                        obs_erai_ts,
                                                        year,
                                                        month)]

##DT_final = DT_small[!(is.na(ecmwf_ts_bar) | is.na(mf_ts_bar) | is.na(ukmo_ts_bar) | is.na(obs_erai_ts))]
DT_final = DT_small[!(is.na(ecmwf_ts_bar)  | is.na(ukmo_ts_bar) | is.na(obs_erai_ts))]


DT_final[,climatology := ff_mean(obs_erai_ts),.(grid_id,month)]
DT_final[,ecmwf_climatology := ff_mean(ecmwf_ts_bar),.(grid_id,month)]
##DT_final[,mf_climatology := ff_mean(mf_ts_bar),.(grid_id,month)]
DT_final[,ukmo_climatology := ff_mean(ukmo_ts_bar),.(grid_id,month)]

DT_final[,obs_anamoly := obs_erai_ts - climatology]
DT_final[,ecmwf_anamoly := ecmwf_ts_bar - ecmwf_climatology]
##DT_final[,mf_anamoly := mf_ts_bar - mf_climatology]
DT_final[,ukmo_anamoly := ukmo_ts_bar - ukmo_climatology]
DT_final[,obs_anamoly_1 := shift(obs_anamoly,1,NA,"lag"),.(grid_id,month)]

save(DT_final, file = "./FcNov2018/ts_hindcast_slimmed.RData")
