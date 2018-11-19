rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

DT = readr::read_rds("./FcNov2018/prect_hindcast.rds")

setkey(DT,"lon","lat","month","year")

Lon_bounds = c(-10,28)
Lat_bounds = c(31,80)
DT[,"Lon" := lon]
DT[, "Lat" := lat]
DT_obs = DT[is.finite(obs_erai_prect) & between(Lon,Lon_bounds[1], Lon_bounds[2]) & between(Lat, Lat_bounds[1], Lat_bounds[2])]
DT_small = DT_obs[,.(Lon,Lat,norcpm_prect_bar,
                     ecmwf_prect_bar,
                     mf_prect_bar,
                     ukmo_prect_bar,
                     obs_erai_prect,
                     year,
                     month)]
nms_adjust = c("norcpm_prect_bar","ecmwf_prect_bar","mf_prect_bar","ukmo_prect_bar")
ff = function(x){return(60 * 60 * 24 * 30 * 1e3 * x)}
for(j in 1:length(nms_adjust)){
    DT_small[,eval(nms_adjust[j]) := ff(get(nms_adjust[j]))]
}
DT_final = DT_small[!(is.na(ecmwf_prect_bar))]
DT_final[,grid_id:=.GRP,.(Lon,Lat)]
DT_final[,climatology := ff_mean(obs_erai_prect),.(grid_id,month)]
DT_final[,norcpm_climatology := ff_mean(norcpm_prect_bar),.(grid_id,month)]
DT_final[,ecmwf_climatology := ff_mean(ecmwf_prect_bar),.(grid_id,month)]
DT_final[,mf_climatology := ff_mean(mf_prect_bar),.(grid_id,month)]
DT_final[,ukmo_climatology := ff_mean(ukmo_prect_bar),.(grid_id,month)]

DT_final[,obs_anamoly := obs_erai_prect - climatology]
DT_final[,norcpm_anamoly := norcpm_prect_bar - norcpm_climatology]
DT_final[,ecmwf_anamoly := ecmwf_prect_bar - ecmwf_climatology]
DT_final[,mf_anamoly := mf_prect_bar - mf_climatology]
DT_final[,ukmo_anamoly := ukmo_prect_bar - ukmo_climatology]
DT_final[,obs_anamoly_1 := shift(obs_anamoly,1,NA,"lag"),.(grid_id,month)]

save(DT_final, file = "./FcNov2018/prect_hindcast_slimmed.RData")

