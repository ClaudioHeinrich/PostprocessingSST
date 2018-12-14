rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

DT_ukmo = readr::read_rds("./FcNov2018/ukmo_prect_hindcast.rds")
DT_ukmo[,ukmo_prect_bar:= rowMeans(.SD),.SDcols = paste0("prect",1:28)]
DT_mf = readr::read_rds("./FcNov2018/mf_prect_hindcast.rds")
DT_mf[, mf_prect_bar := rowMeans(.SD), .SDcols = paste0("prect",1:15)]
DT_train_leftout = merge(DT_ukmo[,.(year, month, lat, lon, ukmo_prect_bar)],
                         DT_mf[,.(year, month, lat, lon, mf_prect_bar)],
                         by = c("year", "month", "lat", "lon"),
                         all = FALSE)

DT = readr::read_rds("./FcNov2018/wfc_prect_hindcast_vintage-11.rds")
setkey(DT,"lon","lat","month","year")
DT = DT[!is.na(lon)]
DT[,grid_id:=.GRP,.(lon,lat)]
DT[,"Lon" := lon]
DT[,"Lat" := lat]
Lon_bounds = c(-10,40)
DT_train_most = DT[is.finite(obs_erai_prect) & between(Lon, Lon_bounds[1], Lon_bounds[2]) & is.finite(ecmwf_prect_bar),
                   .(lon, lat, year, month, grid_id,
                     obs_erai_prect, norcpm_prect_bar, ecmwf_prect_bar,
                     norcpm_vintageyear, norcpm_vintagemonth,
                     cds_vintageyear, cds_vintagemonth)]

DT_train =  merge(DT_train_most,
                  DT_train_leftout,
                  by = c("year", "month", "lat", "lon"),
                  all = FALSE)

nms_adjust = c("ecmwf_prect_bar","ukmo_prect_bar", "norcpm_prect_bar", "mf_prect_bar")
ff = function(x){return(60 * 60 * 24 * 30 * 1e3 * x)}
for(j in 1:length(nms_adjust)){
    DT_train[,eval(nms_adjust[j]) := ff(get(nms_adjust[j]))]
}

 
setkey(DT_train,"lon","lat","month","year")

save(DT_train, file = "./FcNov2018/prect_processed.RData")

