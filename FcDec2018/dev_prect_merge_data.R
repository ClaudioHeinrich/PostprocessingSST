rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

load("./Fc_201812/Obs_precip.RData")
ff = function(x){return(60 * 60 * 24 * 30 * 1e3 * x)}

DT_ukmo = readr::read_rds("./Fc_201812/ukmo_prect_hindcast_12.rds")
DT_ukmo[, Lon := lon]
DT_ukmo[, Lat := lat]
DT_ukmo[, c("lon","lat") := NULL]
setkey(DT_ukmo, "Lon","Lat", "year","month")

DT_merge1 = merge(DT_obs_precip,
                  DT_ukmo[,.(Lon, Lat, year, month, vintagemonth, vintageyear, "ukmo_precip" = ff(prect_bar))],
                  by = c("year","month","Lon","Lat"),
                  all.x = FALSE,
                  all.y = TRUE)

DT_mf = readr::read_rds("./Fc_201812/meteo_france_prect_hindcast_12.rds")
DT_mf[,Lon := lon]
DT_mf[,Lat := lat]
DT_mf[, c("lon","lat") := NULL]
setkey(DT_mf, "Lon","Lat", "year","month")

DT_merge2 = merge(DT_mf[,.(Lon,Lat,year,month,"mf_precip" = ff(prect_bar))],
                 DT_merge1,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)

DT_cmcc = readr::read_rds("./Fc_201812/cmcc_prect_hindcast_12.rds")
DT_cmcc[, Lon := lon]
DT_cmcc[, Lat := lat]
DT_cmcc[, c("lon","lat") := NULL]
setkey(DT_cmcc, "Lon","Lat", "year","month")

DT_merge3 = merge(DT_cmcc[,.(Lon,Lat,year,month,"cmcc_precip" = ff(prect_bar))],
                 DT_merge2,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)

DT_dwd = readr::read_rds("./Fc_201812/dwd_prect_hindcast_12.rds")
DT_dwd[, Lon := lon]
DT_dwd[, Lat := lat]
DT_dwd[, c("lon", "lat") := NULL]
setkey(DT_dwd, "Lon","Lat", "year", "month")

DT_merge4 = merge(DT_dwd[,.(Lon, Lat, year, month, "dwd_precip" = ff(prect_bar))],
                 DT_merge3,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)

DT_ecmwf = readr::read_rds("./Fc_201812/ecmwf_prect_hindcast_12.rds")
DT_ecmwf[, Lon := lon]
DT_ecmwf[, Lat := lat]
DT_ecmwf[, c("lon", "lat") := NULL]
DT_ecmwf[,"prect_bar" := rowMeans(.SD), .SDcols = paste0("prect", 1:25)]
setkey(DT_ecmwf, "Lon","Lat", "year","month")

DT_merge5 = merge(DT_ecmwf[,.(Lon,Lat,year,month,"ecmwf_precip" = ff(prect_bar))],
                 DT_merge4,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)

DT_merge5[Lon > 180, Lon := Lon - 360]

DT_precip = DT_merge5

setkey(DT_precip, "Lat", "Lon", "month")


save(DT_precip, file = "./Fc_201812/DT_merged_precip.RData")
