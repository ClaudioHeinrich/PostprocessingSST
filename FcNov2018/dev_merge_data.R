rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

load("./Fc_201812/Obs_2tm.RData")

DT_ukmo = readr::read_rds("./Fc_201812/ukmo_ts_hindcast_12.rds")
DT_ukmo[,Lon := lon]
DT_ukmo[,Lat := lat]
DT_ukmo[, c("lon","lat") := NULL]
setkey(DT_ukmo, "Lon","Lat", "year","month")

DT_merge1 = merge(DT_obs_t2m,
                  DT_ukmo[,.(Lon,Lat,year,month,"ukmo_ts" = ts_bar)],
                  by = c("year","month","Lon","Lat"),
                  all.x = FALSE,
                  all.y = TRUE)

DT_mf = readr::read_rds("./Fc_201812/meteo_france_ts_hindcast_12.rds")
DT_mf[,Lon := lon]
DT_mf[,Lat := lat]
DT_mf[, c("lon","lat") := NULL]
setkey(DT_mf, "Lon","Lat", "year","month")

DT_merge2 = merge(DT_mf[,.(Lon,Lat,year,month,"mf_ts" = ts_bar)],
                 DT_merge1,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)

DT_cmcc = readr::read_rds("./Fc_201812/cmcc_ts_hindcast_12.rds")
DT_cmcc[,Lon := lon]
DT_cmcc[,Lat := lat]
DT_cmcc[, c("lon","lat") := NULL]
setkey(DT_cmcc, "Lon","Lat", "year","month")

DT_merge3 = merge(DT_cmcc[,.(Lon,Lat,year,month,"cmcc_ts" = ts_bar)],
                 DT_merge2,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)

DT_dwd = readr::read_rds("./Fc_201812/dwd_ts_hindcast_12.rds")
DT_dwd[,Lon := lon]
DT_dwd[,Lat := lat]
DT_dwd[, c("lon","lat") := NULL]
setkey(DT_dwd, "Lon","Lat", "year","month")

DT_merge4 = merge(DT_dwd[,.(Lon,Lat,year,month,"dwd_ts" = ts_bar)],
                 DT_merge3,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)

DT_ecmwf = readr::read_rds("./Fc_201812/ecmwf_ts_hindcast_12.rds")
DT_ecmwf[,Lon := lon]
DT_ecmwf[,Lat := lat]
DT_ecmwf[, c("lon","lat") := NULL]
DT_ecmwf[,"ts_bar" := rowMeans(.SD), .SDcols = paste0("ts",1:25)]
setkey(DT_ecmwf, "Lon","Lat", "year","month")

DT_merge5 = merge(DT_ecmwf[,.(Lon,Lat,year,month,"ecmwf_t2m" = ts_bar)],
                 DT_merge4,
                 by = c("year","month","Lon","Lat"),
                 all.y = TRUE)
DT_merge5[Lon > 180, Lon := Lon - 360]

DT_t2m = DT_merge5

save(DT_t2m, file = "./Fc_201812/DT_merged.RData")
