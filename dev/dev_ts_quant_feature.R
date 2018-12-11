rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)

setwd("~/PostClimDataNoBackup/SFE/")

load("./FcNov2018/ts_hindcast_cleaned.RData")

DT_long_list = list()
for(j in 1:28){
    nms = c("lon","lat","month","year",paste0("ukmo_ts",j))
    DT_long_list[[j]] = DT_clean[,.SD, .SDcols = nms]
    names(DT_long_list[[j]])[5] = "ukmo"
    DT_long_list[[j]][,ensemble:=j]
}

DT_long = rbindlist(DT_long_list)

DT_q = DT_long[,.("ukmo_median" = quantile(ukmo,0.5, na.rm=TRUE)), .(lon,lat,month)]

DT_clean[,obs_climatology := median(obs_erai_ts, na.rm = TRUE),
         .(lon, lat, month)]
DT_clean = merge(DT_clean, DT_q, by = c("lon","lat","month"), all.x = TRUE, all.y = FALSE)

DT_clean[, obs_anomaly := obs_erai_ts - obs_climatology]
for(j in 1:28){
    nm = paste0("ukmo_anamoly_",j)
    DT_clean[,eval(nm) := get(paste0("ukmo_ts",j)) - ukmo_median]
}

save(DT_clean, file = "./FcNov2018/ts_hindcast_cleaned_quantiled.RData")




