## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./FcNov2018/prect_processed.RData")
DT_train_prect = DT_train
load("./FcNov2018/ts_processed.RData")
DT_train_ts = DT_train

DT_train = merge(DT_train_ts,
                 DT_train_prect,
                 by = c("lon","lat","year","month",
                        "grid_id","norcpm_vintageyear","norcpm_vintagemonth",
                        "cds_vintageyear", "cds_vintagemonth"),
                 all = TRUE)
vintage_years = DT_train[,unique(cds_vintageyear)]

months = c(12,1,2,3,4)
mod_list = list()
mod = list()
DT_y_all = list()

DT_train_y = DT_train

DT_train_y[, N := .N, .(lat, lon, month)]
DT_train_y[, ukmo_prect_climatology := mean(ukmo_prect_bar), .(lat, lon, month)]
DT_train_y[, mf_prect_climatology := mean(mf_prect_bar), .(lat,lon,month)]
DT_train_y[, ecmwf_prect_climatology := mean(ecmwf_prect_bar), .(lat,lon,month)]
DT_train_y[, norcpm_prect_climatology := mean(norcpm_prect_bar), .(lat,lon,month)]
DT_train_y[, obs_prect_climatology := mean(obs_erai_prect), .(lat,lon,month)]

DT_train_y[, ukmo_ts_climatology := mean(ukmo_ts_bar), .(lat, lon, month)]
DT_train_y[, mf_ts_climatology := mean(mf_ts_bar), .(lat,lon,month)]
DT_train_y[, ecmwf_ts_climatology := mean(ecmwf_ts_bar), .(lat,lon,month)]
DT_train_y[, norcpm_ts_climatology := mean(norcpm_ts_bar), .(lat,lon,month)]
DT_train_y[, obs_ts_climatology := mean(obs_erai_ts), .(lat,lon,month)]

DT_train_y[, ukmo_ts_anamoly := ukmo_ts_bar - ukmo_ts_climatology, .(lat, lon, month)]
DT_train_y[, mf_ts_anamoly := mf_ts_bar - mf_ts_climatology, .(lat, lon, month)]
DT_train_y[, ecmwf_ts_anamoly := ecmwf_ts_bar - ecmwf_ts_climatology, .(lat, lon, month)]
DT_train_y[, norcpm_ts_anamoly := norcpm_ts_bar - norcpm_ts_climatology, .(lat, lon, month)]
DT_train_y[, obs_ts_anamoly := obs_erai_ts - obs_ts_climatology, .(lat, lon, month)]

DT_train_y[, ukmo_prect_anamoly := ukmo_prect_bar - ukmo_prect_climatology, .(lat, lon, month)]
DT_train_y[, mf_prect_anamoly := mf_prect_bar - mf_prect_climatology, .(lat, lon, month)]
DT_train_y[, ecmwf_prect_anamoly := ecmwf_prect_bar - ecmwf_prect_climatology, .(lat, lon, month)]
DT_train_y[, norcpm_prect_anamoly := norcpm_prect_bar - norcpm_prect_climatology, .(lat, lon, month)]
DT_train_y[, obs_prect_anamoly := obs_erai_prect - obs_prect_climatology, .(lat, lon, month)]


for(i in 1:length(vintage_years)){
    y = vintage_years[i]

    DT_train_y = DT_train[cds_vintageyear != y]
    setkey(DT_train_y,lon,lat,month,year)
    DT_train_y[, N := .N, .(lat, lon, month)]
    DT_train_y[, ukmo_prect_climatology := mean(ukmo_prect_bar), .(lat, lon, month)]
    DT_train_y[, mf_prect_climatology := mean(mf_prect_bar), .(lat,lon,month)]
    DT_train_y[, ecmwf_prect_climatology := mean(ecmwf_prect_bar), .(lat,lon,month)]
    DT_train_y[, norcpm_prect_climatology := mean(norcpm_prect_bar), .(lat,lon,month)]
    DT_train_y[, obs_prect_climatology := mean(obs_erai_prect), .(lat,lon,month)]

    nms = c("lat", "lon", "month", paste0(c("ukmo","mf","ecmwf","norcpm","obs"),"_prect_climatology"))
    DT_y = unique(merge(DT_train[cds_vintageyear == y],
                        DT_train_y[cds_vintageyear == y - 1,.SD,.SDcols = nms],
                        by = c("lat","lon","month"),
                        all = FALSE))

    DT_train_y[, ukmo_prect_anamoly := ukmo_prect_bar - ukmo_prect_climatology, .(lat, lon, month)]
    DT_train_y[, mf_prect_anamoly := mf_prect_bar - mf_prect_climatology, .(lat, lon, month)]
    DT_train_y[, ecmwf_prect_anamoly := ecmwf_prect_bar - ecmwf_prect_climatology, .(lat, lon, month)]
    DT_train_y[, norcpm_prect_anamoly := norcpm_prect_bar - norcpm_prect_climatology, .(lat, lon, month)]
    DT_train_y[, obs_prect_anamoly := obs_erai_prect - obs_prect_climatology, .(lat, lon, month)]

    DT_y[, ukmo_prect_anamoly := ukmo_prect_bar - ukmo_prect_climatology, .(lat, lon, month)]
    DT_y[, mf_prect_anamoly := mf_prect_bar - mf_prect_climatology, .(lat, lon, month)]
    DT_y[, ecmwf_prect_anamoly := ecmwf_prect_bar - ecmwf_prect_climatology, .(lat, lon, month)]
    DT_y[, norcpm_prect_anamoly := norcpm_prect_bar - norcpm_prect_climatology, .(lat, lon, month)]


    mod1 = lm(obs_prect_anamoly ~ ukmo_prect_anamoly + mf_prect_anamoly + norcpm_prect_anamoly + ecmwf_prect_anamoly,
              data = DT_train_y)
    mod2 = lm(obs_prect_anamoly ~ ukmo_prect_anamoly,
              data = DT_train_y)
    mod3 = lm(obs_prect_anamoly ~ mf_prect_anamoly,
              data = DT_train_y)
    mod4 = lm(obs_prect_anamoly ~ norcpm_prect_anamoly,
              data = DT_train_y)
    mod5 = lm(obs_prect_anamoly ~ ecmwf_prect_anamoly,
              data = DT_train_y)
    mod6 = lm(obs_prect_anamoly ~ ukmo_prect_anamoly + mf_prect_anamoly,
              data = DT_train_y)
    for(j in 1:6){
        DT_y[, eval(paste0("model", j)) := obs_prect_climatology + predict(get(paste0("mod",j)), newdata = DT_y)]
    }

    for(m in 1:length(months)){
        DT_train_my = DT_train_y[month == months[m]]
        mod1 = lm(obs_prect_anamoly ~ ukmo_prect_anamoly + mf_prect_anamoly + norcpm_prect_anamoly + ecmwf_prect_anamoly,
                  data = DT_train_my)
        mod2 = lm(obs_prect_anamoly ~ ukmo_prect_anamoly,
                  data = DT_train_my)
        mod3 = lm(obs_prect_anamoly ~ mf_prect_anamoly,
                  data = DT_train_my)
        mod4 = lm(obs_prect_anamoly ~ norcpm_prect_anamoly,
                  data = DT_train_my)
        mod5 = lm(obs_prect_anamoly ~ ecmwf_prect_anamoly,
                  data = DT_train_my)
        mod6 = lm(obs_prect_anamoly ~ ukmo_prect_anamoly + mf_prect_anamoly,
                  data = DT_train_my)
        for(j in 1:6){
            DT_y[month == months[m], eval(paste0("model_monthly_", j)) := obs_prect_climatology + predict(get(paste0("mod",j)), newdata = DT_y)]
        }
    }
    DT_y_all[[i]] = DT_y
}

gg = function(x,y){return(sqrt(mean( (x - y)^2)))}

DT_y = rbindlist(DT_y_all)
mods = c("obs_prect_climatology", paste0("model",1:6), paste0("model_monthly_",1:6))
Scores = DT_y[,lapply(.SD,gg,obs_erai_prect),.SDcols = mods ,month]
