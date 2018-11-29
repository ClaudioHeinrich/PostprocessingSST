## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

load("./FcNov2018/prect_processed.RData")

vintage_years = DT_train[,unique(cds_vintageyear)]

months = c(12,1,2,3,4)
mod_list = list()
mod = list()
DT_y_all = list()
for(i in 1:length(vintage_years)){
    y = vintage_years[i]
    print(paste(y,Sys.time()))
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
                        DT_train_y[,.SD,.SDcols = nms],
                        by = c("lat","lon","month"),
                        all = FALSE))

##    DT_train_y[, ukmo_ts_climatology := (N * ukmo_ts_climatology - ukmo_ts_bar) / (N - 1), .(lat, lon, month)]
##    DT_train_y[, mf_ts_climatology := (N * mf_ts_climatology - mf_ts_bar) / (N - 1), .(lat, lon, month)]
##    DT_train_y[, ecmwf_ts_climatology := (N * ecmwf_ts_climatology - ecmwf_ts_bar) / (N - 1), .(lat, lon, month)]
##    DT_train_y[, norcpm_ts_climatology := (N * norcpm_ts_climatology - norcpm_ts_bar) / (N - 1), .(lat, lon, month)]
##    DT_train_y[, obs_ts_climatology := (N * obs_ts_climatology - obs_erai_ts) / (N - 1), .(lat, lon, month)]

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
    mod7 = party::mob(obs_prect_anamoly ~ norcpm_prect_anamoly | lon + lat,
                      data = DT_train_y)
    mod5 = lm(obs_prect_anamoly ~ ecmwf_prect_anamoly,
              data = DT_train_y)
    mod6 = lm(obs_prect_anamoly ~ ukmo_prect_anamoly + mf_prect_anamoly,
              data = DT_train_y)
    for(j in 1:7){
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
            DT_y[month == months[m], eval(paste0("model_monthly_", j)) := obs_prect_climatology + predict(get(paste0("mod",j)), newdata = DT_y[month == months[m]])]
        }
    }
    DT_y_all[[i]] = DT_y
}

gg = function(x,y){return(sqrt(mean( (x - y)^2)))}
Lat_nordic = c(58,63)
Lon_nordic = c(5,11)

DT_y = rbindlist(DT_y_all)
mods = c("obs_prect_climatology", paste0("model",1:7), paste0("model_monthly_",1:6))
Scores = DT_y[between(lat, Lat_nordic[1], Lat_nordic[2]) & between(lon, Lon_nordic[1], Lon_nordic[2]),
              lapply(.SD,gg,obs_erai_prect),
              .SDcols = mods ,
              month]
