rm(list = ls())

library(PostProcessing)
library(data.table)
library(qgam)

##--- Parameters ----
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"

##---- Setup Training Data ---
load("./FcNov2018/ts_processed.RData")
months = c(12,1,2,3,4)
DT_train[, Lon := lon]
DT_train[, Lat := lat]
DT_train[, obs_ts_climatology := mean(obs_erai_ts), .(lat, lon, month)]
DT_train[, obs_ts_anamoly := obs_erai_ts - obs_ts_climatology]

DT_train[, ukmo_ts_climatology := mean(ukmo_ts_bar), .(lat, lon, month)]
DT_train[, mf_ts_climatology := mean(mf_ts_bar), .(lat,lon,month)]
DT_train[, ecmwf_ts_climatology := mean(ecmwf_ts_bar), .(lat,lon,month)]
DT_train[, obs_ts_climatology := mean(obs_erai_ts), .(lat,lon,month)]

DT_train[, ukmo_ts_anamoly := ukmo_ts_bar - ukmo_ts_climatology, .(lat, lon, month)]
DT_train[, mf_ts_anamoly := mf_ts_bar - mf_ts_climatology, .(lat, lon, month)]
DT_train[, ecmwf_ts_anamoly := ecmwf_ts_bar - ecmwf_ts_climatology, .(lat, lon, month)]
DT_train[, obs_ts_anamoly := obs_erai_ts - obs_ts_climatology, .(lat, lon, month)]

q_print = c(10,25,33,50,67,75,90)
for(j in 1:length(q_print)){
    DT_train[,eval(paste0("q_",q_print[j])) := quantile(obs_erai_ts, q_print[j] / 1e2), .(grid_id,month)]
}

## ---- Train Quantile Models ------------
mod_q = list()
for(i in 1:length(months)){
    m = months[i]
    mod_q[[i]] = list()
    for(s in 1:99){
        print(paste(m,s,Sys.time()))
        mod = qgam::qgam(obs_ts_anamoly ~ te(Lon,Lat) + ecmwf_ts_anamoly, qu = s / 1e2, data = DT_train[month == m])
        mod_q[[i]][[s]] = mod
        DT_train[month == m, eval(paste0("q_hat_",s)) := obs_ts_climatology + predict(mod, newdata = DT_train[month == m])]
    }
}

## --- Train Mean Model(s) --------
mod_mean = lm(obs_ts_anamoly ~ ecmwf_ts_anamoly,
              data = DT_train)
DT_train[, model_mean := obs_ts_climatology + predict(mod_mean, newdata = DT_train)]

##--- Write ------
save(DT_train, mod_q, mod_mean, file = "./FcNov2018/models_ts.RData")

