rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
library(qgam)

##--- Parameters ----
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
##---------------------

##---- Setup Training Data ---
load("./FcNov2018/prect_processed.RData")
months = c(12,1,2,3,4)
DT_train[, Lon := lon]
DT_train[, Lat := lat]
DT_train[, obs_prect_climatology := mean(obs_erai_prect), .(lat, lon, month)]
DT_train[, obs_prect_anamoly := obs_erai_prect - obs_prect_climatology]

DT_train[, ukmo_prect_climatology := mean(ukmo_prect_bar), .(lat, lon, month)]
DT_train[, mf_prect_climatology := mean(mf_prect_bar), .(lat,lon,month)]
DT_train[, ecmwf_prect_climatology := mean(ecmwf_prect_bar), .(lat,lon,month)]
DT_train[, obs_prect_climatology := mean(obs_erai_prect), .(lat,lon,month)]

DT_train[, ukmo_prect_anamoly := ukmo_prect_bar - ukmo_prect_climatology, .(lat, lon, month)]
DT_train[, mf_prect_anamoly := mf_prect_bar - mf_prect_climatology, .(lat, lon, month)]
DT_train[, ecmwf_prect_anamoly := ecmwf_prect_bar - ecmwf_prect_climatology, .(lat, lon, month)]
DT_train[, obs_prect_anamoly := obs_erai_prect - obs_prect_climatology, .(lat, lon, month)]

q_print = c(10,25,33,50,67,75,90)
for(j in 1:length(q_print)){
    DT_train[,eval(paste0("q_",q_print[j])) := quantile(obs_erai_prect, q_print[j] / 1e2), .(grid_id,month)]
}
##-------------------------------

## ---- Train Models ------------
mod_q = list()
for(i in 1:length(months)){
    m = months[i]
    mod_q[[i]] = list()
    for(s in 1:99){
        print(paste(m,s,Sys.time()))
        mod = qgam::qgam(obs_prect_anamoly ~ te(Lon,Lat) + ecmwf_prect_anamoly, qu = s / 1e2, data = DT_train[month == m])
        mod_q[[i]][[s]] = mod
        DT_train[month == m, eval(paste0("q_hat_",s)) := obs_prect_climatology + predict(mod, newdata = DT_train[month == m])]
    }
}
mod_mean = lm(obs_prect_anamoly ~ ecmwf_prect_anamoly,
              data = DT_train)
DT_train[, model_mean := obs_prect_climatology + predict(mod_mean, newdata = DT_train)]
##-------------------------------

save(DT_train, mod_q, mod_mean, file = "./FcNov2018/models_prect.RData")
