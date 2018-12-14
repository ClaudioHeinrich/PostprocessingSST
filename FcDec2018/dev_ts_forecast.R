rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
library(qgam)

setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./FcNov2018/models_ts.RData")
DT_forecast = readr::read_rds("./FcNov2018/wfc_ts_forecast2.rds")

months = c(12,1,2,3,4)


nms = c("lat", "lon", "Lat", "Lon", "month", paste0(c("ukmo","mf","ecmwf","obs"),"_ts_climatology"), paste0("q_", q_print))
DT_forecast = unique(merge(DT_forecast,
                           DT_train[,.SD,.SDcols = nms],
                           by = c("lat","lon","month"),
                           all = FALSE))

DT_forecast[, ukmo_ts_anamoly := ukmo_ts_bar - ukmo_ts_climatology, .(lat, lon, month)]
DT_forecast[, mf_ts_anamoly := mf_ts_bar - mf_ts_climatology, .(lat, lon, month)]
DT_forecast[, ecmwf_ts_anamoly := ecmwf_ts_bar - ecmwf_ts_climatology, .(lat, lon, month)]

DT_forecast[,model_mean := obs_ts_climatology + predict(mod_mean, newdata = DT_forecast)]
for(i in 1:length(months)){
    m = months[i]
    for(s in 1:99){
        mod = mod_q[[i]][[s]]
        DT_forecast[month == m, eval(paste0("q_hat_",s)) := obs_ts_climatology + predict(mod, newdata = DT_forecast[month == m])]
    }
}

rank1 = function(x){return(rank(x)[1])}
for(i in 1:length(q_print)){
    q = q_print[i]
    nm = paste0("p_below_",i)
    DT_forecast[,eval(nm) := apply(.SD,1,"rank1") / 1e2,
                .SDcols = c(paste0("q_",q), paste0("q_hat_",1:99))]    
}

save(DT_forecast,file = "./FcNov2018/Forecast_ts.RData")

