rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
library(qgam)

setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./FcNov2018/models_prect.RData")
DT_forecast = readr::read_rds("./FcNov2018/wfc_prect_forecast.rds")

months = c(12,1,2,3,4)
q_print = c(10,25,33,50,67,75,90)

nms = c("lat", "lon", "Lat", "Lon", "month", paste0(c("ukmo","mf","ecmwf","obs"),"_prect_climatology"), paste0("q_", q_print))
DT_forecast = unique(merge(DT_forecast,
                           DT_train[,.SD,.SDcols = nms],
                           by = c("lat","lon","month"),
                           all = FALSE))

DT_forecast[, ukmo_prect_anamoly := ukmo_prect_bar - ukmo_prect_climatology, .(lat, lon, month)]
DT_forecast[, mf_prect_anamoly := mf_prect_bar - mf_prect_climatology, .(lat, lon, month)]
DT_forecast[, ecmwf_prect_anamoly := ecmwf_prect_bar - ecmwf_prect_climatology, .(lat, lon, month)]

DT_forecast[,model_mean := obs_prect_climatology + predict(mod_mean, newdata = DT_forecast)]
for(i in 1:length(months)){
    m = months[i]
    for(s in 1:9){
        mod = mod_q[[i]][[s]]
        DT_forecast[month == m, eval(paste0("q_hat_",s)) := obs_prect_climatology + predict(mod, newdata = DT_forecast[month == m])]
    }
}

rank1 = function(x){return(rank(x)[1])}
for(i in 1:length(q_print)){
    q = q_print[i]
    nm = paste0("p_below_",i)
    DT_forecast[,eval(nm) := apply(.SD,1,"rank1") / 10,
                .SDcols = c(paste0("q_",q), paste0("q_hat_",1:9))]    
}

save(DT_forecast,file = "./FcNov2018/Forecast_prect.RData")
