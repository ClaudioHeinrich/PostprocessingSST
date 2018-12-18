rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)


setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./Fc_201812/DT_merged.RData")
setkey(DT_t2m, "Lon", "Lat", "month")
months = DT_t2m[,unique(month)]

model_names = c("ukmo",
                "ecmwf",
                "dwd",
                "cmcc",
                "mf")

for(j in 1:length(model_names)){
    nms_v = paste0(model_names[j], "_climatology")
    nms_mean = paste0(model_names[j], "_t2m")
    DT_t2m[, eval(nms_v) := median(get(nms_mean), na.rm = TRUE),
           .(Lon, Lat, month)]
}

DT_t2m[, obs_climatology := median(obs_t2m), .(Lon, Lat, month)]
DT_t2m[, obs_anomaly := obs_t2m - obs_climatology, .(Lon, Lat, month)]
DT_t2m[, anamoly_sd := sd(obs_anomaly), .(Lon, Lat, month)]

for(j in 1:length(model_names)){
    nms_v = paste0(model_names[j], "_anomaly")
    nms_v1 = paste0(model_names[j], "_t2m")
    nms_v2 = paste0(model_names[j], "_climatology")
    DT_t2m[, eval(nms_v) := get(nms_v1) - get(nms_v2)]
}

ff = paste0("obs_anomaly ~ ", paste0(model_names,"_anomaly", collapse = " +  "))
mod = lm(as.formula(ff), data = DT_t2m)

DT_t2m[,anomaly_hat := predict(mod, data = DT_t2m)]
DT_t2m[,residual_hat := obs_anomaly - anomaly_hat]
DT_t2m[, sd_residual := sd(residual_hat), .(Lon, Lat, month)]

load("./Fc_201812/Forecast_ts_processed.RData")

nms_aux = c("Lon", "Lat", "month", paste0(model_names, "_climatology"), "sd_residual")
DT_aux = unique(DT_t2m[,.SD, .SDcols = nms_aux])

DT_forecast_featured = merge(DT_forecast, DT_aux, by = c("Lon", "Lat", "month"))

model_names = c("ukmo",
                "ecmwf",
                "dwd",
                "cmcc",
                "meteo_france")

for(j in 1:length(model_names)){
    nms_v = paste0(model_names[j], "_anomaly")
    nms_v1 = paste0(model_names[j], "_bar")
    nms_v2 = paste0(model_names[j], "_climatology")
    DT_forecast_featured[, eval(nms_v) := get(nms_v1) - get(nms_v2)]
}

DT_forecast_featured[, obs_anomaly_hat := predict(mod, newdata = DT_forecast_featured)]




