rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)

setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./Fc_201812/DT_merged_precip.RData")
setkey(DT_precip, "Lon", "Lat", "month")
months = DT_precip[,unique(month)]

model_names = c("ukmo",
                "ecmwf",
                "dwd",
                "cmcc",
                "mf")

for(j in 1:length(model_names)){
    nms_v = paste0(model_names[j], "_climatology")
    nms_mean = paste0(model_names[j], "_precip")
    DT_precip[, eval(nms_v) := median(get(nms_mean), na.rm = TRUE),
           .(Lon, Lat, month)]
}

DT_precip[, obs_climatology := median(obs_precip), .(Lon, Lat, month)]
DT_precip[, obs_anomaly := obs_precip - obs_climatology, .(Lon, Lat, month)]
DT_precip[, anamoly_sd := sd(obs_anomaly), .(Lon, Lat, month)]

for(j in 1:length(model_names)){
    nms_v = paste0(model_names[j], "_anomaly")
    nms_v1 = paste0(model_names[j], "_precip")
    nms_v2 = paste0(model_names[j], "_climatology")
    DT_precip[, eval(nms_v) := get(nms_v1) - get(nms_v2)]
}

ff = paste0("obs_anomaly ~ ", paste0(model_names,"_anomaly", collapse = " +  "))
mod = qgam::qgam(as.formula(ff), data = DT_precip[month == 1], qu = 0.5)

DT_precip[,anomaly_hat := predict(mod, data = DT_t2m)]
DT_precip[,residual_hat := obs_anomaly - anomaly_hat]
DT_precip[, sd_residual := sd(residual_hat), .(Lon, Lat, month)]

load("./Fc_201812/Forecast_precip_processed.RData")

DT_forecast[, mf_bar := meteo_france_bar]
DT_forecast[, meteo_france_bar := NULL]
nms_aux = c("Lon", "Lat", "month", paste0(model_names, "_climatology"), "sd_residual","obs_climatology")
DT_aux = unique(DT_precip[,.SD, .SDcols = nms_aux])

DT_forecast_featured = merge(DT_forecast, DT_aux, by = c("Lon", "Lat", "month"))

for(j in 1:length(model_names)){
    nms_v = paste0(model_names[j], "_anomaly")
    nms_v1 = paste0(model_names[j], "_bar")
    nms_v2 = paste0(model_names[j], "_climatology")
    DT_forecast_featured[, eval(nms_v) := get(nms_v1) - get(nms_v2)]
}

DT_forecast_featured[, obs_anomaly_hat := predict(mod, newdata = DT_forecast_featured)]
DT_forecast_featured[, forecast_precip := obs_climatology + obs_anomaly_hat]

save(DT_forecast_featured, file = "./Fc_201812/DT_precip_forecast.RData")



