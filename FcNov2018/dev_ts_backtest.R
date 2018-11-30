## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./FcNov2018/ts_processed.RData")

vintage_years = DT_train[,unique(cds_vintageyear)]

months = c(12,1,2,3,4)
mod_list = list()
mod = list()
DT_y_all = list()
for(i in 1:length(vintage_years)){
    y = vintage_years[i]
    DT_train_y = DT_train[cds_vintageyear != y]

    DT_train_y[, N := .N, .(lat, lon, month)]
    DT_train_y[, ukmo_ts_climatology := mean(ukmo_ts_bar), .(lat, lon, month)]
    DT_train_y[, mf_ts_climatology := mean(mf_ts_bar), .(lat,lon,month)]
    DT_train_y[, ecmwf_ts_climatology := mean(ecmwf_ts_bar), .(lat,lon,month)]
    DT_train_y[, norcpm_ts_climatology := mean(norcpm_ts_bar), .(lat,lon,month)]
    DT_train_y[, obs_ts_climatology := mean(obs_erai_ts), .(lat,lon,month)]

    nms = c("lat", "lon", "month", paste0(c("ukmo","mf","ecmwf","norcpm","obs"),"_ts_climatology"))
    DT_y = unique(merge(DT_train[cds_vintageyear == y],
                        DT_train_y[,.SD,.SDcols = nms],
                        by = c("lat","lon","month"),
                        all = FALSE))

    DT_train_y[, ukmo_ts_anamoly := ukmo_ts_bar - ukmo_ts_climatology, .(lat, lon, month)]
    DT_train_y[, mf_ts_anamoly := mf_ts_bar - mf_ts_climatology, .(lat, lon, month)]
    DT_train_y[, ecmwf_ts_anamoly := ecmwf_ts_bar - ecmwf_ts_climatology, .(lat, lon, month)]
    DT_train_y[, norcpm_ts_anamoly := norcpm_ts_bar - norcpm_ts_climatology, .(lat, lon, month)]
    DT_train_y[, obs_ts_anamoly := obs_erai_ts - obs_ts_climatology, .(lat, lon, month)]

    DT_y[, ukmo_ts_anamoly := ukmo_ts_bar - ukmo_ts_climatology, .(lat, lon, month)]
    DT_y[, mf_ts_anamoly := mf_ts_bar - mf_ts_climatology, .(lat, lon, month)]
    DT_y[, ecmwf_ts_anamoly := ecmwf_ts_bar - ecmwf_ts_climatology, .(lat, lon, month)]
    DT_y[, norcpm_ts_anamoly := norcpm_ts_bar - norcpm_ts_climatology, .(lat, lon, month)]
    DT_y[, obs_ts_anamoly := obs_erai_ts - obs_ts_climatology, .(lat, lon, month)]

    mod1 = lm(obs_ts_anamoly ~ ukmo_ts_anamoly + mf_ts_anamoly + norcpm_ts_anamoly + ecmwf_ts_anamoly,
              data = DT_train_y)
    mod2 = lm(obs_ts_anamoly ~ ecmwf_ts_anamoly, data = DT_train_y)
    for(j in 1:2){
        DT_y[, eval(paste0("model", j)) := obs_ts_climatology + predict(get(paste0("mod",j)), newdata = DT_y)]
    }
    for(m in months){
        mod_month = lm(obs_ts_anamoly ~ ukmo_ts_anamoly + mf_ts_anamoly + norcpm_ts_anamoly + ecmwf_ts_anamoly,
                       data = DT_train_y[month == m])
        mod_month2 = lm(obs_ts_anamoly ~ ecmwf_ts_anamoly,
                       data = DT_train_y[month == m])
        j = 1
        DT_y[month == m, eval(paste0("model_month", j)) := obs_ts_climatology + predict(mod_month,newdata = DT_y[month == m])]
        DT_y[month == m, eval(paste0("model_month", 2)) := obs_ts_climatology + predict(mod_month2,newdata = DT_y[month == m])]
    }
    nms_ensemble = c("ecmwf_ts_bar","norcpm_ts_bar","ukmo_ts_bar","mf_ts_bar")
    DT_y[, equal_weighted := rowMeans(.SD), .SDcols = nms_ensemble]
    nms_anamoly = c("ecmwf_ts_anamoly","norcpm_ts_anamoly","ukmo_ts_anamoly","mf_ts_anamoly")
    DT_y[, equal_weighted_anamoly := obs_ts_climatology + rowMeans(.SD), .SDcols = nms_anamoly]
    DT_y_all[[i]] = DT_y
}

gg = function(x,y){return(sqrt(mean( (x - y)^2)))}
Lat_nordic = c(58,63)
Lon_nordic = c(5,11)

mod_names = c("obs_ts_climatology", "equal_weighted", "equal_weighted_anamoly", "model1", "model_month1", "model2", "model_month2")
DT_y = rbindlist(DT_y_all)
Scores_ts = DT_y[between(lat, Lat_nordic[1], Lat_nordic[2]) & between(lon, Lon_nordic[1], Lon_nordic[2]),
              lapply(.SD,gg,obs_erai_ts),
              .SDcols = mod_names,
              month]

Scores_grid = DT_y[,
              lapply(.SD,gg,obs_erai_ts),
              .SDcols = mod_names,
              .(month, lon,lat)]

Scores_grid[, Rel := model1 / obs_ts_climatology]

Scores_grid[,Lon:=lon]
Scores_grid[,Lat:=lat]

X11();plot_smooth(Scores_grid[month == 12,.(Lon,Lat,Rel)], exclude_ocean = TRUE)
X11();plot_smooth(Scores_grid[month == 12,.(Lon,Lat,model_month1)],exclude_ocean = TRUE)

