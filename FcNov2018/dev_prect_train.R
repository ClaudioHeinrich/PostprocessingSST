## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)

setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./FcNov2018/prect_processed.RData")

## First run a backtest
ens_names = c("ecmwf_prect_bar")
DT_train[,equal_weighted := rowMeans(.SD), .SDcol = ens_names]
DT_train[,persistence := obs_anamoly_1]

months = c(11,12,1,2,3,4)
mod_list = list()
stat_model_names = c("ECMWF","Multi")
mod = list()
for(y in 2003:2015){
    print(y)
    DT_y = DT_train[year < y]
    mod[[1]] = mgcv::gam(obs_anamoly ~  ecmwf_anamoly, data = DT_y)
##    mod[[2]] = mgcv::gam(obs_anamoly ~  ukmo_anamoly, data = DT_y)
    mod[[2]] = mgcv::gam(obs_anamoly ~  ecmwf_anamoly, data = DT_y)
    for(j in 1:length(mod)){
        DT_train[year == y, eval(paste0(stat_model_names[j])) := climatology +  predict(mod[[j]],newdata = DT_train[year == y])]
    }
    for(m in months){
        DT_my = DT_train[year < y & month == m]
        mod[[1]] = mgcv::gam(obs_anamoly ~  ecmwf_anamoly, data = DT_my)
##        mod[[2]] = mgcv::gam(obs_anamoly ~  ukmo_anamoly, data = DT_my)
        mod[[2]] = mgcv::gam(obs_anamoly ~  ecmwf_anamoly, data = DT_my)
        for(j in 1:length(mod)){
            DT_train[year == y & month == m,
                     eval(paste0(stat_model_names[j],"_monthly")) := climatology + predict(mod[[j]],newdata = DT_train[year == y & month == m])]
        }
    }
}


gg = function(a,b){return( sqrt(mean( (a - b)^2)))}

Lat_nordic = c(58,63)
Lon_nordic = c(5,11)
##Lat_nordic = c(50,70)
##Lon_nordic = c(-10,11)
cols_all = c("climatology", "persistence", "equal_weighted", stat_model_names, paste0(stat_model_names,"_monthly"))
Score_prect = DT_train[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(Multi),
                       lapply(.SD,gg,obs_erai_prect),
                       month,
                       .SDcols = cols_all]


mod_prect = lm(obs_anamoly ~  ecmwf_anamoly, data = DT_train)


DT_fit_prect = DT_forecast[, .(Lon, Lat, month, year, ecmwf_anamoly, ukmo_anamoly, q_low, q_high, climatology)]
DT_fit_prect[,model_mean := climatology + predict(mod_prect, newdata = DT_fit_prect)]

for(s in 1:99){
    print(s)
    mod = rq(obs_anamoly ~ ecmwf_anamoly, data = DT_train, tau = s/100)
    DT_fit_prect[,eval(paste0("q_",s)) := climatology + predict(mod, newdata = DT_fit_prect)]
}

save(Score_prect, DT_fit_prect,file = "./FcNov2018/Forecast_ts.RData")
