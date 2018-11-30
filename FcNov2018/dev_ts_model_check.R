rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
library(qgam)

setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

rho_tau = function(x, tau){
    return( (tau - 1) * x * (x < 0) + tau * x * (x > 0))
}

qs_loss = function(Y,mu,tau){
    return(sum(rho_tau(Y - mu, tau)))
}

load("./FcNov2018/ts_processed.RData")

DT_train[, Lon := lon]
DT_train[, Lat := lat]
DT_train = DT_train[month == 12]

years = DT_train[, unique(year)]
DT_test = list()
for(j in 1:length(years)){
    y = years[j]
    print(y)

    DT_train_y = DT_train[year != y]

    DT_train_y[ ,obs_q_75 := quantile(obs_erai_ts,.75), .(lon,lat,month)]
    DT_train_y[ ,obs_q_50 := quantile(obs_erai_ts,.5), .(lon,lat,month)]
    DT_train_y[ ,obs_q_25 := quantile(obs_erai_ts,.25), .(lon,lat,month)]
    DT_train_y[ ,obs_IQR := obs_q_75 - obs_q_25]
    DT_train_y[ ,obs_erai_ts_standard := (obs_erai_ts - obs_q_50) / obs_IQR]

    DT_train_y[,ukmo_q_75 := quantile(ukmo_ts_bar,.75), .(lon,lat,month)]
    DT_train_y[,ukmo_q_50 := quantile(ukmo_ts_bar,.5), .(lon,lat,month)]
    DT_train_y[,ukmo_q_25 := quantile(ukmo_ts_bar,.25), .(lon,lat,month)]
    DT_train_y[,ukmo_IQR := ukmo_q_75 - ukmo_q_25]
    DT_train_y[,ukmo_ts_standard := (ukmo_ts_bar - ukmo_q_50) / ukmo_IQR]

    DT_train_y[,mf_q_75 := quantile(mf_ts_bar,.75), .(lon,lat,month)]
    DT_train_y[,mf_q_50 := quantile(mf_ts_bar,.5), .(lon,lat,month)]
    DT_train_y[,mf_q_25 := quantile(mf_ts_bar,.25), .(lon,lat,month)]
    DT_train_y[,mf_IQR := mf_q_75 - mf_q_25]
    DT_train_y[,mf_ts_standard := (mf_ts_bar - mf_q_50) / mf_IQR]

    DT_test[[j]] = unique(merge(DT_train[year == y, .(lon,lat,month,year,ukmo_ts_bar,mf_ts_bar,obs_erai_ts)],
                           DT_train_y[,.(lon,lat,month,obs_q_50, obs_IQR,ukmo_q_50,ukmo_IQR, mf_q_50, mf_IQR)],
                           by = c("lon","lat","month"),
                           all = TRUE))
    mod = qgam(obs_erai_ts_standard ~ ukmo_ts_standard + mf_ts_standard, data = DT_train_y, qu = 0.5)
    DT_test[[j]][, ukmo_ts_standard := (ukmo_ts_bar - ukmo_q_50) / ukmo_IQR]
    DT_test[[j]][, mf_ts_standard := (mf_ts_bar - mf_q_50) / mf_IQR]
    DT_test[[j]][, model := predict(mod, newdata = DT_test[[j]]) * obs_IQR + obs_q_50]
}


DT_test = rbindlist(DT_test)
Lat_nordic = c(58,63)
Lon_nordic = c(5,11)

DT_test[between(lat, Lat_nordic[1], Lat_nordic[2]) & between(lon, Lon_nordic[1], Lon_nordic[2]),
        .(mean(abs(obs_erai_ts - model)),
          mean(abs(obs_erai_ts - obs_q_50)))]
    
plot(DT_train[month == 12, .(ukmo_ts_standard, obs_erai_ts_standard)], pch = ".", cex = 2)

Global_mean = DT_train[,mean(obs_
q_print = c(10,25,33,50,67,75,90)

A = DT_train[month == 12,.("q_25_coverage" = mean(obs_erai_ts < q_hat_25)),.(Lon,Lat)]
plot_smooth(A, exclude_ocean = TRUE, rr = c(0,.5))

X = DT_train[month == 12]

DT_train[,q_25_anomaly := q_25 - obs_ts_climatology]
DT_train[,q_50_anomaly := q_50 - obs_ts_climatology]
a = qgam(obs_ts_anamoly ~ q_25_anomaly + ecmwf_ts_anamoly, data = DT_train[month == 12], qu = 0.25)
a = qgam(obs_ts_anamoly ~ s(q_25_anomaly) + s(q_50_anomaly) + ecmwf_ts_anamoly, data = DT_train[month == 12], qu = 0.1)
DT_train[month == 12, q_10_hat := predict(a, newdata = DT_train[month == 12])]
plot(DT_train[month == 12,.(q_50_anomaly,q_10_hat)])




plot_smooth(DT_train[month == 12 & year == 2005,.(Lon,Lat,q_25_hat)], exclude_ocean = TRUE)
