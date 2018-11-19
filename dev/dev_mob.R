rm(list = ls())

library(data.table)
library(party)

setwd("~/PostClimDataNoBackup/SFE/")

load("./FcNov2018/ts_hindcast_slimmed.RData")

ctl = mob_control(verbose = FALSE)
for(y in 2006:2017){
    print(y)
    DT_train = DT_final[year < y]
    mod = mob(obs_anamoly ~ ecmwf_anamoly + ukmo_anamoly | Lon + Lat,
              data = DT_train,
              control = ctl)
    mod2 = lm(obs_anamoly ~  ecmwf_anamoly, data = DT_train)
    DT_final[year == y, pred:=climatology + predict(mod, newdata = DT_final[year == y])]
    DT_final[year == y, pred2:=climatology + predict(mod2, newdata = DT_final[year == y])]
}

Lat_nordic = c(55,80)
Lon_nordic = c(5,30)
gg = function(a,b){return( sqrt(mean( (a - b)^2)))}
Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(pred2),
                 lapply(.SD,gg,obs_erai_ts),
                 month,
                 .SDcols = c("climatology", "pred","pred2")]



