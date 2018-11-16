rm(list = ls())

library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

load("./FcNov2018/ts_hindcast_slimmed.RData")

DT_final[,ecmwf_anamoly_1:=shift(ecmwf_anamoly,1,NA,"lag")]
for(y  in 2003:2017){

    print(y)
    DT_train = DT_final[year < y & year > 1993]
    mod1 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
    mod2 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    mod3 = lm(obs_anamoly ~ climatology + obs_anamoly_1 + ecmwf_anamoly + ecmwf_anamoly_1, data = DT_train)
    
    for(j in 1:3){
        DT_final[year == y,
                 eval((paste0("pred",j))) := climatology + predict(get(paste0("mod",j)),
                                                                   newdata = DT_final[year == y])]
    }
}

gg = function(a,b){return( mean( (a - b)^2))}

Score = DT_final[year > 2002 & !is.na(pred3),
                 lapply(.SD,gg,obs_erai_ts),
                 .SDcols = c("climatology", paste0("pred",1:3)),month]

Score = DT_final[year > 2002 & !is.na(pred3),
                 lapply(.SD,gg,obs_erai_ts),
                 .SDcols = c("climatology", paste0("pred",1:3)),.(lon,lat)]
