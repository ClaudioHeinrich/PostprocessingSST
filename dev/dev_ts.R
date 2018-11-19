rm(list = ls())

library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

load("./FcNov2018/ts_hindcast_slimmed.RData")

for(y in 2006:2017){

    print(y)
    DT_train = DT_final[year < y]
    
    
    mod0 = lm(obs_anamoly ~  climatology, data = DT_train)
    mod1 = lm(obs_anamoly ~  obs_anamoly_1, data = DT_train)
    mod2 = lm(obs_anamoly ~  climatology + obs_anamoly_1, data = DT_train)
    mod3 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly, data = DT_train)
    mod4 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
    mod5 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + mf_anamoly, data = DT_train)
    mod6 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ukmo_anamoly, data = DT_train)
    mod7 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    mod8 = lm(obs_anamoly ~  ecmwf_anamoly, data = DT_train)
    
    for(j in 0:8){
        DT_final[year == y,
                 eval((paste0("pred",j))) := climatology + predict(get(paste0("mod",j)),
                                                                   newdata = DT_final[year == y])]
    }
}

gg = function(a,b){return( sqrt(mean( (a - b)^2)))}

Lat_nordic = c(55,80)
Lon_nordic = c(5,30)
Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(pred7),
                 lapply(.SD,gg,obs_erai_ts),
                 month,
                 .SDcols = c("climatology", paste0("pred",0:8))]

round(Score[month %in% c(11,12,1,2), .(month,climatology,"ECMWF"=pred4,"All"=pred7,"Just E"=pred8)],3)

Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(pred7),
                 lapply(.SD,gg,obs_erai_ts),
                 .(grid_id,Lon,Lat,month),
                 .SDcols = c("climatology", paste0("pred",0:8))]


DT_531 = DT_final[grid_id == 531 & month == 1]

plot(DT_531[,.(year, obs_erai_ts)], type = "l", xlab = "Year", ylab = "Month")
points(DT_531[,.(year, obs_erai_ts)], pch = 20)
lines(DT_531[,.(year, climatology)], lty = 2)
lines(DT_531[,.(year, ecmwf_ts_bar)], type = "l", xlab = "Year", ylab = "Month", col="red")
points(DT_531[,.(year, ecmwf_ts_bar)], pch = 20, col="red")
lines(DT_531[,.(year, ecmwf_climatology)], lty = 2, col="red")


rr = range(DT_531[,.(obs_anamoly, ecmwf_anamoly)])
plot(DT_531[,.(year, obs_anamoly)], type = "l", xlab = "Year", ylab = "Month", ylim=rr)
points(DT_531[,.(year, obs_anamoly)], pch = 20)
lines(DT_531[,.(year, ecmwf_anamoly)], type = "l", xlab = "Year", ylab = "Month", col="red")
points(DT_531[,.(year, ecmwf_anamoly)], pch = 20, col="red")
lines(DT_531[,.(year, norcpm_anamoly)], type = "l", xlab = "Year", ylab = "Month", col="blue")
points(DT_531[,.(year, norcpm_anamoly)], pch = 20, col="blue")
abline(h = 0, lty = 2, col="grey30")


