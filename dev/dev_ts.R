rm(list = ls())

library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

<<<<<<< HEAD
load("./FcNov2018/ts_hindcast_slimmed.RData")
=======
ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

DT = readr::read_rds("./FcNov2018/ts_hindcast.rds")

setkey(DT,"lon","lat","month","year")

DT_obs = DT[is.finite(obs_erai_ts)]

DT_small = DT_obs[,.(lon,lat,norcpm_ts_bar,
                     ecmwf_ts_bar,
                     mf_ts_bar,
                     ukmo_ts_bar,
                     obs_erai_ts,
                     year,
                     month)]

DT_final = DT_small[!(is.na(ecmwf_ts_bar) | is.na(mf_ts_bar) | is.na(ukmo_ts_bar))]

DT_final[,grid_id:=.GRP,.(lon,lat)]
DT_final[,N:=.N,grid_id]
DT_final = DT_final[N ==268]
DT_final[,climatology := ff_mean(obs_erai_ts),.(grid_id,month)]
DT_final[,norcpm_climatology := ff_mean(norcpm_ts_bar),.(grid_id,month)]
DT_final[,ecmwf_climatology := ff_mean(ecmwf_ts_bar),.(grid_id,month)]
DT_final[,mf_climatology := ff_mean(mf_ts_bar),.(grid_id,month)]
DT_final[,ukmo_climatology := ff_mean(ukmo_ts_bar),.(grid_id,month)]

DT_final[,obs_anamoly := obs_erai_ts - climatology]
DT_final[,norcpm_anamoly := norcpm_ts_bar - norcpm_climatology]
DT_final[,ecmwf_anamoly := ecmwf_ts_bar - ecmwf_climatology]
DT_final[,mf_anamoly := mf_ts_bar - mf_climatology]
DT_final[,ukmo_anamoly := ukmo_ts_bar - ukmo_climatology]
DT_final[,obs_anamoly_1 := shift(obs_anamoly,1,NA,"lag"),.(grid_id,month)]
>>>>>>> 17a98ee186fa79e7244767102712e22e421c6b83

for(y in 2006:2017){

    print(y)
    DT_train = DT_final[year < y]
    
    
    mod0 = lm(obs_anamoly ~  climatology, data = DT_train)
    mod1 = lm(obs_anamoly ~  obs_anamoly_1, data = DT_train)
    mod2 = lm(obs_anamoly ~  climatology + obs_anamoly_1, data = DT_train)
<<<<<<< HEAD
    mod3 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly, data = DT_train)
    mod4 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
    mod5 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + mf_anamoly, data = DT_train)
    mod6 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ukmo_anamoly, data = DT_train)
    mod7 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    mod8 = lm(obs_anamoly ~  ecmwf_anamoly, data = DT_train)
=======
    mod4 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly, data = DT_train)
    mod5 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
    mod6 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + mf_anamoly, data = DT_train)
    mod7 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ukmo_anamoly, data = DT_train)
    mod8 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    mod8 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
>>>>>>> 17a98ee186fa79e7244767102712e22e421c6b83
    
    for(j in 0:8){
        DT_final[year == y,
                 eval((paste0("pred",j))) := climatology + predict(get(paste0("mod",j)),
                                                                   newdata = DT_final[year == y])]
    }
}

<<<<<<< HEAD
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


=======
gg = function(a,b){return( mean( (a - b)^2))}

Score = DT_final[year > 2005 & !is.na(pred8),
                 lapply(.SD,gg,obs_erai_ts),
                 .SDcols = c("climatology", paste0("pred",0:8))]
>>>>>>> 17a98ee186fa79e7244767102712e22e421c6b83
