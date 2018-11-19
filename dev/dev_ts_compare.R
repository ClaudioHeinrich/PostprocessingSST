## Where I compare different statistical models
rm(list = ls())

library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "/home/alex/NR/SFE/Presentations/20181120_TLT/fig_alex/"
print_figs = FALSE

load("./FcNov2018/ts_hindcast_slimmed.RData")

ens_names = c("norcpm_ts_bar", "ecmwf_ts_bar","mf_ts_bar","ukmo_ts_bar")
for(y in 2006:2017){

    print(y)
    DT_train = DT_final[year < y]
    mod1 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly, data = DT_train)
    mod2 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
    mod3 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + mf_anamoly, data = DT_train)
    mod4 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ukmo_anamoly, data = DT_train)
    mod5 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    for(j in 1:5){
        DT_final[year == y,
                 eval((paste0("pred",j))) := climatology + predict(get(paste0("mod",j)),
                                                                   newdata = DT_final[year == y])]
    }
}

gg = function(a,b){return( sqrt(mean( (a - b)^2)))}

Lat_nordic = c(55,80)
Lon_nordic = c(5,30)
Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(pred5),
                 lapply(.SD,gg,obs_erai_ts),
                 .SDcols = paste0("pred",1:5)]

names(Score) = c("NorCPM","ECMWF","MF","UKMO","Multi")


Score_month = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(pred5),
                 lapply(.SD,gg,obs_erai_ts),month,
                 .SDcols = paste0("pred",1:5)]
names(Score_month) = c("month","NorCPM","ECMWF","MF","UKMO","Multi")


Score_year = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(pred5),
                 lapply(.SD,gg,obs_erai_ts),year,
                 .SDcols = paste0("pred",1:5)]
names(Score_year) = c("year","NorCPM","ECMWF","MF","UKMO","Multi")


if(print_figs){pdf(paste0(path_out,"/Models_compare.pdf"))}else{X11()}
rr = range(Score_month[,-1])
plot(1:12,Score_month[[2]], type  = "l", xlab = "Month", ylab = "RMSE", ylim = rr)
points(1:12,Score_month[[2]], pch = 20)
for(j in 3:5){
    lines(1:12,Score_month[[j]], col=j - 1)
    points(1:12,Score_month[[j]], pch = 20, col = j - 1)
}
legend("top", lty = 1, pch = 20,col = 1:4, legend = c("NorCPM","ECMWF","MF","UKMO"))
if(print_figs)dev.off()
