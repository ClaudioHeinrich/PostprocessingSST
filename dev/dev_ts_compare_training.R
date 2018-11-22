## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./FcNov2018/ts_hindcast_slimmed.RData")

ens_names = c("ecmwf_ts_bar","ukmo_ts_bar")
DT_final[,equal_weighted := rowMeans(.SD), .SDcol = ens_names]
DT_final[,persistence := obs_anamoly_1]
months = c(11,12,1,2,3,4)
for(y in 2003:2017){

    print(y)
    DT_train = DT_final[year < y]
    mod = mgcv::gam(obs_anamoly ~  te(Lon,Lat) + climatology + obs_anamoly_1 + ecmwf_anamoly + ukmo_anamoly, data = DT_train)
    mod2 = mgcv::gam(obs_anamoly ~  te(Lon,Lat) + climatology + obs_anamoly_1 + ukmo_anamoly, data = DT_train)
    mod3 = mgcv::gam(obs_anamoly ~  te(Lon,Lat) + climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
    DT_final[year == y,
             multi_model := climatology + predict(mod,newdata = DT_final[year == y])]
    DT_final[year == y,
             ukmo := climatology + predict(mod2,newdata = DT_final[year == y])]
    DT_final[year == y,
             ecmwf := climatology + predict(mod3,newdata = DT_final[year == y])]
}

gg = function(a,b){return( sqrt(mean( (a - b)^2)))}

Lat_nordic = c(58,63)
Lon_nordic = c(5,11)

Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(multi_model),
                 lapply(.SD,gg,obs_erai_ts),
                 month,
                 .SDcols = c("climatology","persistence","equal_weighted","multi_model","ukmo","ecmwf")]


if(print_figs){pdf(paste0(path_out,"/Model_to_Basic.pdf"))}else{X11()}
rr = range(Score[,-1])
plot(Score[[2]], type  = "l", xlab = "Month", ylab = "RMSE", axes = FALSE, ylim = rr)
axis(1, at = 1:dim(Score)[1],labels = Score[[1]])
axis(2)
points(Score[[2]], pch = 20)
for(j in 3:5){
    lines(Score[[j]], col=j - 1)
    points(Score[[j]], pch = 20, col = j - 1)
}
legend("bottomleft", lty = 1, pch = 20, legend = c("Climatology","Equal Weighted", "Trained All", "Trained Month"), col=1:4)
if(print_figs)dev.off()
DT_final = DT_final[month == 12]
DT_final[,residual:=obs_erai_ts - multi_model]
setkey(DT_final, "grid_id","month","year")
DT_final[year >=2006,running_mean:=cumsum(residual)/1:10,.(grid_id,month)]
DT_final[,mean_ar:=shift(running_mean,1,NA,"lag")]
DT_final[,model_ar:=multi_model + mean_ar]

Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(model_ar),
                 lapply(.SD,gg,obs_erai_ts),
                 year,
                 .SDcols = c("climatology","persistence","equal_weighted","multi_model","model_ar")]
