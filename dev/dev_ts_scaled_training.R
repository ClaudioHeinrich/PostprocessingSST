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
DT_final[,year_0:= year - min(year)]
months = c(11,12,1,2,3,4)

for(y in 2006:2017){

    print(y)
    DT_train = DT_final[year < y]
    DT_train[,"sd_obs" := sd(obs_anamoly, na.rm=TRUE),.(Lon,Lat,month)]
    DT_train[,"sd_ecmwf" := sd(ecmwf_anamoly, na.rm=TRUE),.(Lon,Lat,month)]
    DT_train[,"sd_ukmo" := sd(ukmo_anamoly, na.rm=TRUE),.(Lon,Lat,month)]
    DT_train[,"obs_scaled":=obs_anamoly / sd_obs]
    DT_train[,"ecmwf_scaled":=ecmwf_anamoly / sd_ecmwf]
    DT_train[,"ukmo_scaled":=ukmo_anamoly / sd_ukmo]
    mod = mgcv::gam(obs_anamoly ~  climatology + obs_anamoly_1 + te(Lon,Lat) + ecmwf_anamoly + ukmo_anamoly, data = DT_train)
##    mod2 = mgcv::gam(obs_anamoly ~  te(Lon,Lat) + climatology + obs_anamoly_1 + ukmo_anamoly, data = DT_train)
##    mod3 = mgcv::gam(obs_anamoly ~  te(Lon,Lat) + climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
    DT_test = DT_final[year == y]
    DT_test = merge(DT_test,unique(DT_train[,.(Lon,Lat,month,sd_obs,sd_ecmwf,sd_ukmo)]), by = c("Lon","Lat","month"), all.y = FALSE)
    DT_test[,"ecmwf_scaled" := ecmwf_anamoly / sd_ecmwf]
    DT_test[,"ukmo_scaled" := ukmo_anamoly / sd_ukmo]
    DT_final[year == y,
             multi_model := climatology + DT_test[,sd_obs] * predict(mod,newdata = DT_test)]
}

gg = function(a,b){return( sqrt(mean( (a - b)^2)))}

Lat_nordic = c(58,63)
Lon_nordic = c(5,11)

Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(multi_model),
                 lapply(.SD,gg,obs_erai_ts),
                 month,
                 .SDcols = c("climatology","persistence","equal_weighted","multi_model")]


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

Score = DT_final[year > 2005 & !is.na(multi_model) & month == 12,
                 lapply(.SD,gg,obs_erai_ts),
                 .(Lon,Lat),
                 .SDcols = c("climatology","multi_model")]

png("~/rmse_climatology.png")
plot_smooth(Score[,.(Lon,Lat,climatology)],exclude_ocean = TRUE)
plot_smooth(Score[,.(Lon,Lat,multi_model)],exclude_ocean = TRUE)
points(Score[,.(Lon,Lat)], pch = 20, cex = 1.2)
dev.off()

A = DT_final[grid_id == 36694 & month == 12]


DT_final[,residual:=obs_erai_ts - multi_model]
DT_sd = DT_final[month == 12,sd(residual,na.rm=TRUE), .(Lon,Lat)]

r_m = max(abs(range(DT_final[,residual], na.rm= TRUE)))
rr = c(-r_m,r_m)
for(y in 2006:2015){
    X11()
    plot_smooth(DT_final[year == y & month == 12,residual,.(Lon,Lat)], exclude_ocean = TRUE,mn = y, rr = rr)
}
