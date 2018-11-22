## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./FcNov2018/ts_hindcast_slimmed.RData")

ens_names = c("ecmwf_ts_bar","ukmo_ts_bar")
y = 2015
DT_train = DT_final[year < y & month == 12]
DT_train[,year_0:=year - min(year) + 1]
mod = lm(obs_anamoly ~ year_0, data = DT_train)

mod2 = mgcv::gam(obs_anamoly ~ te(Lon,Lat), data = DT_train)
DT_test = DT_final[year == y & month == 12]
mm = predict(mod2, newdata = DT_test)
DT_test[,trend:=mm]
plot_smooth(DT_train[year == 2007,.(Lon,Lat,climatology)], exclude_ocean = TRUE)

DT_ab = DT_train[year == 2005 | year == 2015]

png("climatology_delta.png")
plot_smooth(DT_ab[,(tail(climatology,1) - head(climatology,1)),.(Lon,Lat)], exclude_ocean = TRUE, rr = c(-1.6,1.6),mn = "2015 Climatology - 2005 Climatology")
points(DT_ab[,.(Lon,Lat)], pch = 20)
dev.off()


mod = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ecmwf_anamoly + ukmo_anamoly, data = DT_train)
mod2 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ukmo_anamoly, data = DT_train)
mod3 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + ecmwf_anamoly, data = DT_train)
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

DT_europe = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005]
Score = DT_europe[month == 12,
                 lapply(.SD,gg,obs_erai_ts),
                 .(Lon,Lat),
                 .SDcols = c("climatology","equal_weighted","trained_model","trained_model2")]

png("~/rmse_climatology.png")
plot_smooth(Score[,.(Lon,Lat,climatology)],exclude_ocean = TRUE)
points(Score[,.(Lon,Lat)], pch = 20, cex = 1.2)
dev.off()

A = DT_europe[grid_id == 36694]

png("~/climatology_delta.png")
plot_smooth(DT_ab[,(tail(climatology,1) - head(climatology,1)),.(Lon,Lat)], exclude_ocean = TRUE, rr = c(-1.6,1.6),mn = "2015 Climatology - 2005 Climatology")
points(DT_ab[,.(Lon,Lat)], pch = 20, cex = 0.5)
dev.off()
