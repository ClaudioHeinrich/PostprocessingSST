rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")
load("./FcNov2018/ts_hindcast_cleaned.RData")
print_figs = TRUE

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

Lon_bounds = c(-10,40)
DT_clean = DT_clean[between(lon,Lon_bounds[1],Lon_bounds[2])]
DT_clean[,grid_id:=.GRP,.(lon,lat)]
DT_clean[,"Lon":= lon]
DT_clean[,"Lat":= lat]
DT_clean[,climatology := ff_mean(obs_erai_ts),.(grid_id,month)]
DT_clean[,ecmwf_climatology := ff_mean(ecmwf_ts_bar),.(grid_id,month)]
DT_clean[,ukmo_climatology := ff_mean(ukmo_ts_bar),.(grid_id,month)]



DT_clean[,ukmo_min:= apply(.SD,1,"min"),.SDcols = paste0("ukmo_ts",1:28)]
DT_clean[,ukmo_max:= apply(.SD,1,"max"),.SDcols = paste0("ukmo_ts",1:28)]
DT_clean[,ukmo_min_climatology := ff_mean(ukmo_min),.(grid_id,month)]
DT_clean[,ukmo_max_climatology := ff_mean(ukmo_max),.(grid_id,month)]
DT_clean[,ukmo_min_anamoly := ukmo_min - ukmo_min_climatology,.(grid_id,month)]
DT_clean[,ukmo_max_anamoly := ukmo_max - ukmo_max_climatology,.(grid_id,month)]

X11()
hist(DT_clean[month == 12 & year == 2010,obs_anamoly], breaks = 1e2, col="grey30", border = "grey30", main = "Anamoly Dec 2010", xlab = "Observed Anamoly")

X11()
plot_smooth(DT_clean[month == 12 & year == 2010, .(Lon,Lat,obs_erai_ts)], exclude_ocean = TRUE)


rr_1 = max(abs(range(DT_clean[month == 12 & year > 2005, obs_anamoly])))
rr_2 = max(abs(range(DT_clean[month == 12 & year > 2005, ukmo_anamoly], na.rm=TRUE)))
rr_3 = max(abs(range(DT_clean[month == 12 & year > 2005, ukmo_min_anamoly], na.rm=TRUE)))
rr_4 = max(abs(range(DT_clean[month == 12 & year > 2005, ukmo_max_anamoly], na.rm= TRUE)))
pdf("~/Dec_temps.pdf")
for(y in 2006:2015){
##    X11()
    par(mfrow = c(2,2))
    plot_smooth(DT_clean[month == 12 & year == y, .(Lon,Lat,obs_anamoly)], exclude_ocean = TRUE, mn = paste0("Obs Anamoly ",y), rr = c(-rr_1,rr_1))
    plot_smooth(DT_clean[month == 12 & year == y, .(Lon,Lat,ukmo_anamoly)], exclude_ocean = TRUE, rr = c(-rr_2,rr_2), mn = "UKMO Anamoly")
    plot_smooth(DT_clean[month == 12 & year == y, .(Lon,Lat,ukmo_min_anamoly)], exclude_ocean = TRUE, mn="UKMO Min Anamoly", rr = c(-rr_3,rr_3))
    plot_smooth(DT_clean[month == 12 & year == y, .(Lon,Lat,ukmo_max_anamoly)], exclude_ocean = TRUE, mn="UKMO Max Anamoly", rr = c(-rr_4,rr_4))
}
dev.off()

mod = mgcv::gam(ukmo_anamoly ~ te(Lon,Lat), data = DT_clean[month == 12])
mod2 = mgcv::gam(obs_anamoly ~ te(Lon,Lat), data = DT_clean[month == 12])
mod3 = mgcv::gam(ecmwf_anamoly ~ te(Lon,Lat), data = DT_clean[month == 12])
DT_clean[month == 12,ukmo_climate_trend:=predict(mod, newdata = DT_clean[month==12])]
DT_clean[month == 12,obs_climate_trend:=predict(mod2, newdata = DT_clean[month==12])]
DT_clean[month == 12,ecmwf_climate_trend:=predict(mod3, newdata = DT_clean[month==12])]

X11();plot_smooth(DT_clean[month ==12 & year == 2010, .(Lon,Lat,obs_climate_trend)], rr = c(-1,1),exclude_ocean = TRUE)
X11();plot_smooth(DT_clean[month ==12 & year == 2010, .(Lon,Lat,ukmo_climate_trend)], rr = c(-1,1),exclude_ocean = TRUE)
X11();plot_smooth(DT_clean[month ==12 & year == 2010, .(Lon,Lat,ecmwf_climate_trend)], rr = c(-1,1),exclude_ocean = TRUE)

X11()
par(mfrow = c(2,2))
plot_smooth(DT_clean[month == 12 & year == 2011, .(Lon,Lat,obs_anamoly)], exclude_ocean = TRUE, mn = "Obs Anamoly 2011")
plot_smooth(DT_clean[month == 12 & year == 2011, .(Lon,Lat,ukmo_anamoly)], exclude_ocean = TRUE)
plot_smooth(DT_clean[month == 12 & year == 2011, .(Lon,Lat,ukmo_min_anamoly)], exclude_ocean = TRUE, rr = c(-4.5,4.5), mn="UKMO Min Anamoly")
plot_smooth(DT_clean[month == 12 & year == 2011, .(Lon,Lat,ukmo_max_anamoly)], exclude_ocean = TRUE, rr = c(-4.5,4.5), mn="UKMO Max Anamoly")

X11()
plot_smooth(DT_clean[month == 12 & year == 2011, .(Lon,Lat,obs_anamoly)], exclude_ocean = TRUE, mn = "Obs Anamoly 2011")

X11()
plot_smooth(DT_clean[month == 12 & year == 2010, .(Lon,Lat,ecmwf_anamoly)], exclude_ocean = TRUE)



X11()
plot_smooth(DT_clean[month == 12 & year == 2011, .(Lon,Lat,ukmo_min_anamoly)], exclude_ocean = TRUE, rr = c(-4.5,4.5))


X11()
plot(DT_clean[month == 12, mean(ukmo_max - ukmo_min),year], pch = 20, 
lines(DT_clean[month == 12, mean(ukmo_max - ukmo_min),year], lty = 1



)


plot(DT_clean[grid_id == 208 & year > 1997][,.(ukmo_anamoly, obs_anamoly)], xlim = c(-8,8), ylim = c(-8,8))
mod = lm(obs_anamoly ~ ukmo_anamoly + ukmo_min_anamoly + ecmwf_anamoly, data = DT_clean[grid_id == 208 & year > 1997])


Lat_nordic = c(58,63)
Lon_nordic = c(5,11)
DT_nordic = DT_clean[between(Lat, Lat_nordic[1], Lat_nordic[2]) & between(Lon,Lon_nordic[1], Lon_nordic[2])]

png("~/Dec_obs.png")
plot(DT_nordic[month == 12,mean(obs_erai_ts),year], xlab = "Year",ylab = "Mean Temperature (ERA-I)", main = "Nordic", type="l", xlim = c(1997,2015))
points(DT_nordic[month == 12,mean(obs_erai_ts),year], pch = 20)
lines(DT_nordic[month == 12,mean(climatology),year], lty = 2, col="grey20")
legend("bottomleft",lty = 2, legend = "Running Climatology")
dev.off()
