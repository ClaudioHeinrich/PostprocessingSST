rm(list = ls())

##----- Set up -------
library(SeasonalForecasting)
library(geosphere)
setwd("~/NR/SFE/")
options(max.print = 5e3)
print_figs = FALSE
##--------------------

##- Load Merged Data --
dt = load_combined_wide()
##---------------------

##--- Some restrictions ---
nyc.loc = c(-74.0060,40.7128)
bergen.loc = c(5.3221, 60.3913)
##---------------------------

##------- Calculate path from NYC to bergen ---
d_full = distHaversine(bergen.loc, nyc.loc)
N_points = 1e3
delta_bearing = d_full/N_points
p_curr = bergen.loc
p_route = matrix(NA,N_points,2)
for(j in 1:N_points)
{
  b = bearing(p_curr, nyc.loc)
  p_route[j,] = destPoint(p_curr, b, delta_bearing)
  p_curr = p_route[j,]
}
p_route = rbind(bergen.loc,p_route,nyc.loc)
##---------------------------------------------

##----- Subset dt -----------
lon.box = range(p_route[,1])
lat.box = range(p_route[,2])
dt_box = dt[ (Lon > lon.box[1]) & (Lon < lon.box[2]) & (Lat > lat.box[1]) & (Lat < lat.box[2])]
grid_map = dt_box[YM==23821, .(Lon, Lat, grid_id)]
##---------------------------

##----- Find grid path ------
w_grid = rep(NA,dim(p_route)[1])
for(j in 1:dim(p_route)[1])
{
  d = distHaversine(p_route[j,],grid_map[,.(Lon,Lat)])
  w_grid[j] = grid_map[order(d)[1], grid_id]
}
grid_use = unique(w_grid)
##----------------------------

##----- Restrict path --------
dt_route = dt_box[grid_id %in% grid_use]
##----------------------------

##---- Plot route -----------
if(print_figs){pdf("./figures/bergen_nyc_route.pdf")}else{X11()}
plot(grid_map[,.(Lon,Lat)],pch=".", xlim = lon.box, ylim = lat.box)
map("world",add = TRUE)
points(rbind(nyc.loc,bergen.loc), col="blue",pch=19, cex=2)
segments(head(p_route[,1],-1), head(p_route[,2],-1),
          tail(p_route[,1],-1), tail(p_route[,2],-1),lty = 2, col="blue")
points(dt_route[,.(Lon,Lat)], pch = 20, cex = .25,col="red")
if(print_figs)dev.off()
##----------------------------

##----- Run Local Bias Correction --
setkey(dt_route,"grid_id", "month", "year")
dt_route[,Resid:=SST_bar- Ens_bar]
dt_route[,"MeanResid" := (cumsum(Resid) - Resid) / (year - min(year)+1),.(grid_id, month)]
#setkey(dt,"grid_id","month","year", "Ens")
for(j in 1:9)
{
  dt_route[,paste0("EnsHat",j) := get(paste0("Ens",j)) + MeanResid]
}
##------------------------------------

##----------------------------------
setkey(dt_route, "YM","grid_id")
##----------------------------

##----- Find mins along path ---
dt_min = dt_route[,lapply(.SD,"min",na.rm=TRUE),
                  YM,
                  .SDcols = c("year","month",paste0("SST", 1:10),paste0("EnsHat",1:9))]
dt_min[,SST_bar := rowMeans(dt_min[ , .SD, .SDcols = paste0("SST",1:10)])]
dt_min[,Ens_bar:=rowMeans(dt_min[, .SD, .SDcols = paste0("EnsHat",1:9)])]
dt_min[,ResidMin:=SST_bar - Ens_bar]
dt_min[,MeanResidMin:= (cumsum(ResidMin) - ResidMin) / (year - min(year) + 1), .(month)]
dt_min[,EnsMinCorrect:=Ens_bar + MeanResidMin]
dt_min[,FinalResid := SST_bar - EnsMinCorrect]
##--------------------------------

DD = dt_min[,.("year" = year,
               "one"= 1)]
DD[,sum:=cumsum(one)]
xlab = DD[,min(sum),year]
##------ Plot Series of Mins ----
rr= range(dt_min[,.(SST_bar,Ens_bar)])

if(print_figs){pdf("./figures/bergen_nyc_bias.pdf")}else{X11()}
plot(dt_min[,SST_bar], axes = FALSE, xlab = "date", ylab="Min SST Across Route",ylim = rr, type="l")
axis(1, at = xlab[,V1],label=xlab[,year])
axis(2)
lines(dt_min[,Ens_bar], col="red")
legend("topleft",lty = 1, col=c("blue","red"), legend = c("SST","Ensemble"))
if(print_figs)dev.off()
##-------------------------------

if(print_figs){pdf("./figures/bergen_nyc_bias.pdf")}else{X11()}
plot(dt_min[month==8, year],dt_min[month == 8,ResidMin], xlab = "date", ylab="Residual",ylim = rr, type="l")
lines(dt_min[month==2, year],dt_min[month == 2,ResidMin])
axis(1, at = xlab[,V1],label=xlab[,year])
axis(2)
lines(dt_min[,Ens_bar], col="red")
legend("topleft",lty = 1, col=c("blue","red"), legend = c("SST","Ensemble"))
if(print_figs)dev.off()

if(print_figs){pdf("./figures/bias_by_month.pdf")}else{X11()}
A = dt_min[,.(mean(ResidMin),mean(FinalResid)),month]
plot(A[,month], A[,V1], xlab="Month",ylab="Mean Bias", ylim = range(A[,.(V1,V2)]), type="l")
lines(A[,month], A[,V2], col="red")
legend("topleft",lty = 1, col=c("black","red"), legend = c("Ensemble","Corrected"))
if(print_figs)dev.off()
  
if(print_figs){pdf("./figures/bias_by_year.pdf")}else{X11()}
A = dt_min[,.(mean(ResidMin),mean(FinalResid)),year]
plot(A[,year], A[,V1], xlab="year",ylab="Mean Bias", ylim = range(A[,.(V1,V2)]), type="l")
lines(A[,year], A[,V2], col="red")
abline(h = 0, col="grey35", lty = 2)
legend("topleft",lty = 1, col=c("black","red"), legend = c("Ensemble","Corrected"))
if(print_figs)dev.off()

if(print_figs){pdf("./figures/rmse_by_month.pdf")}else{X11()}
A = dt_min[,.(sqrt(mean(ResidMin^2)),sqrt(mean(FinalResid^2))),month]
plot(A[,month], A[,V1], xlab="Month",ylab="RMSE", ylim = range(A[,.(V1,V2)]), type="l")
lines(A[,month], A[,V2], col="red")
legend("topleft",lty = 1, col=c("black","red"), legend = c("Ensemble","Corrected"))
if(print_figs)dev.off()
  
if(print_figs){pdf("./figures/rmse_by_year.pdf")}else{X11()}
A = dt_min[,.(sqrt(mean(ResidMin^2)),sqrt(mean(FinalResid^2))),year]

plot(A[,year], A[,V1], xlab="year",ylab="RMSE", ylim = range(A[,.(V1,V2)]), type="l")
lines(A[,year], A[,V2], col="red")
abline(h = 0, col="grey35", lty = 2)
legend("topleft",lty = 1, col=c("black","red"), legend = c("Ensemble","Corrected"))
if(print_figs)dev.off()
