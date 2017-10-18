rm(list = ls())

##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)
print_figs = FALSE
##------------------------


##------- Specify desired vintages ---------


vintage.vec=c("mr")



num_vin = length(vintage.vec)

vintage_label = function (vin){
  if (vin == "mr") return("most recent vintage")
  if (vin == "2r") return("2nd recent vintage")
  if (vin == "3r") return("3rd recent vintage")
  if (vin == "4r") return("4th recent vintage")
  if (vin == "Jan") return("January vintage")
  if (vin == "Apr") return("April vintage")
  if (vin == "Jul") return("July vintage")
  if (vin == "Oct") return("October vintage")
  if (vin == "test") return("test vintage")
}
  


##------ Set up------
dt=list()
for(k in 1:num_vin)
{
  print(vintage.vec[k])
    dt_vin = load_combined(vintage = vintage.vec[k])
    dt_vin[,vintage := vintage.vec[k]]
    dt[[k]]=dt_vin
    }
dt = rbindlist(dt)
dt[,residual := SST_bar - Forecast]

##---------------------

##--- Calculate Running Bias -------
setkey(dt,"vintage","grid_id","Ens","month","year")
dt[,"MeanResid" := mean(residual),.(vintage, grid_id, month,year)]
#setkey(dt,"grid_id","Ens","month","year")
dt[,"MeanResidCum" := (cumsum(MeanResid) - MeanResid) / (year - min(year)+1),.(vintage, grid_id, Ens, month)]
#setkey(dt,"grid_id","month","year", "Ens")
dt[,"SST_hat_local":=Forecast + MeanResidCum]
dt[,"Obs_above_orig":=SST_bar > Forecast]
dt[,"Obs_above_local":=SST_bar > SST_hat_local]
dt[,"Rank_orig":=cumsum(Obs_above_orig),.(vintage, grid_id,year,month)]
dt[,"Rank_local":=cumsum(Obs_above_local),.(vintage, grid_id,year,month)]
dt[,"residual_local":=SST_bar - SST_hat_local]
##----------------------------------

dt_rank = dt[,.("Q_Orig"=max(Rank_orig) + 1,"Q_local"=max(Rank_local) + 1),.(vintage,grid_id,year,month)]

prop_orig = dt_rank[,.N,by=.(Q_Orig, vintage)]
prop_orig = prop_orig[!is.na(Q_Orig)]
prop_orig[,"Prop":=N/sum(N), by = "vintage"]

prop_local = dt_rank[,.N,by=.(Q_local, vintage)]
prop_local = prop_local[!is.na(Q_local)]
prop_local[,"Prop":=N/sum(N), by = "vintage"]

##---- Rank Histograms -------------------
col_vec = c("blue","darkgreen","darkmagenta","darkorange2")
col_vec <- col_vec[1:num_vin]

if(print_figs){pdf(paste0("./figures/calibration_",vintage.vec[num_vin],"_vin.pdf"))}else{X11()}
plot(prop_orig[vintage == vintage.vec[1],.(Q_Orig -.2 ,Prop)],type="h",lwd = 3, 
     xlab = "Rank of Observation in Ensemble", ylab = "Proportion", 
     xlim = c(.5 , max(prop_orig[,Q_Orig])+.5), main = "Rank histogram" )
leg = c("Raw, most recent vintage")
for(l in 1:num_vin){
    lines(prop_local[vintage == vintage.vec[l],.(Q_local + l*.1 - .2, Prop)], type="h", lwd = 3, col=col_vec[l])
    leg <- c(leg,vintage_label(vintage.vec[l]))
}
legend("topleft", lty = 1, col = c("black",col_vec),legend = leg) 
if(print_figs)dev.off()

##----------------------------------------

##----- Reduce Observation Space ---------
setkey(dt,"vintage","grid_id","YM","Ens")
dt_reduced = dt[,.("Lon" = head(Lon,1),
                   "Lat" = head(Lat,1),
                   "SST_bar" = head(SST_bar,1),
                   "MeanOrig" = mean(Forecast),
                   "MedianOrig" = median(Forecast),
                   "MeanLocal" = mean(SST_hat_local),
                   "MedianLocal" = median(SST_hat_local),
                   "year" = head(year,1),
                   "month" = head(month,1)),
                .(vintage, grid_id, YM)]
##------------------------------------------

##-------- Get Results -------------
results_ym = dt_reduced[,.("Bias"= mean(SST_bar - MeanOrig,na.rm=TRUE),
                   "RMSE" = sqrt(mean( (SST_bar - MeanOrig)^2, na.rm=TRUE)),
                   "MAE" = mean(abs(SST_bar - MedianOrig), na.rm=TRUE),
                   "Bias_Local"= mean(SST_bar - MeanLocal,na.rm=TRUE),
                   "RMSE_Local" = sqrt(mean( (SST_bar - MeanLocal)^2, na.rm=TRUE)),
                   "MAE_Local" = mean(abs(SST_bar - MedianLocal), na.rm=TRUE),
                   "Year"= year,
                   "Month" = min(month)),
                keyby = .(vintage,YM)]
minYear = results_ym[,.(min_y = min(Year)),by = vintage]
results_ym = results_ym[Year > max(minYear[,min_y])]
yy_all = results_ym[,range(Year)]
##-----------------------------------

##----- Plot -------------

for(k in 1:length(vintage.vec))
{
  print(k)
  if(print_figs){pdf(paste0("./figures/bias_ym_",vintage.vec[k],".pdf"))}else{X11()}
  rr = range(results_ym[,.(Bias,Bias_Local)])
  plot(results_ym[vintage == vintage.vec[1],.(YM,Bias)], type="l", axes = FALSE,
     xlab="Year",ylab="Global Mean Bias (C)", ylim = rr*c(1,1.1))
     leg = c("Raw, most recent vintage")
  for(l in 1:k)
  {
    print(l)
    lines(results_ym[vintage == vintage.vec[l],.(YM,Bias_Local)],col=col_vec[l])  
    leg = c(leg,paste0("Local, ",vintage_label(vintage.vec[l])))
  }
  axis(2)
  axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
  legend("topright", lty = 1, col = c("black",col_vec), legend = leg)
if(print_figs)dev.off()
}


for(k in 1:length(vintage.vec))
{
  print(k)
  if(print_figs){pdf(paste0("./figures/RMSE_ym_",vintage.vec[k],".pdf"))}else{X11()}
  rr = range(results_ym[,.(RMSE,RMSE_Local)])
  plot(results_ym[vintage == vintage.vec[1],.(YM,RMSE)], type="l", axes = FALSE,
       xlab="Year",ylab="Global Mean RMSE (C)", ylim = rr*c(1,1.1))
  leg = c("Raw, most recent vintage")
  for(l in 1:k)
  {
    print(l)
    lines(results_ym[vintage == vintage.vec[l],.(YM,RMSE_Local)],col=col_vec[l])  
    abline(h = 1.5, lty = "dashed")
    leg = c(leg,paste0("Local, ",vintage_label(vintage.vec[l])))
  }
  axis(2)
  axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
  legend("topright", lty = 1, col = c("black",col_vec), legend = leg)
  if(print_figs)dev.off()
}

for(k in 1:length(vintage.vec))
{
  print(k)
  if(print_figs){pdf(paste0("./figures/MAE_ym_",vintage.vec[k],".pdf"))}else{X11()}
  rr = range(results_ym[,.(MAE,MAE_Local)])
  plot(results_ym[vintage == vintage.vec[1],.(YM,MAE)], type="l", axes = FALSE,
       xlab="Year",ylab="Global Mean MAE (C)", ylim = rr*c(1,1.1))
  leg = c("Raw, most recent vintage")
  for(l in 1:k)
  {
    print(l)
    lines(results_ym[vintage == vintage.vec[l],.(YM,MAE_Local)],col=col_vec[l])  
    leg = c(leg,paste0("Local, ",vintage_label(vintage.vec[l])))
  }
  axis(2)
  axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
  legend("topright", lty = 1, col = c("black",col_vec), legend = leg)
  if(print_figs)dev.off()
}


##-----------------------------------

YM_all = sort(dt_reduced[,unique(YM)])
YM_all = YM_all[-(1:12)]
dt_reduced[,residual:=SST_bar - MeanOrig]
dt_reduced[,residual_local:=SST_bar - MeanLocal]

##----- Plot a few systems -----------
ym_ex = YM_all[8]
rr = range(dt_reduced[YM == ym_ex, .(SST_bar,MeanOrig)],na.rm=TRUE)
rr_resid = range(dt_reduced[YM == ym_ex, .(residual,residual_local)],na.rm=TRUE)
plot_system(dt_reduced,ym_ex,var_plot="SST_bar",rr =rr, mn_add = " Observation", file_out = "./figures/obs_example")
plot_system(dt_reduced,ym_ex,var_plot="MeanOrig",rr =rr, mn_add = " Ensemble", file_out = "./figures/obs_ensemble")
plot_system(dt_reduced,ym_ex,var_plot="residual", rr = rr_resid, mn_add = " Residual", file_out = "./figures/obs_residual")
plot_system(dt_reduced,ym_ex,var_plot="residual_local", rr = rr_resid, mn_add = " Local Residual", file_out = "./figures/obs_local")
##------------------------------------

##----- Make a movie ----------------
rr = range(dt_reduced[,.(residual,residual_local)], na.rm=TRUE)
setkey(dt_reduced,"YM","grid_id")
for(j in 1:length(YM_all))
{
  print(j)
  if(print_figs){png(paste0("./figures/system_joined_",YM_all[j],".png"),width=1000,height=500)}else{X11()}
  par(mfrow = c(1,2))
  plot_system(dt_reduced,YM_all[j],
              var_plot = "residual", rr = rr,
              mn_add = " Ensemble", outside_control=TRUE)
  plot_system(dt_reduced,YM_all[j],
              var_plot = "residual_local", rr = rr,
              mn_add = " Local", outside_control = TRUE)
  if(print_figs)dev.off()
}
##------------------------------------


A = dt_reduced[,.(Lon = head(Lon,1),
                  Lat = head(Lat,1),
                  month = 12,
                  year= 2010,
                  YM = 1,
                  "SigmaOrig" = var(residual),
                  "SigmaLocal" = var(residual_local,na.rm=TRUE)),.(grid_id)]

rr = range(A[,SigmaLocal],na.rm=TRUE)

##----- HACK, fix -------
##Load globals
dt_obs = load_observations(2000,1)
lon_all = sort(dt_obs[,unique(Lon)])
lat_all = sort(dt_obs[,unique(Lat)])
n_lon = length(lon_all)
n_lat = length(lat_all)
##------------------
  
##------- Setup --------
dt_sub = A
mn = "Variance of Local Residuals"
##----------------------

  ##----- Bias --------------
var_plot = "SigmaLocal"
A_plot = matrix(NA, n_lon, n_lat)
A_plot[dt_sub[, grid_id]] = dt_sub[, eval(parse(text = var_plot))]
  ##-------------------------

##------- Scaling ----------
brk = seq(0,rr[2],length = 500)
brk.ind = round(seq(1,length(brk),length = 10))
brk.at = brk[brk.ind]
brk.lab = round(brk[brk.ind],2)
color <- designer.colors(n=length(brk)-1, col = c("white", "darkred"))
##--------------------------

##------- Scope ------------
if(is.null(lons))lons = range(dt_sub[,Lon])
if(is.null(lats))lats = range(dt_sub[,Lat])
##--------------------------

##------- Plot -----------
if(print_figs){ pdf("./figures/local_variance.pdf")}else{X11()}
image.plot(lon_all,lat_all,A_plot,
           main=mn,xlab="Longitude",ylab="Latitude",
           zlim=rr,
           xlim = lons,
           ylim = lats,
           breaks=brk,
           col=color,
           cex.main=1.8,cex.lab=1.4,
           cex.axis=1,
           axis.args=list(cex.axis=1,
                          at = brk.at,
                          label = brk.lab))
map("world", add = TRUE)
if(print_figs) dev.off()


