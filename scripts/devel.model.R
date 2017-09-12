rm(list = ls())

##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)
print_figs = TRUE
##------------------------

##------ Set up -------
dt = load_combined()
dt[,residual := SST_bar - Forecast]
##---------------------

##--- Calculate Running Bias -------
setkey(dt,"grid_id","Ens","month","year")
dt[,"MeanResid" := mean(residual),.(grid_id, month,year)]
setkey(dt,"grid_id","Ens","month","year")
dt[,"MeanResidCum" := (cumsum(MeanResid) - MeanResid) / (year - min(year)),.(grid_id, Ens, month)]
setkey(dt,"grid_id","month","year", "Ens")
dt[,"SST_hat_local":=Forecast + MeanResidCum]
dt[,"Obs_above_orig":=SST_bar > Forecast]
dt[,"Obs_above_local":=SST_bar > SST_hat_local]
dt[,"Rank_orig":=cumsum(Obs_above_orig),.(grid_id,year,month)]
dt[,"Rank_local":=cumsum(Obs_above_local),.(grid_id,year,month)]
dt[,"residual_local":=SST_bar - SST_hat_local]
##----------------------------------

dt_rank = dt[,.("Q_Orig"=max(Rank_orig) + 1,"Q_local"=max(Rank_local) + 1),.(grid_id,year,month)]

prop_orig = dt_rank[,.N,by="Q_Orig"]
prop_orig = prop_orig[!is.na(Q_Orig)]
prop_orig[,"Prop":=N/sum(N)]

prop_local = dt_rank[,.N,by="Q_local"]
prop_local = prop_local[!is.na(Q_local)]
prop_local[,"Prop":=N/sum(N)]

##---- Rank Histograms -------------------
if(print_figs){pdf("./figures/calibration.pdf")}else{X11()}
plot(prop_orig[,.(Q_Orig,Prop)],type="h",lwd = 3, xlab = "Rank of Observation in Ensemble", ylab = "Proportion")
lines(prop_local[,.(Q_local + .1, Prop)], type="h", lwd = 3, col="blue")
legend("topleft", lty =1, col=c("black","blue"),legend = c("Raw Ensemble","Local Seasonally Adjusted")) 
if(print_figs)dev.off()
##----------------------------------------

##----- Reduce Observation Space ---------
setkey(dt,"grid_id","YM","Ens")
dt_reduced = dt[,.("Lon" = head(Lon,1),
                   "Lat" = head(Lat,1),
                   "SST_bar" = head(SST_bar,1),
                   "MeanOrig" = mean(Forecast),
                   "MedianOrig" = median(Forecast),
                   "MeanLocal" = mean(SST_hat_local),
                   "MedianLocal" = median(SST_hat_local),
                   "year" = head(year,1),
                   "month" = head(month,1)),
                .(grid_id, YM)]
##------------------------------------------

##-------- Get Results -------------
results_ym = dt_reduced[,.("Bias"= mean(SST_bar - MeanOrig,na.rm=TRUE),
                   "RMSE" = sqrt(mean( (SST_bar - MeanOrig)^2, na.rm=TRUE)),
                   "MAE" = mean(abs(SST_bar - MedianOrig), na.rm=TRUE),
                   "Bias_Local"= mean(SST_bar - MeanLocal,na.rm=TRUE),
                   "RMSE_Local" = sqrt(mean( (SST_bar - MeanLocal)^2, na.rm=TRUE)),
                   "MAE_Local" = mean(abs(SST_bar - MedianLocal), na.rm=TRUE),
                   "Year"= min(year),
                   "Month" = min(month)),
                keyby = YM]
results_ym = results_ym[!(Year == 1985)]
yy_all = results_ym[,range(Year)]
##-----------------------------------

##----- Plot -------------
if(print_figs){pdf("./figures/bias_ym.pdf")}else{X11()}
rr = range(results_ym[,.(Bias,Bias_Local)])
plot(results_ym[,.(YM,Bias)], type="l", axes = FALSE,
     xlab="Year",ylab="Global Mean Bias (C)", ylim = rr)
lines(results_ym[,.(YM,Bias_Local)],col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
legend("bottomright", lty = 1, col = c("black","blue"), legend = c("Raw","Local"))
if(print_figs)dev.off()

if(print_figs){pdf("./figures/rmse_ym.pdf")}else{X11()}
rr = range(results_ym[,.(RMSE,RMSE_Local)])
plot(results_ym[,.(YM,RMSE)], type="l", axes = FALSE,
     xlab="Year",ylab="RMSE", ylim = c(0,rr[2]))
lines(results_ym[,.(YM,RMSE_Local)],col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
legend("bottomright", lty = 1, col = c("black","blue"), legend = c("Raw","Local"))
if(print_figs)dev.off()

if(print_figs){pdf("./figures/mae_ym.pdf")}else{X11()}
rr = range(results_ym[,.(MAE,MAE_Local)])
plot(results_ym[,.(YM,MAE)], type="l", axes = FALSE,
     xlab="Year",ylab="MAE", ylim = c(0,rr[2]))
lines(results_ym[,.(YM,MAE_Local)],col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
legend("bottomright", lty = 1, col = c("black","blue"), legend = c("Raw","Local"))
if(print_figs)dev.off()
##-----------------------------------

YM_all = sort(dt_reduced[,unique(YM)])
YM_all = YM_all[-(1:12)]
dt_reduced[,residual:=SST_bar - MeanOrig]
dt_reduced[,residual_local:=SST_bar - MeanLocal]
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



