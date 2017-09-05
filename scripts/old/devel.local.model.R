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
dt[,residual := SST_bar - SST_hat_grid]
##---------------------

##--- Calculate Running Bias -------
setkey(dt,"grid_id","month","year")
dt[,meanresid := (cumsum(residual) - residual) / (year - min(year)),.(month,grid_id)]
dt[,SST_hat_local:=SST_hat_grid + meanresid]
dt[,residual_local:=SST_bar  - SST_hat_local]
##----------------------------------

##--- Now do Global Model --------
setkey(dt,"month","year")
bias_ym = dt[,.("Bias" = mean(residual, na.rm=TRUE)),keyby=.(month,year)]
bias_ym[,Trailing_Bias:=(cumsum(Bias) - Bias)/(year - min(year)),month]
bias_ym[!is.finite(Trailing_Bias), Trailing_Bias:=NA]
bias_ym[,Bias:=NULL]
dt = dt[bias_ym,on=.(month,year)]
dt[,SST_hat_global:=SST_hat_grid + Trailing_Bias]
dt[,residual_global:=SST_bar - SST_hat_global]
##----------------------------------

##-------- Get Results -------------
results_ym = dt[,.("Bias"= mean(residual,na.rm=TRUE),
                   "RMSE" = sqrt(mean(residual^2, na.rm=TRUE)),
                   "MAE" = mean(abs(residual), na.rm=TRUE),
                   "Bias_Global"= mean(residual_global,na.rm=TRUE),
                   "RMSE_Global" = sqrt(mean(residual_global^2, na.rm=TRUE)),
                   "MAE_Global" = mean(abs(residual_global), na.rm=TRUE),
                   "Bias_Local"= mean(residual_local,na.rm=TRUE),
                   "RMSE_Local" = sqrt(mean(residual_local^2, na.rm=TRUE)),
                   "MAE_Local" = mean(abs(residual_local), na.rm=TRUE),
                   "Year"= min(year),
                   "Month" = min(month)),
                keyby = YM]
results_ym = results_ym[!(Year == 1985)]
yy_all = results_ym[,range(Year)]
##-----------------------------------

##----- Plot -------------
if(print_figs){pdf("./figures/bias_ym.pdf")}else{X11()}
rr = range(results_ym[,.(Bias,Bias_Global,Bias_Local)])
plot(results_ym[,.(YM,Bias)], type="l", axes = FALSE,
     xlab="Year",ylab="Global Mean Bias (C)", ylim = rr)
lines(results_ym[,.(YM,Bias_Global)],col="red")
lines(results_ym[,.(YM,Bias_Local)],col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
legend("bottomright", lty = 1, col = c("black","red","blue"), legend = c("Raw","Global","Local"))
if(print_figs)dev.off()

if(print_figs){pdf("./figures/rmse_ym.pdf")}else{X11()}
rr = range(results_ym[,.(RMSE,RMSE_Global,RMSE_Local)])
plot(results_ym[,.(YM,RMSE)], type="l", axes = FALSE,
     xlab="Year",ylab="RMSE", ylim = c(0,rr[2]))
lines(results_ym[,.(YM,RMSE_Global)],col="red")
lines(results_ym[,.(YM,RMSE_Local)],col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
legend("bottomright", lty = 1, col = c("black","red","blue"), legend = c("Raw","Global","Local"))
if(print_figs)dev.off()

if(print_figs){pdf("./figures/mae_ym.pdf")}else{X11()}
rr = range(results_ym[,.(MAE,MAE_Global,MAE_Local)])
plot(results_ym[,.(YM,MAE)], type="l", axes = FALSE,
     xlab="Year",ylab="MAE", ylim = c(0,rr[2]))
lines(results_ym[,.(YM,MAE_Global)],col="red")
lines(results_ym[,.(YM,MAE_Local)],col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
legend("bottomright", lty = 1, col = c("black","red","blue"), legend = c("Raw","Global","Local"))
if(print_figs)dev.off()
##-----------------------------------
