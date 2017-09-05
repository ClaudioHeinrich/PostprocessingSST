rm(list = ls())

##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)
print_figs = TRUE
##------------------------

##------ Set up -------
dt_combined_all = load_combined()
dt = load_combined()
##---------------------

##----- Calculate Errors ---
dt[,bias := SST_bar - SST_hat_grid]
##--------------------------

##------ Let's look at overall bias --
ym_index = dt[, .("Year" = min(year),"Month" = min(month)),keyby=YM]
bias_ym = dt[,.("Bias"= mean(bias,na.rm=TRUE),
                "RMSE" = sqrt(mean(bias^2, na.rm=TRUE)),
                "MAE" = mean(abs(bias), na.rm=TRUE),
                "Year"= min(year),
                "Month" = min(month)),
             keyby = YM]
yy_all = dt[,range(year)]
##------------------------------------

##--------------------------------------

##------  Simple Model ------------
year_all = bias_ym[,unique(Year)]
bias_ym[,bias_coef:=0.0]
mod_save = list()
for(j in 2:length(year_all))
{
  mod = lm("Bias ~ as.factor(Month) - 1", data = bias_ym[Year < year_all[2]])
  f_mod = approxfun(1:12,mod$coeff,method = "constant",rule=2)
  mod_save[[j]] = mod
  bias_ym[Year == year_all[j],bias_coef:=f_mod(Month)]
}
##---------------------------------

##------ Propogate simple model forwards -----
f_bias_ym = approxfun(bias_ym$YM,bias_ym$bias_coef, method="constant", rule=2)
dt[,bias_coef:= f_bias_ym(YM)]
dt[,SST_hat_grid_global_bias:=SST_hat_grid + bias_coef]
dt[,bias_global:= SST_bar - SST_hat_grid_global_bias]
##---------------------------------------------

##------- New results -------------------------
bias_ym = dt[,.("Bias"= mean(bias,na.rm=TRUE),
                "RMSE" = sqrt(mean(bias^2, na.rm=TRUE)),
                "MAE" = mean(abs(bias), na.rm=TRUE),
                "Bias_Global"= mean(bias_global,na.rm=TRUE),
                "RMSE_Global" = sqrt(mean(bias_global^2, na.rm=TRUE)),
                "MAE_Global" = mean(abs(bias_global), na.rm=TRUE),
                "Year"= min(year),
                "Month" = min(month)),
             keyby = YM]
##----------------------------------------------


yy_all = range(dt[,year])

if(print_figs){pdf("./figures/bias_ym.pdf")}else{X11()}
plot(bias_ym[,.(YM,Bias)], type="l", axes = FALSE, xlab="Year",ylab="Global Mean Bias (C)", ylim = c(0, bias_ym[,max(Bias)]))
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
if(print_figs)dev.off()

if(print_figs){pdf("./figures/rmse_ym.pdf")}else{X11()}
plot(bias_ym[,.(YM,RMSE)], type="l", axes = FALSE, xlab="Year",ylab="Global RMSE")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
if(print_figs)dev.off()

if(print_figs){pdf("./figures/mae_ym.pdf")}else{X11()}
plot(bias_ym[,.(YM,MAE)],type="l", axes = FALSE,xlab="Year",ylab="Global MAE")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
if(print_figs)dev.off()

##---------- Plot --------------------------
if(print_figs){pdf("./figures/bias_global_ym.pdf")}else{X11()}
plot(bias_ym[,.(YM,Bias)], type="l", axes = FALSE, xlab="Year",ylab="Global Mean Bias (C)",
     ylim = c(min(0,min(bias_ym[,.(Bias, Bias_Global)])), max(0,max(bias_ym[,.(Bias, Bias_Global)]))))
lines(bias_ym[,.(YM,Bias_Global)],col = "blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
if(print_figs)dev.off()

if(print_figs){pdf("./figures/rmse_global_ym.pdf")}else{X11()}
plot(bias_ym[,.(YM,RMSE)], type="l", axes = FALSE, xlab="Year",ylab="Global RMSE",
     ylim = range(bias_ym[,.(RMSE,RMSE_Global)]))
lines(bias_ym[,.(YM,RMSE_Global)], col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
if(print_figs)dev.off()

if(print_figs){pdf("./figures/mae_global_ym.pdf")}else{X11()}
plot(bias_ym[,.(YM,MAE)],type="l", axes = FALSE,xlab="Year",ylab="Global MAE",
     ylim = range(bias_ym[,.(MAE,MAE_Global)]))
lines(bias_ym[,.(YM,MAE_Global)], col="blue")
axis(2)
axis(1, at = 0:(diff(yy_all)) * 12 + dt[,min(YM)], label = yy_all[1]:yy_all[2])
if(print_figs)dev.off()
##--------------------------------------

YM_j = 23949
dt_sub = dt[YM == YM_j]
n_lon = length(dt_sub[, unique(Lon)])
n_lat = length(dt_sub[, unique(Lat)])
mn = paste0(dt_sub[, min(month)], "/", dt_sub[, min(year)])

A_bias_global = matrix(NA, n_lon, n_lat)
A_bias_global[dt_sub[grid_id < length(A_bias_global), grid_id]] =
  dt_sub[grid_id < length(A_bias_global), bias_global]
if(print_figs){pdf(paste0("./figures/system_",YM_j,"_global_bias.pdf"))}else{X11()}
image(A_bias_global, main = paste0(mn, " global bias"))
if (print_figs) dev.off()
