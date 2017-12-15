rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)


load(file = "~/PostClimDataNoBackup/SFE/Derived/senorge2_gcfs1_combined_bias.RData")

num.years = 35

# for fitting the linear model we need to shift index in climatology/look at in-sample-climatology: 
# the regression parameters for the 2000 forecast (e.g.) are trained on temperature data up to 1999 

dt_senorge[,"oos_clim" := mean(temp,na.rm = TRUE) - temp/num.years, by = .(month,grid_id)]
dt_senorge[,"oos_clim_past" := (cumsum(temp) - temp)/(year - min(year) ),
           by = .(month,grid_id)]

dt_senorge[,"fc_clim" := mean(Ens_bar,na.rm = TRUE) - Ens_bar/num.years, by = .(month,grid_id)]
dt_senorge[,"fc_clim_past" := (cumsum(Ens_bar) - Ens_bar)/(year - min(year) ),
           by = .(month,grid_id)]
dt_senorge[year == min(year),"fc_clim_past" := 0]





#--- do linear regression 'manually' for faster updating ---

# minimizing the squared error of the regression eq. x_i = a y_i + b z_i + eps gives explicit solution in terms of product sums, so we compute those first
# the regression coefficients in year y are the one available for the forecast in year y, i.e. they are trained on the past up to y 
# here we use y_i to be the past climatology and z_i the forecast anomaly


cumsum.na.rm = function(x){
  miss = is.na(x)
  x[miss] = 0
  cs = cumsum(x)
  cs[miss] = NA
  return(cs)
}
 

dt_senorge = merge(dt_senorge, A, by = c("month","year"), all = TRUE)

dt_senorge[,"prev_clim_diff" = shift()]

regression_coeff = dt_senorge[,.(year, grid_id,
                                 "xy_sum" = cumsum.na.rm(temp * oos_clim_past) - temp * oos_clim_past,
                                 "xz_sum" = cumsum.na.rm(temp * (Ens_bar - fc_clim_past)) - temp * (Ens_bar - fc_clim_past),
                                 "yy_sum" = cumsum.na.rm(oos_clim_past^2) - oos_clim_past^2,
                                 "yz_sum" = cumsum.na.rm(oos_clim_past * (Ens_bar - fc_clim_past)) - oos_clim_past * (Ens_bar - fc_clim_past),
                                 "zz_sum" = cumsum.na.rm((Ens_bar - fc_clim_past)^2) - (Ens_bar - fc_clim_past)^2),
                              by = .(month)]

regression_coeff_2 = regression_coeff_2[!is.na(xy_sum)][!is.na(zz_sum)][grid_id == max(grid_id)]
setkey(regression_coeff_2, "year","month")

regression_coeff_2[year > min(year) + 1, "a" := (xy_sum * zz_sum - xz_sum * yz_sum)/(yy_sum * zz_sum - yz_sum^2) ]
regression_coeff_2[year > min(year) + 1, "b" := (xz_sum - a * yz_sum)/ zz_sum]

regression_coeff_2[,c("grid_id","xy_sum","xz_sum","yy_sum","yz_sum","zz_sum") := NULL]


dt_senorge = merge(dt_senorge, regression_coeff_2, by = c("month","year"), all = TRUE)



dt_senorge[,"temp_hat_lr_2" := a * temp_hat + b * bias]
dt_senorge[,c("a","b") := NULL]



#----- do bias correction based on bias of last month --- 

# use linear regression model x_i = a y_i + b z_i + eps where y_i is the forecast corrected by the average bias of the month, and z_i is fc - obs of the previous month

k=1

dt_senorge[,"bias":= shift(Ens_bar - temp,k), by = grid_id]


regression_coeff_3 = dt_senorge[,.(year, 
                                 "xy_sum" = cumsum.na.rm(temp * temp_hat) - temp * temp_hat,
                                 "xz_sum" = cumsum.na.rm(temp * bias) - temp * bias,
                                 "yy_sum" = cumsum.na.rm(temp_hat^2) - temp_hat^2,
                                 "yz_sum" = cumsum.na.rm(temp_hat * bias) - temp_hat * bias,
                                 "zz_sum" = cumsum.na.rm(bias^2) - bias^2),
                              by = .(month,grid_id)]

regression_coeff_3[year > min(year) + 1, "a" := (xy_sum * zz_sum - xz_sum * yz_sum)/(yy_sum * zz_sum - yz_sum^2) ]
regression_coeff_3[year > min(year) + 1, "b" := (xz_sum - a * yz_sum)/ zz_sum  ]

regression_coeff_3 = regression_coeff_3[order(year,month,grid_id)]

dt_senorge[,"temp_hat_lr_3" := regression_coeff_3[,a] * temp_hat + regression_coeff_3[,b] * bias]

#------ not by grid id -----

regression_coeff_4 = dt_senorge[,.(year, grid_id, 
                                   "xy_sum" = cumsum.na.rm(temp * temp_hat) - temp * temp_hat,
                                   "xz_sum" = cumsum.na.rm(temp * bias) - temp * bias,
                                   "yy_sum" = cumsum.na.rm(temp_hat^2) - temp_hat^2,
                                   "yz_sum" = cumsum.na.rm(temp_hat * bias) - temp_hat * bias,
                                   "zz_sum" = cumsum.na.rm(bias^2) - bias^2),
                                by = .(month)]

regression_coeff_4 = regression_coeff_4[!is.na(xy_sum)][!is.na(zz_sum)][grid_id == max(grid_id)]
setkey(regression_coeff_4, "year","month")

regression_coeff_4[year > min(year) + 1, "a" := (xy_sum * zz_sum - xz_sum * yz_sum)/(yy_sum * zz_sum - yz_sum^2) ]
regression_coeff_4[year > min(year) + 1, "b" := (xz_sum - a * yz_sum)/ zz_sum]

regression_coeff_4[,c("grid_id","xy_sum","xz_sum","yy_sum","yz_sum","zz_sum") := NULL]


dt_senorge = merge(dt_senorge, regression_coeff_4, by = c("month","year"), all = TRUE)



dt_senorge[,"temp_hat_lr_4" := a * temp_hat + b * bias]
dt_senorge[,c("a","b") := NULL]

#---- plot parameters ---

pdf("./figures/senorge/Regression_Coeff.pdf")

rr = range(na.omit(regression_coeff_4[,.(a,b)]))

means = regression_coeff_4[,"a_mean":= mean(a,na.rm = TRUE),by = month]
means[,"b_mean":= mean(b,na.rm = TRUE), keyby = month]

plot( means[year == 2000,month], 
      means[year == 2000,a_mean],
      main = "Regression parameters by month",
      ylab = "value",
      xlab = "month",
      
      type = 'b',
      ylim = rr, 
      col = "blue") 

lines(means[year == 2000,month], 
      means[year == 2000,b_mean], type = 'b', ylog = TRUE, ylim = rr, col = "darkred")

legend("topright",legend = c("a", "b"), col = c("blue","darkred"),lty = 1)

dev.off()



#---- plot mean scores by month ---

RMSE_scores_by_month = dt_senorge[year>1999,.( "RMSE_Ens_bar" = sqrt(mean((Ens_bar-temp)^2,na.rm = TRUE)),
                                               "RMSE_bc" = sqrt(mean((temp_hat-temp)^2,na.rm = TRUE)),
                                               "RMSE_clim" = sqrt(mean((temp-oos_clim_past)^2,na.rm = TRUE)),
                                               "RMSE_lr" = sqrt(mean((temp-temp_hat_lr)^2,na.rm = TRUE)),
                                               #"RMSE_lr_2" = sqrt(mean((temp-temp_hat_lr_2)^2,na.rm = TRUE)),
                                               "RMSE_lr_3" = sqrt(mean((temp-temp_hat_lr_3)^2,na.rm = TRUE)),
                                               "RMSE_lr_4" = sqrt(mean((temp-temp_hat_lr_4)^2,na.rm = TRUE))),
                                  by = month]



rr = range(na.omit(RMSE_scores_by_month[,.(RMSE_Ens_bar,RMSE_bc,RMSE_clim)]))


pdf("./figures/senorge/RMSE_by_month_2_2000.pdf")

plot( RMSE_scores_by_month[,month], 
      RMSE_scores_by_month[,RMSE_Ens_bar],
      main = "RMSE spatially averaged by month, 2000 - 2015",
      ylab = "mean RMSE",
      xlab = "month",
      
      type = 'b',
      ylim = rr, 
      col = "blue") 

lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_bc], type = 'b', ylog = TRUE, ylim = rr, col = "darkred")
lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_clim], type = 'b', ylog = TRUE, ylim = rr, col = "darkgreen")
lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_lr_3], type = 'b', ylog = TRUE, ylim = rr, col = "black")
lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_lr_4], type = 'b', ylog = TRUE, ylim = rr, col = "yellow2")
#lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_lr_3], type = 'b', ylog = TRUE, ylim = rr, col = "orange2")

legend("topright",legend = c("raw ensemble", "bias corrected", "climatology", "lr by gridpoint", "lr globally"), col = c("blue","darkred","darkgreen","black","yellow2"),lty = 1)

dev.off()


#--- January 2005+ 

RMSE_scores_jan = dt_senorge[year>2004 & month ==1,.( "RMSE_Ens_bar" = sqrt(mean((Ens_bar-temp)^2,na.rm = TRUE)),
                                              "RMSE_bc" = sqrt(mean((temp_hat-temp)^2,na.rm = TRUE)),
                                              "RMSE_clim" = sqrt(mean((temp-oos_clim_past)^2,na.rm = TRUE)),
                                              "RMSE_lr" = sqrt(mean((temp-temp_hat_lr)^2,na.rm = TRUE)),
                                              #"RMSE_lr_2" = sqrt(mean((temp-temp_hat_lr_2)^2,na.rm = TRUE)),
                                              "RMSE_lr_3" = sqrt(mean((temp-temp_hat_lr_3)^2,na.rm = TRUE)),
                                              "RMSE_lr_4" = sqrt(mean((temp-temp_hat_lr_4)^2,na.rm = TRUE))),
                                 by = year]



rr = range(na.omit(RMSE_scores_jan[,.(RMSE_Ens_bar,RMSE_bc,RMSE_clim)]))


pdf("./figures/senorge/RMSE_by_year_jan2005.pdf")

plot( RMSE_scores_jan[,year], 
      RMSE_scores_jan[,RMSE_Ens_bar],
      main = "RMSE spatially averaged by year, January 2005 +",
      ylab = "mean RMSE",
      xlab = "year",
      
      type = 'b',
      ylim = rr, 
      col = "blue") 

lines(RMSE_scores_jan[,year], RMSE_scores_jan[,RMSE_bc], type = 'b', ylog = TRUE, ylim = rr, col = "darkred")
lines(RMSE_scores_jan[,year], RMSE_scores_jan[,RMSE_clim], type = 'b', ylog = TRUE, ylim = rr, col = "darkgreen")
lines(RMSE_scores_jan[,year], RMSE_scores_jan[,RMSE_lr_3], type = 'b', ylog = TRUE, ylim = rr, col = "black")
lines(RMSE_scores_jan[,year], RMSE_scores_jan[,RMSE_lr_4], type = 'b', ylog = TRUE, ylim = rr, col = "yellow2")

legend("topright",legend = c("raw ensemble", "bias corrected", "climatology", "lr. by gridpoint", "lr. globally"), col = c("blue","darkred","darkgreen","black", "yellow2"),lty = 1)

dev.off()



#---- plot mean scores by year ---

RMSE_scores_by_year = dt_senorge[year>1985,.( "RMSE_Ens_bar" = sqrt(mean((Ens_bar-temp)^2,na.rm = TRUE)),
                                              "RMSE_bc" = sqrt(mean((temp_hat-temp)^2,na.rm = TRUE)),
                                              "RMSE_clim" = sqrt(mean((temp-oos_clim_past)^2,na.rm = TRUE)),
                                              "RMSE_lr" = sqrt(mean((temp-temp_hat_lr)^2,na.rm = TRUE)),
                                              #"RMSE_lr_2" = sqrt(mean((temp-temp_hat_lr_2)^2,na.rm = TRUE)),
                                              "RMSE_lr_3" = sqrt(mean((temp-temp_hat_lr_3)^2,na.rm = TRUE)),
                                              "RMSE_lr_4" = sqrt(mean((temp-temp_hat_lr_4)^2,na.rm = TRUE))),
                                 by = year]



rr = range(na.omit(RMSE_scores_by_year[,.(RMSE_Ens_bar,RMSE_bc,RMSE_clim)]))


pdf("./figures/senorge/RMSE_by_year.pdf")

plot( RMSE_scores_by_year[,year], 
      RMSE_scores_by_year[,RMSE_Ens_bar],
      main = "RMSE spatially averaged by year",
      ylab = "mean RMSE",
      xlab = "year",
      
      type = 'b',
      ylim = rr, 
      col = "blue") 

lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_bc], type = 'b', ylog = TRUE, ylim = rr, col = "darkred")
lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_clim], type = 'b', ylog = TRUE, ylim = rr, col = "darkgreen")
lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_lr], type = 'b', ylog = TRUE, ylim = rr, col = "black")
lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_lr_2], type = 'b', ylog = TRUE, ylim = rr, col = "yellow2")

legend("topright",legend = c("raw ensemble", "bias corrected", "climatology", "lr. by gridpoint", "lr. globally"), col = c("blue","darkred","darkgreen","black", "yellow2"),lty = 1)

dev.off()

#-----------------------------------------------------------------------------------------------------------------------

###############################


#-------------- Do the same for SST --------------------

load(file = "~/PostClimDataNoBackup/SFE/Derived/dt_combine_wide_bias.RData")

num.years = 26

dt[,"oos_clim" := (num.years * mean(SST_bar,na.rm = TRUE) - SST_bar)/(num.years-1), by = .(month,grid_id)]
dt[,"oos_clim_past" := (cumsum(SST_bar) - SST_bar)/(year - min(year) ),
           by = .(month,grid_id)]

dt[,"fc_clim" := mean(Ens_bar,na.rm = TRUE) - Ens_bar/num.years, by = .(month,grid_id)]
dt[,"fc_clim_past" := (cumsum(Ens_bar) - Ens_bar)/(year - min(year) ),
           by = .(month,grid_id)]
dt[year == min(year),"fc_clim_past" := 0]





#--- do linear regression 'manually' for faster updating ---

# minimizing the squared error of the regression eq. x_i = a y_i + b z_i + eps gives explicit solution in terms of product sums, so we compute those first
# the regression coefficients in year y are the one available for the forecast in year y, i.e. they are trained on the past up to y 
# here we use y_i to be the past climatology and z_i the forecast anomaly


cumsum.na.rm = function(x){
  miss = is.na(x)
  x[miss] = 0
  cs = cumsum(x)
  cs[miss] = NA
  return(cs)
}

regression_coeff = dt[,.(year, 
                                 "xy_sum" = cumsum.na.rm(SST_bar * oos_clim_past) - SST_bar * oos_clim_past,
                                 "xz_sum" = cumsum.na.rm(SST_bar * (Ens_bar - fc_clim_past)) - SST_bar * (Ens_bar - fc_clim_past),
                                 "yy_sum" = cumsum.na.rm(oos_clim_past^2) - oos_clim_past^2,
                                 "yz_sum" = cumsum.na.rm(oos_clim_past * (Ens_bar - fc_clim_past)) - oos_clim_past * (Ens_bar - fc_clim_past),
                                 "zz_sum" = cumsum.na.rm((Ens_bar - fc_clim_past)^2) - (Ens_bar - fc_clim_past)^2),
                              by = .(month,grid_id)]

regression_coeff[year > min(year) + 1, "a" := (xy_sum * zz_sum - xz_sum * yz_sum)/(yy_sum * zz_sum - yz_sum^2) ]
regression_coeff[year > min(year) + 1, "b" := (xz_sum - a * yz_sum)/ zz_sum  ]

regression_coeff = regression_coeff[order(year,month,grid_id)]

dt[,"SST_hat_lr" := regression_coeff[,a] * oos_clim_past + regression_coeff[,b] * (Ens_bar - fc_clim_past)]


#----- not by grid_id ------

regression_coeff_2 = dt[,.(year, grid_id,
                                   "xy_sum" = cumsum.na.rm(SST_bar * oos_clim_past) - SST_bar * oos_clim_past,
                                   "xz_sum" = cumsum.na.rm(SST_bar * (Ens_bar - fc_clim_past)) - SST_bar * (Ens_bar - fc_clim_past),
                                   "yy_sum" = cumsum.na.rm(oos_clim_past^2) - oos_clim_past^2,
                                   "yz_sum" = cumsum.na.rm(oos_clim_past * (Ens_bar - fc_clim_past)) - oos_clim_past * (Ens_bar - fc_clim_past),
                                   "zz_sum" = cumsum.na.rm((Ens_bar - fc_clim_past)^2) - (Ens_bar - fc_clim_past)^2),
                                by = .(month)]

regression_coeff_2 = regression_coeff_2[!is.na(xy_sum)][!is.na(zz_sum)][grid_id == max(grid_id)]
setkey(regression_coeff_2, "year","month")

regression_coeff_2[year > min(year) + 1, "a" := (xy_sum * zz_sum - xz_sum * yz_sum)/(yy_sum * zz_sum - yz_sum^2) ]
regression_coeff_2[year > min(year) + 1, "b" := (xz_sum - a * yz_sum)/ zz_sum]

regression_coeff_2[,c("grid_id","xy_sum","xz_sum","yy_sum","yz_sum","zz_sum") := NULL]


dt = merge(dt, regression_coeff_2, by = c("month","year"), all = TRUE)



dt[,"SST_hat_lr_2" := a * oos_clim_past + b * (Ens_bar - fc_clim_past)]
dt[,c("a","b") := NULL]



#----- do bias correction based on bias of last month --- 

# use linear regression model x_i = a y_i + b z_i + eps where y_i is the forecast corrected by the average bias of the month, and z_i is fc - obs of the previous month

k=1

dt[,"bias":= shift(Ens_bar - SST_bar,k), by = grid_id]


regression_coeff_3 = dt[,.(year, 
                                   "xy_sum" = cumsum.na.rm(SST_bar * SST_hat) - SST_bar * SST_hat,
                                   "xz_sum" = cumsum.na.rm(SST_bar * bias) - SST_bar * bias,
                                   "yy_sum" = cumsum.na.rm(SST_hat^2) - SST_hat^2,
                                   "yz_sum" = cumsum.na.rm(SST_hat * bias) - SST_hat * bias,
                                   "zz_sum" = cumsum.na.rm(bias^2) - bias^2),
                                by = .(month,grid_id)]

regression_coeff_3[year > min(year) + 1, "a" := (xy_sum * zz_sum - xz_sum * yz_sum)/(yy_sum * zz_sum - yz_sum^2) ]
regression_coeff_3[year > min(year) + 1, "b" := (xz_sum - a * yz_sum)/ zz_sum  ]

regression_coeff_3 = regression_coeff_3[order(year,month,grid_id)]

dt[,"SST_hat_lr_3" := regression_coeff_3[,a] * SST_hat + regression_coeff_3[,b] * bias]

#------ not by grid id -----

regression_coeff_4 = dt[,.(year, grid_id, 
                                   "xy_sum" = cumsum.na.rm(SST_bar * SST_hat) - SST_bar * SST_hat,
                                   "xz_sum" = cumsum.na.rm(SST_bar * bias) - SST_bar * bias,
                                   "yy_sum" = cumsum.na.rm(SST_hat^2) - SST_hat^2,
                                   "yz_sum" = cumsum.na.rm(SST_hat * bias) - SST_hat * bias,
                                   "zz_sum" = cumsum.na.rm(bias^2) - bias^2),
                                by = .(month)]

regression_coeff_4 = regression_coeff_4[!is.na(xy_sum)][!is.na(zz_sum)][grid_id == max(grid_id)]
setkey(regression_coeff_4, "year","month")

regression_coeff_4[year > min(year) + 1, "a" := (xy_sum * zz_sum - xz_sum * yz_sum)/(yy_sum * zz_sum - yz_sum^2) ]
regression_coeff_4[year > min(year) + 1, "b" := (xz_sum - a * yz_sum)/ zz_sum]

regression_coeff_4[,c("grid_id","xy_sum","xz_sum","yy_sum","yz_sum","zz_sum") := NULL]


dt = merge(dt, regression_coeff_4, by = c("month","year"), all = TRUE)



dt[,"SST_hat_lr_4" := a * SST_hat + b * bias]
dt[,c("a","b") := NULL]



#---- plot mean scores by month ---

RMSE_scores_by_month = dt[year>1980,.( "RMSE_Ens_bar" = sqrt(mean((Ens_bar-SST_bar)^2,na.rm = TRUE)),
                                               "RMSE_bc" = sqrt(mean((SST_hat-SST_bar)^2,na.rm = TRUE)),
                                               "RMSE_clim" = sqrt(mean((SST_bar-oos_clim)^2,na.rm = TRUE)),
                                                "RMSE_clim_past" = sqrt(mean((SST_bar-oos_clim_past)^2,na.rm = TRUE)),
                                                "RMSE_lr" = sqrt(mean((SST_bar-SST_hat_lr)^2,na.rm = TRUE)),
                                                "RMSE_lr_2" = sqrt(mean((SST_bar-SST_hat_lr_2)^2,na.rm = TRUE)),
                                               "RMSE_lr_3" = sqrt(mean((SST_bar-SST_hat_lr_3)^2,na.rm = TRUE)),
                                               "RMSE_lr_4" = sqrt(mean((SST_bar-SST_hat_lr_4)^2,na.rm = TRUE))),
                                  by = month]



rr = range(na.omit(RMSE_scores_by_month[,.(RMSE_Ens_bar,RMSE_bc,RMSE_clim)]))


pdf("./figures/RMSE_by_month_2.pdf")

plot( RMSE_scores_by_month[,month], 
      RMSE_scores_by_month[,RMSE_Ens_bar],
      main = "RMSE spatially averaged by month",
      ylab = "mean RMSE",
      xlab = "month",
      
      type = 'b',
      ylim = rr, 
      col = "blue") 

lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_clim], type = 'b', ylog = TRUE, ylim = rr, col = "darkgreen")
lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_bc], type = 'b', ylog = TRUE, ylim = rr, col = "darkred")
#lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_lr], type = 'b', ylog = TRUE, ylim = rr, col = "black")
lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_clim_past], type = 'b', ylog = TRUE, ylim = rr, col = "yellow2")
#lines(RMSE_scores_by_month[,month], RMSE_scores_by_month[,RMSE_lr_3], type = 'b', ylog = TRUE, ylim = rr, col = "orange2")

legend("topleft",legend = c("raw ensemble", "bias corrected ens", "climatology", "past climatology"), col = c("blue","darkred","darkgreen","yellow2"),lty = 1)

dev.off()


#---- plot mean scores by year ---

RMSE_scores_by_year = dt[year>1985,.( "RMSE_Ens_bar" = sqrt(mean((Ens_bar-SST_bar)^2,na.rm = TRUE)),
                                              "RMSE_bc" = sqrt(mean((SST_hat-SST_bar)^2,na.rm = TRUE)),
                                              "RMSE_clim" = sqrt(mean((SST_bar-oos_clim_past)^2,na.rm = TRUE)),
                                              "RMSE_lr" = sqrt(mean((SST_bar-SST_hat_lr)^2,na.rm = TRUE)),
                                              "RMSE_lr_2" = sqrt(mean((SST_bar-SST_hat_lr_2)^2,na.rm = TRUE)),
                                              "RMSE_lr_3" = sqrt(mean((SST_bar-SST_hat_lr_3)^2,na.rm = TRUE)),
                                              "RMSE_lr_4" = sqrt(mean((SST_bar-SST_hat_lr_4)^2,na.rm = TRUE))),
                                 by = year]



rr = range(na.omit(RMSE_scores_by_year[,.(RMSE_Ens_bar,RMSE_bc,RMSE_clim)]))


pdf("./figures/senorge/RMSE_by_year.pdf")

plot( RMSE_scores_by_year[,year], 
      RMSE_scores_by_year[,RMSE_Ens_bar],
      main = "RMSE spatially averaged by year",
      ylab = "mean RMSE",
      xlab = "year",
      
      type = 'b',
      ylim = rr, 
      col = "blue") 

lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_bc], type = 'b', ylog = TRUE, ylim = rr, col = "darkred")
lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_clim], type = 'b', ylog = TRUE, ylim = rr, col = "darkgreen")
lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_lr], type = 'b', ylog = TRUE, ylim = rr, col = "black")
lines(RMSE_scores_by_year[,year], RMSE_scores_by_year[,RMSE_lr_2], type = 'b', ylog = TRUE, ylim = rr, col = "yellow2")

legend("topright",legend = c("raw ensemble", "bias corrected", "climatology", "lr. by gridpoint", "lr. globally"), col = c("blue","darkred","darkgreen","black", "yellow2"),lty = 1)

dev.off()





