## This produces the plots for the PIT marginal stuff

rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)

#------- Check marginal calibration of PCA --------------------

month = 1:12
y = 1985:2010



dist_fun_tn = function(value, mean, sd, truncate = -1.619995){ # distribution function for a truncated normal distribution evaluated at value; value, mean and sd need to be vectors of equal length
  a=rep(0,times = length(value))
  for(i in  1:length(value))
  if(is.na(value[i])) { a[i] = NA
  }else if(value[i] > truncate) {a[i] = pnorm(value[i], mean = mean[i], sd = sd[i])
  } else a[i] =  runif(1, max = pnorm(truncate, mean = mean[i], sd = sd[i]))
  return(a)
}

setup_PCA(m = month, y = y, oriented = TRUE)

vec = c(30,50)
for(k in vec){
  print(paste0("k = ",k))
  
  calib = forecast_PCA(y = y, m = month, PCA_depth = k, output_opts = "mar_sd" , saveorgo = FALSE, truncate = FALSE)
  print("forecast_PCA finished")
  
  calib[,"PIT" := dist_fun_tn(SST_bar, mean = SST_hat, sd = forecast)]
  print("PIT complete")
  
  for(moment in 1:2){
    
    if(moment == 1) calib[,"moment":= mean(PIT,na.rm = TRUE),by = grid_id]
    if(moment == 2) calib[,"moment":= var(PIT,na.rm = TRUE),by = grid_id]
    
    
    save.dir="./Data/PostClim/SFE/Derived/PCA"
    save(calib, file = paste0(save.dir,"/cal",k,"_mom",moment,".RData"))
    print("moment complete")
  }
  
  plot_system(type = "cal", moment = 1, depth = k, plot_title = paste0("PIT mean, ",k, " prin.com."))
  plot_system(type = "cal", moment = 2, depth = k, plot_title = paste0("PIT SD, ",k, " prin.com."))
  print("plots complete")
  
}


#-------- Calibration of Ensemble --------


dt = load_combined_wide(bias = TRUE)


calib = dt[,.(month, year, Lon, Lat,grid_id, SST_bar,SST_hat,Ens_bar,Ens_sd)]
calib[,"PIT" := pnorm(SST_bar, mean = SST_hat, sd = Ens_sd), key = c("month","year","Lat")]
  
for(moment in 1:2){
    
  if(moment == 1) calib[,"moment":= mean(PIT,na.rm = TRUE), by = grid_id]
  if(moment == 2) calib[,"moment":= var(PIT,na.rm = TRUE), by = grid_id]
    
  save.dir="./Data/PostClim/SFE/Derived/PCA"
  save(calib, file = paste0(save.dir,"/cal_ens_mom",moment,".RData"))
  }
  
  plot_system(type = "cal", moment = 1, depth = 0, plot_title = paste0("PIT mean, bias corrected ensemble forecast"))
  plot_system(type = "cal", moment = 2, depth = 0, plot_title = paste0("PIT SD, raw ensemble forecast"))
  
  
  #-------- Calibration of (past) climatology --------
  
  
  dt = load_combined_wide(bias = TRUE)
  
  dt[,"oos_clim" := mean(SST_bar,na.rm = TRUE) - SST_bar/num.years, by = .(month,grid_id)]
  dt[,"oos_clim_past" := (cumsum(SST_bar) - SST_bar)/(year - min(year) ),
     by = .(month,grid_id)]
  calib = dt[,.(month, year, Lon, Lat,grid_id, SST_bar,oos_clim_past)]
  
  past_sd = function (x){
    a = c() 
    for(k in 1:length(x)) a = c(a, sd(x[1:k], na.rm = TRUE))
    return(a)
  }
  calib[,"clim_sd" := past_sd(oos_clim_past), by = .(grid_id, month)]
  calib[,"PIT" := pnorm(SST_bar, mean = oos_clim_past, sd = clim_sd), key = c("month","year","Lat")]
  
  for(moment in 1:2){
    
    if(moment == 1) calib[year>1985,"moment":= mean(PIT,na.rm = TRUE), by = grid_id]
    if(moment == 2) calib[year>1985,"moment":= var(PIT,na.rm = TRUE), by = grid_id]
    
    save.dir="./Data/PostClim/SFE/Derived/PCA"
    save(calib, file = paste0(save.dir,"/cal_clim_mom",moment,".RData"))
  }
  
  plot_system(type = "cal", moment = 1, depth = -1, plot_title = paste0("PIT mean, past climatology"))
  plot_system(type = "cal", moment = 2, depth = -1, plot_title = paste0("PIT variance, past climatology"))
  
  
  


