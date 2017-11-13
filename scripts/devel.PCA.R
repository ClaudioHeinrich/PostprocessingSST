rm(list = ls())

##------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)
##------------------------


forecast_PCA = function(dt_reduced = NULL, 
                        y = 1999, 
                        m = 1, 
                        PCA_depth = 5, 
                        save.dir="./Data/PostClim/SFE/Derived/PCA",
                        cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                        data.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                        saveorgo = TRUE,
                        truncate = TRUE){
  
  if(is.null(dt_reduced))
  {
    load(file = paste0(data.dir,"dt_reduced_PCA.RData"))
  }
  
  #find land grid ids:
  
  land_grid_id = dt_reduced[variable == "SST1" & 
                            year == min(y) & 
                            month == min(m) & 
                            is.na(SST_hat_local),
                            .(Lon,Lat,grid_id,month,year)]
  
  fc <- na.omit( dt_reduced[variable == "SST1" 
                            & year %in% y 
                            & month %in% m,
                            .(Lon,Lat,grid_id,month,year,SST_hat_local)])
  
  
  #----- generate noise ------------
  if(PCA_depth == 0){
    no = rep(0,length(fc[,SST_hat_local]))
  }else {
    no <- c()
    for(mon in m){
      load(file = paste0(cov.dir,"CovOU_mon",mon,".RData"))  
      load(file = paste0(cov.dir,"Cov_FU_mon",mon,".RData"))  
      
      A = matrix(c(cov_for_unc,cov_obs_unc), 
                 nrow = dim(cov_obs_unc)[1])
      
      PCA <- irlba(crossprod(A),nv = PCA_depth)                       
      
      for(year in y)
      {
        a <- A %*% PCA$v  %*% rnorm(PCA_depth)
        no <- c(no, a)
      }
    }
  }
  fc = fc[,"noise":= no]
  fc=fc[,"forecast":=SST_hat_local+noise]
  fc=fc[,c("SST_hat_local","noise"):=NULL]
  
  #-------- add land --------------
  
  fc_land <- list()
  fc_land[[1]] = fc
  fc_land[[2]] = land_grid_id
  fc_land = rbindlist(fc_land, fill = TRUE)
  fc_land = fc_land[order(fc_land[["grid_id"]])]
  
  #-------- truncate negative temperatures
  
  if(truncate) fc_land[forecast < -1.62, forecast := -1.62]
  
  
  #-------- save -------
  
  if(saveorgo){ 
    save(fc_land, file = paste0(save.dir,"/fc_",PCA_depth,"pc_",y,"_",m,".RData"))
  }else  return(fc_land)
}

    
load(file = paste0("~/PostClimDataNoBackup/SFE/Derived/dt_reduced_PCA.RData"))




vec1 = c(0,5,10,15,25,50,100,200)
for(d in vec1){
  print(d)
  forecast_PCA(dt_reduced = dt_reduced,y=2000, m=1, PCA_depth = d)
  plot_forecast(YM_j=2000*12+1,depth = d)
}
  
