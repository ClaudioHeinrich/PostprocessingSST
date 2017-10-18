rm(list = ls())



##-------- Setup ---------
library(SeasonalForecasting)
library(irlba)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/SFE/PCACov/"
options(max.print = 1e3)

#------setup PCA

dt = load_combined_wide()
obs.num = 10

#---- bias correction -----

dt[,"Loc_Bias_Est" := (cumsum(SST_bar-Ens_bar) - (SST_bar-Ens_bar)) / (year - min(year)),.(grid_id, month)]
dt[,"SST_hat_local":=Ens_bar + Loc_Bias_Est]

#--- bring out the trash, reshape and save ---

dt_reduced <- dt[,c(paste0("Ens", 1:9),"Ens_bar","SST_sd","Ens_sd"):=NULL,]
dt_reduced = dt_reduced[,"obs_mean_y" := mean( SST_bar, na.rm = TRUE), by = .(month, grid_id)]
dt_reduced = dt_reduced[,"SST_bar":=NULL,]
dt_reduced = melt(dt_reduced, measure.vars = paste0("SST",1:obs.num))
dt_reduced = dt_reduced[,"cov_vec" :=  value - obs_mean_y , by = .(month, grid_id)]
save(dt_reduced, file = paste0(data.dir,"../Derived/dt_reduced_PCA.RData"))








#---- compute empirical covariance matrices

dt_reduced = na.omit(load(file = paste0(data.dir,"../Derived/dt_reduced_PCA.RData")))
y_range= c(min(dt_reduced[,year]),max(dt_reduced[,year]))
num_y = y_range[2]-y_range[1]+1
cov_size_row = length(unique(dt_reduced[,grid_id]))
cov_size_col = num_y * obs.num

for(mon in 1:12){
  A=matrix(dt_reduced[month == mon, cov_vec], nrow = cov_size_row, ncol = cov_size_col, byrow = FALSE)/sqrt(cov_size_col-1)
  # the empirical covariance matrix is A %*% t(A), but does not need to be computed
  save(A, file = paste0(data.dir,"/Cov_",mon,".RData"))
}




#--- PCA ---
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
load(file = paste0("~/PostClimDataNoBackup/SFE/Derived/dt_reduced_PCA.RData"))

forecast_PCA = function(dt_reduced, 
                        y = 1999, 
                        m = 11, 
                        PCA_depth = 5, 
                        save.dir="./Data/PostClim/SFE/Derived/PCA",
                        data.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                        saveorgo = TRUE){
  
  land_grid_id = dt_reduced[variable == "SST1" & year %in% y & month %in% m & is.na(SST_hat_local),.(Lon,Lat,grid_id,month,year)]
  
  fc <- na.omit(dt_reduced[variable == "SST1" & year %in% y & month %in% m,.(Lon,Lat,grid_id,month,year,SST_hat_local)])
  
  
  #----- generate noise ------------
  if(PCA_depth == 0){ no = rep(0,length(fc[,SST_hat_local]))
  }else {
    no <- c()
    for(mon in m){
      load(file = paste0(data.dir,"Cov_",mon,".RData"))  
      PCA <- irlba(crossprod(A),nv = PCA_depth)                       
    
      for(year in y){
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
  
  if(saveorgo){ 
    save(fc_land, file = paste0(save.dir,"/fc_",PCA_depth,"pc_",y,"_",m,".RData"))
  }else  return(fc_land)
}




