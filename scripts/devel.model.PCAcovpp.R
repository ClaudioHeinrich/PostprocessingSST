rm(list = ls())



##-------- Setup ---------
library(SeasonalForecasting)
library(irlba)
setwd("~/NR/SFE/")
options(max.print = 1e3)


  
dt = load_combined_wide()

ens.num = 9

#--- bring out the trash and reshape ---

dt_reduced <- dt[,c("Lon","Lat",paste0("SST", 1:10),"SST_bar","SST_sd","Ens_sd"):=NULL,]
dt_reduced = dt_reduced[,"obs_mean" := mean( Ens_bar, na.rm = TRUE), by = .(month, grid_id)]
dt_reduced = dt_reduced[,"Ens_bar":=NULL,]
dt_reduced = melt(dt_reduced, measure.vars = paste0("Ens",1:ens.num))
dt_reduced = na.omit(dt_reduced[,"cov_vec" :=  value - obs_mean , by = .(month, grid_id)])

#---- compute empirical covariance matrix

y_range= c(min(dt_reduced[,year]),max(dt_reduced[,year]))
cov_size = length(unique(dt_reduced[,grid_id]))

for(mon in 1:12){
  print(paste0("Month = ",mon))
  emp_cov = mat.or.vec(cov_size,cov_size)
  gc()
  
  for (y in y_range[1]:y_range[2]){
    print(paste0("year = ",y))
    for(ens in 1:ens.num){
      print(paste0("ensemble = ",ens))
      emp_cov = emp_cov + 
                tcrossprod(dt_reduced[month == mon & year == y & variable == paste0("Ens",ens),cov_vec],dt_reduced[month == mon & year == y & variable == paste0("Ens",ens),cov_vec])
      gc()
    }
  }

emp_cov = emp_cov/((y_range[2]-y_range[1])*ens.num-1)
data.dir = "~/PostClimDataNoBackup/SFE/PCACov"
save(emp_cov,
     file = paste0(data.dir,"Cov_",mon,".RData"))

}


temp <- dt_reduced[,cov_vec]

mon = 1

Cov_Mat = function(mon) 




corr_fac = 9*y_range / (9*y_range - 1)

dt_reduced = dt_reduced[, "emp_cov_vec" :=  - 9*y_range/(9*y_range-1)*obs_mean, by = .(month, grid_id)]
emp_cov[obs_mean == mean_obs,.N]

emp_cov=melt(emp_cov, measure.vars = paste0("Ens",1:9))


emp_cov = dt[,sum .SD - obs_mean, .SD = paste0("Ens_dev_",1:9)]



dt_cov = dt[,c(.SD,"SST_bar","SST_sd") := NULL,.SD=paste0("SST",1:10)]

dt[,paste0("SST_", 1:10,"plus_1"):=.SD +1.0,.SD=paste0("SST",1:10)]



emp_cov = dt[, "sum_obs" := rowSums(.SD), .SD = paste0("Ens",1:9)]
emp_cov = dt[, "mean_obs" := mean(sum_obs, na.rm = TRUE)/9, by = .(month, grid_id