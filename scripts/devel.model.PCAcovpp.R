rm(list = ls())



##-------- Setup ---------
library(SeasonalForecasting)
library(irlba)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/SFE/PCACov/"
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
save(emp_cov,
     file = paste0(data.dir,"Cov_",mon,".RData"))

}


#--- PCA ---

mon = 1
PCA_depth = 5
load(file = paste0(data.dir,"Cov_",mon,".RData"))
PCA <- irlba(emp_cov,nv = PCA_depth, right_only = TRUE)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

#cov_reduced <- tcrossprod(PCA$v %*% Diagonal(PCA_depth,PCA$d), PCA$v)

noise <- PCA$v %*% sqrt(Diagonal(PCA_depth,PCA$d)) %*% rnorm(PCA_depth)


