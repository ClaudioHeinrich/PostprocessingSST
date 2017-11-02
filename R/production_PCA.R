

for_unc_cov = function(dt = NULL,
                       Y = 1985:2010,
                       M = 1:12,
                       save.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       obs.num = 10,
                       ens.num = 9
                       ){
  
  #------setup PCA
  
  if(is.null(dt)) dt = load_combined_wide()
  
  #---- bias correction -----
  
  dt = bias_correct(dt = dt)
  
  year.num = max(Y)-min(Y)+1
  
  for(mon in M){
    print(paste0("month =",mon))  
    
    for (obs in 1:obs.num){
      
      print(paste0("obs =",obs))
      
      dt_PCA = dt[month == mon,]
      trash = c(paste0("SST", 1:obs.num),"SST_sd","SST_bar","Ens_sd")
      trash = trash[! trash == paste0("SST",obs)]
      dt_PCA = dt_PCA[,(trash):=NULL,]
      
      for(ens in  1:ens.num){
        dt_PCA = dt_PCA[,paste0("Res",ens):= eval(parse(text = paste0("Ens",ens)))+Bias_Est-eval(parse(text = paste0("SST",obs)))]
        dt_PCA = dt_PCA[, paste0("Ens",ens):=NULL]
      }
      dt_PCA = dt_PCA[,"res_mean" := mean(Ens_bar +Bias_Est- eval(parse(text = paste0("SST",obs)))), by =  grid_id]
      
      sqrt_cov_mat = c()
      for( ens in 1:ens.num){
        sqrt_cov_mat = c(sqrt_cov_mat, 
                         na.omit(dt_PCA[,eval(parse(text = paste0("Res",ens))) - res_mean]))
      }
      sqrt_cov_mat = matrix(sqrt_cov_mat,
                            nrow = length(na.omit(dt_PCA[,res_mean]))/year.num,
                            ncol = year.num*ens.num) 
      sqrt_cov_mat = sqrt_cov_mat/(sqrt(year.num*ens.num -1))
      
      save(sqrt_cov_mat, file = paste0(save.dir,"CovFU_mon",mon,"_obs",obs,".RData"))
      
      rm(dt_PCA)
    }
  }
}


for_unc_cov_combined = function(M = 1:12,
                                cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                                obs.num = 10){
  
  for(mon in M){
    load(file = paste0(cov.dir,"CovFU_mon",mon,"_obs",1,".RData"))
    
    cov_for_unc = sqrt_cov_mat
    
    if (obs.num > 1){
      for(obs in 2:obs.num){
        load(file = paste0(cov.dir,"CovFU_mon",mon,"_obs",obs,".RData"))
        cov_for_unc = matrix(c(cov_for_unc,sqrt_cov_mat),
                             nrow = dim(cov_for_unc)[1],
                             ncol = dim(cov_for_unc)[2] + dim(sqrt_cov_mat)[2]
        )
        
      }
    }
    
    cov_for_unc = cov_for_unc/sqrt(obs.num)
    
    save(cov_for_unc, file = paste0(cov.dir,"Cov_FU_mon",mon,".RData"))
  }
}





obs_unc_cov = function(M = 1:12,
                       save.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       obs.num = 10,
                       ens.num = 9,
                       load_reduced_data = TRUE,
                       data.dir = "~/PostClimDataNoBackup/SFE/PCACov/"
){
  
  #--- bring out the trash, reshape and save ---
  
  if(! load_reduced_data){
    
    dt = load_combined_wide()
    
    dt = bias_correct(dt = dt)
    
    dt_reduced <- dt[,c(paste0("Ens", 1:9),"Ens_bar","SST_sd","Ens_sd"):=NULL,]
    dt_reduced = dt_reduced[,"obs_mean_y" := mean( SST_bar, na.rm = TRUE), by = .(month, grid_id)]
    dt_reduced = dt_reduced[,"SST_bar":=NULL,]
    dt_reduced = melt(dt_reduced, measure.vars = paste0("SST",1:obs.num))
    dt_reduced = dt_reduced[,"cov_vec" :=  value - obs_mean_y , by = .(month, grid_id)]
    save(dt_reduced, file = paste0(data.dir,"../Derived/dt_reduced_PCA.RData"))
    print("reduced data saved")
  }
  
  #---- compute empirical covariance matrices
  
  load(file = paste0(data.dir,"../Derived/dt_reduced_PCA.RData"))
  dt_reduced_water = na.omit(dt_reduced)
  y_range = range(dt_reduced_water[,year])
  num_y = y_range[2]-y_range[1]+1
  cov_size_row = length(unique(dt_reduced_water[,grid_id]))
  cov_size_col = num_y * obs.num
  
  
  for(mon in M){
    print(paste0("month = ",mon))
    cov_obs_unc=matrix(na.omit(dt_reduced_water[month == mon, cov_vec]), 
                       nrow = cov_size_row, 
                       ncol = cov_size_col, 
                       byrow = FALSE)/sqrt(cov_size_col-1)
    
    # the empirical covariance matrix is cov_obs_unc %*% t(cov_obs_unc), but does not need to be computed
    save(cov_obs_unc, file = paste0(data.dir,"/CovOU_mon",mon,".RData"))
  }
}
