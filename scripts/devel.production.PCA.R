for_res_cov = function(dt = NULL,
                       Y = 1985:2010,
                       M = 1:12,
                       save.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       obs.num = 10,
                       ens.num = 9
){  
  
  if(is.null(dt)) dt = load_combined_wide()
  
  #---- bias correction -----
  
  dt = bias_correct(dt = dt)
  
  year.num = max(Y)-min(Y)+1
  
  for(mon in M){
    print(paste0("month =",mon))  
    
      dt_PCA = copy(dt)
      trash = c(paste0("SST", 1:obs.num),"SST_sd","Ens_sd")
      dt_PCA = dt_PCA[,(trash):=NULL,]
      dt_PCA = dt[month == mon,]
      
      for(ens in  1:ens.num){
        print(paste0("ensemble =",ens))  
          dt_PCA = dt_PCA[,paste0("Res",ens):= eval(parse(text = paste0("Ens",ens)))+Bias_Est-SST_bar]
        }
      dt_PCA = dt_PCA[,"res_mean" := mean(Ens_bar + Bias_Est- SST_bar), by =  grid_id]
      
      sqrt_cov_mat = c()
      for( ens in 1:ens.num){
        sqrt_cov_mat = c(sqrt_cov_mat, 
                         na.omit(dt_PCA[,eval(parse(text = paste0("Res",ens))) - res_mean]))
      }
      sqrt_cov_mat = matrix(sqrt_cov_mat,
                            nrow = length(na.omit(dt_PCA[,res_mean]))/year.num)
      
      res_cov = sqrt_cov_mat/(sqrt(year.num*ens.num -1))
      
      save(res_cov, file = paste0(save.dir,"CovRes_mon",mon,".RData"))
      
      rm(dt_PCA)
    }
  
}