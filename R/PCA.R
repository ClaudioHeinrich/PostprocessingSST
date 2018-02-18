#' Creates and saves the covariance matrix for the residuals trained on a specific set of years
#'
#' @param dt The wide data.table
#' @param Y training years.  So this will use all years in order to form the covariance matrix
#' @param M training months the PCA should be computed for
#' @param save.dir
#' @param obs.num The number of observations in the observational product
#' @param ens.num The number of ensemble members
#'
#' @export
for_res_cov = function(dt = NULL,
                       Y = 1985:2010,
                       M = 1:12,
                       save.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       obs.num = 10,
                       ens.num = 9
){  
  
  if(is.null(dt)){
    print("loading data")
    dt = load_combined_wide(bias = TRUE)
  }
  
  
  year.num = max(Y)-min(Y)+1
  
  for(mon in M){
    print(paste0("month =",mon))  
    
    dt_PCA = copy(dt)
    trash = c(paste0("SST", 1:obs.num),"SST_sd","Ens_sd")
    dt_PCA = dt_PCA[,(trash):=NULL]
    dt_PCA = dt_PCA[month == mon & year %in% Y,]
    
    for(ens in  1:ens.num){
      print(paste0("ensemble =",ens))  
      dt_PCA = dt_PCA[,paste0("Res",ens):= eval(parse(text = paste0("Ens",ens)))+Bias_Est-SST_bar]
    }
    dt_PCA = dt_PCA[,"res_mean" := mean(SST_hat- SST_bar), by =  grid_id]
    
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


#' Sets up the PCA
#'
#' @param dt The wide data table
#' @param m the month you want the PCA to be performed for
#' @param y the years you want to perform the PCA for
#' @param max_PCA_depth The total number of potential PCs
#' @param cov.dir The directory with the covariances to be saved
#' @param data.dir Save directory
#' @param oriented tries to orient the eigenvectors in the same direction (e.g. for visualisation of PCs) if TRUE. This is done
#' by multiplying the left  singular vectors by -1 if this 
#' decreases the Euclidean distance to the 1st left singular vector
#' @export
setup_PCA = function(dt=NULL,
                     m=7,
                     y = 1999,
                     max_PCA_depth = 100,
                     cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                     data.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                     oriented = TRUE
){
  
  if(is.null(dt)) { print("load and prepare data")
    dt = load_combined_wide(bias = TRUE)
    trash = c(paste0("SST",1:10),paste0("Ens",1:9))
    dt[, (trash):=NULL]
    dt = dt[year %in% y & month %in% m,]
    
  }
  
  print("loading and data reduction complete")
  
  
  #find land grid ids:
  
  land_grid_id <<- dt[year %in% y & month %in% m & (is.na(Ens_bar) | is.na(SST_bar)),
                      .(Lon,Lat,grid_id,month,year)]
  print("finding land complete")
  
  
  fc <<- na.omit( dt[,.(Lon,Lat,grid_id,month,year,SST_hat,SST_bar)])
  
  
  print("creating fc complete")
  
  
                                        #get covariance matrices
  
  for(mon in m){
    
    load(file = paste0(cov.dir,"CovRes_mon",mon,".RData"))
    
    assign(paste0("A",mon),
           matrix(res_cov, 
                  nrow = dim(res_cov)[1]),
           envir = globalenv())
  }
  
  PCA = irlba(eval(parse(text = paste0 ("A",mon))), nv = max_PCA_depth)
  
  if(oriented){
    for(d in 2:max_PCA_depth){
      if(sum((PCA$u[,1]-PCA$u[,d])^2) > sum((PCA$u[,1]+PCA$u[,d])^2)) PCA$u[,d] = -PCA$u[,d]     
    }
  }
  
  assign(paste0("PCA",mon), PCA, envir = globalenv())
  print(paste0("Month ",mon," complete"))
}



#' Forecast using PCA
#' @export
forecast_PCA = function(y = 1999, 
                        m = 7, 
                        PCA_depth = 5,  #accepts 0, then the mean of the bias corrected 
                        #ensemble is returned. If saveorgo = TRUE, PCA_depth can be
                        #a vector
                        save.dir="./Data/PostClim/SFE/Derived/PCA",
                        cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                        data.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                        output_opts = "forecast", # also takes "mar_sd", then the approximate marginal 
                        # variance is returned, 
                        # or "PC" where the dth Eigenvector
                        # is returned (d=PCA_depth)
                        # or "PCsum" where the sum over the first d PCs is returned
                        saveorgo = TRUE,
                        truncate = TRUE,
                        max_PCA_depth = 100) {
  
  # Check whether setup_PCA has been run, beware that setup_PCA has been called with the right parameters
  
  if(! exists("land_grid_id")) setup_PCA(m = m, max_PCA_depth = max_PCA_depth)
  print("setup complete")
  
  fc = fc[year %in% y][month %in% m]
  land_grid_id = land_grid_id[year %in% y][month %in% m]
  
  #----- generate noise ------------
  
  for (d in PCA_depth){
    no <- c()
    
    for(mon in m){
      A <- eval(parse(text = paste0("A",mon)))
      PCA <- eval(parse(text = paste0("PCA",mon)))
      eigen_vectors <- PCA$u[,1:d]
      sing_values <- diag(x = PCA$d[1:d], nrow = length(PCA$d[1:d]))
      
      for(year in y){
        
        if(output_opts == "forecast") a <- eigen_vectors  %*% sing_values %*% rnorm(d)
        if(output_opts == "forecast" & d == 0) a <- rep(0, times = dim(PCA$u)[1])
        if(output_opts == "mar_sd") a <- sqrt( (eigen_vectors %*% sing_values)^2  %*% rep(1,d))
        if(output_opts == "PC") {
          if(d == 1) a <-  PCA$d[d]*eigen_vectors else a <-  PCA$d[d]*eigen_vectors[,d]
        }
        if(output_opts == "PCsum") a <-  eigen_vectors %*% sing_values  %*% rep(1,d)
        
        no <- c(no, a)
      }
      
    }
    
    if(output_opts == "forecast")  {
      fc = fc[,"noise":= no]
      fc=fc[,"forecast":=SST_hat + noise]
    }
    
    if(output_opts == "mar_sd" | output_opts == "PC" | output_opts == "PCsum"){
      fc[, "forecast" := no]
      #because the plotting functions plot the forecast
    }
    
    #-------- truncate negative temperatures
    
    if(truncate & output_opts %in% c("forecast","PC")) fc[forecast < -1.619995, forecast := -1.619995]
    if(truncate & output_opts %in% c("mar_sd","PCsum")) fc[SST_hat < -1.619995, forecast := 0]
    
    #-------- add land --------------
    
    fc_land <- list()
    fc_land[[1]] = fc
    fc_land[[2]] = land_grid_id
    fc_land = rbindlist(fc_land, fill = TRUE)
    fc_land = fc_land[ order(year,month,grid_id)]
    
    
    #-------- save -------
    
    if(saveorgo & output_opts == "forecast"){ 
      save(fc_land, file = paste0(save.dir,"/fc_",d,"pc_",y,"_",m,".RData"))
    }
    if(saveorgo & output_opts == "mar_sd"){ 
      save(fc_land, file = paste0(save.dir,"/PCA_mar_sd",d,"_m",m,".RData"))
    }
    
    if(saveorgo & output_opts == "PC"){ 
      save(fc_land, file = paste0(save.dir,"/PCA_",d,"PC.RData"))
    }
    
    if(saveorgo & output_opts == "PCsum"){ 
      save(fc_land, file = paste0(save.dir,"/PCA_",d,"sum.RData"))
    }
    
    
  }
  
  return(fc_land)
}




