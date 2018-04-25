#' Creates and saves the covariance matrix for the residuals trained on a specific set of years
#'
#' @param dt The wide data.table
#' @param Y training years.  So this will use all years in order to form the covariance matrix
#' @param M training months the PCA should be computed for.
#' @param save_dir Where we save.
#' @param ens_size The number of ensemble members.
#'
#' @author Claudio Heinrich
#' @examples 
#' \dontrun{for_res_cov()}
#' 
#' @export
for_res_cov = function(dt = NULL,
                       Y = 1985:2010,
                       M = 1:12,
                       save_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       ens_size = 9
){  
  
  if(is.null(dt))
  {
    print("loading data")
    dt = load_combined_wide(bias = TRUE)
  }
  
  year.num = max(Y)-min(Y)+1
  
  for(mon in M){
    print(paste0("month =",mon))  
    
    dt_PCA = copy(dt)
    dt_PCA = dt_PCA[month == mon & year %in% Y,]
    
    for(ens in  1:ens_size){
      print(paste0("ensemble =",ens))  
      dt_PCA = dt_PCA[,paste0("Res",ens):= eval(parse(text = paste0("Ens",ens)))+Bias_Est-SST_bar]
    }
    dt_PCA = dt_PCA[,"res_mean" := mean(SST_hat- SST_bar), by =  grid_id]
    
    sqrt_cov_mat = c()
    for( ens in 1:ens_size){
      sqrt_cov_mat = c(sqrt_cov_mat, 
                       na.omit(dt_PCA[,eval(parse(text = paste0("Res",ens))) - res_mean]))
    }
    sqrt_cov_mat = matrix(sqrt_cov_mat,
                          nrow = length(na.omit(dt_PCA[,res_mean]))/year.num)
    
    res_cov = sqrt_cov_mat/(sqrt(year.num*ens_size -1))
    
    save(res_cov, file = paste0(save_dir,"CovRes_mon",mon,".RData"))
    
    rm(dt_PCA)
  }
  
}


#' Helping function
#' 
#' @description Sets up for_unc_cov_combined.
#'
#' @param dt The wide data table.
#' @param M The months you want the PCA to be performed for.
#' @param Y The years you want to perform the PCA for.
#' @param save_dir The directory in which the covariances are saved.
#' @param obs_size The number of observations in the observational product.
#' @param ens_size The number of ensemble members.
#' 
#' @examples \dontrun{for_unc_cov(Y = 1985:1986)}
#'
#' @author Claudio Heinrich
#' 
#' @export


for_unc_cov = function(dt = NULL,
                       Y = 1985:2010,
                       M = 1:12,
                       save_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       obs_size = 10,
                       ens_size = 9
                       ){
  
  if(is.null(dt)) dt = load_combined_wide(bias = TRUE)
  
  year.num = max(Y)-min(Y)+1
  
  for(mon in M){
    print(paste0("month =",mon))  
    
    for (obs in 1:obs_size){
      
      print(paste0("obs =",obs))
      
      dt_PCA = dt[month == mon,]
      trash = c(paste0("SST", 1:obs_size),"SST_sd","SST_bar","Ens_sd")
      trash = trash[! trash == paste0("SST",obs)]
      dt_PCA = dt_PCA[,(trash):=NULL,]
      
      for(ens in  1:ens_size){
        dt_PCA = dt_PCA[,paste0("Res",ens):= eval(parse(text = paste0("Ens",ens)))+Bias_Est-eval(parse(text = paste0("SST",obs)))]
        dt_PCA = dt_PCA[, paste0("Ens",ens):=NULL]
      }
      dt_PCA = dt_PCA[,"res_mean" := mean(Ens_bar +Bias_Est- eval(parse(text = paste0("SST",obs)))), by =  grid_id]
      
      sqrt_cov_mat = c()
      for( ens in 1:ens_size){
        sqrt_cov_mat = c(sqrt_cov_mat, 
                         na.omit(dt_PCA[,eval(parse(text = paste0("Res",ens))) - res_mean]))
      }
      sqrt_cov_mat = matrix(sqrt_cov_mat,
                            nrow = length(na.omit(dt_PCA[,res_mean]))/year.num,
                            ncol = year.num*ens_size) 
      sqrt_cov_mat = sqrt_cov_mat/(sqrt(year.num*ens_size -1))
      
      save(sqrt_cov_mat, file = paste0(save_dir,"CovFU_mon",mon,"_obs",obs,".RData"))
      
      rm(dt_PCA)
    }
  }
}


#' Creates and saves the covariance matrix for the residuals taking observation uncertainty into account.
#'
#' @param M Training months the PCA should be computed for.
#' @param cov_dir Where the covariances are saved, should match save_dir of for_unc_cov.
#' @param obs_size The number of observations in the observational product.
#'
#' @author Claudio Heinrich
#' 
#' @examples \dontrun{for_unc_cov_combined()}
#' 
#' @export

for_unc_cov_combined = function(M = 1:12,
                                cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                                obs_size = 10){
  
  for(mon in M){
    load(file = paste0(cov_dir,"CovFU_mon",mon,"_obs",1,".RData"))
    
    cov_for_unc = sqrt_cov_mat
    
    if (obs_size > 1){
      for(obs in 2:obs_size){
        load(file = paste0(cov_dir,"CovFU_mon",mon,"_obs",obs,".RData"))
        cov_for_unc = matrix(c(cov_for_unc,sqrt_cov_mat),
                             nrow = dim(cov_for_unc)[1],
                             ncol = dim(cov_for_unc)[2] + dim(sqrt_cov_mat)[2]
        )
        
      }
    }
    
    cov_for_unc = cov_for_unc/sqrt(obs_size)
    
    save(cov_for_unc, file = paste0(cov_dir,"Cov_FU_mon",mon,".RData"))
  }
}


#' Creates and saves the covariance matrix for the observation for modelling only observation uncertainty.
#'
#' @param M Training months the PCA should be computed for.
#' @param save_dir The directory in which the covariances are saved.
#' @param obs_size The number of observations in the observational product.
#' @param ens_size The number of ensemble members.
#' @param load_reduced_data If FALSE, the reduced data table dt_reduced is computed and saved in data_dir, which takes time. If TRUE, the data table is loaded.
#' @param data_dir Where everything is saved.
#'
#' @examples \dontrun{obs_unc_cov()}
#'
#' @author Claudio Heinrich
#' 
#' @export

obs_unc_cov = function(M = 1:12,
                       save_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       obs_size = 10,
                       ens_size = 9,
                       load_reduced_data = TRUE,
                       data_dir = "~/PostClimDataNoBackup/SFE/PCACov/"
){
  
  #--- bring out the trash, reshape and save ---
  
  if(! load_reduced_data){
    
    dt = load_combined_wide()
    
    dt = bias_correct(dt = dt)
    
    dt_reduced <- dt[,c(paste0("Ens", 1:9),"Ens_bar","SST_sd","Ens_sd"):=NULL,]
    dt_reduced = dt_reduced[,"obs_mean_y" := mean( SST_bar, na.rm = TRUE), by = .(month, grid_id)]
    dt_reduced = dt_reduced[,"SST_bar":=NULL,]
    dt_reduced = melt(dt_reduced, measure.vars = paste0("SST",1:obs_size))
    dt_reduced = dt_reduced[,"cov_vec" :=  value - obs_mean_y , by = .(month, grid_id)]
    save(dt_reduced, file = paste0(data_dir,"../Derived/dt_reduced_PCA.RData"))
    print("reduced data saved")
  }
  
  #---- compute empirical covariance matrices
  
  load(file = paste0(data_dir,"../Derived/dt_reduced_PCA.RData"))
  dt_reduced_water = na.omit(dt_reduced)
  y_range = range(dt_reduced_water[,year])
  num_y = y_range[2]-y_range[1]+1
  cov_size_row = length(unique(dt_reduced_water[,grid_id]))
  cov_size_col = num_y * obs_size
  
  
  for(mon in M){
    print(paste0("month = ",mon))
    cov_obs_unc=matrix(na.omit(dt_reduced_water[month == mon, cov_vec]), 
                       nrow = cov_size_row, 
                       ncol = cov_size_col, 
                       byrow = FALSE)/sqrt(cov_size_col-1)
    
    # the empirical covariance matrix is cov_obs_unc %*% t(cov_obs_unc), but does not need to be computed
    save(cov_obs_unc, file = paste0(data_dir,"/CovOU_mon",mon,".RData"))
  }
}



#-----------------------------

#----  functions for issuing PCA forecasts ----

#' Organize data to be run in a PCA analysis
#'
#' @description A preprossesing script that creates the example
#' dataset from the original raw data
#'
#' @importFrom irlba irlba
#' 
#' @param dt The data.table of ensemble forecasts and associated observations including bias estimate, is loaded if NULL.
#' @param m Integer vector. Which month(s) to process.
#' @param y Integer vector. Which year(s) to process.
#' @param max_PCA_depth Integer. How many principal components should we consider at a maximum.
#' @param cov_dir String.  Where should we store the Covariance matrices?
#' @param data_dir String. Where should the data be stored?
#' 
#' @examples \dontrun{setup_PCA(m = 7)}
#' 
#' @author Claudio Heinrich
#' @export
setup_PCA = function(dt=NULL,
                     m=1:12,
                     y = 2001:2010,
                     max_PCA_depth = 50,
                     cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/")
{
  
  if(is.null(dt))
  {
    print("load and prepare data")
    dt = load_combined_wide(bias = TRUE)
    trash = c(paste0("SST",1:10),paste0("Ens",1:9))
    dt[, (trash):=NULL]
    dt = dt[year %in% y & month %in% m,]
  }
  
    #find land grid ids:
  
  land_grid_id <<- dt[year %in% y & month %in% m & (is.na(Ens_bar) | is.na(SST_bar)),
                      .(Lon,Lat,grid_id,month,year)]
  print("finding land complete")
  
  
  fc <<- na.omit( dt[,.(Lon,Lat,grid_id,month,year,SST_hat,SST_bar)])
  
  
  print("creating fc complete")
  
  
                                        #get covariance matrices
  
  for(mon in m){
    load(file = paste0(cov_dir,"CovRes_mon",mon,".RData"))
    
    assign(paste0("A",mon),
           matrix(res_cov, 
                  nrow = dim(res_cov)[1]),
           envir = globalenv())
  
  
  PCA = irlba::irlba(eval(parse(text = paste0 ("A",mon))), nv = max_PCA_depth)
  
  
  assign(paste0("PCA",mon), PCA, envir = globalenv())
  print(paste0("Month ",mon," complete"))

  }
  print("setup complete")
}



#' Forecast using PCA
#' 
#' @description Generates forecasts and related stuff for the PCA post-processing approach, \cr
#'      requires setup_PCA to be run with the corresponding year(s), month(s), and with max_PCA_depth >= PCA_depth
#' 
#' @param y,m Integer vectors containing year(s) and month(s).
#' @param output_opts,PCA_depth One of the following: \cr
#'      "forecast", generates a forecast for (y,m) perturbed by PCA-noise generated by d = PCA_depth principal components, d=0 returns the (unperturbed) bias corrected ensemble mean as forecast\cr
#'      "mar_sd", returns the approximate marginal standard deviation for d PCs, \cr 
#'      "PC", returns the dth principal component (upscaled eigenvector), \cr
#'      "PCsum", returns the sum over the first d PCs, \cr
#'           if saveorgo = TRUE, PCA_depth accepts integer vectors.
#' @param saveorgo Logical, whether we save or not. 
#' @param save_dir The directory to save in.
#' @param truncate logical, whether the forecasted temperature is truncated at freezing temperature
#' 
#' @return data table containing a column with output_opts as name, if output_opts = "forecast" it contains additionally a column with the noise.
#' 
#' @examples \dontrun{forecast_PCA(y = 1999, m = 7)}
#' 
#' @author Claudio Heinrich        
#' 
#' @export
forecast_PCA = function(y, m, 
                        output_opts = "forecast", 
                        PCA_depth = 15,  #accepts 0, then the mean of the bias corrected 
                        #ensemble is returned. If saveorgo = TRUE, PCA_depth can be
                        #a vector
                        saveorgo = TRUE,
                        save_dir="./Data/PostClim/SFE/Derived/PCA",
                        truncate = TRUE) {
  
  # Check whether setup_PCA has been run, beware that setup_PCA has been called with the right parameters
  
  if(! exists("land_grid_id")) stop("You need to run setup_PCA() first.")
    
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
      fc[, eval(quote(output_opts)) := no]
    }
    
    #-------- truncate negative temperatures
    
    if(truncate){
        tr.value = fc[,min(SST_bar)]
        if( output_opts == "forecast") fc[forecast < tr.value, forecast := tr.value]
        if( output_opts == "PC") fc[PC < tr.value, PC := tr.value]
        if(output_opts %in% c("mar_sd","PCsum")) fc[SST_hat < tr.value, eval(quote(output_opts)) := 0]
    }
    #-------- add land --------------
    
    fc_land <- list()
    fc_land[[1]] = fc
    fc_land[[2]] = land_grid_id
    fc_land = rbindlist(fc_land, fill = TRUE)
    fc_land = fc_land[ order(year,month,grid_id)]
    
    
    #-------- save -------
    
    if(saveorgo & output_opts == "forecast"){ 
      save(fc_land, file = paste0(save_dir,"/fc_",d,"pc_",y,"_",m,".RData"))
    }
    if(saveorgo & output_opts == "mar_sd"){ 
      save(fc_land, file = paste0(save_dir,"/PCA_mar_sd",d,"_m",m,".RData"))
    }
    
    if(saveorgo & output_opts == "PC"){ 
      save(fc_land, file = paste0(save_dir,"/PCA_",d,"PC_month",m,".RData"))
    }
    
    if(saveorgo & output_opts == "PCsum"){ 
      save(fc_land, file = paste0(save_dir,"/PCA_",d,"sum.RData"))
    }
    
    
  }
  
  return(fc_land)
}




#' New forecasts using PCA
#' 
#' @description Generates forecasts for the PCA post-processing approach, \cr
#'      does not require setup_PCA to be run first.
#'      Creates output samples, see n below, and can use ensemble members plus noise as forecast. 
#' 
#' @param dt the data table, if NULL, it is loaded.
#' @param y,m Integer vectors containing year(s) and month(s).
#' @param n Integer. Size of the desired forecast ensemble. 
#' @param PCA_depth Integer vector containing the numbers of considered principal components. Is allowed to contain 0 which refers to the forecast without noise.
#' @param ens_member Logical. If FALSE, the noise is added to the ensemble mean. 
#' @param ens_size Integer. Only needed if ens_member == TRUE. Size of the NWP ensemble.
#' @param saveorgo Logical, whether we save or not. 
#' @param save_dir The directory to save in.
#' @param cov_dir Where the PCA covariance matrices are stored.
#' @param truncate logical, whether the forecasted temperature is truncated at freezing temperature
#' 
#' @return data table containing n columns with noise and n columns with forecasts.
#' 
#' @examples \dontrun{forecast_PCA_new(y = 1999, m = 7)}
#' 
#' @author Claudio Heinrich        
#' 
#' @export

forecast_PCA_new = function(dt = NULL,
                            y,
                            m,
                            n = 1, 
                            PCA_depth = 15,  
                            ens_member = TRUE,
                            ens_size = 9,
                            saveorgo = TRUE,
                            save_dir = "./Data/PostClim/SFE/Derived/PCA/",
                            file_name = "forecastPCA",
                            cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                            truncate = TRUE) {
  
  
  if(is.null(dt))
  {
    print("load and prepare data")
    dt = load_combined_wide(bias = TRUE)
    trash = c(paste0("SST",1:10))
    dt[, (trash):=NULL]
    dt = dt[year %in% y & month %in% m,]
  }
  
  #find land grid ids:
  
  land_grid_id <- dt[year %in% y & month %in% m & (is.na(Ens_bar) | is.na(SST_bar)),
                      .(Lon,Lat,grid_id,month,year)]
  SD_cols = c("Lon","Lat","grid_id","month","year","YM",
              "SST_hat","SST_bar",paste0("Ens",1:ens_size),"Ens_bar","Bias_Est","SD_hat")
  fc <- na.omit( dt[year %in% y & month %in% m ,.SD,.SDcols = SD_cols])
  
  
  
  #get covariance matrices
  
  print("data preparation complete - getting covariance matrices next:")
  
  for(mon in m){
      load(file = paste0(cov_dir,"CovRes_mon",mon,".RData"))
      
      assign(paste0("A",mon),
             matrix(res_cov, nrow = dim(res_cov)[1]))
      
      
      PCA = irlba::irlba(eval(parse(text = paste0 ("A",mon))), nv = max(PCA_depth))
      
      
      assign(paste0("PCA",mon), PCA)
      print(paste0("Month ",mon," complete"))
      
    }
  
  
  print("setting up PCA complete - moving to forecasting")
  
  
  #----- generate noise ------------
  
  for (d in PCA_depth){
    
    print(paste0("Depth of PCA: ",d))
    
    no <- c()
    
    for(mon in m){
      A <- eval(parse(text = paste0("A",mon)))
      PCA <- eval(parse(text = paste0("PCA",mon)))
      eigen_vectors <- PCA$u[,1:d]
      sing_values <- diag(x = PCA$d[1:d], nrow = length(PCA$d[1:d]))
      
      for(year in y){
        
        if(d == 0)  {
          a <- matrix(0,nrow =  dim(PCA$u)[1],ncol = n)
        }else{
          a = eigen_vectors  %*% sing_values %*% matrix(rnorm(d * n),nrow = d, ncol = n)
        }
            
        no <- rbind(no, a)
      }
      
    }
    
    for(i in 1:n) {
      
      fc = fc[,  paste0("no",i,"PC",d):= no[,i]]
      
      if(ens_member){
      ens_mem = sample.int(ens_size,1) 
      
      fc=fc[,paste0("fc",i,"PC",d):= rowSums(.SD), .SDcols = c(paste0("Ens",ens_mem),"Bias_Est",paste0("no",i,"PC",d))]
      }else{
        fc=fc[,paste0("fc",i,"PC",d):= rowSums(.SD), 
              .SDcols = c("Ens_bar","Bias_Est",paste0("no",i,"PC",d))]
      }
      
      if(truncate){
        trc = function (x){ 
          truncation.value = -1.769995
          x = truncation.value * (x < truncation.value) + x * (x >= truncation.value)
          return(x)}
        
        fc[, paste0("fc",i,"PC",d) := lapply(.SD,trc), .SDcols = paste0("fc",i,"PC",d) ]
        }  
      }
  }

    
    #-------- add land --------------
    
    fc_land <- list()
    fc_land[[1]] = fc
    fc_land[[2]] = land_grid_id
    fc_land = rbindlist(fc_land, fill = TRUE)
    fc_land = fc_land[ order(year,month,grid_id)]
    fc = fc_land
    
    #-------- save -------
    
    if(saveorgo){ 
      save(fc, file = paste0(save_dir,file_name,".RData"))
    }
  return(fc)
}

