
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

for_res_cov = function(dt = NULL,
                       Y = 1985:2000,
                       M = 1:12,
                       save_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       ens_size = 9,
                       centering = "b"){  
  
  if(is.null(dt))
  {
    print("loading data")
    dt = load_combined_wide(var = TRUE)
  }
  
  
  
  for(mon in M){
    print(paste0("month =",mon))  
    
    dt_PCA = copy(dt[month == mon & year %in% Y,])
    
    dt_PCA[,paste0("Res",1:ens_size):= .SD-SST_bar,.SDcols = paste0("Ens",1:ens_size)]
    
    
    if(centering == "b")
    {
      dt_PCA[,"center" := - Bias_Est]
      decor_factor = 1 / sqrt( length(Y) * ens_size )
    }
    if(centering == "frm")
    {
      dt_PCA[,"center" := mean(rowSums(.SD)/ens_size),.SDcols = paste0("Res",1:9), by =.(grid_id, month)]
      decor_factor = 1 / sqrt(length(Y) * ens_size - 1)
    }
    if(centering == "rmby")
    {
      dt_PCA[,"center" := rowSums(.SD)/ens_size,.SDcols = paste0("Res",1:9)]  
      decor_factor = 1 / sqrt(length(Y)*(ens_size -1))
    }
    
    dt_PCA[,c(paste0("Res",1:ens_size,"new")) := .SD - center, .SDcols = paste0("Res",1:9)]
    
    sqrt_cov_mat = as.matrix(na.omit(dt_PCA[,.SD - center,.SDcols = paste0("Res",1:9)]))
    
    # sqrt_cov_mat = c()
    # for( ens in 1:ens_size){
    #   sqrt_cov_mat = c(sqrt_cov_mat, 
    #                    na.omit(dt_PCA[,eval(parse(text = paste0("Res",ens,"new")))]))
    # }
    # sqrt_cov_mat = matrix(sqrt_cov_mat,
    #                       ncol = ens_size * length(Y))
    
    res_cov = decor_factor * matrix(sqrt_cov_mat,ncol = length(Y) * ens_size)
    
    save(res_cov, file = paste0(save_dir,"CovRes_mon",mon,".RData"))
    
    rm(dt_PCA)
  }
}



for_res_cov_wrtm = function(dt = NULL,
                       Y = 1985:2000,
                       M = 1:12,
                       save_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                       centering = "b",
                       ens_size = 9){  
  
  if(is.null(dt))
  {
    print("loading data")
    dt = load_combined_wide(var = TRUE)
  }
  
  
  
  for(mon in M){
    print(paste0("month =",mon))  
    
    dt_PCA = copy(dt[month == mon & year %in% Y,])
    
    if(centering == "b")
    {
      dt_PCA[,"Res":= SST_bar - trc(Ens_bar + Bias_Est)]  
      cor_factor = 1 / sqrt( length(Y))
    }
    
    
    sqrt_cov_mat = as.matrix(na.omit(dt_PCA[,Res]))
    
    # sqrt_cov_mat = c()
    # for( ens in 1:ens_size){
    #   sqrt_cov_mat = c(sqrt_cov_mat, 
    #                    na.omit(dt_PCA[,eval(parse(text = paste0("Res",ens,"new")))]))
    # }
    # sqrt_cov_mat = matrix(sqrt_cov_mat,
    #                       ncol = ens_size * length(Y))
    
    res_cov = cor_factor * matrix(sqrt_cov_mat,ncol = length(Y) )
    
    save(res_cov, file = paste0(save_dir,"CovRes_mon",mon,".RData"))
    
    rm(dt_PCA)
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
    dt = load_combined_wide(var = TRUE)
    trash = c(paste0("SST",1:10),paste0("Ens",1:9))
    dt[, (trash):=NULL]
    dt = dt[year %in% y & month %in% m,]
  }
  
    #find land grid ids:
  
  land_grid_id <<- dt[year %in% y & month %in% m & (is.na(Ens_bar) | is.na(SST_bar)),
                      .(Lon,Lat,grid_id,month,year)]
  
  
  fc <<- na.omit( dt[,.(Lon,Lat,grid_id,month,year,SST_hat,SST_bar)])
  
  
                                        #get covariance matrices
  
  for(mon in m){
    load(file = paste0(cov_dir,"CovRes_mon",mon,".RData"))
    
    assign(paste0("A",mon),
           matrix(res_cov, 
                  nrow = dim(res_cov)[1]),
           envir = globalenv())
  
  
  PCA = irlba::irlba(eval(parse(text = paste0 ("A",mon))), nv = max_PCA_depth)
  
  
  assign(paste0("PCA",mon), PCA, envir = globalenv())
  }
  print("setup for PCA done.")
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
                            n = 10, 
                            PCA_depth = 15,  
                            ens_member = FALSE,
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
                      .(Lon,Lat,grid_id,month,year,YM)]
  SD_cols = c("Lon","Lat","grid_id","month","year","YM",
              "SST_hat","SST_bar",paste0("Ens",1:ens_size),"Ens_bar","Bias_Est","var_bar","SD_hat")
  fc <- na.omit( dt[year %in% y & month %in% m ,.SD,.SDcols = SD_cols])
  
  
  #get covariance matrices
  
  print("data preparation complete - getting covariance matrices next:")
  
  for(mon in m){
      load(file = paste0(cov_dir,"CovRes_mon",mon,".RData"))
      
      assign(paste0("A",mon),
             matrix(res_cov, nrow = dim(res_cov)[1]))
      
      
      PCA = irlba::irlba(eval(parse(text = paste0 ("A",mon))), nv = max(PCA_depth))
      
      
      assign(paste0("PCA",mon), PCA)
      
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
      mar_sds <- sqrt(rowSums((eigen_vectors %*% sing_values)^2))
      crit_ind  = which(mar_sds < 1e-20)
      
      for(y_0 in y){
        
        if(d == 0)  {
          a <- matrix(0,nrow =  dim(PCA$u)[1],ncol = n)
        }else{ 
          var_correct_vec = fc[year == y_0,][month == mon,SD_hat] / mar_sds
          var_correct_vec[crit_ind] = 1
          var_correct = diag(var_correct_vec)
          a = var_correct %*% eigen_vectors  %*% sing_values %*% matrix(rnorm(d * n),nrow = d, ncol = n) / sqrt(ens_size)
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


#' getting PCs for visualization/visual analysis
#' 
#' 
#' @param dt the data table, if NULL, it is loaded.
#' @param y,m Integer vectors containing year(s) and month(s).
#' @param PCA_depth the output contains the columns PC1, ..., PC(PCA_depth)
#' @param saveorgo Logical, whether we save or not. 
#' @param save_dir The directory to save in.
#' @param cov_dir Where the PCA covariance matrices are stored.
#' 
#' @return data table containing PCA_depth columns containing the PCs
#' 
#' @examples \dontrun{get_PCs(y = 1999, m = 7)}
#' 
#' @author Claudio Heinrich        
#' 
#' @export

get_PCs = function(dt = NULL,
                   y,
                   m,
                   PCA_depth = 4,
                   mar_var_correct = TRUE,
                   saveorgo = FALSE,
                   save_dir = "./Data/PostClim/SFE/Derived/PCA/",
                   file_name = "PCs.RData",
                   cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/")
{
  
  
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
                     .(Lon,Lat,grid_id,month,year,YM)]
  SD_cols = c("Lon","Lat","grid_id","month","year","YM",
              "SST_hat","SST_bar",paste0("Ens",1:ens_size),"Ens_bar","Bias_Est","var_bar","SD_hat")
  
  fc <- na.omit( dt[year %in% y & month %in% m ,.SD,.SDcols = SD_cols])
  
  
  #get covariance matrices
  
  print("data preparation complete - getting covariance matrices next:")
  
  for(mon in m)
  {
    load(file = paste0(cov_dir,"CovRes_mon",mon,".RData"))
    
    PCA = irlba::irlba(res_cov, nv = PCA_depth)
    
    for(y_0 in y)
    {
      for (d in 1:PCA_depth)
        {
          fc[month == mon & year == y_0, paste0("PC",d) := PCA$d[d] * PCA$u[,d]]
        }
      
      if(mar_var_correct)
      {
        eigen_vectors <- PCA$u[,1:PCA_depth]
        sing_values <- diag(x = PCA$d[1:PCA_depth], nrow = length(PCA$d[1:PCA_depth]))
        mar_sds <- sqrt(rowSums(( eigen_vectors %*% sing_values)^2))
        
        for(d in 1:PCA_depth)
        {
          fc[month == mon & year == y_0, paste0("PC_marcor_",d) := SD_hat * PCA$d[d] * PCA$u[,d] / mar_sds]  
        }
      }
    }
  }
  
 
  fc_land <- list()
  fc_land[[1]] = fc
  fc_land[[2]] = land_grid_id
  fc_land = rbindlist(fc_land, fill = TRUE)
  fc_land = fc_land[ order(year,month,grid_id)]
  fc = fc_land
  
  #-------- save -------
  
  if(saveorgo)
  { 
    save(fc, file = paste0(save_dir,file_name))
  }
  
  return(fc)
}


