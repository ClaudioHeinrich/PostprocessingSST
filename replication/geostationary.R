
#######################################################################
##  Estimate a monthly Exponential semi-variogram function.
#######################################################################




#' preparation for geostationary forecasts
#' 
#' @description Estimates the variogram for the residuals assuming anexponential covariance model and saves the results.
#' 
#' 
#' @param dt the data table, if NULL, it is loaded.
#' @param training_years,m Integer vectors containing the year(s) and month(s) of the training dataset.
#' @param ens_size Integer. Size of the NWP ensemble.
#' @param saveorgo Logical, whether we save or not. 
#' @param save_dir,file_name The directory to save in. The name of the file of the saved variogram for a given month is \code{file_name <month> .RData}.
#' @param n_intv Integer. How many distance bins are considered for the empirical variogram.
#' @param truncate Logical. The empirical variogram oftentimes is far from the fitted variogram for the 10% largest distances considered. If truncate == TRUE, those are ignored leading to a visually much better fit of the variogram.
#' 
#' @return data table containing n columns with noise and n columns with forecasts.
#' 
#' @examples \dontrun{ DT = load_combined_wide(data_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO", output_name = "dt_combine_NAO_wide_bc.RData") \cr
#'                     geostationary_training(dt = DT)}
#' 
#' @author Claudio Heinrich        
#' 
#' @importFrom spacetime STFDF
#' @importFrom sp SpatialPoints spDists SpatialPointsDataFrame
#' @importFrom gstat variogramST fit.variogram
#' 
#' @export


geostationary_training = function (dt = NULL,
                                   training_years = 1986:2000,
                                   m = 1:12,
                                   ens_size = 9,
                                   saveorgo = TRUE,
                                   save_dir = "./Data/PostClim/SFE/Derived/GeoStat/",
                                   file_name = "variogram_exp_m",
                                   nintv = 100,
                                   truncate = TRUE){

  
  if(is.null(dt))
  {
    print("load and prepare data")
    dt = load_combined_wide(bias = TRUE)
    dt = dt[year %in% training_years & month %in% m,]
  }
  
  Lon_min  = dt[,range(Lon)][1]
  Lon_max  = dt[,range(Lon)][2]
  Lat_min  = dt[,range(Lat)][1]
  Lat_max  = dt[,range(Lat)][2]
  
  
  for(i in 1:ens_size){
      dt = dt[,  paste0("Res",i) := .SD + Bias_Est - SST_bar,.SDcols = paste0("Ens",i)]
  }
  
  
  for(mon in m){
    print(paste0("month = ",mon))
    
    DT = dt[month == mon, .SD,.SDcols = c("Lon","Lat","year","month","YM", paste0("Res",1:ens_size))]
    
    sp <- SpatialPoints(cbind(x=DT[YM == min(YM), Lon],y=DT[YM == min(YM), Lat]), proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    # Convert YM into date
    
    time_convert = function(YM){
      M = YM %% 12
      M[M == 0] = 12
      Y = (YM - M)/ 12
      M[M < 10]  = paste0(0,M[M < 10])
      as.Date(paste0(Y,"-",M,"-15"))
    }
    
    time = as.POSIXct( time_convert(unique(DT[,YM])), tz = "GMT")
    
    setkey(DT,YM,Lon,Lat) # for creating STFDFs the data should be ordered such that the 'spatial index is moving fastest'
    random_pick = sample(ens_size,1)
    data = DT[,.SD,.SDcols = paste0("Res",random_pick)] 
    setnames(data,"Res")
    
    
    stfdf = STFDF(sp, time, data)
    
    
    #### Calculate the empirical semi-variogram
    #### ----------------------------------------------
    
    ## calculate the distance matrix [km], full symmetric matrix
    Dist <- spDists(sp, longlat = TRUE)
    
    ## set the intervals
    up_Dist <- Dist[upper.tri(Dist, diag = FALSE)] 
    sort_up <- sort(up_Dist)
    
    bound_id = seq(1,length(up_Dist),length.out = nintv + 1)
    boundaries <- sort_up[bound_id]
    
    
    
    empVgm <- variogramST( Res ~ 1, stfdf, tlags=0, boundaries = boundaries,assumeRegular = TRUE, na.omit = TRUE) #"gstat"
  
    #### Fit to the Exponential semi-variogram function
    #### -------------------------------------------------- 
  
    
    ## set the cutoff - usually at least the last 10% of the variogram fit look quite bad
    if(truncate) 
    {
      cutoff_ind  = ceiling(0.9*nintv)
      cutoff = empVgm$dist[cutoff_ind]
    }
    
  
    ## prepare empirical variograms for fitting
    spEmpVgm <- empVgm[empVgm$dist<=cutoff,]  
    spEmpVgm <- spEmpVgm[spEmpVgm$timelag==0,]
    sSpEmpVgm <- spEmpVgm[spEmpVgm$np!=0,] 
    spEmpVgm <- sSpEmpVgm[,1:3] 
    class(spEmpVgm) <- c("gstatVariogram", "data.frame")
    spEmpVgm$dir.hor <- 0
    spEmpVgm$dir.ver <- 0
    
    
    ## Exponential semi-variogram function with nugget, see the 1st argument. 
    ## Fixing "psill", fit "nugget" and "range", the 2nd arg.
    ## Using fit.method = 7, the 3rd arg.
    Mod <- fit.variogram(spEmpVgm, vgm("Exp"), fit.sills = c(T,F), fit.method = 7)
    if(saveorgo)
    {
      save(Mod,Dist,file = paste0(save_dir,file_name,mon,".RData"))
    }
  }
}


#' geostationary forecasts 
#' 
#' @description Generates forecasts using a geostationary model with exponential covariance function for post-processing
#'      Creates output samples, see n below, and uses ensemble members plus noise as forecast. 
#' 
#' @param dt the data table, if NULL, it is loaded.
#' @param y,m Integer vectors containing year(s) and month(s).
#' @param n Integer. Size of the desired forecast ensemble. 
#' @param ens_size Integer. Size of the NWP ensemble.
#' @param truncate logical, whether the forecasted temperature is truncated at freezing temperature.
#' @param saveorgo Logical, whether we save or not. 
#' @param save_dir,file_name The directory to save in and the name of the saved file.
#' @param data_dir,var_file_names Where the fitted variograms are stored.
#' 
#' @return data table containing n columns with noise and n columns with forecasts.
#' 
#' @examples \dontrun{ DT = load_combined_wide(data_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO", output_name = "dt_combine_NAO_wide_bc.RData") \cr
#'                     forecast_geostat(dt = DT)}
#' 
#' @author Claudio Heinrich        
#' 
#' @importFrom MASS mvrnorm
#' 
#' @export

forecast_geostat = function(dt = NULL,
                            n=100,
                            y = 2001:2010,
                            m = 1:12,
                            ens_size = 9,
                            truncate = TRUE,
                            saveorgo = TRUE,
                            save_dir = "./Data/PostClim/SFE/Derived/GeoStat/",
                            file_name = "geostat_fc.RData",
                            data_dir = "./Data/PostClim/SFE/Derived/GeoStat/",
                            var_file_names = paste0("variogram_exp_m",m,".RData")){
    
    if(is.null(dt))
    {
      print("load and prepare data")
      dt = load_combined_wide(bias = TRUE)
    }
    
    SD_cols = c("Lon","Lat","grid_id","month","year","YM",
                "SST_hat","SST_bar",paste0("Ens",1:ens_size),"Ens_bar","Bias_Est")
    
    fc = dt[year %in% y & month %in% m ,.SD,.SDcols = SD_cols]
    
    
    for(mon in m){
    
      # --- get variogram ---
      
      load(paste0(data_dir,var_file_names[which(mon == m)]))
      
      psill <- Mod$psill[2]
      range <- Mod$range[2]
      nugget <- Mod$psill[1]
      
      Sigma <- psill*exp(-Dist/range)
      sills <- diag(Sigma) + nugget
      diag(Sigma) <- sills
      
      ns <- length(sp)
      
      print(paste0("generating noise for month ",mon))
      no <- MASS::mvrnorm(n=n*length(y), mu=rep(0,ns), Sigma=Sigma)
      for (year in y){
        for (i in 1:n){
            y_ind = which(year == y)
            yn_ind = (y_ind - 1)*n + i
            fc[year == year & month == mon,  paste0("no",i):= no[yn_ind,]]
            fc[year == year & month == mon & is.na(SST_hat) | is.na(SST_bar),paste0("no",i) := NA]
            ens_mem = sample.int(ens_size,1) 
            fc[year == year & month == mon,paste0("fc",i):= rowSums(.SD), .SDcols = c(paste0("Ens",ens_mem),"Bias_Est",paste0("no",i))]
        }
      }
      
      if(truncate){
        for(i in 1:n){
          fc[, paste0("fc",i) := lapply(.SD,trc), .SDcols = paste0("fc",i) ]
        }
      }  
      
      if(saveorgo){
        save(fc,file = paste0(save_dir,file_name))
      }
    
    }
return(fc)
}
