make_grid_mapping = function(out.file = "./Data/PostClim/SFE/Derived")
{
  dt_map = construct_grid_map()
  save(dt_map, file = paste0(out.file,"dt_map.RData"))
}

make_combined_dataset = function(y_start = 1985,
                                 y_stop = 2010,
                                 vintage = "mr",
                                 data.dir = "~/PostClimDataNoBackup/",
                                 grid_mapping_loc = "./Data/PostClim/SFE/Derived/")
{

  ##----- Load Grid Mapping ---
  ff = paste0(grid_mapping_loc,"dt_map.RData")
  if(file.exists(ff))
  {
    load(ff)
    names(dt_map)= c("Lon_Obs","Lat_Obs","Lon_Ens","Lat_Ens") ## Get rid of this eventually.
  }else{
    stop("Could not find grid mapping info")
  }
  ##--------------------------

 ##------ Loop ----------
  dt_combine_all = list()
  k = 1
  for(y in y_start:y_stop)
  {
    for(m in 1:12)
    {
      print(c(y,m))
      dt_ens = load_ensemble(y,m,vintage)
      dt_obs = load_observations(y,m)
      dt_combine_all[[k]] = combine_data(dt_ens, dt_obs, dt_map)
      dt_combine_all[[k]][,year:=y]
      dt_combine_all[[k]][,month:=m]
      k = k + 1
    }
  }
  ##------------------------

  ##--------- Combine -----
  dt = rbindlist(dt_combine_all)
  dt[, YM := year * 12 + month]
  setkey(dt, "YM", "Lon", "Lat")
  ##------------------------
  
  ##----- Should I save or should I go? -----
  if(is.null(data.dir))
  {
    return(dt)
  }else{
    save(dt,
        file = paste0(data.dir,"./SFE/Derived/dt_combine_",vintage,".RData"))
    return(1)
  }
  ##-------------------------------------------

  
}

make_combined_wide_dataset = function(y_start = 1985,
                                      y_stop = 2010,
                                      vintage = "mr",
                                      data.dir = "~/PostClimDataNoBackup/",
                                      grid_mapping_loc = "./Data/PostClim/SFE/Derived/")
{

  ##----- Load Grid Mapping ---
  ff = paste0(grid_mapping_loc,"dt_map.RData")
  if(file.exists(ff))
  {
    load(ff)
    names(dt_map)= c("Lon_Obs","Lat_Obs","Lon_Ens","Lat_Ens") ## Get rid of this eventually.
  }else{
    stop("Could not find grid mapping info")
  }
  ##--------------------------

 ##------ Loop ----------
  dt_combine_all = list()
  k = 1
  for(y in y_start:y_stop)
  {
    for(m in 1:12)
    {
      print(c(y,m))
      dt_ens = load_ensemble(y,m,vintage)
      dt_obs = load_observations(y,m)
      dt_combine_all[[k]] = combine_data_wide(dt_ens, dt_obs, dt_map)
      dt_combine_all[[k]][,year:=y]
      dt_combine_all[[k]][,month:=m]
      k = k + 1
    }
  }
  ##------------------------

  ##--------- Combine -----
  dt = rbindlist(dt_combine_all)
  dt[, YM := year * 12 + month]
  setkey(dt, "YM", "Lon", "Lat")
  ##------------------------
  
  ##----- Should I save or should I go? -----
  if(is.null(data.dir))
  {
    return(dt)
  }else{
    save(dt,
        file = paste0(data.dir,"./SFE/Derived/dt_combine_",vintage,"_wide.RData"))
    return(1)
  }
  ##-------------------------------------------

  
}


#---- for bias correction with exponential moving averages ----

exp_mov_av_bc = function (vec,ratio){
  a = c(0,vec[1])   
  for(k in 3:length(vec)){
    temp=ratio * vec[k-1] + (1-ratio)*a[length(a)]
    a = c(a,temp)
  }
  return(a)
}




bias_correct = function(dt = NULL,
                        Y = 1985:2010,
                        M = 1:12,
                        method = "ema",  # "gwa" for 'gliding window average', i.e. simple moving average, 
                                         # "ema" for 'exponential moving average'
                        par_1 = .2,   # if method = gwa, par_1 is the length of gliding window
                                      # if method = ema, par_1 is the ratio of the exp. mov. av.
                        global_mean_scores = FALSE,
                        reduced_output = FALSE,
                        saveorgo = FALSE,
                        data.dir = "~/PostClimDataNoBackup/"
){
  
  if(is.null(dt)) dt = load_combined_wide()
  
  dt = dt[year %in% Y & month %in% M,]
  dt_new = copy(dt)
  
  if(method == "gwa"){
    dt_new = dt_new[!is.na(Ens_bar) & !is.na(SST_bar) ,
                    "Bias_Est" := par_1/(par_1-1) * SMA(SST_bar - Ens_bar,n = par_1) - (SST_bar - Ens_bar)/(par_1-1),
                    by = .(grid_id, month)]
    
    dt_new = dt_new[year < min(Y) + par_1, 
                    Bias_Est := (cumsum(SST_bar-Ens_bar) - (SST_bar-Ens_bar)) / (year - min(year)+1),
                    by = .(grid_id, month)]
  }
  if (method == "ema"){
    dt_new = dt_new[!is.na(Ens_bar) & !is.na(SST_bar) ,
                    "Bias_Est" := exp_mov_av_bc(SST_bar - Ens_bar, ratio = par_1),
                    by = .(grid_id, month)]
  }
  
  
  dt_new[,"SST_hat" := Ens_bar + Bias_Est]
  
  if(saveorgo){dt = dt_new
  save(dt, file = paste0(data.dir,"SFE/Derived/dt_combine_wide_bias.RData"))
  }
  
  if(global_mean_scores){
    glob_mean_sc = dt_new[,.( "RMSE" = 
                                sqrt(mean( (SST_bar - SST_hat)^2, 
                                           na.rm=TRUE)),
                              "MAE" = mean(abs(SST_bar - SST_hat),
                                           na.rm=TRUE)),
                          keyby = YM]
    return(glob_mean_sc)
  }else if(reduced_output)
  {result = dt_new[,.(grid_id,year,month,YM,SST_hat)]
  return(result)
  } else return(dt_new)
  
}


                                 
