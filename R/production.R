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
                                      grid_mapping_loc = "./Data/PostClim/SFE/Derived/",
                                      include_germans = TRUE)
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
                        model = "NorESM",
                        method = "ema",  # "gwa" for 'gliding window average', i.e. simple moving average, 
                                         # "ema" for 'exponential moving average'
                        par_1 = .2,   # if method = gwa, par_1 is the length of gliding window
                                      # if method = ema, par_1 is the ratio of the exp. mov. av.
                        global_mean_scores = FALSE,
                        reduced_output = FALSE,
                        saveorgo = FALSE,
                        data.dir = "~/PostClimDataNoBackup/"
){
  
  if(model == "NorESM"){
    if(is.null(dt)) dt = load_combined_wide() 
    
  
  Y = dt[,min(year)]:dt[,max(year)]
  M = 1:12
  
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
  
  if(saveorgo & model == "NorESM"){dt = dt_new
  save(dt, file = paste0(data.dir,"SFE/Derived/dt_combine_wide_bias.RData"))
  }
  if(saveorgo & model == "senorge"){dt_senorge = dt_new
  save(dt_senorge, file = paste0(data.dir,"SFE/Derived/senorge2_gfsc1_combined_bias.RData"))
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
  
  
  if(model == "senorge"){
    
    if(is.null(dt)) {
      dt = load_combined_wide(model = "senorge") 
      setnames(dt,"senorge_grid_id","grid_id")
    }
    
    Y = dt[,min(year)]:dt[,max(year)]
    M = 1:12
    
    dt_new = copy(dt)
    
    if(method == "gwa"){
      dt_new = dt_new[!is.na(Ens_bar) & !is.na(temp) ,
                      "Bias_Est" := par_1/(par_1-1) * SMA(temp - Ens_bar,n = par_1) - (temp - Ens_bar)/(par_1-1),
                      by = .(grid_id, month)]
      
      dt_new = dt_new[year < min(Y) + par_1, 
                      Bias_Est := (cumsum(temp-Ens_bar) - (temp-Ens_bar)) / (year - min(year)+1),
                      by = .(grid_id, month)]
    }
    if (method == "ema"){
      dt_new = dt_new[!is.na(Ens_bar) & !is.na(temp) ,
                      "Bias_Est" := exp_mov_av_bc(temp - Ens_bar, ratio = par_1),
                      by = .(grid_id, month)]
    }
    
    
    dt_new[,"temp_hat" := Ens_bar + Bias_Est]
    
    if(saveorgo & model == "NorESM"){dt = dt_new
    save(dt, file = paste0(data.dir,"SFE/Derived/dt_combine_wide_bias.RData"))
    }
    if(saveorgo & model == "senorge"){dt_senorge = dt_new
    save(dt_senorge, file = paste0(data.dir,"SFE/Derived/senorge2_gfsc1_combined_bias.RData"))
    }
    
    
    
    if(global_mean_scores){
      glob_mean_sc = dt_new[,.( "RMSE" = 
                                  sqrt(mean( (temp - temp_hat)^2, 
                                             na.rm=TRUE)),
                                "MAE" = mean(abs(temp - temp_hat),
                                             na.rm=TRUE)),
                            keyby = YM]
      return(glob_mean_sc)
    }else if(reduced_output)
    {result = dt_new[,.(grid_id,year,month,YM,temp_hat)]
    return(result)
    } else return(dt_new)
    
    
    
    }
  
}


make_GCFS1_wide = function(in.dir = "~/PostClimDataNoBackup/SFE/GCFS1",
                           out.dir = "~PostClimDataNoBackup/SFE/Derived/",
                           N_ens = 15,
                           verbose = TRUE)
{

  temp2_ens_all = list()

  ##------ Loop through Ensemble members ----
  for(j in 1:N_ens)
  {
    if(verbose) {print(paste0("Loading ",j," of ",N_ens))}

    ##---- Stupid ---
    if(j < 10)
    {
      ss = paste0("0",j)
    }else{
      ss = j
    }
    ##---------------

    ##----- Extact Temp ------
    ncens = nc_open(paste0(in.dir,"/Mem", ss, "_GCFS1_temp2_mm_1981-2015_smon11.nc"))
    temp2_ens_all[[j]] =  ncvar_get(ncens, "temp2") - 273
    ##------------------------

    ##---- Bookkeeping -------  
    if(j == 1)
    {
      grid_lon_ens = ncvar_get(ncens, "lon")
      grid_lon_ens = grid_lon_ens - 180.0
      grid_lat_ens = ncvar_get(ncens, "lat")
      tt = ncvar_get(ncens, "time")
      tt_s = as.Date(tt, origin = "1978-10-31 23:50:00")
      YM_all = data.table(year = as.numeric(format(tt_s,"%Y")),
                          month = as.numeric(format(tt_s,"%m")))
      N_lon = length(grid_lon_ens)
      N_lat = length(grid_lat_ens)
    }
    ##------------------------
  }
  ##---- End Loop through ensemble - ---------

  ##------ Now construct wide tables -------
  dt_ens_all = list()
  for(t in 1:YM_all[,.N])
  {
    ##----- Dates --------
    y = YM_all[t,year]
    m = YM_all[t,month]
    YM = y * 12 + m
    ##--------------------
    
    ##---- Init DT -------
    dt_ens_all[[t]] = data.table(YM = YM, year = y, month = m,
                                 Lon = rep(grid_lon_ens, times = N_lat),
                                 Lat = rep(grid_lat_ens, each = N_lon),
                                 GCFS1_id = 1:(N_lat * N_lon))
    ##--------------------
    
    ##---- Now add Ensemble info ------
    for(j in 1:N_ens)
    {
      dt_ens_all[[t]][,paste0("Ens_",j) := as.vector(temp2_ens_all[[j]][,,t])]
    }
    ##----------------------------------
  }
  ##------ Finish writing wide table ------------------

  ##----- Write ---------------
  dt_ens = rbindlist(dt_ens_all)
  dt_ens[ , Ens_bar := rowMeans(.SD),.SDcols = paste0("Ens_",1:15)]
  save(dt_ens, file = paste0(out.dir, "/GCFS1_wide.RData"))
  ##----------------------------

}

make_GCFS1_wide_sst = function(in.dir = "~/PostClimDataNoBackup/SFE/GCFS1",
                           out.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                           N_ens = 15,
                           verbose = TRUE)

{

  variable = "sst"
  gcfs1_ens_all = list()

  ##------ Loop through Ensemble members ----
  for(j in 1:N_ens)
  {
    if(verbose) {print(paste0("Loading ",j," of ",N_ens))}

    ##---- Stupid ---
    if(j < 10)
    {
      ss = paste0("0",j)
    }else{
      ss = j
    }
    ##---------------

    ##----- Extact Temp ------
    ncens = nc_open(paste0(in.dir,"/Mem", ss, "_GCFS1_",variable,"_mm_1981-2015_smon11.nc"))
    gcfs1_ens_all[[j]] =  ncvar_get(ncens, variable) - 273
    ##------------------------

    ##---- Bookkeeping -------  
    if(j == 1)
    {
      grid_lon_ens = ncvar_get(ncens, "lon")
      grid_lon_ens = grid_lon_ens
      grid_lat_ens = ncvar_get(ncens, "lat")
      tt = ncvar_get(ncens, "time")
      tt_s = as.Date(tt, origin = "1981-11-30 22:48:00")
      YM_all = data.table(year = as.numeric(format(tt_s,"%Y")),
                          month = as.numeric(format(tt_s,"%m")))
      N_lon = length(grid_lon_ens)
      N_lat = length(grid_lat_ens)
    }
    ##------------------------
  }
  ##---- End Loop through ensemble - ---------

  ##------ Now construct wide tables -------
  dt_ens_all = list()
  for(t in 1:YM_all[,.N])
  {
    ##----- Dates --------
    y = YM_all[t,year]
    m = YM_all[t,month]
    YM = y * 12 + m
    ##--------------------
    
    ##---- Init DT -------
    dt_ens_all[[t]] = data.table(YM = YM, year = y, month = m,
                                 Lon = as.vector(grid_lon_ens),
                                 Lat = as.vector(grid_lat_ens),
                                 GCFS1_id = 1:N_lon)
    ##--------------------
    
    ##---- Now add Ensemble info ------
    for(j in 1:N_ens)
    {
      dt_ens_all[[t]][,paste0("GCSF1_",variable,"_",j) := as.vector(gcfs1_ens_all[[j]][,,t])]
    }
    ##----------------------------------
  }
  ##------ Finish writing wide table ------------------

  ##----- Write ---------------
  dt_ens = rbindlist(dt_ens_all)
  dt_ens[ , Ens_bar := rowMeans(.SD),.SDcols = paste0("GCSF1_",variable,"_",1:15)]
  save(dt_ens, file = paste0(out.dir, "/GCFS1_",variable,"_wide.RData"))
  ##----------------------------

}
                                 

make_senorge_data = function(in.dir = "~/PostClimDataNoBackup/seNorge/",
                             in.file = "seNorge2_TEMP1d_dataset_1991_2015.nc",
                             out.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                             out.file = "dt_senorge_1991_2015.RData",
                             lon_bound = c(4.5,7.8),
                             lat_bound = c(57.95, 61.04),
                             verbose = TRUE)
{
  ##------ Load Observations -------------------
  ff = paste0(in.dir,in.file)
  ncobs = nc_open(ff)
  ##--------------------------------------------

  ##------- Get high-level info -------
  x_coord = ncvar_get(ncobs,"X")
  y_coord = ncvar_get(ncobs,"Y")
  days = ncvar_get(ncobs,"time")
  DD = as.Date(days,origin = "1900-01-01")
  year = as.numeric(format(DD,"%Y"))
  month = as.numeric(format(DD,"%m"))
  YM = year * 12 + month
  YM_all = unique(YM)
  ##------------------------------------

  ##------ Organize grid ---------
  N_grid = length(x_coord) * length(y_coord)
  Coord_long = cbind(rep(x_coord,length = N_grid),rep(y_coord,each = length(x_coord)))
  SP <- SpatialPoints(Coord_long, proj4string=CRS("+proj=utm +zone=33") )
  SP_lon_lat = spTransform(SP, CRS("+proj=longlat"))
  ##-------------------------------

  ##------ Now construct sub grid-level YM averages---
  dt_list = list()
  for(i in 1:length(YM_all))
  {
    print(paste0("On ",i," of ",length(YM_all)))
    w_in = which(YM == YM_all[i])
    dt_list[[i]] = data.table(year = year[w_in[1]],
                              month = month[w_in[1]],
                              YM= YM_all[i],
                              Lon = SP_lon_lat$coords.x1,
                              Lat = SP_lon_lat$coords.x2,
                              senorge_grid_id = 1:length(SP_lon_lat$coords.x1),
                              temp=0.0)
    n_in = length(w_in)
    for(j in 1:length(w_in))
    {
      temp_j = as.vector(ncvar_get(ncobs, "mean_temperature",start = c(1,1,w_in[j]), count = c(-1,-1,1)))
      dt_list[[i]][,temp := temp + temp_j/n_in]
    }
    
    dt_list[[i]] = dt_list[[i]][ (Lon > lon_bound[1]) &
                                 (Lon < lon_bound[2]) &
                                 (Lat > lat_bound[1]) &
                                 (Lat < lat_bound[2]) ]
  }
  ##-------------------------------------------------------

  ##------ Organize grid ---------
  N_grid = length(x_coord) * length(y_coord)
  Coord_long = cbind(rep(x_coord,length = N_grid),rep(y_coord,each = length(x_coord)))
  SP <- SpatialPoints(Coord_long, proj4string=CRS("+proj=utm +zone=33") )
  SP_lon_lat = spTransform(SP, CRS("+proj=longlat"))
  ##-------------------------------
  
  dt_senorge = rbindlist(dt_list)
  
  save(dt_senorge,file = paste0(out.dir, out.file))

}
                             

construct_NorCPM_GCSF1_map = function(in.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                                      out.dir = "~/PostClimDataNoBackup/SFE/Derived/")
{
  
  ##----- Load Observation --------
  dt_obs = load_observations(1985,1)
  dt_obs_grid = unique(dt_obs[,.(Lon,Lat)])
  dt_obs_grid[,"Obs_grid_id":=1:.N]
  ##-------------------------------

  ##----- Now get Ensemble together ---
  load(paste0(in.dir,"GCFS1_sst_wide.RData"))
  dt_ens_grid = unique(dt_ens[YM == head(YM,1),.(Lon,Lat,GCFS1_id)])
  rm(dt_ens);gc()
  ##------------------------------------
  
  ##---- Line things up -----------
  point_match = unlist(mclapply(1:dim(dt_obs_grid)[1], "closest_point_helper",
                                dt_obs_grid, dt_ens_grid,
                                mc.cores = 20, mc.silent = FALSE))
  ##-------------------------------
  
  ##------ Make final object ---
  dt_map_NorCPM_GCSF1 = data.table(dt_obs_grid[,.(Lon,Lat)], dt_ens_grid[point_match,.(Lon,Lat)])
  names(dt_map_NorCPM_GCSF1) = c("Lon_Obs","Lat_Obs","Lon_GCSF1", "Lat_GCSF1")
  ##-----------------------------

  save(dt_map_NorCPM_GCSF1, file = paste0(out.dir,"NorCPM_GCSF1_map.RData"))

  return(TRUE)
}

make_combined_wide_dataset_add_germans = function(y_start = 1985,
                                                  y_stop = 2010,
                                                  vintage = "mr",
                                                  data_dir = "~/PostClimDataNoBackup/",
                                                  out_dir = "~/PostClimDataNoBackup/SFE/Derived/")
{

  
  dt_wide = load_combined_wide()

  ##----- Load Grid Mapping ---
  ff = paste0(out_dir,"NorCPM_GCSF1_map.RData")
  if(file.exists(ff))
  {
    load(ff)
  }else{
    stop("Could not find grid mapping info")
  }
  ##--------------------------

  ##------ Load the germans ---
  load(paste0(out_dir, "GCFS1_sst_wide.RData"))
  ##---------------------------

  ##----- Set up Germans -------
  setkey(dt_ens, "Lon","Lat")
  setkey(dt_map_NorCPM_GCSF1, "Lon_GCSF1","Lat_GCSF1")
  dt_ens = merge(dt_ens,dt_map_NorCPM_GCSF1,by.x = c("Lon","Lat"), by.y = c("Lon_GCSF1","Lat_GCSF1"), all.x = FALSE,all.y = FALSE, allow.cartesian = TRUE)
 ##------ Loop ----------

  ##---- Now get organized ------
  setkey(dt_ens,YM,Lon_Obs,Lat_Obs)
  setkey(dt_wide,YM, Lon, Lat)
  ##------------------------------

  ##------ MERGE!!! ------
  dt_wide = merge(dt_wide, dt_ens,
                  by.x = c("Lon","Lat","YM"),
                  by.y = c("Lon_Obs","Lat_Obs","YM"),
                  all.x=TRUE,
                  all.y=FALSE)
  ##----------------------

  ##------- WRITE!!! --------
  save(dt_wide,
        file = paste0(data.dir,"/dt_combine_both.RData"))
  ##-------------------------------------------

  
}
