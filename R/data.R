contruct_grid_map = function(dt_ens = load_ensemble(1985,1),
                             dt_obs = load_observations(1985,1))
{

  dt_obs_grid = unique(dt_obs[,.(Lon,Lat)])
  point_match = NULL
  for(j in 1:dim(dt_obs_grid)[1])
  {
    if(j %% 1e2 == 0)print(j)
    ##---------------------------------------
    a = geosphere::distHaversine(as.vector(dt_obs_grid[j,.(Lon,Lat)]),as.matrix(dt_ens_grid[,.(Lon,Lat)]))
    point_match[j] = order(a)[1]
    ##--------------------------------------
  }

  dt_map = data.table(dt_obs_grid[,.(Lon,Lat)], dt_ens_grid[point_match,.(Lon,Lat)])

  return(dt_map)
}

load_ensemble = function(year,
                         month,
                         data.dir = "~/PostClimDataNoBackup")
{

  ##------ Setup ----------
  if(month < 10)month = paste0("0",month)
  ##-----------------------

  ##---- Collect Grid Info ------------------
  gridname <- "./Data/grid.nc"
  ncgrid <- nc_open(gridname)
  grid_lon_ens <- ncvar_get(ncgrid, "plon")
  grid_lat_ens <- ncvar_get(ncgrid, "plat")
  nc_close(ncgrid)
  ##-----------------------------------------

  ##-------- Find target run ----------------
  filedir <- paste0(data.dir,"/SFE/NorCPM_Ocean/")
  ff_all = system(paste0("ls ",filedir,"*_mem01.*.",year,"-",month,"*"), intern = TRUE)
  ff_use = tail(ff_all,1)
  ##-----------------------------------------
  
  ##----- Get Raw Ensemble ------------------
  sst_ensemble = array(NA,dim = c(dim(grid_lon_ens)[1], dim(grid_lat_ens)[2],9))
  for(j in 1:9)
  {
    ff_new = gsub("mem01",paste0("mem0",j),ff_use)
    ncin <- nc_open(ff_new)
    sst_ensemble[,,j] <- ncvar_get(ncin, "sst")
    nc_close(ncin)
  }
  ##-----------------------------------

  ## --- Make a Data Table ----------------
  dt_ensemble_list = list()
  for(j in 1:9)
  {
    dt_ensemble_list[[j]] = data.table(Lon = as.vector(grid_lon_ens),
                                       Lat = as.vector(grid_lat_ens),
                                       Ens = j,
                                       Forecast = as.vector(sst_ensemble[,,j]))
  }
  dt_ensemble = rbindlist(dt_ensemble_list)
  setkey(dt_ensemble,"Lon","Lat","Ens")
  ##-----------------------------------------

  return(dt_ensemble)
}

load_observations = function(year, month,
                             obsdir = "./Data/HadiSST2/")
{
  ##------ Load Observations -------------------
  obsdir = "./Data/HadiSST2/"
  if(month < 10)month = paste0("0",month)
  obsname = paste0("SST_ens_", year,"_", month,".nc")
  ncnameobs = paste0(obsdir, obsname)
  ncobs = nc_open(ncnameobs)
  ##--------------------------------------------
  
  ##------- Extract Observations ---------------
  sst_obs = ncvar_get(ncobs, "sst")
  sst_obs = sst_obs - 273 ## Convert from Kelvin
  lon_obs = ncvar_get(ncobs, "longitude")
  lat_obs = ncvar_get(ncobs, "latitude")
  n_lon = length(lon_obs)
  n_lat = length(lat_obs)
  grid_lon_obs = matrix(rep(lon_obs,length(lat_obs)),n_lon,n_lat)
  grid_lat_obs = matrix(rep(lat_obs,length(lon_obs)),n_lon,n_lat, byrow=TRUE)
  dt_obs_list = list()
  for(j in 1:dim(sst_obs)[3])
  {
    dt_obs_list[[j]] = data.table(Lon = as.vector(grid_lon_obs),
                                  Lat = as.vector(grid_lat_obs),
                                  Obs = j,
                                  SST = as.vector(sst_obs[,,j]))
  }
  dt_obs = rbindlist(dt_obs_list)
  setkey(dt_obs, "Lon", "Lat")
  ##-------------------------------------------

  nc_close(ncobs)
  
  return(dt_obs)
}

combine_data = function(dt_ens, dt_obs, dt_map)
{

  ##------- Collapse Observations -------
  dt_obs_mean = dt_obs[,.("SST_bar" = mean(SST)),.(Lon,Lat)]
  setkey(dt_obs_mean, "Lon","Lat")
  ##-------------------------------------
  
  ##------- First fill out ensemble with obs lookup ---
  setkey(dt_map, "Lon_Ens","Lat_Ens")
  setkey(dt_ens, "Lon", "Lat")
  dt_ens = dt_ens[dt_map]
  ##-----------------------------------------------
  
  ##------- Now fill in Observations --------
  setkey(dt_ens, "Lon_Obs", "Lat_Obs")
  dt_combine = dt_obs_mean[dt_ens]
  dt_combine[,i.Lon := NULL]
  dt_combine[,i.Lat := NULL]
  ##-----------------------------------------
  
  ##------ Form a convenient key -------
  lon_all = sort(dt_combine[,unique(Lon)])
  n_lon = length(lon_all)
  cutoff_lon = c(-Inf,head(lon_all,-1) + diff(lon_all)/2)
  f_lon = approxfun(cutoff_lon, 1:n_lon, method="constant", rule = 2)
  
  lat_all = sort(dt_combine[,unique(Lat)])
  n_lat = length(lat_all)
  cutoff_lat = c(-Inf, head(lat_all,-1) + diff(lat_all)/2)
  f_lat = approxfun(cutoff_lat, 0:(n_lat - 1) * n_lon, method="constant", rule = 2)
  dt_combine[,grid_id := f_lon(Lon) + f_lat(Lat)]
  ##--------------------------------------
  
  return(dt_combine)
  
}


load_combined = function(data.dir="~/PostClimDataNoBackup/")
{
  file = paste0(data.dir,"/SFE/Derived/dt_combine.RData")
  load(file)
  return(dt)
}
