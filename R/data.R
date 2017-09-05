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
  setkey(dt_ensemble,"Lon","Lat")
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

combine_data = function(dt_ens, dt_obs)
{
  ##------- Collapse Observations -------
  dt_obs_mean = dt_obs[,.("SST_bar" = mean(SST)),.(Lon,Lat)]
  setkey(dt_obs_mean, "Lat","Lon")
  ##-------------------------------------
  
  ##------ First Organize the Obs -------
  lon_all = sort(dt_obs_mean[,unique(Lon)])
  n_lon = length(lon_all)
  cutoff_lon = c(-Inf,head(lon_all,-1) + diff(lon_all)/2)
  f_lon = approxfun(cutoff_lon, 1:n_lon, method="constant", rule = 2)
  
  lat_all = sort(dt_obs_mean[,unique(Lat)])
  n_lat = length(lat_all)
  cutoff_lat = c(-Inf, head(lat_all,-1) + diff(lat_all)/2)
  f_lat = approxfun(cutoff_lat, 0:(n_lat - 1) * n_lon, method="constant", rule = 2)
  dt_obs_mean[,grid_id := f_lon(Lon) + f_lat(Lat)]
  ##--------------------------------------
  
  ##------ Now merge the obs and ens data ---
  dt_ens_mean = dt_ens[,.("SST_hat" = mean(Forecast)),.(Lon,Lat)]
  setkey(dt_ens_mean, "Lat","Lon")
  dt_ens_mean[,grid_id := f_lon(Lon) + f_lat(Lat)]
  setkey(dt_ens_mean,grid_id)
  dt_ens_grid_mean = dt_ens_mean[,.("Lon_bar" = mean(Lon),
                                    "Lat_bar" = mean(Lat),
                                    "SST_hat_grid" = mean(SST_hat, na.rm=TRUE)),grid_id]
  dt_combine = dt_obs_mean[dt_ens_grid_mean,on="grid_id"]
  ##------------------------------------------
  
  
  return(dt_combine)

}

load_combined = function(data.dir="~/PostClimDataNoBackup/")
{
  file = paste0(data.dir,"/SFE/Derived/dt_combine.RData")
  load(file)
  return(dt)
}
