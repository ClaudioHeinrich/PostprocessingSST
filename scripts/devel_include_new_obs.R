
obsdir = "./Data/HadiSST2/"
if(month < 10)month = paste0("0",month)
obsname = paste0("HadISST_sst.nc")
ncnameobs = paste0(obsdir, obsname)
ncobs = nc_open(ncnameobs)


##------- Extract Observations ---------------
sst_obs = ncvar_get(ncobs, "sst")

sst_date = as.Date(ncvar_get(ncobs,"time"),origin = "1870-1-1")

lon_obs = ncvar_get(ncobs, "longitude")
lat_obs = ncvar_get(ncobs, "latitude")
n_lon = length(lon_obs)
n_lat = length(lat_obs)
grid_lon_obs = matrix(rep(lon_obs,length(lat_obs)),n_lon,n_lat)
grid_lat_obs = matrix(rep(lat_obs,length(lon_obs)),n_lon,n_lat, byrow=TRUE)

YM = 12*as.numeric(format(sst_date,"%Y"))+as.numeric(format(sst_date,"%m"))

grid_YM = t(replicate(n_lon*n_lat,YM))

  dt_obs = data.table(Lon = as.vector(grid_lon_obs),
                                Lat = as.vector(grid_lat_obs),
                                YM = as.vector(grid_YM),
                                SST = sst_obs)
setkey(dt_obs, "YM","Lon", "Lat")

dt_obs = dt_obs[SST == -1000, SST := NA]

dt = load_combined_wide()
setkey(dt, "YM","Lon", "Lat")

fill_in_YM_min_max = c(2011*12+1,2017*12+10)

SST_new = dt_obs[YM >= fill_in_YM_min_max[1],SST]

dt[YM >= fill_in_YM_min_max[1] & YM <= fill_in_YM_min_max[2], SST_bar := SST_new]

mean(abs(dt[YM %in% 23821:24000,SST_bar] - dt_obs[YM %in% 23821:24000,SST]),na.rm = TRUE)
