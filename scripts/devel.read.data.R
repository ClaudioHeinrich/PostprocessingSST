rm(list = ls())

##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)
##------------------------

##*** Extract Ensemble Forecasts *************
##---- Collect Grid Info ------------------
gridname <- "./Data/grid.nc"
ncgrid <- nc_open(gridname)
grid_lon_ens <- ncvar_get(ncgrid, "plon")
grid_lat_ens <- ncvar_get(ncgrid, "plat")
##-----------------------------------------

sst_ensemble = array(NA,dim = c(dim(grid_lon_ens)[1], dim(grid_lat_ens)[2],9))
for(j in 1:9)
{
  ##----- Location, location --------
  filedir <- paste0(data.dir,"/SFE/NorCPM_Ocean/")
  filename <- paste0("ana_19800115_me_20080715_mem0",j,".micom.hm.2009-06.nc") ## 9 members
  ncname <- paste0(filedir, filename)
  ##--------------------------------

  ##----- Load File ---------
  ncin <- nc_open(ncname)
  ##------------------------
  
  ##------ Extract Covariates ---
  sst_ensemble[,,j] <- ncvar_get(ncin, "sst")
  ##----------------------------------
}
##*** End Extracting Ensemble Memebers *******


##******** Get Observations ******************
##------ Load Observations -------------------
obsdir = "./Data/HadiSST2/"
obsname = "SST_ens_2009_06.nc"
ncnameobs = paste0(obsdir, obsname)
ncobs = nc_open(ncnameobs)
##--------------------------------------------

##------- Extract Observations ---------------
sst_obs = ncvar_get(ncobs, "sst")
lon_obs = ncvar_get(ncobs, "longitude")
lat_obs = ncvar_get(ncobs, "latitude")
##-------------------------------------------
##***** End Observations *********************

##******* Scale sst_ensemble to grid of obs_sst *********
sst_ensemble_rescaled = array(NA,dim = c(dim(sst_obs)[1:2], dim(sst_ensemble)[3]))
for(i in 1:length(lon_obs))
{
  print(i)
  ##------- Get Lon Bounds ---------
  if(i == 1){
    lon_lower = -Inf
  }else{
    lon_lower = (lon_obs[i] + lon_obs[i - 1]) / 2
  }
  if(i == length(lon_obs))
  {
    lon_upper = Inf
  }else{
    lon_upper = (lon_obs[i + 1] + lon_obs[i]) / 2
  }
  ##--------------------------------
  
  for(j in 1:length(lat_obs))
  {

    ##------ Get Lat Bounds ----------
    if(j == 1)
    {
      lat_upper = Inf
    }else{
      lat_upper = (lat_obs[j - 1] + lat_obs[j]) / 2
    }
    if(j == length(lat_obs))
    {
      lat_lower == -Inf
    }else{
      lat_lower = (lat_obs[j] + lat_obs[j + 1])/2
    }
    ##---------------------------------

    ##----- Determine Indicies --------
    w_lon = which( (lon_lower <= grid_lon_ens) & (lon_upper > grid_lon_ens) )
    w_lat = which( (lat_lower <= grid_lat_ens) & (lat_upper > grid_lat_ens) )
    w_use = intersect(w_lon,w_lat)
    ##---------------------------------

    ##------ Load In Ensemble Values ------
    if(length(w_use) == 0)
    {
      sst_ensemble[i,j,] = NA
    }else{
      for(k in 1:dim(sst_ensemble_rescaled)[3])
      {
        m = sst_ensemble[,,k][w_use]
        if(all(is.na(m)))
        {
          sst_ensemble_rescaled[i,j,k] = NA
        }else{
          sst_ensemble_rescaled[i,j,k] = mean(m, na.rm=TRUE)
        }
      }
    }
    ##--------------------------------------
  }
}
## *********** End SST Rescaling *****************************
