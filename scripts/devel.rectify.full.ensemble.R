rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)

y = 1985
m = 1

##------ Loop ----------
dt_ens = load_ensemble(y,m)
dt_obs = load_observations(y,m)
##------------------------

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
##------------------------------------------
  
  
  return(dt_combine)

}

dt_combine_all[[k]][,year:=y]
    dt_combine_all[[k]][,month:=m]
    k = k + 1
  }
}
##------------------------

  ##--------- Combine -----
  dt_combined_all = load_combined()
  dt = rbindlist(dt_combined_all)
  dt[, YM := year * 12 + month]
  setkey(dt, "YM", "Lon", "Lat")
  ##------------------------

