#' A utility function that makes and stores the grid mapping
#'
#' @description This is essentially a wrapper for \code{construct_grid_map}
#'
#' @param out.file Location to store the output, which is then saved at \code{/dt_map.RData}
#'
#' @export
make_grid_mapping = function(out.file = "./Data/PostClim/SFE/Derived")
{
  dt_map = construct_grid_map()
  save(dt_map, file = paste0(out.file,"dt_map.RData"))
}


#' Make the combined dataset for all years and store in a wide format
#'
#' @param y_start The first year for which we have obs/ens data
#' @param y_stop The final year for which we have obs/ens data
#' @param vintage Which vintage are we using, often "mr" for most recent
#' @param data.dir The root directory that stores our data
#' @param grid_mapping_loc The location where the precomputed grid alignment object is stored.
#'
#' @export
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

#' @export
make_NorCPM_wide_again_precip_2mtemp = function(save.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                                                data.names = c("dt_combine_mr_wide_new.RData","dt_prect_NorCPM_wide.RData","dt_2mtemp_NorCPM_wide.RData")){
  
  data = fread("~/PostClimDataNoBackup/SFE/NorCPM_2mTemp/wide_data")
  
  old_names = paste0("wide_data_2.",c("ym",
                                      "lon",
                                      "lat",
                                      paste0("sst",1:10),
                                      "sst_bar",
                                      "sst_sd",
                                      paste0("ens",1:9),
                                      "ens_bar",
                                      "ens_sd",
                                      paste0("ts",1:9),
                                      "ts_bar",
                                      "ts_sd",
                                      paste0("prect",1:9),
                                      "prect_bar",
                                      "prect_sd",
                                      "sst_noaa",
                                      "year",
                                      "month"))
  
  new_names = c("YM",
                "Lon",
                "Lat",
                paste0("SST",1:10),
                "SST_bar",
                "SST_sd",
                paste0("Ens",1:9),
                "Ens_bar",
                "Ens_sd",
                paste0("ts",1:9),
                "ts_bar",
                "ts_sd",
                paste0("prect",1:9),
                "prect_bar",
                "prect_sd",
                "sst_noaa",
                "year",
                "month")
  
  
  setnames(data,old_names,new_names)
  
  # --- change "NULL" into NA and convert strings into doubles ---
  
  is.na(data) <- data == "NULL"
  
  num_names = c(paste0("SST",1:10),
                "SST_bar",
                "SST_sd",
                paste0("Ens",1:9),
                "Ens_bar",
                "Ens_sd",
                paste0("ts",1:9),
                "ts_bar",
                "ts_sd",
                paste0("prect",1:9),
                "prect_bar",
                "prect_sd",
                "sst_noaa")
  
  data[,(num_names) := lapply(.SD,as.numeric),.SDcols = num_names]
  
  ##------ include grid_id -------
  lon_all = sort(data[,unique(Lon)])
  n_lon = length(lon_all)
  cutoff_lon = c(-Inf,head(lon_all,-1) + diff(lon_all)/2)
  f_lon = approxfun(cutoff_lon, 1:n_lon, method="constant", rule = 2)
  
  lat_all = sort(data[,unique(Lat)])
  n_lat = length(lat_all)
  cutoff_lat = c(-Inf, head(lat_all,-1) + diff(lat_all)/2)
  f_lat = approxfun(cutoff_lat, 0:(n_lat - 1) * n_lon, method="constant", rule = 2)
  data[,grid_id := f_lon(Lon) + f_lat(Lat)]
  ##--------------------------------------
  
  #----break into smaller chunks and save
  
  dt = data[,.SD,.SDcols = c("Lon","Lat", paste0("SST",1:10), "SST_bar","SST_sd", paste0("Ens",1:9),"Ens_bar","Ens_sd","grid_id","year","month","YM","sst_noaa")]
  dt[,,.SDcols =  c(paste0("SST",1:10), "SST_bar","SST_sd","Ens_bar","Ens_sd")]
  save(dt, file = paste0(save.dir,data.names[1]))
  
  dt_prect = data[,.SD,.SDcols = c("Lon","Lat",paste0("prect",1:9),"prect_bar","prect_sd","grid_id","year","month","YM")]
  save(dt_prect, file = paste0(save.dir,data.names[2]))
  
  dt_2mtemp = data[,.SD,.SDcols = c("Lon","Lat",paste0("ts",1:9),"ts_bar","ts_sd","grid_id","year","month","YM")]
  save(dt_2mtemp, file = paste0(save.dir,data.names[3]))
  
}

