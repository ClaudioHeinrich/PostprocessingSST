make_grid_mapping = function(out.file = "./Data/PostClim/SFE/Derived")
{
  dt_map = construct_grid_map()
  save(dt_map, file = paste0(out.file,"dt_map.RData"))
}

make_combined_dataset = function(y_start = 1985,
                                 y_stop = 2010,
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
      dt_ens = load_ensemble(y,m)
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
         file = paste0(data.dir,"./SFE/Derived/dt_combine.RData"))
    return(1)
  }
  ##-------------------------------------------

  
}

                                 
