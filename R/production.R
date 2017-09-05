make_combined_dataset = function(y_start = 1985,
                                 y_stop = 2010,
                                 data.dir = "~/PostClimDataNoBackup/")
{

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
      dt_combine_all[[k]] = combine_data(dt_ens, dt_obs)
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

                                 
