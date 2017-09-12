make_combined_dataset = function(y_start = 1985,
                                 y_stop = 2010, vintage = "mr",
                                 data.dir = "~/PostClimDataNoBackup/")
{ ##------ Loop ----------
  dt_combine_all = list()
  k = 1
  for(y in y_start:y_stop)
  {
    for(m in 1:12)
    {
      print(c(y,m))
      dt_ens = load_ensemble(y,m,vintage)
      dt_obs = load_observations(y,m)
      dt_combine_all[[k]] = combine_data(dt_ens, dt_obs)
      dt_combine_all[[k]][,year:=y]
      dt_combine_all[[k]][,month:=m]
      k = k + 1
    }
  }
  ##------------------------

  ##------Create data.table-----
   dt = rbindlist(dt_combine_all)
   dt[, YM := year * 12 + month]
   setkey(dt, "YM", "Lon", "Lat")
  # ##------------------------
  
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

                                 
