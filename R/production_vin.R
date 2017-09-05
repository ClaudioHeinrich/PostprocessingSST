
# vintage accepts "mr" for most recent, and "Jan","Apr","Jul","Oct" which correspond to 
# the prediction that was started in the most recent corresponding month

make_combined_dataset_vin = function(y_start = 1986, y_stop = 2010, vintage = "mr",
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
      dt_ens = load_ensemble_vin(y,m,vintage)
      dt_obs = load_observations(y,m)
      dt_combine_all[[k]] = combine_data(dt_ens, dt_obs)
      dt_combine_all[[k]][,year:=y]
      dt_combine_all[[k]][,month:=m]
      k = k + 1
    }
  }
  ##------------------------
  
  ##----- Should I save or should I go? -----
  if(is.null(data.dir))
  {
    return(dt_combine_all)
  }else{
    save(dt_combine_all,
         file = paste0(data.dir,"./SFE/Derived/dt_combine_", vintage, ".RData"))
    return(1)
  }
  ##-------------------------------------------
  
}

load_combined_vin = function(vintage = "mr",data.dir="~/PostClimDataNoBackup/")
{
  file = paste0(data.dir,"/SFE/Derived/dt_combine_",vintage,".RData")
  load(file)
  return(dt_combine_all)
}
                                 
