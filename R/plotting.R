plot_system = function(dt,
                       YM_j,
                       file_out = paste0("./figures/system_",YM_j))
{

  ##------- Setup --------
  dt_sub = dt[YM == YM_j]
  n_lon = length(dt_sub[,unique(Lon)])
  n_lat = length(dt_sub[,unique(Lat)])
  mn = paste0(dt_sub[,min(month)], "/",dt_sub[,min(year)])
  ##----------------------

  ##----- Observation -------
  A_obs = matrix(NA, n_lon,n_lat)
  A_obs[ dt_sub[grid_id < length(A_obs), grid_id] ] = dt_sub[grid_id < length(A_obs),SST_bar] ## FIX
  if(print_figs){pdf(paste0(file_out, "_obs.pdf"))}else{X11()}
  image(A_obs, main = paste0(mn," observation"))
  if(print_figs)dev.off()
  ##-------------------------

  ##----- Ensemble -----------
  A_ens = matrix(NA, n_lon,n_lat)
  A_ens[ dt_sub[grid_id < length(A_ens), grid_id] ] = dt_sub[grid_id < length(A_ens),SST_hat_grid] ## FIX
  if(print_figs){pdf(paste0(file_out, "_ens.pdf"))}else{X11()}
  image(A_ens,main = paste0(mn, " ensemble"))
  if(print_figs)dev.off()
  ##-------------------------

  ##----- Bias --------------
  A_bias = A_obs - A_ens
  if(print_figs){pdf(paste0(file_out, "_bias.pdf"))}else{X11()}
  image(A_bias,main = paste0(mn," bias"))
  if(print_figs)dev.off()

  if(print_figs)
  {
    png(paste0(file_out, "_bias.png"))
    image(A_bias,main = paste0(mn," bias"))
    dev.off()
  }
  ##------------------------
  
}

plot_animation = function(dt,
                          file_out = "./figures/system_animation.pdf")
{

  YM_all = dt[,unique(YM)]
  if(print_figs){pdf(file_out)}else{X11()}

  
  for(j in 1:length(YM_all))
  {
    print(j)
    YM_j = YM_all[j]
    ##------- Setup --------
    dt_sub = dt[YM == YM_j]
    n_lon = length(dt_sub[,unique(Lon)])
    n_lat = length(dt_sub[,unique(Lat)])
    mn = paste0(dt_sub[,min(month)], "/",dt_sub[,min(year)])
    ##----------------------
    
    ##----- Observation -------
    A_obs = matrix(NA, n_lon,n_lat)
    A_obs[ dt_sub[grid_id < length(A_obs), grid_id] ] = dt_sub[grid_id < length(A_obs),SST_bar] ## FIX
    ##-------------------------

    ##----- Ensemble -----------
    A_ens = matrix(NA, n_lon,n_lat)
    A_ens[ dt_sub[grid_id < length(A_ens), grid_id] ] = dt_sub[grid_id < length(A_ens),SST_hat_grid] ## FIX
    ##-------------------------

    ##----- Bias --------------
    A_bias = A_obs - A_ens
    image(A_bias,main = paste0(mn," bias"))
  }

  if(print_figs) dev.off()
  
}
               
