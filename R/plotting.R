plot_system = function(dt,
                       YM_j,
                       file_out = paste0("./figures/system_",YM_j),
                       lons = NULL,
                       lats = NULL,
                       var_plot = "residual",
                       rr = NULL,
                       mn_add = "",
                       outside_control = FALSE,
                       vintage = "mr")
{
  
  if(vintage != "mr") file_out <- paste0(file_out,vintage)
  
  ##----- HACK, fix -------
  ##Load globals
  dt_obs = load_observations(2000,1)
  lon_all = sort(dt_obs[,unique(Lon)])
  lat_all = sort(dt_obs[,unique(Lat)])
  n_lon = length(lon_all)
  n_lat = length(lat_all)
  ##------------------
  
  ##------- Setup --------
  dt_sub = dt[YM == YM_j]
  mn = paste0(dt_sub[,month][1], "/",dt_sub[,year][1],mn_add)
  if(is.null(rr)) rr = range(dt[,eval(parse(text = var_plot))],na.rm=TRUE)
  ##----------------------

  ##----- Bias --------------
  A_bias = matrix(NA, n_lon, n_lat)
  A_bias[dt_sub[, grid_id]] = dt_sub[, eval(parse(text = var_plot))]
  ##-------------------------

  ##------- Scaling ----------
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 10))
  brk.at = brk[brk.ind]
  brk.lab = round(brk[brk.ind],2)
  color <- designer.colors(n=length(brk)-1)
  ##--------------------------

  ##------- Scope ------------
  if(is.null(lons))lons = range(dt_sub[,Lon])
  if(is.null(lats))lats = range(dt_sub[,Lat])
  ##--------------------------

  ##------- Plot -----------
  if(print_figs & !outside_control)
  {
    pdf(paste0(file_out, "_bias.pdf"))
  }
  else{
    if(!outside_control) X11()
  }

  image.plot(lon_all,lat_all,A_bias,
             main=mn,xlab="Longitude",ylab="Latitude",
             zlim=rr,
             xlim = lons,
             ylim = lats,
             breaks=brk,
             col=color,
             cex.main=1.8,cex.lab=1.4,
             cex.axis=1,
             axis.args=list(cex.axis=1,
                            at = brk.at,
                            label = brk.lab))
  map("world", add = TRUE)
  if(print_figs & !outside_control) dev.off()

  if(print_figs & !outside_control)
  {
    png(paste0(file_out, "_bias.png"))
    image.plot(lon_all,lat_all,A_bias,
               main=mn,xlab="Longitude",ylab="Latitude",
               xlim = lons,
               ylim = lats,
               zlim=rr,
               breaks=brk,
               col=color,
               cex.main=1.8,cex.lab=1.4,
               cex.axis=1,
               axis.args=list(cex.axis=1,
                              at = brk.at,
                              label = brk.lab))
    map("world", add = TRUE)
    dev.off()
  }
  ##------------------------
  
}


plot_obs= function(obs_num = 1,
                   YM_j,
                   file_out = paste0("./figures/obs_",YM_j),
                   lons = NULL,
                   lats = NULL,
                   rr = NULL,
                   mn_add = "",
                   outside_control = FALSE)
{
  
  ##----- HACK, fix -------
  ##Load globals
  Y = floor((YM_j-1)/12)
  M = YM_j %% 12
  if(M == 0) M = 12
  dt_obs = load_observations(Y,M)
  dt_obs = dt_obs[Obs == obs_num,]
  lon_all = sort(dt_obs[,unique(Lon)])
  lat_all = sort(dt_obs[,unique(Lat)])
  n_lon = length(lon_all)
  n_lat = length(lat_all)
  ##------------------
  
  ##------- Setup --------
  mn = paste0("Observation ",M, "/",Y,mn_add)
  if(is.null(rr)) rr = range(dt_obs[,SST],na.rm=TRUE)
  ##----------------------
  
  ##----- SST --------------
  A_SST = matrix(dt_obs[, SST],  n_lon, n_lat, byrow = TRUE)
  ##-------------------------
  
  ##------- Scaling ----------
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 10))
  brk.at = brk[brk.ind]
  brk.lab = round(brk[brk.ind],2)
  color <- designer.colors(n=length(brk)-1)
  ##--------------------------
  
  ##------- Scope ------------
  if(is.null(lons))lons = range(dt_obs[,Lon])
  if(is.null(lats))lats = range(dt_obs[,Lat])
  ##--------------------------
  
  ##------- Plot -----------
  if(print_figs & !outside_control)
  {
    pdf(paste0(file_out, "_obs",obs_num,".pdf"))
  }else{
    if(!outside_control) X11()
  }
  
  image.plot(lon_all,lat_all,A_SST,
             main=mn,xlab="Longitude",ylab="Latitude",
             zlim=rr,
             xlim = lons,
             ylim = lats,
             breaks=brk,
             col=color,
             cex.main=1.8,cex.lab=1.4,
             cex.axis=1,
             axis.args=list(cex.axis=1,
                            at = brk.at,
                            label = brk.lab))
  map("world", add = TRUE)
  if(print_figs & !outside_control) dev.off()
  
  if(print_figs & !outside_control)
  {
    png(paste0(file_out, "_bias.png"))
    image.plot(lon_all,lat_all,A_SST,
               main=mn,xlab="Longitude",ylab="Latitude",
               xlim = lons,
               ylim = lats,
               zlim=rr,
               breaks=brk,
               col=color,
               cex.main=1.8,cex.lab=1.4,
               cex.axis=1,
               axis.args=list(cex.axis=1,
                              at = brk.at,
                              label = brk.lab))
    map("world", add = TRUE)
    dev.off()
  }
  ##------------------------
  
}





plot_forecast = function(YM_j,
                         file_out = paste0("./figures/for_",YM_j),
                         lons = NULL,
                         lats = NULL,
                         rr = NULL,
                         method = "PCA",
                         depth=5,
                         outside_control = FALSE,
                         print_figs = TRUE)
{
  
  
  ##------Load globals-----
  
  Y = floor((YM_j-1)/12)
  M = YM_j %% 12
  if(M == 0) M = 12
  
  
  if(method == "PCA"){
    dir.name = "./Data/PostClim/SFE/Derived/PCA"
    file.name = paste0("/fc_",depth,"pc_",Y,"_",M,".RData")
    load(paste0(dir.name,file.name))
    dt_for = fc_land
    if(depth == 0){mn_add = ""
    } else{ mn_add = paste0(", ",depth," principal components")}
    file_out = paste0(file_out,"PCA",depth)
    }
   
  lon_all = sort(dt_for[,unique(Lon)])
  lat_all = sort(dt_for[,unique(Lat)])
  n_lon = length(lon_all)
  n_lat = length(lat_all)
  ##------------------
  
    ##------- Setup --------
  mn = paste0("FC ",M, "/",Y,mn_add)
  if(is.null(rr)) rr = range(dt_for[,forecast],na.rm=TRUE)
  ##----------------------
  
  ##----- SST --------------
  A_SST = matrix(dt_for[, forecast],  n_lon, n_lat)
  ##-------------------------
  
  ##------- Scaling ----------
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 10))
  brk.at = brk[brk.ind]
  brk.lab = round(brk[brk.ind],2)
  color <- designer.colors(n=length(brk)-1)
  ##--------------------------
  
  ##------- Scope ------------
  if(is.null(lons))lons = range(dt_for[,Lon])
  if(is.null(lats))lats = range(dt_for[,Lat])
  ##--------------------------
  
  ##------- Plot -----------
  if(print_figs & !outside_control)
  {
    pdf(paste0(file_out,".pdf"))
  }else{
    if(!outside_control) X11()
  }
  
  image.plot(lon_all,lat_all,A_SST,
             main=mn,xlab="Longitude",ylab="Latitude",
             zlim=rr,
             xlim = lons,
             ylim = lats,
             breaks=brk,
             col=color,
             cex.main=1.8,cex.lab=1.4,
             cex.axis=1,
             axis.args=list(cex.axis=1,
                            at = brk.at,
                            label = brk.lab))
  map("world", add = TRUE)
  if(print_figs & !outside_control) dev.off()
  
  if(print_figs & !outside_control)
  {
    png(paste0(file_out,".png"))
    image.plot(lon_all,lat_all,A_SST,
               main=mn,xlab="Longitude",ylab="Latitude",
               xlim = lons,
               ylim = lats,
               zlim=rr,
               breaks=brk,
               col=color,
               cex.main=1.8,cex.lab=1.4,
               cex.axis=1,
               axis.args=list(cex.axis=1,
                              at = brk.at,
                              label = brk.lab))
    map("world", add = TRUE)
    dev.off()
  }

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
               
