
#---- for coloring NA values

image.plot.na <- function(x,y,z,zlim,  col, na.color='gray', breaks, ...)
{
  newz.na <- zlim[2]+(zlim[2]-zlim[1])/length(col) # new z for NA
  
  z[which(is.na(z))] <- newz.na # we affect newz.outside
  
  zlim[2] <- newz.na # we finally extend the z limits to include the new value 
  
  col <- c(col, na.color) # we construct the new color range by including: na.color and outside.color
  
  breaks = c(breaks,zlim[2])
  
  image.plot(x,y,z,zlim=zlim, col=col, breaks = breaks,...) 
}




plot_system = function(Y = 1999,
                          M = 7,
                          type = "res",    #'res' plots residuals, 
                                           #'obs' observed SST, 
                                           #'for' forecasted SST using PCA generated noise,
                                           #'PC' the dth principal component (upscaled eigenvector) where d=depth
                                           #'PCsum' sum over the first d PCs
                                           #'mar_sd' marginal standard deviation computed for the first d PCs
                          obs_num = "mean",     # takes numbers from 1 to 9 or "mean", only used for 
                                                # type = 'obs' or 'res'
                          depth = 10,  
                          file_dir = "./figures/",
                          data.dir = "./Data/PostClim/SFE/Derived/",
                          lons = NULL,
                          lats = NULL,
                          rr = NULL,
                          outside_control = FALSE,
                          print_figs = TRUE,
                          png_out = FALSE)
{
  #------- file name -----------
  
  if(type == "obs")file_out = paste0(file_dir,type,"_y",Y,"_m",M)
  if(type == "for")file_out = paste0(file_dir,type,"_y",Y,"_m",M)
  if(type == "res")file_out = paste0(file_dir,type,"_y",Y,"_m",M)
  if(type == "mar_sd") file_out = paste0("./figures/",type,depth,"_m",M)
  if(type == "PC") file_out = paste0("./figures/",type,depth,"_m",M)
  if(type == "PCsum") file_out = paste0("./figures/",type,depth,"_m",M)
  
  ## ----- load observation ------
  if(type == "obs" | type == "res"){
  
    dt_obs = load_observations(Y,M)
    mn_obs = paste0(" Obs ",obs_num)
    file_out = paste0(file_out,"_",obs_num,"_")
  
    if(obs_num == "mean") {
      dt_obs[,"mean" := mean(SST), by = .(Lon, Lat)]
      dt_obs = dt_obs[Obs == 1,]
      dt_obs = dt_obs[,SST := mean]
      obs_num = 1
    }
  
    dt_obs = dt_obs[Obs == obs_num,]
    dt_obs = dt_obs[order(Lat,Lon)]
    
  }
  ##------------------
  
  #---- load forecast ------
  
  if(type == "res" | type == "for"){
  
      dir.name = "./Data/PostClim/SFE/Derived/PCA"
      file.name = paste0("/fc_",depth,"pc_",Y,"_",M,".RData")
      load(paste0(dir.name,file.name))
      dt_for = fc_land
      
      mn_for = paste0(" FC depth ",depth)
  }
  
  
  if(type == "res"){
    dt_res = copy(dt_obs)
    dt_res[,"Res":= dt_for[,forecast] - SST]
    }
  
  #---- load rest ------
  
  if(type == "mar_sd") {
    file.name = paste0("PCA/PCA_mar_sd",depth,".RData")
    load(paste0(data.dir,file.name))
    dt_for = fc_land
  }
  
  if(type == "PC") {
    file.name = paste0("PCA/PCA_",depth,"PC.RData")
    load(paste0(data.dir,file.name))
    dt_for = fc_land
  }
  
  if(type == "PCsum") {
    file.name = paste0("PCA/PCA_",depth,"sum.RData")
    load(paste0(data.dir,file.name))
    dt_for = fc_land
  }
  
  
  #----- get number of lons and lats ----
  
  if(type %in% c("res","obs")) {
    lon_all = sort(dt_obs[,unique(Lon)])
    lat_all = sort(dt_obs[,unique(Lat)])}
  
  if(type %in% c("for","mar_sd","PC","PCsum")) {
    lon_all = sort(dt_for[,unique(Lon)])
    lat_all = sort(dt_for[,unique(Lat)])}
  
  n_lon = length(lon_all)
  n_lat = length(lat_all)
  
  #------- adjusting range -------
  
  if(is.null(rr)) {
    if(type == "obs")   rr = range(dt_obs[,SST],na.rm=TRUE)
    if(type == "for")   rr = range(dt_for[,forecast],na.rm=TRUE)
    if(type == "res")   rr = range(dt_res[,Res],na.rm=TRUE)
    if(type == "mar_sd")rr = range(dt_for[,forecast],na.rm=TRUE)
    if(type == "PC")    rr = range(dt_for[,forecast],na.rm=TRUE)
    if(type == "PCsum") rr = range(dt_for[,forecast],na.rm=TRUE)
  } 
  
  
  #------- titles for plot -------
  
  if(type == "res")     mn = paste0("Residual ",M, "/",Y,mn_obs,mn_for)
  if(type == "obs")     mn = paste0(mn_obs," ",M, "/",Y)
  if(type == "for")     mn = paste0(mn_for,M, "/",Y)
  if(type == "mar_sd")  mn =  paste0("marginal standard deviation for ",depth," PCs")
  if(type == "PC")      mn =  paste0(depth,". principal component")
  if(type == "PCsum")   mn =  paste0("sum over first ",depth," principal components")
  
  ##----- get plotting data --------------
  
  if(type == "res")   A = matrix(dt_res[, Res],  n_lon, n_lat)
  if(type == "obs")   A = matrix(dt_obs[, SST],  n_lon, n_lat)
  if(type == "for")   A = matrix(dt_for[, forecast],  n_lon, n_lat)
  if(type == "mar_sd")A = matrix(dt_for[, forecast],  n_lon, n_lat)
  if(type == "PC")    A = matrix(dt_for[, forecast],  n_lon, n_lat)
  if(type == "PCsum") A = matrix(dt_for[, forecast],  n_lon, n_lat)
    
  
  
  ##------- Scaling and colors----------
  if (type %in% c("obs","for")){
    brk = seq(rr[1],rr[2],length = 500)
    brk.ind = round(seq(1,length(brk),length = 10))
    brk.lab = round(brk[brk.ind],2)
    brk.at = brk[brk.ind]
    color <- designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"))
  }
  
  if (type %in% c("res","PC","PCsum")){
    brk = seq(rr[1],rr[2],length = 500)
    brk.ind = round(seq(1,length(brk),length = 10))
    zero.ind = min(which(brk > 0))/length(brk)
    brk.at = brk[brk.ind]
    brk.lab = round(brk[brk.ind],2)
    color <- designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
  }
  
  if (type %in% c("mar_sd")){
    brk = seq(rr[1],rr[2],length = 500)
    brk.ind = round(seq(1,length(brk),length = 10))
    zero.ind = min(which(brk > 0))/length(brk)
    brk.at = brk[brk.ind]
    brk.lab = round(brk[brk.ind],2)
    color <- rgb(139,0,0, # is the rgb specification for darkred, relative to maxColor = 255
                 alpha = seq(0,255,,length(brk)-1),
                 maxColorValue = 255)
  }
  
  
  
  
  
   
  ##--------------------------
  
  ##------- Scope ------------
  if(is.null(lons)){
    if(type == "res")   lons = range(dt_res[,Lon])
    if(type == "obs")   lons = range(dt_obs[,Lon])
    if(type == "for")   lons = range(dt_for[,Lon])
    if(type == "mar_sd")lons = range(dt_for[,Lon])
    if(type == "PC")    lons = range(dt_for[,Lon])
    if(type == "PCsum") lons = range(dt_for[,Lon])
  }
  
  if(is.null(lats)){
    if(type == "res")   lats = range(dt_res[,Lat])
    if(type == "obs")   lats = range(dt_obs[,Lat])
    if(type == "for")   lats = range(dt_for[,Lat])
    if(type == "mar_sd")lats = range(dt_for[,Lat])
    if(type == "PC")    lats = range(dt_for[,Lat])
    if(type == "PCsum") lats = range(dt_for[,Lat])
  }
    
    
  ##--------------------------
  
  ##------- Plot -----------
  if(print_figs & !outside_control)
  {
    pdf(paste0(file_out, ".pdf"))
  }else{
    if(!outside_control) X11()
  }
  
  image.plot.na(lon_all,lat_all,A,
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
  
  if(png_out){
    if(print_figs & !outside_control)
    {
      png(paste0(file_out, ".png"))
      image.plot.na(lon_all,lat_all,A_res,
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
}


##------------------------



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
   
  