
#' Creates a plot on the globe with a specified color for NA entries
#' 
#' @param na.color The color for NA values.
#' @param x,y,z,zlim,col,breaks,... Parameters passed on to image.plot. 
#' 
#' @return none
#'
#' @export image.plot.na
#' 
#' @author Claudio Heinrich
#' @importFrom fields image.plot
image.plot.na <- function(x,y,z,zlim,  col, na.color='gray', breaks, ...)
{
  newz.na <- zlim[2]+(zlim[2]-zlim[1])/length(col) # new z for NA
  
  z[which(is.na(z))] <- newz.na 
  
  zlim[2] <- newz.na # we extend the z limits to include the new value 
  
  col <- c(col, na.color) # we construct the new color range by including na.color 
  
  breaks = c(breaks,zlim[2])
  
  image.plot(x,y,z,zlim = zlim, col = col, breaks = breaks,...) 
}


#' Diagnostic plotting function
#' 
#' @description Takes a data table of the form .(Lon,Lat,value) or .(Lat,Lon,value) and plots value on the globe
#' 
#' @param dt The data table containing the values for plotting.
#' @param mn Title of the plot.
#' @param save_pdf If TRUE, the plot is saved as pdf.
#' @param save_dir,file_name Directory and file name for saving the plot, only used if save_pdf = TRUE (string).
#' @param lons,lats Vectors with two entries containing min and max longitude and latitude for plotting rectangle.
#'             If NULL, the entire globe is used.
#' @param rr Range of the plot, if not specified the range of value is used.
#' @param set_white Forces the blue-white-red color scheme to center white at the set value if specified.
#' @param col_scheme Either of "bwr" for blue - white - red, "wr" for white - red, or "wb" for white - blue. Specifies the color scheme of the plot. 
#' @param stretch_par Numeric. Only used when save_pdf == TRUE. Stretches the pdf output. Default is NULL, where it is stretched to #lons/#lats.
#'
#' @return none
#'  
#' @export
#' 
#' @author Claudio Heinrich
#' @examples \dontrun{
#' dt = load_combined_wide()
#' plot_diagnostic(dt[year == 1990 & month == 1,.(Lon,Lat,Ens_bar)])
#' }
#' 
#' @importFrom fields designer.colors 
#' @importFrom maps map

plot_diagnostic = function( dt, 
                            mn = "",
                            save_pdf = FALSE,
                            save_dir = "./figures/",
                            file_name = "diag_plot",
                            lons = NULL,
                            lats = NULL,
                            rr = NULL,
                            set_white = NULL,
                            col_scheme = "bwr",
                            stretch_par = NULL)
  {
  #--- get longitudes and latitudes
  
  if(is.null(lons)) Lons = unique(dt[,Lon])
  if(is.null(lats)) Lats = unique(dt[,Lat])
  
  if(!is.null(lons)){
    lon_min = floor(lons[1]+0.5)-0.5
    lon_max = floor(lons[2]+0.5)-0.5
    Lons = seq(lon_min,lon_max,by = 1)  
  }
  
  if(!is.null(lats)){
    lat_min = floor(lats[1]+0.5)-0.5
    lat_max = floor(lats[2]+0.5)-0.5
    Lats = seq(lat_min,lat_max,by = 1)  
  }
  
  #--- data ---
  
  dt = dt[Lat %in% Lats & Lon %in% Lons] [order(Lat,Lon)]
  
  n_lon = length(Lons)
  n_lat = length(Lats)
  
  if(is.null(rr))  rr = range(dt[,3],na.rm=TRUE)
  if(!is.null(rr)){
    dt[dt[[3]]<min(rr),3] = min(rr)
    dt[dt[[3]]> max(rr),3] = max(rr)
  }
  
  A = matrix(dt[[3]],  n_lon, n_lat)
  
  # --- scaling and colors ---
  
  # make color white for the value achieved by perfect calibration
 
  
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 10))
  brk.lab = round(brk[brk.ind],2)
  brk.at = brk[brk.ind]
  
  if(col_scheme == "bwr"){
    if(is.null(set_white)){
    color <- designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"))
    }else{
       zero.ind = min(which(brk > set_white))/length(brk)
       color <- designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
    }
  }
  if(col_scheme == "wr"){
    color <- designer.colors(n=length(brk)-1, col = c("white","darkred"))
  }
  if(col_scheme == "wb"){
    color <- designer.colors(n=length(brk)-1, col = c("white","blue"))
  }
    
  #--- plotting ---
  
  if (is.null(stretch_par)) stretch_par = n_lat/n_lon
  
  if(save_pdf) pdf(paste0(save_dir,file_name,".pdf"),width = 7,height = stretch_par * 7)
  
  image.plot.na(Lons,Lats,A,
                  xlab="Longitude",ylab="Latitude",
                  zlim=rr,
                  xlim = range(Lons),
                  ylim = range(Lats),
                  breaks=brk,
                  col=color,
                  main = mn,
                  cex.main=1.8,cex.lab=1.4,
                  cex.axis=1,
                  axis.args=list(cex.axis=1,
                               at = brk.at,
                               label = brk.lab))
      map("world", add = TRUE)
  
  if(save_pdf) dev.off()
  
}



#' Plots the system with various options 
#'
#' @param type Specifies the type of plot: \cr
#'      'res' plots residuals \cr
#'      'obs' observed SST \cr
#'      'ens' the raw ensemble forecast \cr
#'      'for' forecasted SST using PCA generated noise, \cr
#'      'PC' the dth principal component (upscaled eigenvector) where d=depth \cr
#'      'PCsum' sum over the first d PCs \cr
#'      'mar_sd' marginal standard deviation for the PCA method using the first d PCs \cr
#'      'cal' plots moment estimates of the PIT, uses all years and months.
#' @param obs_num Index of used observation (if any) in the ensemble, also takes "mean".
#' @param Y,M Year and month.
#' @param depth Specifies the numbers of principal components to consider.
#' @param moment Takes 1 or 2, only used for type = "cal", 1 plots mean, 2 plots SD of the PIT.
#' @param file_dir,data_dir Directories for saving the plot and for loading the data.
#' @param lons,lats Vector with two entries containing min and max longitude and latitude for plotting rectangle.
#'             If NULL, the entire globe is used.
#' @param rr Range of the plot. 
#' @param plot_title The title of the plot.
#' @param print_figs  If FALSE, the plot is displayed directly and no file is saved.
#' @param png_out If TRUE, a .png is saved rather than a .pdf.
#' 
#' @export 
#' 
#' @author Claudio Heinrich
#' @examples \dontrun{
#' setwd("~/NR/SFE")
#' plot_system()}
#'
#' 
#' @importFrom fields designer.colors
#' @importFrom maps map

plot_system = function(   type = "res",    
                          obs_num = "mean",     # takes numbers from 1 to 9 or "mean", only used for 
                                                # type = 'obs' or 'res'
                          Y = 1999,
                          M = 7,
                          moment = 1, #only used for type = 'cal'
                          depth = 0,  
                          file_dir = "./figures/",
                          data_dir = "./Data/PostClim/SFE/Derived/",
                          lons = NULL,
                          lats = NULL,
                          rr = NULL,
                          plot_title = NULL,
                          print_figs = TRUE,
                          png_out = FALSE)
{
  #------- file name -----------
  
  if(type == "obs")file_out = paste0(file_dir,type,"_y",Y,"_m",M)
  if(type == "for")file_out = paste0(file_dir,type,"_y",Y,"_m",M)
  if(type == "ens")file_out = paste0(file_dir,type,"_y",Y,"_m",M)
  if(type == "res")file_out = paste0(file_dir,type,"_y",Y,"_m",M)
  if(type == "mar_sd") file_out = paste0("./figures/",type,depth,"_m",M)
  if(type == "PC") file_out = paste0("./figures/",type,depth,"_m",M)
  if(type == "PCsum") file_out = paste0("./figures/",type,depth,"_m",M)
  if(type == "cal" & depth > 0) file_out = paste0("./figures/",type,"_PC",depth,"_mom",moment)
  if(type == "cal" & depth == 0) file_out = paste0("./figures/",type,"ens_mom",moment)
  if(type == "cal" & depth == -1) file_out = paste0("./figures/",type,"_clim_mom",moment)
  
  
  #--- get longitudes and latitudes
  
  if(is.null(lons)) {Lons = seq(-179.5,179.5,by = 1)
                     lons = c(-179.5,179.5)
  }
  if(is.null(lats)) {Lats = seq(-89.5,89.5,by = 1)
                     lats = c(-89.5,89.5)
  }
  
  if(!is.null(lons)){
    lon_min = floor(lons[1]+0.5)-0.5
    lon_max = floor(lons[2]+0.5)-0.5
    Lons = seq(lon_min,lon_max,by = 1)  
  }
  
  if(!is.null(lats)){
    lat_min = floor(lats[1]+0.5)-0.5
    lat_max = floor(lats[2]+0.5)-0.5
    Lats = seq(lat_min,lat_max,by = 1)  
  }
  
  
  ## ----- load observation ------
  if(type == "obs" | type == "res"){
  
    dt_obs = load_observations(Y,M)
    dt_obs = dt_obs[Lon %in% Lons & Lat %in% Lats]
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
      file_name = paste0("/fc_",depth,"pc_",Y,"_",M,".RData")
      load(paste0(dir.name,file_name))
      dt_for = fc_land[Lon %in% Lons & Lat %in% Lats]
      
      mn_for = paste0(" FC depth ",depth)
  }
  
  
  if(type == "res"){
    dt_res = copy(dt_obs)
    dt_res[Lon %in% Lons & Lat %in% Lats,"Res":= dt_for[,forecast] - SST]
    }
  
  #---- load rest ------
  
  if(type == "ens"){
    file_out = paste0(file_out,"_",obs_num)
    
    dt = load_combined_wide()
    mn_ens = paste0(" Ens ",obs_num)
    ens_ind = paste0("Ens",obs_num)
    if(obs_num == "mean") ens_ind = "Ens_bar"
    
   dt_for = dt[Lon %in% Lons & Lat %in% Lats][year %in% Y & month %in% M,.(Lon,Lat, eval(parse(text = ens_ind)))]
   setnames(dt_for,"V3", "Ens")
  }
    
  
  if(type == "mar_sd") {
    file_name = paste0("PCA/PCA_mar_sd",depth,"_m",M,".RData")
    load(paste0(data_dir,file_name))
    dt_for = fc_land[Lon %in% Lons & Lat %in% Lats]
  }
  
  if(type == "PC") {
    file_name = paste0("PCA/PCA_",depth,"PC_month",M,".RData")
    load(paste0(data_dir,file_name))
    dt_for = fc_land[Lon %in% Lons & Lat %in% Lats]
  }
  
  if(type == "PCsum") {
    file_name = paste0("PCA/PCA_",depth,"sum.RData")
    load(paste0(data_dir,file_name))
    dt_for = fc_land[Lon %in% Lons & Lat %in% Lats]
  }
  
  if(type == "cal"){
    file_name = paste0("PCA/cal",depth,"_mom",moment,".RData")
    if(depth == 0) file_name = paste0("PCA/cal_ens_mom",moment,".RData")
    if(depth == -1) file_name = paste0("PCA/cal_clim_mom",moment,".RData")
    load(paste0(data_dir,file_name))
  }
  
  
  #----- get number of lons and lats ----
  
  if(type %in% c("res","obs")) {
    lon_all = sort(dt_obs[,unique(Lon)])
    lat_all = sort(dt_obs[,unique(Lat)])}
  
  if(type %in% c("ens","for","mar_sd","PC","PCsum")) {
    lon_all = sort(dt_for[,unique(Lon)])
    lat_all = sort(dt_for[,unique(Lat)])}
  
  if(type == "cal") {
    lon_all = sort(calib[,unique(Lon)])
    lat_all = sort(calib[,unique(Lat)])
  }
  
  n_lon = length(lon_all)
  n_lat = length(lat_all)
  
  #------- adjusting range -------
  
  if(is.null(rr)) {
    if(type == "obs")   rr = range(dt_obs[,SST],na.rm=TRUE)
    if(type == "for")   rr = range(dt_for[,forecast],na.rm=TRUE)
    if(type == "ens")   rr = range(dt_for[,Ens],na.rm=TRUE)
    if(type == "res")   rr = range(dt_res[,Res],na.rm=TRUE)
    if(type == "mar_sd")rr = range(dt_for[,forecast],na.rm=TRUE)
    if(type == "PC")    rr = range(dt_for[,forecast],na.rm=TRUE)
    if(type == "PCsum") rr = range(dt_for[,forecast],na.rm=TRUE)
    if(type == "cal" & moment ==1) rr = c(0,1)
    if(type == "cal" & moment ==2) rr = c(0,.5)
  } 
  
  
  #------- titles for plot -------
  
  if(!is.null(plot_title)){
    mn = plot_title
  }else{
    if(type == "res")     mn = paste0("Residual ",M, "/",Y,mn_obs,mn_for)
    if(type == "obs")     mn = paste0(mn_obs," ",M, "/",Y)
    if(type == "ens")     mn = paste0(mn_ens," ",M, "/",Y)
    if(type == "for")     mn = paste0(mn_for,M, "/",Y)
    if(type == "mar_sd")  mn =  paste0("marginal standard deviation for ",depth," PCs")
    if(type == "PC")      mn =  paste0(depth,". principal component")
    if(type == "PCsum")   mn =  paste0("sum over first ",depth," principal components")
    if(type == "cal")     mn =  paste0("PIT variance, ",depth," principal components")
  }
    
  ##----- get plotting data --------------
  
  if(type == "res")   A = matrix(dt_res[, Res],  n_lon, n_lat)
  if(type == "obs")   A = matrix(dt_obs[, SST],  n_lon, n_lat)
  if(type == "ens")   A = matrix(dt_for[, Ens],  n_lon, n_lat,byrow = TRUE)
  if(type == "for")   A = matrix(dt_for[, forecast],  n_lon, n_lat)
  if(type == "mar_sd")A = matrix(dt_for[, forecast],  n_lon, n_lat)
  if(type == "PC")    A = matrix(dt_for[, forecast],  n_lon, n_lat)
  if(type == "PCsum") A = matrix(dt_for[, forecast],  n_lon, n_lat)
  if(type == "cal" & moment == 1)   A = matrix(calib[month == min(month) & year == min(year), moment],
                                 n_lon, n_lat)
  if(type == "cal" & moment == 2)   A = matrix(calib[month == min(month) & year == min(year), 
                                                     sqrt(moment)], n_lon, n_lat)
  
    
  
  
  ##------- Scaling and colors----------
  
  if (type %in% c("obs","for","ens")){
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
  
  
  if (type %in% c("cal")){
    brk = seq(rr[1],rr[2],length = 500)
    brk.ind = round(seq(1,length(brk),length = 10))
    
    # make color white for the value achieved by perfect calibration
    if(moment == 1) zero.ind = min(which(brk > 0.5))/length(brk)
    if(moment == 2) zero.ind = min(which(brk > .2887))/length(brk)
    
    brk.at = brk[brk.ind]
    brk.lab = round(brk[brk.ind],2)
    color <- designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
  }
  
 
  
   
  ##--------------------------
 
  
  ##------- Plot -----------
  if(print_figs)  pdf(paste0(file_out, ".pdf"))
  
  
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
  if(print_figs ) dev.off()
  
  if(png_out){
    if(print_figs)
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

   
  
