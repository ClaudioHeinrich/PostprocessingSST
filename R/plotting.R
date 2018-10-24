
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
  
  fields::image.plot(x,y,z,zlim = zlim, col = col, breaks = breaks, font.main = 2,...) 
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
#' @param xlab,ylab labeling of the plot.
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
                            lons = NULL, lats = NULL,
                            xlab = NULL, ylab = NULL,
                            rr = NULL,
                            cex = 1,
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
  
  Lons = sort(Lons)
  Lats = sort(Lats)
  
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
  
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 10))
  brk.lab = round(brk[brk.ind],2)
  brk.at = brk[brk.ind]
  
  if(col_scheme == "bwr"){
    if(is.null(set_white)){
    color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"))
    }else{
       zero.ind = min(which(brk > set_white))/length(brk)
       color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
    }
  }
  if(col_scheme == "wr"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","darkred"))
  }
  if(col_scheme == "wb"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","blue"))
  }
    
  #--- plotting ---
  
  if (is.null(stretch_par)) stretch_par = n_lat/n_lon
  
  if(save_pdf) pdf(paste0(save_dir,file_name,".pdf"),width = 7,height = stretch_par * 7)
  
  if(is.null(xlab))
  {
    xlab = "Longitude"
  }
  
  if(is.null(ylab))
  {
    ylab = "Latitude"
  }
  
  
  image.plot.na(Lons,Lats,A,
                  xlab=xlab,ylab=ylab,
                  zlim=rr,
                  xlim = range(Lons),
                  ylim = range(Lats),
                  breaks=brk,
                  col=color,
                  main = mn,
                  cex.main=cex,cex.lab=1.4,
                  cex.axis=1,
                  axis.args=list(cex.axis=1,
                               at = brk.at,
                               label = brk.lab))
      maps::map("world", add = TRUE)
  
  if(save_pdf) dev.off()
  
}


#' smooth plotting function
#' 
#' @description Takes a data table of the form .(Lon,Lat,value) or .(Lat,Lon,value) and plots values on the globe after applying a kernel smoothing
#' 
#' @param dt The data table containing the values for plotting.
#' @param theta parameter for the Gaussian smoothing kernel
#' @param pixels Resolution of the plot
#' @param mn Title of the plot.
#' @param save_pdf If TRUE, the plot is saved as pdf.
#' @param save_dir,file_name Directory and file name for saving the plot, only used if save_pdf = TRUE (string).
#' @param lons,lats Vectors with two entries containing min and max longitude and latitude for plotting rectangle.
#'             If NULL, the entire globe is used.
#' @param xlab,ylab labeling of the plot.
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
#' @importFrom fields as.image designer.colors image.plot image.smooth
#' @importFrom maps map
#' @importFrom sp SpatialPoints


plot_smooth = function( dt, 
                        theta = 0.5,
                        pixels = 256,
                        mn = "",
                        save_pdf = FALSE,
                        save_dir = "./figures/",
                        file_name = "diag_plot",
                        xlab = NULL, ylab = NULL,
                        rr = NULL,
                        cex = 1,
                        set_white = NULL,
                        col_scheme = "bwr",
                        stretch_par = NULL)
{
  
  data(wrld_simpl) 
  
  #--- get longitudes and latitudes
  
  Lons = sort(unique(dt[,Lon]))
  Lats = sort(unique(dt[,Lat]))
  
  
  #--- data ---
  
  dt = dt[Lat %in% Lats & Lon %in% Lons] [order(Lat,Lon)]
  
  
  n_lon = length(Lons)
  n_lat = length(Lats)
  
  x = dt[,.(Lon,Lat)]
  setnames(x,c("Lon","Lat"), c("lon","lat"))
  
  A = matrix(dt[[3]],  n_lon, n_lat)
  
  im_0 = fields::image.smooth(fields::as.image(A,x = x,nx = pixels,ny = pixels),theta = theta)
  
  all_loc = expand.grid(lat = im_0$x,lon = im_0$y)
  
  pts <- sp::SpatialPoints(all_loc, proj4string=CRS(proj4string(wrld_simpl)))
  
  ## Find which points fall over land
  ii <- !is.na(over(pts, wrld_simpl)$FIPS)
  
  im_0$z[ii] = NA
  
  if(is.null(rr))  rr = range(im_0$z,na.rm=TRUE)
  if(!is.null(rr)){
    im_0$z[im_0$z< min(rr)] = min(rr)
    im_0$z[im_0$z> max(rr)] = max(rr)
  }
  
  # --- scaling and colors ---
  
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 10))
  brk.lab = round(brk[brk.ind],2)
  brk.at = brk[brk.ind]
  
  if(col_scheme == "bwr"){
    if(is.null(set_white)){
      color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"))
    }else{
      zero.ind = min(which(brk > set_white))/length(brk)
      color <- fields::designer.colors(n=length(brk)-1, col = c("darkblue","white","darkred"), x = c(0,zero.ind,1))
    }
  }
  if(col_scheme == "wr"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","darkred"))
  }
  if(col_scheme == "wb"){
    color <- fields::designer.colors(n=length(brk)-1, col = c("white","blue"))
  }
  
  
  
  # color NAs grey
  
  newz.na <- rr[2]+(rr[2]-rr[1])/length(color) # new z for NA
  
  im_0$z[which(is.na(im_0$z))] <- newz.na 
  
  rr[2] <- newz.na # we extend the z limits to include the new value 
  
  color <- c(color, 'gray') # we construct the new color range by including na.color 
  
  brk = c(brk,rr[2])
  
  
  
  #--- plotting ---
  
  if (is.null(stretch_par)) stretch_par = n_lat/n_lon
  
  if(save_pdf) pdf(paste0(save_dir,file_name,".pdf"),width = 7,height = stretch_par * 7)
  
  if(is.null(xlab))
  {
    xlab = "Longitude"
  }
  
  if(is.null(ylab))
  {
    ylab = "Latitude"
  }
  
  
  fields::image.plot(im_0,
             xlab=xlab,ylab=ylab,
             zlim=rr,
             xlim = range(Lons),
             ylim = range(Lats),
             breaks=brk,
             col=color,
             main = mn,
             cex.main=cex,cex.lab=1.4,
             cex.axis=1,
             axis.args=list(cex.axis=1,
                            at = brk.at,
                            label = brk.lab))
  maps::map("world", add = TRUE)
  
  if(save_pdf) dev.off()
  
}


