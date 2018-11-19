## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(data.table)
library(qgam)

plot_smooth2 = function( dt, var = colnames(dt)[3], mn = var, rr = NULL,
                        theta = 0.5, pixels = 256,
                        col_scheme = "bwr", set_white = NULL,
                        xlab = "", ylab = "",
                        save_pdf = FALSE, save_dir = "./figures/", file_name = "diag_plot", stretch_par = NULL)
{
  # prepare data table
  
    dt = dt[,.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]

  
  
  #--- create image ---
  
  x = dt[,.(Lon,Lat)]
  setnames(x,c("Lon","Lat"), c("lon","lat"))
  
  Lons = unique(dt[,Lon])
  Lats = unique(dt[,Lat])
  
  n_lon = length(Lons)
  n_lat = length(Lats)
  
  A = matrix(dt[[3]],  n_lon, n_lat)
  
  im_0 = fields::image.smooth(fields::as.image(A,x = x,nx = pixels,ny = pixels),theta = theta)
  
  ## Find the points that fall over land
  
  if(!exists("wrld_simpl")) data(wrld_simpl, package = 'maptools') 
  
  all_loc = expand.grid(lat = im_0$x,lon = im_0$y)
  pts <- sp::SpatialPoints(all_loc, proj4string=sp::CRS(proj4string(wrld_simpl)))
  ii <- !is.na(over(pts, wrld_simpl)$FIPS)
    im_0$z[-which(ii)] = NA
  
  # --- fix range of plot and fill in values for points out of range ---
  
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
  rr[2] <- newz.na # extend the range to include the new value 
  color <- c(color, 'gray') # extend the color range by gray
  brk = c(brk,rr[2]) # extend the vector of breaks
  
  #--- plotting ---
  
  
  
  if(save_pdf) 
    {
    if (is.null(stretch_par)) stretch_par = n_lat/n_lon
    
    par_0 = par() # allow to set par manually before calling the function
    
    pdf(paste0(save_dir,file_name,".pdf"),width = 7,height = stretch_par * 7)
    
    suppressWarnings(par(par_0))
    }
  
  par(mar = c(2,2,2,2))
  
  fields::image.plot(im_0,
                     zlim=rr, main = mn,
                     xlim = range(Lons), xlab=xlab,
                     ylim = range(Lats), ylab=ylab,
                     breaks=brk,
                     col=color,
                     axis.args=list(at = brk.at,
                                    label = brk.lab))
  
  # add world map
  
  maps::map("world", add = TRUE)
  
  if(save_pdf) dev.off()
  
}

square_image = function(dt, var = names(dt)[3])
{
    lons = dt[,unique(Lon)]
    lats = dt[,unique(Lat)]

    DT_c = data.table(Lon = rep(lons,each = length(lats)), Lat  = rep(lats, times = length(lons)))
    DT_plot = merge(DT_c,
                    dt[,.SD,.SDcols = c("Lon","Lat",var)],
                    by = c("Lon","Lat"),
                    all.x = TRUE,
                    all.y = FALSE)

    return(DT_plot)
    

}

setwd("~/PostClimDataNoBackup/SFE/")
out_path = "/home/alex/NR/SFE/Presentations/20181120_TLT/fig_alex/"
print_figs = FALSE

load("./FcNov2018/ts_hindcast_slimmed.RData")

DT_train_list = list()
DT_final[,residual := obs_erai_ts  - climatology]
nms_run = c("var_global","var_grid","var_month","var_month_grid")
for(y in 2006:2017){
    print(y)
    mod = qgam(obs_anamoly ~ te(Lon,Lat, k = 10) + ecmwf_anamoly, data = DT_final[year < y], qu = .05)
    p_hat = predict(mod, newdata = DT_final[year == y])
    DT_final[year == y, Q_pred := p_hat]

}

DT_plot1 = square_image(DT_final[year == 2010 &  month == 7,.(Lon,Lat,climatology + Q_pred)])
DT_plot2 = square_image(DT_final[year == 2012 &  month == 7,.(Lon,Lat,climatology + Q_pred)])
X11();plot_smooth2(DT_plot1, rr = c(6.5,27))##, rr = c(0,.1),set_white = .05)
X11();plot_smooth2(DT_plot2, rr = c(6.5,27))##, rr = c(0,.1),set_white = .05)

A = DT_final[year > 2010,mean(obs_erai_ts < climatology + Q_pred, na.rm = TRUE),.(Lon,Lat)]
DT_plot3 = square_image(A)
plot_smooth2(DT_plot3)

A_fake = A
A_fake[,V1:=predict(mod, newdata = data.table(Lon = Lon, Lat = Lat, ecmwf_anamoly = 0))]
DT_plot4 = square_image(A_fake)
plot_smooth2(DT_plot4)
