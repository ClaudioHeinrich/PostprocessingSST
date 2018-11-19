rm(list = ls())

library(PostProcessing)
library(data.table)
library(party)
print_figs = FALSE

setwd("~/PostClimDataNoBackup/SFE/")

load("./FcNov2018/ts_hindcast_slimmed.RData")
out_path = "/home/alex/NR/SFE/Presentations/20181120_TLT/fig_alex/"
theta = 0.5
pixels = 256
col_scheme = "bwr"
set_white = NULL

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

ctl = mob_control(verbose = FALSE)
y = 2010
print(y)
DT_train = DT_final[year < y]
mod = mob(obs_anamoly ~ ecmwf_anamoly| Lon + Lat,
          data = DT_train,
          control = ctl)

DT_test = DT_final[year == y & month == 7]
DT_test[, ecmwf_anamoly := 0]

DT_test[,pred_alpha:=predict(mod, newdata = DT_test)]

DT_test[, ecmwf_anamoly := 1.0]

DT_test[,pred_beta:=predict(mod, newdata = DT_test) - pred_alpha]

DT_plot = square_image(DT_test,var = "pred_alpha")
if(print_figs){pdf(paste0(out_path,"mob_alpha.pdf"))}else{X11()}
plot_smooth2(DT_plot, col_scheme = "bwr", rr = c(-.33,.33), mn = "Intercept")
if(print_figs) dev.off()

DT_plot = square_image(DT_test,var = "pred_beta")
if(print_figs){pdf(paste0(out_path,"mob_beta.pdf"))}else{X11()}
plot_smooth2(DT_plot, col_scheme = "bwr", rr = c(.8,1.2), mn = "ECMWF Slope")
if(print_figs) dev.off()
