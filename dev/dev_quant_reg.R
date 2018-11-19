## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")
out_path = "/home/alex/NR/SFE/Presentations/20181120_TLT/fig_alex/"
print_figs = FALSE

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

load("./FcNov2018/ts_hindcast_slimmed.RData")

ens_names = c("norcpm_ts_bar", "ecmwf_ts_bar","mf_ts_bar","ukmo_ts_bar")
DT_final[,equal_weighted := rowMeans(.SD), .SDcol = ens_names]
DT_final[,persistence := obs_anamoly_1]
for(y in 2006:2017){

    print(y)
    DT_train = DT_final[year < y]
    mod = rq(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    DT_final[year == y,
             trained_model := climatology + predict(mod,newdata = DT_final[year == y])]
}

DT_final[,residual:=obs_erai_ts - trained_model]
DT_train_list = list()
nms_run = c("var_global","var_grid","var_month","var_month_grid")
for(y in 2010:2017){

    print(y)
    mod = rq("residual ~ as.factor(month) + trained_model + climatology + Lon + Lat", data = DT_final[year < y], tau = .05)
    mod2 = rq("residual ~ as.factor(month) + trained_model + climatology", data = DT_final[year < y], tau = .05)
    mod3 = rq("residual ~ Lon + Lat", data = DT_final[year < y], tau = .05)
    DT_train = DT_final[year < y]
    DT_final[year == y,Q_reg:=predict(mod, newdata = DT_final[year == y])]
    DT_final[year == y,Q_reg2:=predict(mod2, newdata = DT_final[year == y])]
    DT_final[year == y,Q_reg3:=predict(mod3, newdata = DT_final[year == y])]
}

rr_all = c(0,.1)

X11()
DT_plot = square_image(DT_final[,mean(residual < Q_reg, na.rm=TRUE),.(Lon,Lat)])
plot_smooth2(DT_plot, rr = rr_all, set_white = .05, mn = "Monthly Local Variance")
X11()
DT_plot = square_image(DT_final[,mean(residual < Q_reg2, na.rm=TRUE),.(Lon,Lat)])
plot_smooth2(DT_plot, rr = rr_all, set_white = .05, mn = "Monthly Local Variance")
X11()
DT_plot = square_image(DT_final[,mean(residual < Q_reg3, na.rm=TRUE),.(Lon,Lat)])
plot_smooth2(DT_plot, rr = rr_all, set_white = .05, mn = "Monthly Local Variance")



if(print_figs) dev.off()

DT_plot = square_image(DT_q[,mean(obs_erai_ts < Q_var_month),.(Lon,Lat)])
if(print_figs){pdf(paste0(out_path,"q_valid_month_global.pdf"))}else{X11()}
plot_smooth2(DT_plot, rr = c(0,.1), set_white = .05, mn = "Monthly Global Variance")
if(print_figs) dev.off()

DT_plot = square_image(DT_q[,mean(obs_erai_ts < Q_var_global),.(Lon,Lat)])
if(print_figs){pdf(paste0(out_path,"q_valid_annual_global.pdf"))}else{X11()}
plot_smooth2(DT_plot, rr = c(0,.1), set_white = .05, mn = "Annual Global Variance")
if(print_figs) dev.off()

DT_plot = square_image(DT_q[,mean(obs_erai_ts < Q_var_grid),.(Lon,Lat)])
if(print_figs){pdf(paste0(out_path,"q_valid_annual_grid.pdf"))}else{X11()}
plot_smooth2(DT_plot, rr = c(0,.1), set_white = .05, mn = "Annual Local Variance")
if(print_figs) dev.off()

