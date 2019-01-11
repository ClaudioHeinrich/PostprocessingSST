

##################################################################
### This script generates the plots for section 4 of the paper ###
##################################################################


rm(list = ls())

time_s4 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


########### set parameters for plots #######################

par('cex' = 0.75, 'cex.lab' = 0.6,'cex.axis' = 0.6)

plot_dir0 = './figures/paper/'

dir.create(plot_dir0,showWarnings = FALSE)


Lat_res = c(-75,80) # Latitude restrictions for area plots in order to exclude the polar regions


###############################################

############## Sec 4.0 ########################

###############################################




##################################################

######## Plots for univariate calibration ########

##################################################

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


# reset plot directory

plot_dir = plot_dir0




######################################################

####### Plot MSE by weighting parameter ##############

#### plotting ####

# get data for last year
y = max(validation_years)

row_sma = msc_sma[year == y,-1,with = FALSE]
row_ema = msc_ema[year == y,-1,with = FALSE]

# reduce window length for plot to make it better readable

wl = win_length[5:length(win_length)]

values = as.vector(row_sma[,-c('min_MSE','min_l'),with = FALSE])

values = values[,5:length(win_length)]

y_range = range(list(row_sma[,-'min_l',with = FALSE][,5:length(win_length)],row_ema[,-'min_a',with = FALSE]))  




pdf(paste0(plot_dir,"MSE_by_par.pdf"),width = 15)
  
  par('mfrow' = c(1,2))
  par('cex' = 1.25, 'cex.lab' = 0.9,'cex.axis' = 0.75)
  
  ## plot for sma ##
  
  plot(x = wl, 
       y = values,
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("MSE for bias correction by SMA"),
       xlab = "window length",
       ylab = "MSE"
  )
  
  # highlight minimum and add minimum reference line 
  abline(h = row_ema[,min_MSE], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  
  points(x = row_sma[,min_l],
         y = row_sma[,min_MSE],
         col = "blue",
         bg = "blue",
         pch = 21)
  
  
  ## plot for ema ##
  
  plot(x = par_vec, 
       y = row_ema[,-c('min_MSE','min_a'),with = FALSE],
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("MSE for bias correction by EMA"),
       xlab = "scale parameter",
       ylab = "MSE"
  )
  
  # highlight minimum and add minimum reference line 
  abline(h = row_ema[,min_MSE], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  
  points(x = row_ema[,min_a],
         y = row_ema[,min_MSE],
         col = "blue",
         bg = "blue",
         pch = 21)

dev.off()


##################################################

############ Sec 4.1 #############################

##################################################



##################################################

################## PIT plots #####################


######## get the distribution fct. of a censored normal distribution #############

# value, mean and sd need to be vectors of equal length.
# returns F(value), where F is dist. fct. of a normal with parameters mean and sd, censored at trc_value

dist_fun_tn = function(value, mean, sd, trc_value = -1.79){ 
  a=rep(0,times = length(value))
  na_loc = which(is.na(value) | is.na(sd) | sd == 0)
  trc_loc = which(value <= trc_value & sd > 0)
  nor_loc = which(value > trc_value & sd > 0)
  
  a[na_loc] = NA
  a[trc_loc] = runif(length(trc_loc), max = pnorm(trc_value, mean = mean[trc_loc], sd = sd[trc_loc]))
  a[nor_loc] = pnorm(value[nor_loc], mean = mean[nor_loc], sd = sd[nor_loc])
  
  return(a)
}

########### get PITs ###############

DT_calib = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]

DT_calib[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

########### plot mean and standard deviation ################

plot_diagnostic(DT_calib[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_mean")

unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0, 2*unif_sd), set_white = unif_sd,
                mn = paste0("PIT standard deviation"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_sd")


###### get PIT for estimated mean by linear regression ######

DT_calib_3 = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,T_hat_lr_both,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_3[,"PIT_mc" := dist_fun_tn(SST_bar, mean = T_hat_lr_both, sd = SD_hat)]

DT_calib_3[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_3[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean, linear regression"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_mean_lr")


# the following plots the PIT standard deviation for mean estimation by linear regression. This is not in the paper but nevertheless interesting.

unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0,2*unif_sd),set_white = unif_sd,
                mn = paste0("PIT mean, linear regression"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_sd_lr")


######################################

########### Sec. 4.2. ################

######################################

name_abbr = "NAO/lv/2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


# reset plot directory

plot_dir = plot_dir0


####################################################################

##################### plot example residual ########################


# for showing fewer breaks than usual we define a variation on the plotting function 

library(fields,sp)

plot_smooth = function( dt, var = colnames(dt)[3], mn = var, rr = NULL,...,
                        theta = 0.5, pixels = 256,
                        col_scheme = "bwr", set_white = NULL,
                        xlab = "", ylab = "",
                        save_pdf = FALSE, save_dir = "./figures/", file_name = "diag_plot", stretch_par = NULL)
{
  # prepare data table
  
  if("year" %in% colnames(dt))
  {
    if("month" %in% colnames(dt))
    {
      dt = dt[year == min(year) & month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    } else {
      dt = dt[month == min(month),.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
    }
  } else {
    dt = dt[,.SD,.SDcols = c('Lon','Lat',var)][order(Lat,Lon)]
  }
  
  
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
  pts <- sp::SpatialPoints(all_loc, proj4string=sp::CRS(sp::proj4string(wrld_simpl)))
  ii <- !is.na(sp::over(pts, wrld_simpl)$FIPS)
  im_0$z[ii] = NA
  
  # --- fix range of plot and fill in values for points out of range ---
  
  if(is.null(rr))  rr = range(im_0$z,na.rm=TRUE)
  if(!is.null(rr)){
    im_0$z[im_0$z< min(rr)] = min(rr)
    im_0$z[im_0$z> max(rr)] = max(rr)
  }
  
  # --- scaling and colors ---
  
  brk = seq(rr[1],rr[2],length = 500)
  brk.ind = round(seq(1,length(brk),length = 5))
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
                     zlim=rr, main = mn,...,
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



m = 6
y = 2016

#in order to add route, we need to set the stretch parameter manually:

temp = DT[year == y & month == m,.(Lon,Lat,SST_bar - SST_hat)]

Lons = unique(temp[,Lon])
Lats = unique(temp[,Lat])

n_lon = length(Lons)
n_lat = length(Lats)

save_cex = par('cex')


pdf(paste0(plot_dir,'Example_res.pdf'),width = 7,height = 7 * n_lat/n_lon)

par('cex' = save_cex)
plot_smooth(DT[year == y & month == m,.(Lon,Lat,SST_bar - SST_hat)], 
            mn = latex2exp::TeX('observed forecast residual June 2016'), 
            rr = c(-3,3),
            pixels = 512,
            save_pdf = FALSE,
            save_dir = plot_dir,
            file_name = 'Example_res'
)  

# add Bordeaux and Norfolk and shipping route

Bordeaux = c(-0.57,44.8)
Norfolk = c(-76.3,36.9)

p1 = data.table(Lon = Bordeaux[1], Lat = Bordeaux[2], Loc = 'Bordeaux') 
p2 = data.table(Lon = Norfolk[1], Lat = Norfolk[2], Loc = "Norfolk") 

cities = rbindlist(list(p1,p2))

points(cities[,Lon],cities[,Lat],col="black", cex=2, pch=20)

par('cex' = save_cex)

# Connection between Bordeaux and Norfolk
inter <- geosphere::gcIntermediate(Bordeaux, Norfolk, n=100, addStartEnd=TRUE, breakAtDateLine=F)             
lines(inter, col="black", lwd=2)

dev.off()


##############################################

### multivariate rank histograms for route ###

# range of probability the probability plots
rr = c(0,100)

### average rank histograms ###

pdf(paste0(plot_dir,'avg_rhs_route.pdf'))

par(oma = c(1,1,1,1), mfrow=c(2,2), mar=c(2,1,2,1) )
par('cex' = 1.4,'cex.axis' = 0.75)

rhist_dt(rks_PCA_mc_route[,.(YM,av.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, $\\widehat{\\Sigma}^{mc}$')
)

rhist_dt(rks_PCA_ac_route[,.(YM,av.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, $\\widehat{\\Sigma}^{ac}$')
)

rhist_dt(rks_GS_route[,.(YM,av.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, GS')
)

rhist_dt(rks_ECC_route[,.(YM,av.rk.obs)], 
         max_rk = ens_size +1,
         breaks = breaks, 
         hist_xlab = "average",
         hist_ylim = rr,
         mn = latex2exp::TeX('average rank, ECC')
)



dev.off()

### band depth rank ###

### average rank histograms ###

pdf(paste0(plot_dir,'bd_rhs_route.pdf'))

par(oma = c(1,1,1,1), mfrow=c(2,2), mar=c(2,1,2,1) )
par('cex' = 1.4,'cex.axis' = 0.75)

rhist_dt(rks_PCA_mc_route[,.(YM,bd.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, $\\widehat{\\Sigma}^{mc}$')
)

rhist_dt(rks_PCA_ac_route[,.(YM,bd.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, $\\widehat{\\Sigma}^{ac}$')
)

rhist_dt(rks_GS_route[,.(YM,bd.rk.obs)], 
         max_rk = fc_ens_size +1,
         breaks = breaks, 
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, GS')
)

rhist_dt(rks_ECC_route[,.(YM,bd.rk.obs)], 
         max_rk = ens_size +1,
         breaks = breaks, 
         hist_ylim = rr,
         mn = latex2exp::TeX('band depth rank, ECC')
)



dev.off()


