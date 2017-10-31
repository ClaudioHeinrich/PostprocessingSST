rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)



#------ plot uncertainty by location ------------


##----Load globals---

dt = bias_correct()
lon_all = sort(dt[,unique(Lon)])
lat_all = sort(dt[,unique(Lat)])
n_lon = length(lon_all)
n_lat = length(lat_all)
##------------------


##------- compute standard deviations by month and location --------

sddt = function(dt){
  return(sd(as.matrix(dt),na.rm = TRUE))
}
dt_sd = dt[,obs_sd_by_loc := lapply(.SD[,paste0("SST",1:10)],sddt), by =.(month,Lon,Lat)]
dt_sd[,ens_sd_by_loc := lapply(.SD[,paste0("Ens",1:9)],sddt), by = .(month,Lon,Lat)]
dt_sd[,sd_by_loc := sqrt(ens_sd_by_loc^2 + obs_sd_by_loc^2)] # perhaps not very informative?

#--- Add residuals (w.r.t. the mean observation) and sd of residuals by month and location for the residuals---

for(k in 1:9){
  ens_lab = paste0("Ens",k)
  dt_sd[, paste0("Res",k) := eval(parse(text = ens_lab)) + Bias_Est - SST_bar]
}

dt_sd[,res_sd_by_loc := lapply(.SD[,paste0("Res",1:9)],sddt), by = .(month,Lon,Lat)]

# This all takes forever, so better save it:
save(dt_sd, file = "~/PostClimDataNoBackup/SFE/Derived/dtcombine_mr_wide_sd.RData")



global_uncertainty_plot = function(dt_sd = NULL,
                                   M=1:12,
                                   type = "comb",  #type of uncertainty takes 'ens','obs', 'comb' or 'res'
                                   file_out = paste0("./figures/sd_by_loc_",type),
                                   lons = NULL,
                                   lats = NULL,
                                   rr = NULL,
                                   outside_control = FALSE,
                                   print_figs = TRUE,
                                   png_out = FALSE){
  

  if(is.null(dt_sd))  load(file = "~/PostClimDataNoBackup/SFE/Derived/dtcombine_mr_wide_sd.RData")

  if(type == "comb") type_id = ""
  if(type == "ens") type_id = "ens_"
  if(type == "obs") type_id = "obs_"
  if(type == "res") type_id = "res_"
  
  type_id = paste0(type_id,"sd_by_loc")
  
  lon_all = sort(dt_sd[,unique(Lon)])
  lat_all = sort(dt_sd[,unique(Lat)])
  n_lon = length(lon_all)
  n_lat = length(lat_all)
  
  for(m in M){
    print(paste0("Month =",m))
  
    mn = paste0("sample variance for ",type,", month",m)
    if(is.null(rr)) rr = range(dt_sd[month == m,eval(parse(text = type_id))],na.rm=TRUE)
  
  
  #-------get sd as matrix---------------


  A_sd = matrix(NA, n_lon, n_lat)
  A_sd[dt_sd[month == m & year == 1985, grid_id]] = 
    dt_sd[month == m & year == 1985, eval(parse(text = type_id))]
##-------------------------

##------- Scaling ----------
brk = seq(rr[1],rr[2],length = 500)
brk.ind = round(seq(1,length(brk),length = 10))
brk.at = brk[brk.ind]
brk.lab = round(brk[brk.ind],2)
color <- designer.colors(n=length(brk)-1)
##--------------------------

##------- Scope ------------
if(is.null(lons))lons = range(dt_sd[,Lon])
if(is.null(lats))lats = range(dt_sd[,Lat])
##--------------------------

##------- Plot -----------
if(print_figs & !outside_control)
{
  pdf(paste0(file_out,"_m",m, ".pdf"))
}
else{
  if(!outside_control) X11()
}

image.plot(lon_all,lat_all,A_sd,
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
    png(paste0(file_out,"_m",m, ".png"))
    image.plot(lon_all,lat_all,A_sd,
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
##------------------------

}

}
