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
                                   type = "res",  #type of uncertainty takes 'ens','obs', 'comb' or 'res'
                                   file_out = paste0("./figures/sd_by_loc_",type),
                                   lons = NULL,
                                   lats = NULL,
                                   rr = NULL,
                                   outside_control = FALSE,
                                   print_figs = TRUE){
  

  if(is.null(dt_sd))  load(file = "~/PostClimDataNoBackup/SFE/Derived/dtcombine_mr_wide_sd.RData")

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
  
  
  n_lon = length(Lons)
  n_lat = length(Lats)
  
  dt_sd = dt_sd[Lon %in% Lons & Lat %in% Lats]
  
  #-----
  
  
  if(type == "comb") type_id = ""
  if(type == "ens") type_id = "ens_"
  if(type == "obs") type_id = "obs_"
  if(type == "res") type_id = "res_"
  
  type_id = paste0(type_id,"sd_by_loc")
  
  
  for(m in M){
    print(paste0("Month =",m))
  
    mn = paste0("residual SD for month ",m)
    if(is.null(rr)) rr = range(dt_sd[month == m,eval(parse(text = type_id))],na.rm=TRUE)
  
  
  #-------get sd as matrix---------------

  A_sd =  dt_sd[month == m & year == 1985, eval(parse(text = type_id))]
  A_sd = matrix(A_sd, nrow = n_lon, byrow = TRUE)
##-------------------------

##------- Scaling ----------
brk = seq(rr[1],rr[2],length = 500)
brk.ind = round(seq(1,length(brk),length = 10))
brk.at = brk[brk.ind]
brk.lab = round(brk[brk.ind],2)
color <- rgb(139,0,0, # is the rgb specification for darkred, relative to maxColor = 255
             alpha = seq(0,255,,length(brk)-1),
             maxColorValue = 255)
##--------------------------

##------- Plot -----------
if(print_figs ) pdf(paste0(file_out,"_m",m, ".pdf"))

image.plot.na(Lons,Lats,A_sd,
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


}

}

#--------

global_uncertainty_plot()



get_ranges = function(M=1:12,
                      file_out = "~/PostClimDataNoBackup/SFE/Derived/range_sd_res.RData"
                      ){
  # saves the ranges of the standard deviation by month to allow direct comparison with PCA

  if(is.null(dt_sd))  load(file = "~/PostClimDataNoBackup/SFE/Derived/dtcombine_mr_wide_sd.RData")
  
  
  rr = dt_sd[,.(month,res_sd_by_loc)]
  rr = rr[, .(max_sd_res = max(res_sd_by_loc,na.rm = TRUE),
         min_sd_res = min(res_sd_by_loc,na.rm = TRUE)), month]   
  
  save(rr,file = file_out)
  }
  
  