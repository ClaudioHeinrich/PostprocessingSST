
###########################################################################

############## script for issuing the summer forecast 2018 ################

###########################################################################


# This script is a version of the current master.NAO script with some stuff left out and some minor modifications.
# Models are trained on April vintage instead of most recent and have a gap in the data of 7 years prior to the forecast, which we take into account for bias modelling.
# PCA uses all data available until 2010.


rm(list = ls())

library(PostProcessing)
library(data.table)
library(parallel)
library(fExtremes)
library(fields)
library(vegan)

setwd("~/NR/SFE")
options(max.print = 1e3)

# box for analysis, name abbreviation for saved files and directories

lat_box = c(30,70)
lon_box = c(-60,20)

name_abbr = "Apr18" 

save.dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr)
dir.create(save.dir, showWarnings = FALSE)

plot.dir = paste0("./figures/", name_abbr)
dir.create(plot.dir, showWarnings = FALSE)


########################################
### construct and load wide data set ###
########################################

# this function also attaches the forecast issued in March 2018:

make_combined_wide_dataset_SFc2018 = function(y_start = 1985,
                                      y_stop = 2010,
                                      vintage = "Apr",
                                      data.dir = "~/PostClimDataNoBackup/SFE/FcApr2018/",
                                      grid_mapping_loc = "~/PostClimDataNoBackup/SFE/Derived/",
                                      output_loc = save.dir, 
                                      output_name = paste0("dt_combine_",name_abbr,"_wide.RData"),
                                      lat_box = lat_box,
                                      lon_box = lon_box)
{
  
  ##----- Load Grid Mapping ---
  ff = paste0(grid_mapping_loc,"dt_map.RData")
  if(file.exists(ff))
  {
    load(ff)
    names(dt_map)= c("Lon_Obs","Lat_Obs","Lon_Ens","Lat_Ens") ## Get rid of this eventually.
  }else{
    stop("Could not find grid mapping info")
  }
  ##--------------------------
  
  ##------ Loop ----------
  dt_combine_all = list()
  k = 1
  for(y in y_start:y_stop)
  {
    for(m in 1:12)
    {
      print(c(y,m))
      dt_ens = load_ensemble(y,m,vintage,data.dir = data.dir)
      dt_obs = load_observations(y,m)
      dt_combine_all[[k]] = combine_data_wide(dt_ens, dt_obs, dt_map)
      dt_combine_all[[k]][,year:=y]
      dt_combine_all[[k]][,month:=m]
      k = k + 1
    }
  }
  
  # attach 2018 forecast
  
  y = 2018
  for(m in 3:12)
  {
    print(c(y,m))
    dt_ens = load_ensemble(y,m,vintage,data.dir = data.dir)
    dt_obs = load_observations(2010,1) # hacked - we don't have an observation for 2018, so we just attach any observation
    dt_combine_all[[k]] = combine_data_wide(dt_ens, dt_obs, dt_map)
    dt_combine_all[[k]][,year:=y]
    dt_combine_all[[k]][,month:=m]
    k = k + 1
  }
  
  ##------------------------
  
  ##--------- Combine -----
  dt = rbindlist(dt_combine_all)
  dt[, YM := year * 12 + month]
  setkey(dt, "YM", "Lon", "Lat")
  ##------------------------
  
  ##---- Restrict ---
  dt = dt[ (Lon >= lon_box[1]) & (Lon <= lon_box[2]) & (Lat >= lat_box[1]) & (Lat <= lat_box[2])]
  
  ##----- Should I save or should I go? -----
  if(is.null(data.dir))
  {
    return(dt)
  }else{
    if(is.null(output_name))
    {
      output_name = paste0("dt_combine_",vintage,"_wide.RData")
    }
    f_name = paste0(output_loc,"/",output_name)
    save(dt,file = f_name)
    return(1)
  }
  ##-------------------------------------------
  
  
}


make_combined_wide_dataset_SFc2018()

DT = load_combined_wide(data.dir = save.dir, output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))



########################################
####### univariate modelling ###########
########################################


############### bias ###################


### run bias analysis ###

validation_years = 2001:2010

# for simple moving averages

num.years = DT[,range(year)][2] - DT[,range(year)][1] + 1

sc_sma = list()
dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "sma", 
                      par_1 = k,
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE,
                      skip = 7)
  sc_sma[[k]] = temp[,"win_length" := k]
}

sc_sma = mclapply(X = 1:(num.years-1), FUN = dummy_function, mc.cores = 8)
sc_sma = rbindlist(sc_sma)

save(sc_sma, file = paste0(save.dir,"/scores.bc.sma.Rdata"))

# for exponential moving averages

par.vec = seq(0.05,0.4,length.out = 24)

sc_ema = list()

dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "ema",
                      par_1 = par.vec[k], 
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE,
                      skip = 7)
  sc_ema[[k]] = temp[,"a" := par.vec[k]]
}

sc_ema = mclapply(X = 1:length(par.vec), FUN = dummy_function,mc.cores = 8)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save.dir,"/scores.bc.ema.Rdata"))


#### plotting scores for different ways of bias correction ####

load(paste0(save.dir,"/scores.bc.sma.Rdata"))
load(paste0(save.dir,"/scores.bc.ema.Rdata"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma[,sqrt(MSE)],sc_ema[,sqrt(MSE)]))  

## plot for sma ##

pdf(paste0(plot.dir,"/mean_scores_sma.pdf"))
plot(x = sc_sma[,win_length],
     y = sc_sma[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for ",name_abbr," bias correction by SMA"),
     xlab = "window length",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = sc_sma[,min(sqrt(MSE))], lty = "dashed", col = adjustcolor("blue",alpha = .5))

min_loc_RMSE = sc_sma[,which.min(MSE)]

points(x = sc_sma[,win_length][min_loc_RMSE],
       y = sc_sma[,sqrt(MSE)][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)

dev.off()

## plot for ema ##

pdf(paste0(plot.dir,"/mean_scores_ema.pdf"))
plot(x = sc_ema[,a],
     y = sc_ema[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for ",name_abbr," bias correction by EMA"),
     xlab = "weight parameter a",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = sc_ema[,min(sqrt(MSE))], lty = "dashed", col = adjustcolor("blue",alpha = .5))

min_loc_RMSE = sc_ema[,which.min(MSE)]

points(x = sc_ema[,a][min_loc_RMSE],
       y = sc_ema[,sqrt(MSE)][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)

dev.off()

### finding optimal way of bias correcion ###

if(sc_sma[,min(MSE)] < sc_ema[,min(MSE)]){
  print(paste0("optimal bias correction uses simple moving averages with window length of ",sc_sma[,which.min(MSE)], " years, and achieves a RMSE of ",round(sc_sma[,sqrt(min(MSE))],3),
               ". Best correction with exponential moving averages achieves a RMSE of ",round(sc_ema[,sqrt(min(MSE))],3),"."))
  opt_par = c("sma",sc_sma[,which.min(MSE)])
} else{
  print(paste0("optimal bias correction uses exponential moving averages with parameter a = ",round(sc_ema[,a][sc_sma[,which.min(MSE)]],3),
               ", and achieves a RMSE of ",round(sc_ema[,sqrt(min(MSE))],3),". Best correction with simple moving averages achieves a RMSE of ",round(sc_sma[,sqrt(min(MSE))],3),"."))
  opt_par = c("ema",sc_ema[,a][sc_sma[,which.min(MSE)]])
}

### bias correction ###

bias_correct(dt = DT,
             method = opt_par[1],
             par_1 = as.double(opt_par[2]),
             save.dir = paste0(save.dir,"/"),
             file.name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"),
             skip = 7)

DT = load_combined_wide(data.dir = save.dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"))




#################################################
########## multivariate modelling ###############
#################################################


###################  PCA  #######################

training_years = 1985:2010
eval_years = 2018

ens.size = 9

##### compute and save principal components #####

cov.dir = paste0(save.dir,"/PCACov")
dir.create(cov.dir, showWarnings = FALSE)

for_res_cov(Y = training_years,dt = DT, save.dir = cov.dir,ens.size = ens.size)

############## set up PCA #######################

setup_PCA(dt = DT, y = eval_years, cov.dir = cov.dir)

#### a brief example ####

ex_depth = 15
ex_month = 1
ex_year = 2001

result = forecast_PCA(y=ex_year,m=ex_month,PCA_depth = ex_depth, saveorgo = FALSE)

plot_diagnostic(result[,.(Lon,Lat,noise)], mn = paste0("Noise for ",ex_month,"/",ex_year," generated by ",ex_depth," PCs"))




###### Example plots of forecasted SST and anomalies w.r.t climatology #######

opt_num_PCs = 17 # This is an educated guess: for 16 years training the results are best for 11 PCs and very good for 9 to 17 PCs, now we have 26 years of training.

ex_depth = opt_num_PCs
ex_months = 4:9
ex_month_names = c("April","May","June","July","August","September")
ex_year = 2018

clim_years = 1985:2010 # the years to compute the climatology from

MC_sample_size = 30   # number of plots with independently generated noise

ens.size = 9   # size of forecast ensemble

PCs = opt_num_PCs      # number of considered principal components


#compute climatology
climatology = DT[year %in% clim_years, clim := mean(SST_bar),by = .(grid_id,month)][year == min(year) ,.(Lon,Lat,month,clim)]


for(m in ex_months){
  print(paste0("month = ",m))
  
  #generate noise:
  no.dt = list()
  for(i in 1:MC_sample_size){
    no.dt[[i]] = forecast_PCA(m = m, y = ex_year, PCA_depth = PCs, saveorgo = FALSE)[,.(Lon,Lat,noise), keyby = .(Lon,Lat)][,noise]
  }
  no.dt = as.data.table(no.dt)
  setnames(no.dt,paste0("no",1:MC_sample_size))
  
  DT_pca_plot = no.dt[,c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens.size)) := 
                        DT[year == ex_year & month == m, c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens.size)),
                           with = FALSE]]
  DT_pca_plot[,clim := climatology[month == m, clim]]
  
  # choose random ensemble members (REM) and generate forecast as REM + bias + noise
  
  ens_mem = sample.int(ens.size,MC_sample_size,replace = TRUE)
  for(i in 1:MC_sample_size){
    dummy_dt = DT_pca_plot[,.SD,.SDcols = c(paste0("no",i),paste0("Ens",ens_mem[i]),"Bias_Est")]
    forecast = dummy_dt[[1]] + dummy_dt[[2]] + dummy_dt[[3]]
    DT_pca_plot = DT_pca_plot[,paste0("fc",i) := forecast]
  }
  
  #forecast plots:
  
  rr_sst = range(na.omit(DT_pca_plot[,.SD,.SDcols = c(paste0("fc",1:MC_sample_size))]))
  
  for(i in 1:MC_sample_size){
    plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i))))],
                    rr = rr_sst,
                    mn = paste0("SST forecast for ",ex_month_names[which(ex_months == m)]),
                    save.pdf = TRUE, 
                    save.dir = paste0(plot.dir,"/"),
                    file.name = paste0("m",m,"_fc",i),
                    stretch_par = .8)
  }
  
  #anomaly plots
  rr_clim = range(na.omit(DT_pca_plot[,.SD - clim,.SDcols = c(paste0("fc",1:MC_sample_size))]))
  
  for(i in 1:MC_sample_size){
    plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i)))-clim)],
                    rr = rr_clim,
                    mn = paste0("Anomaly forecast for ",ex_month_names[which(ex_months == m)]),
                    save.pdf = TRUE, 
                    save.dir = paste0(plot.dir,"/"),
                    file.name = paste0("m",m,"_afc",i),
                    stretch_par = .8)
  }
  
}



