

#############################################################################################

###### BIG test script for a bunch of models #########

#############################################################################################

# This script tests variants of the PCA method on three different 

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

# choose your favourite area for analysis and give it a name abbreviation

#NAO_2:

lat_box = c(40,50)
lon_box = c(-60,-50)

# lat_box = c(40,70)
# lon_box = c(-60,-30)

name_abbr = "NAO_2" # for northern atlantic ocean

ens_size = 9 # size of forecast ensemble

validation_years = 2001:2010 # all previous years are used for training 
months = 1:12

# create directories

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)

### construct or load wide data set ###

DT = load_combined_wide()[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]


training_year_index = !(DT[,unique(year)] %in% validation_years) 
training_years = DT[,unique(year)][training_year_index]



##########################################################
###### run bias analysis for simple moving averages ######
##########################################################


num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1

sc_sma = list()
dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "sma", 
                      par_1 = k,
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE)
  sc_sma[[k]] = temp[,"win_length" := k]
}

sc_sma = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = 8)
sc_sma = rbindlist(sc_sma)

save(sc_sma, file = paste0(save_dir,"scores.bc.sma.RData"))

###### run bias analysis for exponential moving averages ######

par_vec = seq(0.05,0.4,length.out = 24) 

sc_ema = list()

dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "ema",
                      par_1 = par_vec[k], 
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE)
  sc_ema[[k]] = temp[,"a" := par_vec[k]]
}

sc_ema = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = 8)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save_dir,"scores.bc.ema.RData"))

###### plotting scores for different ways of bias correction ######

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma[,sqrt(MSE)],sc_ema[,sqrt(MSE)]))  

## plot for sma ##

pdf(paste0(plot_dir,"mean_scores_sma.pdf"))
plot(x = sc_sma[,win_length],
     y = sc_sma[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for bias correction by SMA"),
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

pdf(paste0(plot_dir,"mean_scores_ema.pdf"))
plot(x = sc_ema[,a],
     y = sc_ema[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for bias correction by EMA"),
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





###### finding and applying optimal way of bias correcion ######

if(sc_sma[,min(MSE)] < sc_ema[,min(MSE)]){
  opt_par = c("sma",sc_sma[,which.min(MSE)])
} else{
  opt_par = c("ema",sc_ema[,a][sc_sma[,which.min(MSE)]])
}

DT = bias_correct(dt = DT,
                  method = opt_par[1],
                  par_1 = as.double(opt_par[2]),
                  save_dir = save_dir)


DT = bias_correct_training(dt = DT,
                           method = opt_par,
                           training_years = training_years,
                           save_dir = save_dir)

###### compare to linear regression models ######

RMSE_linear_models = bias_lr(DT,validation_years = validation_years)

RMSE_linear_models



#######################################################################################

variance_methods = c("sv","bcf")


for(i in 1:length(variance_methods))
{
  
  DT = ens_sd_est(dt = DT,
                  ens_size = ens_size,
                  save_dir = save_dir,
                  mean_est = variance_methods[i],
                  file_name = paste0("dt_combine_wide_bc_var.RData"))
  
  
  ###### getting scores for simple moving averages ######
  
  num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1
  
  sc_sma_var = list()
  dummy_function = function(k){
    temp = sd_est(dt = DT, 
                  method = "sma", 
                  par_1 = k,
                  scores = TRUE,
                  eval_years = validation_years,
                  saveorgo = FALSE)
    sc_sma_var[[k]] = temp[,"win_length" := k]
  }
  
  sc_sma_var = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = 8)
  sc_sma_var = rbindlist(sc_sma_var)
  
    ###### getting scores for exponential moving averages ######
  
  par_vec = seq(0.01,0.4,length.out = 24)
  
  sc_ema_var = list()
  
  dummy_function = function(k){
    temp = sd_est(dt = DT[year > min(year),], 
                  method = "ema",
                  par_1 = par_vec[k], 
                  scores = TRUE,
                  eval_years = validation_years,
                  saveorgo = FALSE)
    sc_ema_var[[k]] = temp[,"a" := par_vec[k]]
  }
  
  sc_ema_var = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = 8)
  sc_ema_var = rbindlist(sc_ema_var)
  
  
  
  ### finding optimal way of variance modelling ###
  
  if(sc_sma_var[,min(CRPS)] < sc_ema_var[,min(CRPS)]){
    opt_par = c("sma",sc_sma_var[,which.min(CRPS)])
  } else{
    opt_par = c("ema",sc_ema_var[,a][sc_sma_var[,which.min(CRPS)]])
  }
  
  ############# saving #####################
  
  save(sc_sma_var,sc_ema_var,opt_par, file = paste0(save_dir,"scores.bc.sd.",variance_methods[i],".Rdata"))
  
  
}


############ plots #################

#### plotting CRPS for different ways of variance estimation ####

sc_sma = data.table(win_length = 1:(num_years-1))
sc_ema = data.table(a = par_vec)

for(i in 1:length(variance_methods))
{
  var_met = variance_methods[i]
  
  load(file = paste0(save_dir,"scores.bc.sd.",var_met,".Rdata"))
  sc_sma = merge(sc_sma,sc_sma_var)
  setnames(sc_sma,"CRPS",paste0("CRPS_",var_met))
  
  sc_ema = merge(sc_ema,sc_ema_var)
  setnames(sc_ema,"CRPS",paste0("CRPS_",var_met))
  
  # ensure that they are plotted on the the same range
  
}


##### CRPS plots for different versions of marginal variance estimation #####

col_vec = c("blue","darkgreen","darkred","black","yellow2")
leg_vec = variance_methods

y_range = range(sc_sma[,2:( 1+length(variance_methods))], sc_ema[,2:( 1+length(variance_methods))])  


# sma plot

pdf(paste0(plot_dir,"/mean_scores_sd_sma.pdf"))
  plot(sc_sma[,.SD,.SDcols = c("win_length",paste0("CRPS_",variance_methods[1]))],
       ylim = y_range,
       type = "l",
       col = col_vec[1],
       main = paste0("CRPS for ",name_abbr," SD estimation by SMA"),
       xlab = "window length",
       ylab = "variance score"
  )
  
  # highlight minimum 
  min_loc_CRPS = which.min(sc_sma[,eval(parse(text = paste0("CRPS_",variance_methods[1])))])
  points(x = sc_sma[,win_length][min_loc_CRPS],
         y = sc_sma[,2][min_loc_CRPS],
         col = col_vec[1],
         bg = col_vec[1],
         pch = 21)
  
  for(i in 2:length(variance_methods))
  {
    lines(sc_sma[,.SD,.SDcols = c("win_length",paste0("CRPS_",variance_methods[i]))],
          type = "l",col = col_vec[i]) 
    
    min_loc_CRPS = which.min(sc_sma[,eval(parse(text = paste0("CRPS_",variance_methods[i])))])
    points(x = sc_sma[,win_length][min_loc_CRPS],
           y = sc_sma[,1+i, with = FALSE][min_loc_CRPS],
           col = col_vec[i],
           bg = col_vec[i],
           pch = 21)
  }
  
  # highlight minimum and add minimum reference line 
  abline(h = min(sc_sma[,-1, with = FALSE]), lty = "dashed", col = adjustcolor("black",alpha = .5))
  
  legend("topright",legend = leg_vec,lty = 1,col = col_vec[1:length(variance_methods)])
  
  dev.off()
  

# ema plot
  
  
  pdf(paste0(plot_dir,"/mean_scores_sd_ema.pdf"))
  plot(sc_ema[,.SD,.SDcols = c("a",paste0("CRPS_",variance_methods[1]))],
       ylim = y_range,
       type = "l",
       col = col_vec[1],
       main = paste0("CRPS for ",name_abbr," SD estimation by EMA"),
       xlab = "a",
       ylab = "CRPS"
  )
  
  # highlight minimum 
  min_loc_CRPS = which.min(sc_ema[,eval(parse(text = paste0("CRPS_",variance_methods[1])))])
  points(x = sc_ema[,a][min_loc_CRPS],
         y = sc_ema[,2][min_loc_CRPS],
         col = col_vec[1],
         bg = col_vec[1],
         pch = 21)
  
  for(i in 2:length(variance_methods))
  {
    lines(sc_ema[,.SD,.SDcols = c("a",paste0("CRPS_",variance_methods[i]))],
          type = "l",col = col_vec[i]) 
    
    min_loc_CRPS = which.min(sc_ema[,eval(parse(text = paste0("CRPS_",variance_methods[i])))])
    points(x = sc_ema[,a][min_loc_CRPS],
           y = sc_ema[,1+i, with = FALSE][min_loc_CRPS],
           col = col_vec[i],
           bg = col_vec[i],
           pch = 21)
  }
  
  # highlight minimum and add minimum reference line 
  abline(h = min(sc_ema[,-1, with = FALSE]), lty = "dashed", col = adjustcolor("black",alpha = .5))
  
  legend("topright",legend = leg_vec,lty = 1,col = col_vec[1:length(variance_methods)])
  
  dev.off()
  
  
  
