
#######################################################################################

###################  master script part 3 - variance estimation  ######################

#######################################################################################

# This script estimates the marginal variance of the ensemble forecast by 
# a range of different methods and computes average scores for each method.
# It then selects and applies the optimal way of variance estimation and
# complements the data table DT by a new column SD_hat.
#
# 
# Files generated:
#   
# Data files: dt_combine_wide_bc_var.RData, scores.bc.sd.sma.Rdata, scores.bc.sd.ema.Rdata
# Plots: mean_scores_sd_sma.pdf, mean_scores_sd_ema.pdf 
#
# Requires previous run of 02.master.bias.correct 
# with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_small" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))
DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_wide_bc.RData"))


###### getting sample variances of ensemble  ######

DT = ens_sd_est(dt = DT,
                ens_size = ens_size,
                save_dir = save_dir,
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

save(sc_sma_var, file = paste0(save_dir,"/scores.bc.sd.sma.Rdata"))


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

save(sc_ema_var, file = paste0(save_dir,"/scores.bc.sd.ema.Rdata"))


#### plotting CRPS for different ways of variance estimation ####

load(paste0(save_dir,"/scores.bc.sd.sma.Rdata"))
load(paste0(save_dir,"/scores.bc.sd.ema.Rdata"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma_var[,CRPS],sc_ema_var[,CRPS]))  

## plot for sma ##

pdf(paste0(plot_dir,"/mean_scores_sd_sma.pdf"))
  plot(x = sc_sma_var[,win_length],
       y = sc_sma_var[,CRPS],
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("CRPS for ",name_abbr," SD estimation by SMA"),
       xlab = "window length",
       ylab = "variance score"
  )
  
  # highlight minimum and add minimum reference line 
  abline(h = sc_sma_var[,min(CRPS)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  min_loc_CRPS = sc_sma_var[,which.min(CRPS)]
  points(x = sc_sma_var[,win_length][min_loc_CRPS],
         y = sc_sma_var[,CRPS][min_loc_CRPS],
         col = "blue",
         bg = "blue",
         pch = 21)
  
dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"/mean_scores_sd_ema.pdf"))
  plot(x = sc_ema_var[,a],
       y = sc_ema_var[,CRPS],
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("CRPS for ",name_abbr," SD estimation by EMA"),
       xlab = "weight parameter a",
       ylab = "variance score"
  )
  
  # highlight minimum and add minimum reference line 
  abline(h = sc_ema_var[,min(CRPS)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  min_loc_CRPS = sc_ema_var[,which.min(CRPS)]
  points(x = sc_ema_var[,a][min_loc_CRPS],
         y = sc_ema_var[,CRPS][min_loc_CRPS],
         col = "blue",
         bg = "blue",
         pch = 21)

dev.off()

### finding optimal way of variance modelling ###

# if the optimal scores for simple moving averages and exponential moving acerages are essentially the same (less than 1 % difference),
# we pick simple moving averages.

if(sc_sma_var[,min(CRPS)] < 1.01*sc_ema_var[,min(CRPS)]){
  print(paste0("optimal variance estimation uses simple moving average over sample variances with window length of ",sc_sma_var[,which.min(CRPS)], " years, and achieves a CRPS of ",round(sc_sma_var[,min(CRPS)],3),
               ". Best estimation with exponentially weighted sample variance achieves a RMSE of ",round(sc_ema_var[,min(CRPS)],3),"."))
  opt_par = c("sma",sc_sma_var[,which.min(CRPS)])
} else{
  print(paste0("optimal variance estimation uses exponentially weighted sample variance with parameter a = ",round(sc_ema_var[,a][sc_ema_var[,which.min(CRPS)]],3),
               ", and achieves a CRPS of ",round(sc_ema_var[,min(CRPS)],3),". Best estimation with simple moving averages achieves a CRPS of ",round(sc_sma_var[,min(CRPS)],3),"."))
  opt_par = c("ema",sc_ema_var[,a][sc_sma_var[,which.min(CRPS)]])
}

### variance estimation ###

sd_est(dt = DT,
       method = opt_par[1],
       par_1 = as.double(opt_par[2]),
       save_dir =save_dir,
       file_name = paste0("dt_combine_wide_bc_var.RData")
)

