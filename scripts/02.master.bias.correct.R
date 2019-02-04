################################################################################

###################  master script part 2 - bias correct  ######################

################################################################################

# This script computes average scores for a range of different methods 
# of bias correction. It then chooses and applies the optimal way of 
# bias correction. It complements the data table DT by two new 
# columns Bias_Est and SST_hat (which is just Ens_bar + Bias_Est, truncated
# at freezing temperature).
#
# 
# Files generated:
#   
# Data files: "scores.bc.sma.RData", "scores.bc.ema.RData", "dt_combine_wide_bc.RData"
# Plots: "mean_scores_sma.pdf", "mean_scores_ema.pdf"
#
# Requires previous run of 01.master.setup.R with the same value of name_abbr as below.



##### setting up ######

rm(list = ls())

time_s2 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "test" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))




###### run bias analysis for simple moving averages ######

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

sc_sma = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = mc_cores)
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

sc_ema = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = mc_cores)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save_dir,"scores.bc.ema.RData"))


###### plotting scores for different ways of bias correction ######

load(paste0(save_dir,"scores.bc.sma.RData"))
load(paste0(save_dir,"scores.bc.ema.RData"))

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


###### finding and applying optimal way of bias correcion, exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

if(sc_sma[,min(MSE)] < 0.9 * sc_ema[,min(MSE)]){
  opt_par = c("sma",sc_sma[,win_length][,which.min(MSE)])
} else{
  opt_par = c("ema",sc_ema[,a][sc_ema[,which.min(MSE)]])
}

DT = bias_correct(dt = DT,
                  method = opt_par[1],
                  par_1 = as.double(opt_par[2]),
                  save_dir = save_dir)

# For the training years the bias correction considers also the future

DT = bias_correct_training(dt = DT,
                           method = opt_par,
                           training_years = training_years,
                           save_dir = save_dir)



time_s2 = proc.time() - time_s2

#### save ####

save.image(file = paste0(save_dir,"setup.RData"))
