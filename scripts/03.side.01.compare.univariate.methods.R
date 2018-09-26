###############################################################################

#############  side script 3.1 - comparing univariate methods  ################

###############################################################################

# This script compares the univariate modelling of bias and variance by moving averages to the 
# most commonly used methods in NGR: for bias correcting this is modelling the bias as a + b x,
# where x is the estimated bias, for variance estimation it is c^2 + d^2 SD, where SD is the ensemble spread.
# 
# Files generated: univ_scores_comp_NGR.RData
#   
# Requires previous run of 03.master.var.est.R with the same value of name_abbr as below.



##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Atl" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


###### compare to linear regression models ######

RMSE_linear_models = bias_lr(DT,validation_years = validation_years)

# compare to simple moving averages:
load(paste0(save_dir,"scores.bc.sma.RData"))
load(paste0(save_dir,"scores.bc.ema.RData"))

RMSE_sma = sc_sma[,sqrt(min(MSE))]
RMSE_ema = sc_ema[,sqrt(min(MSE))]

RMSE_linear_models[,c("RMSE_sma",'RMSE_ema'):=.(RMSE_sma,RMSE_ema)]

RMSE_linear_models = round(RMSE_linear_models,4)

###### Compare to SD estimation in NGR fashion: #####

DT = var_est_NGR(DT, months = months, validation_years = validation_years)

load(file = paste0(save_dir,"scores.bc.sd.sma.Rdata"))
load(file = paste0(save_dir,"scores.bc.sd.ema.Rdata"))
  
CRPS_sma = sc_sma_var[,min(CRPS)]
CRPS_ema = sc_ema_var[,min(CRPS)]

CRPS_comparison = data.table(mean_CRPS_bm = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bm),na.rm = TRUE)],
                             mean_CRPS_bl = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bl),na.rm = TRUE)],
                             mean_CRPS_bb = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bb),na.rm = TRUE)],
                             mean_CRPS_sma = CRPS_sma,
                             mean_CRPS_ema = CRPS_ema)

CRPS_comparison = round(CRPS_comparison,3)

save(RMSE_linear_models,CRPS_comparison,file = paste0(save_dir,'univ_scores_comp_NGR.RData'))
