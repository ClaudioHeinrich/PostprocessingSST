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

time_s31 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


###### compare to linear regression models ######

# get estimates, grouped by month, location and both
bias_lr_bm(DT,months = months,validation_years = validation_years)
bias_lr_bl(DT,months = months,validation_years = validation_years)
bias_lr_bb(DT,months = months,validation_years = validation_years)

# get RMSEs
RMSE_lr_m = sqrt(DT[year %in% validation_years,mean((T_hat_lr_m-SST_bar)^2, na.rm = TRUE)])
RMSE_lr_loc = sqrt(DT[year %in% validation_years,mean((T_hat_lr_loc-SST_bar)^2, na.rm = TRUE)])
RMSE_lr_both = sqrt(DT[year %in% validation_years,mean((T_hat_lr_both-SST_bar)^2, na.rm = TRUE)])

# compare to simple moving averages:
load(paste0(save_dir,"scores.bc.sma.RData"))
load(paste0(save_dir,"scores.bc.ema.RData"))

RMSE_sma = sc_sma[,sqrt(min(MSE))]
RMSE_ema = sc_ema[,sqrt(min(MSE))]


RMSE_linear_models = data.table(RMSE_lr_m = RMSE_lr_m, 
                                RMSE_lr_loc = RMSE_lr_loc, 
                                RMSE_lr_both = RMSE_lr_both,
                                RMSE_sma = RMSE_sma,
                                RMSE_ema = RMSE_ema)

RMSE_linear_models = round(RMSE_linear_models,4)

###### Compare to SD estimation in NGR fashion: #####

DT = var_est_NGR_bm(DT, months = months, validation_years = validation_years)
mean_CRPS_bm = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bm),na.rm = TRUE)]

DT = var_est_NGR_bl(DT, months = months, validation_years = validation_years)
mean_CRPS_bl = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bl),na.rm = TRUE)]

DT = var_est_NGR_bb(DT, months = months, validation_years = validation_years)
mean_CRPS_bb = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bb),na.rm = TRUE)]



load(file = paste0(save_dir,"scores.bc.sd.sma.Rdata"))
load(file = paste0(save_dir,"scores.bc.sd.ema.Rdata"))
  
CRPS_sma = sc_sma_var[,min(CRPS)]
CRPS_ema = sc_ema_var[,min(CRPS)]

CRPS_comparison = data.table(mean_CRPS_bm = mean_CRPS_bm,
                             mean_CRPS_bl = mean_CRPS_bl,
                             mean_CRPS_bb = mean_CRPS_bb,
                             mean_CRPS_sma = CRPS_sma,
                             mean_CRPS_ema = CRPS_ema)

CRPS_comparison = round(CRPS_comparison,5)



### permutation test for moving average vs lr_bb ###

perm_test_dt = DT[year %in% validation_years & month %in% months,.(year,month,SST_bar,SST_hat,SD_hat,SD_hat_lr_bb)]

# getting CRPSs

perm_test_dt[,CRPS_ma := crps_na_rm(SST_bar,SST_hat,SD_hat)]
perm_test_dt[,CRPS_lr_bb := crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bb)]

# permutation test for CRPS_ma ~ CRPS_lr_bb
pt_CRPS = permutation_test_difference(na.omit(perm_test_dt[,CRPS_ma]),na.omit(perm_test_dt[,CRPS_lr_bb]), N = 500  )



time_s31 = proc.time() - time_s31

save(RMSE_linear_models,CRPS_comparison,file = paste0(save_dir,'univ_scores_comp_NGR.RData'))

save.image(file = paste0(save_dir,"setup.RData"))
