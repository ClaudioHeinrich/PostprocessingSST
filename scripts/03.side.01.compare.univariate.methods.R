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

name_abbr = 'Full/lv' 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


###### compare to linear regression models ######


# get estimates, grouped by month, location and both
DT = bias_lr_bm(DT,months = months,validation_years = validation_years)
DT = suppressWarnings(bias_lr_bl(DT,months = months,validation_years = validation_years)) # linear regressions at grid points with no data causes warning
DT = suppressWarnings(bias_lr_bb(DT,months = months,validation_years = validation_years))

# get RMSEs
RMSE_lr_m = sqrt(DT[year %in% validation_years,mean((T_hat_lr_m-SST_bar)^2, na.rm = TRUE)])
RMSE_lr_loc = sqrt(DT[year %in% validation_years,mean((T_hat_lr_loc-SST_bar)^2, na.rm = TRUE)])
RMSE_lr_both = sqrt(DT[year %in% validation_years,mean((T_hat_lr_both-SST_bar)^2, na.rm = TRUE)])

# compare to simple moving averages:

RMSE_sma = msc_sma[,sqrt(mean(min_MSE))]
RMSE_ema = msc_ema[,sqrt(mean(min_MSE))]


RMSE_linear_models = data.table(RMSE_lr_m = RMSE_lr_m, 
                                RMSE_lr_loc = RMSE_lr_loc, 
                                RMSE_lr_both = RMSE_lr_both,
                                RMSE_sma = RMSE_sma,
                                RMSE_ema = RMSE_ema)

RMSE_linear_models = round(RMSE_linear_models,4)

###### Compare to SD estimation in NGR fashion: #####

DT = var_est_NGR_bm(DT, months = months, validation_years = validation_years)
mean_CRPS_bm = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bm),na.rm = TRUE)]

DT = var_est_NGR_bl(DT, months = months, validation_years = validation_years,mc.cores = mc_cores)
mean_CRPS_bl = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bl),na.rm = TRUE)]

DT = var_est_NGR_bb(DT, months = months, validation_years = validation_years,mc.cores = mc_cores)
mean_CRPS_bb = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bb),na.rm = TRUE)]


CRPS_sma = msc_sd_sma[,mean(min_crps)]
CRPS_ema = msc_sd_ema[,mean(min_crps)]

CRPS_comparison = data.table(mean_CRPS_bm = mean_CRPS_bm,
                             mean_CRPS_bl = mean_CRPS_bl,
                             mean_CRPS_bb = mean_CRPS_bb,
                             mean_CRPS_sma = CRPS_sma,
                             mean_CRPS_ema = CRPS_ema)

CRPS_comparison = round(CRPS_comparison,5)


#####################################################
### permutation tests for moving average vs lr_bb ###
#####################################################

N=500

# CRPS (for variance estimation)

perm_test_dt = DT[year %in% validation_years & month %in% months,.(year,month,SST_bar,SST_hat,SD_hat,SD_hat_lr_bb)]

# getting CRPSs

perm_test_dt[,CRPS_ma := crps_na_rm(SST_bar, SST_hat,SD_hat)]
perm_test_dt[,CRPS_lr_bb := crps_na_rm(SST_bar, SST_hat, SD_hat_lr_bb)]


### permutation test for CRPS_ma ~ CRPS_lr_bb ###

pt_CRPS = permutation_test_difference(na.omit(perm_test_dt[,CRPS_ma]),na.omit(perm_test_dt[,CRPS_lr_bb]), N = N  )

pdf(paste0(plot_dir,'Perm_test_CRPS.pdf'))

rr = max(abs(1.1*pt_CRPS$d_bar),abs(1.1*pt_CRPS$D))
rr = c(-rr,rr)

hist(pt_CRPS$D, xlim = rr,breaks = 20,
     xlab = '', main = 'permutation test CRPS: LR_es vs. MA')

abline(v = pt_CRPS$d_bar,col = 'red')


qq = quantile(pt_CRPS$D,c(0.025,0.975))

abline(v = qq,lty = 2)

dev.off()


# permutation test for CRPS_ma ~ CRPS_lr_bb, averaged over the globe

ptbm = perm_test_dt[,.('CRPS_ma' = mean(CRPS_ma,na.rm = TRUE),'CRPS_lr_bb' = mean(CRPS_lr_bb,na.rm = TRUE)),by = .(year,month)]

pt_CRPS_bm = permutation_test_difference(ptbm[,CRPS_ma],ptbm[,CRPS_lr_bb], N = 5000  )

pdf(paste0(plot_dir,'Perm_test_glob_mean_CRPS.pdf'))

rr = max(abs(1.1*pt_CRPS_bm$d_bar),abs(1.1*pt_CRPS_bm$D))
rr = c(-rr,rr)

hist(pt_CRPS_bm$D, xlim = rr,breaks = 20,
     xlab = '', main = 'permutation test global mean CRPS: LR_es vs. MA')

abline(v = pt_CRPS_bm$d_bar,col = 'red')

qq = quantile(pt_CRPS_bm$D,c(0.025,0.975))

abline(v = qq,lty = 2)


dev.off()




#####################################################
### permutation tests for MSE MA vs lr_bb ###
#####################################################

perm_test_dt = DT[year %in% validation_years & month %in% months,.(year,month,SST_bar,SST_hat,T_hat_lr_both)]

# getting MSEs

perm_test_dt[,MSE_ma := (SST_bar - SST_hat)^2]
perm_test_dt[,MSE_lr_bb := (SST_bar - T_hat_lr_both)^2]


### permutation test for MSE_ma ~ MSE_lr_bb ###

pt_MSE = permutation_test_difference(na.omit(perm_test_dt[,MSE_ma]),na.omit(perm_test_dt[,MSE_lr_bb]), N = N )

pdf(paste0(plot_dir,'Perm_test_MSE.pdf'))

rr = max(abs(1.1*pt_MSE$d_bar),abs(1.1*pt_MSE$D))
rr = c(-rr,rr)

hist(pt_MSE$D, xlim = rr,breaks = 10,
     xlab = '', main = latex2exp::TeX('MSE permutation test for EMA vs. $NGR_{m,s}$'))

abline(v = pt_MSE$d_bar,col = 'red')


qq = quantile(pt_MSE$D,c(0.05))
abline(v = qq,lty = 2)

qq_2 = quantile(pt_MSE$D,c(0.01))
abline(v = qq_2,lty = 3)


dev.off()

# permutation test for MSE_ma vs MSE_lr_bb, averaged over the globe

ptbm = perm_test_dt[,.('MSE_ma' = mean(MSE_ma,na.rm = TRUE),'MSE_lr_bb' = mean(MSE_lr_bb,na.rm = TRUE)),by = .(year,month)]

pt_MSE_bm = permutation_test_difference(ptbm[,MSE_ma],ptbm[,MSE_lr_bb], N = 10000  )

pdf(paste0(plot_dir,'Perm_test_glob_mean_MSE.pdf'))

rr = max(abs(1.1*pt_MSE_bm$d_bar),abs(1.1*pt_MSE_bm$D))
rr = c(-rr,rr)

hist(pt_MSE_bm$D, xlim = rr,breaks = 20,
     xlab = '', main = latex2exp::TeX('globally averaged MSE permutation test for EMA vs. $NGR_{m,s}$'))

abline(v = pt_MSE_bm$d_bar,col = 'red')

qq = quantile(pt_MSE_bm$D,c(0.05))
abline(v = qq,lty = 2)

qq_2 = quantile(pt_MSE_bm$D,c(0.01))
abline(v = qq_2,lty = 3)


dev.off()





#### save stuff ####

time_s31 = proc.time() - time_s31

save(RMSE_linear_models,CRPS_comparison,file = paste0(save_dir,'univ_scores_comp_NGR.RData'))




save.image(file = paste0(save_dir,"setup.RData"))
