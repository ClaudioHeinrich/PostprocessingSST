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

mc_cores = 1

###### compare to linear regression models ######


mean_CRPS_bm = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bm),na.rm = TRUE)]


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







#### save stuff ####

time_s31 = proc.time() - time_s31

save(RMSE_linear_models,CRPS_comparison,file = paste0(save_dir,'univ_scores_comp_NGR.RData'))




save.image(file = paste0(save_dir,"setup.RData"))
