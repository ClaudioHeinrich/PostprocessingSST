
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

name_abbr = "test" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

#start timer:

time_s3is = proc.time()

###### getting sample variances of ensemble  ######

DT = ens_sd_est(dt = DT,
                ens_size = ens_size,
                save_dir = save_dir,
                mean_est = "sv",
                file_name = paste0("dt_combine_wide_bc_var.RData"))


###### getting scores for simple moving averages ######

num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1

dummy_function = function(k){
  temp = sd_est_2(dt = DT, method = "sma", par_1 = k)
  
  temp[,SST_hat := DT[,SST_hat]]
  temp[,SST_bar := DT[,SST_bar]]
  
  CRPS = temp[, crps_na_rm(SST_bar,mean = SST_hat,sd = SD_hat),by = .(year,month)]
  setnames(CRPS,old = 3,new = 'CRPS')
  
  CRPS[,l:=k]
  return(CRPS)
}


sc_sma_var = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = mc_cores)
sc_sma_var = rbindlist(sc_sma_var)


###### getting scores for exponential moving averages ######

par_vec = seq(0.01,0.4,length.out = 24)


dummy_function = function(k){
  temp = sd_est_2(dt = DT, method = "ema", par_1 = k)
  
  temp[,SST_hat := DT[,SST_hat]]
  temp[,SST_bar := DT[,SST_bar]]
  
  CRPS = temp[, crps_na_rm(SST_bar,mean = SST_hat,sd = SD_hat),by = .(year,month)]
  setnames(CRPS,old = 3,new = 'CRPS')
  
  CRPS[,l:=k]
  return(CRPS)
}

sc_ema_var = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = mc_cores)
sc_ema_var = rbindlist(sc_ema_var)



############### get optimal parameter and estimate variance ################

### sma ###

temp = sc_sma_var[,mean(CRPS,na.rm = TRUE),by = l]
setnames(temp,2,'CRPS')
sma_par_var = temp[,which.min(CRPS)]
sma_crps_var = temp[,min(CRPS)]

temp_sma_var = sd_est_2(dt = DT,
                        method = 'sma',
                        par_1 = sma_par_var)
DT[,SD_hat_SMA := temp_sma_var[,SD_hat]]

#twosided training:

DT_sma_tr = sd_est_2(dt = DT,
                         method = 'sma',
                         par_1 = sma_par_var,
                         twosided = TRUE)

DT[year %in% training_years, 'SD_hat_SMA' := DT_sma_tr[year %in% training_years,SD_hat]]

rm(DT_sma_tr)

### ema ###

temp = sc_ema_var[,mean(CRPS,na.rm = TRUE),by = l]
setnames(temp,2,'CRPS')
ema_par_var = temp[,which.min(CRPS)]
ema_crps_var = temp[,min(CRPS)]

temp_ema_var = sd_est_2(dt = DT,
                        method = 'ema',
                        par_1 = ema_par_var)
DT[,SD_hat_EMA := temp_ema_var[,SD_hat]]

#twosided training:

DT_ema_tr = sd_est_2(dt = DT,
                     method = 'ema',
                     par_1 = ema_par_var,
                     twosided = TRUE)

DT[year %in% training_years, 'SD_hat_EMA' := DT_ema_tr[year %in% training_years,SD_hat]]

rm(DT_ema_tr)

###### finding and applying optimal way of bias correcion, exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

if(sma_crps_var < 0.9 * ema_crps_var){
  DT[,SD_hat := SD_hat_SMA]
} else{
  DT[,SD_hat := SD_hat_EMA]
}




#### time, update script counter, save ####

time_s3is = proc.time() - time_s3is

script_counter = 3

save.image(file = paste0(save_dir,"setup.RData"))
