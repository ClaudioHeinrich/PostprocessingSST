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

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

#start timer:

time_s2is = proc.time()

###### run bias analysis for simple moving averages ######

num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1


dummy_function = function(k){
  temp = bias_correct(dt = DT, method = "sma", par_1 = k)
  
  temp[,SST := DT[,SST_bar]]
  
  MSE = temp[, mean((SST_hat - SST)^2,na.rm = TRUE),by = .(year,month)]
  setnames(MSE,old = 3,new = 'MSE')
  
  MSE[,l:=k]
  return(MSE)
}

sc_sma = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = mc_cores)
sc_sma = rbindlist(sc_sma)

save(sc_sma, file = paste0(save_dir,"scores.bc.sma.RData"))

###### run bias analysis for exponential moving averages ######

par_vec = seq(0.05,0.4,length.out = 24) 


dummy_function = function(k){
  temp = bias_correct(dt = DT, method = "ema", par_1 = par_vec[k])
  
  temp[,SST := DT[,SST_bar]]
  
  MSE = temp[,mean((SST_hat - SST)^2,na.rm = TRUE),by = .(year,month)]
  setnames(MSE,old = 3,new = 'MSE')
  
  MSE[,a:=par_vec[k]]
  return(MSE)
}

sc_ema = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = mc_cores)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save_dir,"scores.bc.ema.RData"))


############### get optimal parameter and estimate bias ################

### sma ###

temp = sc_sma[,mean(MSE),by = l]
setnames(temp,2,'MSE')
sma_par = temp[,which.min(MSE)]
sma_mse = temp[,min(MSE)]

temp_sma = bias_correct(dt = DT,
                          method = 'sma',
                          par_1 = sma_par)
DT[,Bias_Est_SMA := temp_sma[,Bias_est]][,SST_hat_sma := temp_sma[,SST_hat]]

#twosided training:

DT_sma_tr = bias_correct(dt = DT,
                           method = 'sma',
                           par_1 = sma_par,
                           twosided = TRUE)

DT[year %in% training_years, c('Bias_Est_SMA','SST_hat_sma') := DT_sma_tr[year %in% training_years,.(Bias_est,SST_hat)]]

rm(DT_sma_tr)

### ema ###

temp = sc_ema[,mean(MSE),by = a]
setnames(temp,2,'MSE')
ema_par = temp[,which.min(MSE)]
ema_mse = temp[,min(MSE)]

temp_ema = bias_correct(dt = DT,
                          method = 'ema',
                          par_1 = ema_par)
DT[,Bias_Est_EMA := temp_ema[,Bias_est]][,SST_hat_ema := temp_ema[,SST_hat]]

#twosided training:

DT_ema_tr = bias_correct(dt = DT,
                           method = 'ema',
                           par_1 = ema_par,
                           twosided = TRUE)

DT[year %in% training_years, c('Bias_Est_EMA','SST_hat_ema') := DT_ema_tr[year %in% training_years,.(Bias_est,SST_hat)]]

rm(DT_ema_tr)


###### finding and applying optimal way of bias correcion, exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

if(sma_mse < 0.9 * ema_mse){
  DT[,Bias_Est := Bias_Est_SMA]
  DT[,SST_hat := SST_hat_sma]
} else{
  DT[,Bias_Est := Bias_Est_EMA]
  DT[,SST_hat := SST_hat_ema]
}




#### time, update script counter, save ####

time_s2is = proc.time() - time_s2is

script_counter = 2

save.image(file = paste0(save_dir,"setup.RData"))
