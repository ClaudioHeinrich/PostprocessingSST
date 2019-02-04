#######################################################################################

###################  side script 3.2 - PIT plots  ######################

#######################################################################################

# This script generates plots of the mean and standard deviation of the marginal probability integral transforms 
# for the marginally calibrated forecast.
#
# 
# Files generated:
#   
# Plots: PIT_mean.pdf, PIT_sd.pdf 
#
# Requires previous run of 03.master.var.est 
# with the same value of name_abbr as below.

############# setup ###############

rm(list = ls())

time_s32 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


######## get the distribution fct. of a censored normal #############

# value, mean and sd need to be vectors of equal length.
# returns F(value), where F is dist. fct. of a normal with parameters mean and sd, censored at trc_value

dist_fun_tn = function(value, mean, sd, trc_value = -1.79){ 
  a=rep(0,times = length(value))
  na_loc = which(is.na(value) | is.na(sd) | sd == 0)
  trc_loc = which(value <= trc_value & sd > 0)
  nor_loc = which(value > trc_value & sd > 0)
  
  a[na_loc] = NA
  a[trc_loc] = runif(length(trc_loc), max = pnorm(trc_value, mean = mean[trc_loc], sd = sd[trc_loc]))
  a[nor_loc] = pnorm(value[nor_loc], mean = mean[nor_loc], sd = sd[nor_loc])
    
  return(a)
}

########### get PITs ###############

DT_calib_1 = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_1[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]

DT_calib_1[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_1[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

########### plot ################
Lat_res = c(-75,80)


par('cex' = 0.75, 'cex.axis' = 0.75)
plot_diagnostic(DT_calib_1[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_mean")


unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib_1[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0,2*unif_sd),set_white = unif_sd,
                mn = paste0("PIT standard deviation"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_sd")




#################### sensitivity to validation period ################################

for(f_y in seq(2000,2016,by = 4))
{print(f_y)
validation_years_2 = 2010:2016


DT_calib_2 = DT[year %in% validation_years_2,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_2[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]

DT_calib_2[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_2[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

########### plot ################


plot_diagnostic(DT_calib_2[year == min(year) & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean ", min(validation_years_2),' - ',max(validation_years_2)),
                save_pdf = TRUE, save_dir = plot_dir, file_name = paste0("PIT_mean_sv"))

}






###### for estimated mean by linear regression ######

DT_calib_3 = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,T_hat_lr_both,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_3[,"PIT_mc" := dist_fun_tn(SST_bar, mean = T_hat_lr_both, sd = SD_hat)]

DT_calib_3[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_3[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

########### plot ################

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean, linear regression"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_mean_lr")

unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0,2*unif_sd),set_white = unif_sd,
                mn = paste0("PIT mean, linear regression"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_sd_lr")


################# look at PIT means over different validation periods month by month ##############################

for(m in 1:12)
{print(m)
  validation_years_2 = 2010:2016
  validation_years_3 = 2000:2009
  
  DT_calib_1 = DT[year %in% validation_years & month == m,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]
  DT_calib_2 = DT[year %in% validation_years_2 & month == m,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]
  DT_calib_3 = DT[year %in% validation_years_3 & month == m,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]
  
  # PIT for marginally corrected forecast:
  DT_calib_1[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]
  DT_calib_2[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]
  DT_calib_3[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]
  
  DT_calib_1[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id][,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]
  DT_calib_2[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id][,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]
  DT_calib_3[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id][,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]
  
  ########### plot ################
  
  
  plot_diagnostic(DT_calib_1[year == min(year) & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                  rr = c(0,1),
                  mn = paste0("PIT mean ", min(validation_years),' - ',max(validation_years),', M ',m),
                  save_pdf = TRUE, save_dir = plot_dir, file_name = paste0("PIT_mean_m",m))
  plot_diagnostic(DT_calib_2[year == min(year) & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                  rr = c(0,1),
                  mn = paste0("PIT mean ", min(validation_years_2),' - ',max(validation_years_2),', M ',m),
                  save_pdf = TRUE, save_dir = plot_dir, file_name = paste0("PIT_mean_sv_m",m))
  
  plot_diagnostic(DT_calib_3[year == min(year) & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                  rr = c(0,1),
                  mn = paste0("PIT mean ", min(validation_years_3),' - ',max(validation_years_3),', M ',m),
                  save_pdf = TRUE, save_dir = plot_dir, file_name = paste0("PIT_mean_ev_m",m))
  
}





time_s32 = proc.time() - time_s32

save.image(file = paste0(save_dir,"setup.RData"))

