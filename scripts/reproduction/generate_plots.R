

#############################################################
### This script generates all the plots used in the paper ###
#############################################################



rm(list = ls())

time_s4 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


########### set parameters for plots #######################

par('cex' = 0.75, 'cex.lab' = 0.6,'cex.axis' = 0.6)

plot_dir0 = './figures/paper/'
dir.create(plot_dir0,showWarnings = FALSE)


Lat_res = c(-75,80) # Latitude restrictions


############################################################


#### plot of the first principal component ####


# get data 

name_abbr = "NAO" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


# reset plot directory

plot_dir = plot_dir0


# year and month

y = 2016
m = 7

# get the singular value decomposition of weight matrix * sample covariance matrix

wm = weight_mat(DT,L = 7000) 


num_loc = DT[year == min(year) & month == min(month)][!(is.na(SST_bar) | is.na(SST_hat)) ,.N]

train_years = DT[month == m][year < y & year > min(year),][,unique(year)]

data_mat = matrix(DT[month == m][!(is.na(SST_bar) | is.na(SST_hat)) & year %in% train_years,
                                 SST_bar - SST_hat],
                  nrow = num_loc)

sam_cov_mat = 1/length(train_years) * data_mat %*% t(data_mat) 


sin_val_dec_1 = svd(wm * sam_cov_mat)

# some data preparation

dt_water = DT[!(is.na(Ens_bar) | is.na(SST_bar))]

SD_cols = c("Lon","Lat","grid_id","month","year","YM",
            "SST_hat","SST_bar","Ens_bar","Bias_Est","var_bar","SD_hat")
SD_cols = SD_cols[which(SD_cols %in% colnames(dt))]
fc_water <- na.omit( dt_water[,.SD,.SDcols = SD_cols])

dt_ym = fc_water[month == m & year == y,]


# complement dt_ym by principal components

for(i in 1:30)
{
  temp =  sin_val_dec_1$u[,i]
  
  dt_ym[,paste0('PC',i):= temp]
  
  dt_test = rbindlist(list(dt_ym,dt[year == y & month ==m][is.na(Ens_bar) | is.na(SST_bar),.SD,.SDcols = SD_cols]), fill = TRUE)
  
}


# generate plot for restricted area

dt_test_new = dt_test[Lon >= -20 & Lat >= 50,]

rr = dt_test_new[,range(PC1,na.rm = TRUE)]
rr = c(-(max(abs(rr))),(max(abs(rr))))

plot_smooth(dt_test_new,paste0('PC1'),mn = '1st PC of covariance matrix for June',save_pdf = TRUE,file_name = '1stPCtapered',save_dir = plot_dir,xlab = '',ylab = '' )  

########################################################


############# Plot normalized SST #######################


name_abbr = "Full" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# reset plot directory

plot_dir = plot_dir0

# year and month
y = 2016
m = 5


DT = DT[Lat %between% Lat_res]

dt_ym = DT[year == y & month == m]

# get sample mean and sample variance

test_dt = DT[month == m & year < y & !is.na(SST_hat),][,.('clim' = mean(SST_bar,na.rm = TRUE),'clim_sd' = sd(SST_bar,na.rm = TRUE)),by = grid_id]

n_years = dt[month == m & year < y & !is.na(SST_hat),][grid_id == min(grid_id),.N]

temp = dt_ym[!is.na(SST_hat),][,'clim' := test_dt[,clim]][,'clim_sd' := test_dt[, clim_sd]]

dt_ym2 = merge(temp,dt_ym,by = colnames(dt_ym),all.y = TRUE)


# bound standard deviations away from 0

dt_ym2[clim_sd < 0.1,clim_sd := 0.1]

# get standardized SST and plot
dt_ym2[ ,'SST_stan' := (SST_bar - clim)/clim_sd]

plot_smooth(dt_ym2,'SST_stan',rr = c(-5,5),
            mn = 'normalized sea surface temperature in May 2016',
            pixels = 512,
            save_pdf = TRUE, save_dir = plot_dir,file_name = 'sst_nor_m05_y2016')



#########################################################################


name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


# reset plot directory

plot_dir = plot_dir0


######## get the distribution fct. of a censored normal distribution #############

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

DT_calib = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,SST_hat,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib[,"PIT_mc" := dist_fun_tn(SST_bar, mean = SST_hat, sd = SD_hat)]

DT_calib[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

########### plot mean and standard deviation ################

plot_diagnostic(DT_calib[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_mean")


unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0, 2*unif_sd), set_white = unif_sd,
                mn = paste0("PIT standard deviation"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_sd")



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


