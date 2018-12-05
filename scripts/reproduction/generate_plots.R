

#############################################################
### This script generates ALL the plots used in the paper ###
#############################################################

# this script is not optimized to be computationally efficient. For any plot, it loads a saved workspace image to retain the corresponding data set.
# It is adviced not to run the entire script if possible, but to use this to reconstruct and or modify single plots by running the corresponding section below.



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


Lat_res = c(-75,80) # Latitude restrictions for area plots in order to exclude the polar regions


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


############# Plot of normalized SST #######################


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



##################################################

######## Plots for univariate calibration ########

##################################################

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


# reset plot directory

plot_dir = plot_dir0

##################################################

################## PIT plots #####################



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

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean, linear regression"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_mean_lr")


# the following plots the PIT standard deviation for mean estimation by linear regression. This is not in the paper but nevertheless interesting.

unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0,2*unif_sd),set_white = unif_sd,
                mn = paste0("PIT mean, linear regression"),
                save_pdf = TRUE, save_dir = plot_dir, file_name = "PIT_sd_lr")


###########################################################

######### Permutation tests for MSE #######################

# permutation test where each location is considered separately
pt_MSE = permutation_test_difference(na.omit(perm_test_dt[,MSE_ma]),na.omit(perm_test_dt[,MSE_lr_bb]), N = 1000 )

# for the spatially averaged permutation test:
ptbm = perm_test_dt[,.('MSE_ma' = mean(MSE_ma,na.rm = TRUE),'MSE_lr_bb' = mean(MSE_lr_bb,na.rm = TRUE)),by = .(year,month)]

pt_MSE_bm = permutation_test_difference(ptbm[,MSE_ma],ptbm[,MSE_lr_bb], N = 10000  )



pdf(file = paste0(plot_dir,'Perm_test_MSE.pdf'),width = 15)

par('mfrow' = c(1,2))
  
  ### first plot: separate locations ###
  #range
  rr = max(abs(1.1*pt_MSE$d_bar),abs(1.1*pt_MSE$D))
  rr = c(-rr,rr)
  
  #histogram
  hist(pt_MSE$D, xlim = rr,breaks = 10,
       xlab = '', main = latex2exp::TeX('MSE permutation test for EMA vs. $NGR_{m,s}$'))
  
  #actual observation
  abline(v = pt_MSE$d_bar,col = 'red')
  
  # 5th percentile
  qq = quantile(pt_MSE$D,c(0.05))
  abline(v = qq,lty = 2)
  
  # 1st percentile
  qq_2 = quantile(pt_MSE$D,c(0.01))
  abline(v = qq_2,lty = 3)

  ### second plot: averaged over the globe ###

  #range
  rr = max(abs(1.1*pt_MSE_bm$d_bar),abs(1.1*pt_MSE_bm$D))
  rr = c(-rr,rr)
  
  #histogram
  hist(pt_MSE_bm$D, xlim = rr,breaks = 20,
       xlab = '', main = latex2exp::TeX('globally averaged MSE permutation test for EMA vs. $NGR_{m,s}$'))
  
  # actual observation
  abline(v = pt_MSE_bm$d_bar,col = 'red')
  
  # 5th percentile
  qq = quantile(pt_MSE_bm$D,c(0.05))
  abline(v = qq,lty = 2)
  
  # 1st percentile
  qq_2 = quantile(pt_MSE_bm$D,c(0.01))
  abline(v = qq_2,lty = 3)

dev.off()


#######################################################

####### Plot RMSE by weighting parameter ##############

#### plotting ####

# get data for last year
y = max(validation_years)

row_sma = msc_sma[year == y,-1,with = FALSE]
row_ema = msc_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))  

y_range = sqrt(y_range)



pdf(paste0(plot_dir,"RMSE_by_par.pdf"),width = 15)

  par('mfrow' = c(1,2))
  
  ## plot for sma ##
  
  plot(x = win_length, 
       y = as.vector(sqrt(row_sma[,-c('min_MSE','min_l'),with = FALSE])),
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("RMSE for bias correction by SMA"),
       xlab = "window length",
       ylab = "RMSE"
  )
  
  # highlight minimum and add minimum reference line 
  abline(h = row_sma[,sqrt(min_MSE)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  
  points(x = row_sma[,min_l],
         y = row_sma[,sqrt(min_MSE)],
         col = "blue",
         bg = "blue",
         pch = 21)
  
  
  ## plot for ema ##
  
  plot(x = par_vec, 
       y = sqrt(row_ema[,-c('min_MSE','min_a'),with = FALSE]),
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("RMSE for bias correction by EMA"),
       xlab = "scale parameter",
       ylab = "RMSE"
  )
  
  # highlight minimum and add minimum reference line 
  abline(h = row_ema[,sqrt(min_MSE)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  
  points(x = row_ema[,min_a],
         y = row_ema[,sqrt(min_MSE)],
         col = "blue",
         bg = "blue",
         pch = 21)
  
dev.off()



#############################################################

################## CRPS by parameter ########################

# get data for last year
y = max(validation_years)


row_sma = msc_sd_sma[year == y,-1,with = FALSE]
row_ema = msc_sd_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))  

pdf(paste0(plot_dir,"/mean_CRPS_sd.pdf"),width = 15)

  par('mfrow' = c(1,2))
  
  ## plot for sma ##
  
  plot(x = win_length, 
       y = as.vector(row_sma[,-c('min_crps','min_l'),with = FALSE]),
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("CRPS for variance estimation by SMA"),
       xlab = "window length",
       ylab = "CRPS")
  
  # highlight minimum and add minimum reference line 
  abline(h = row_sma[,min_crps], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  
  points(x = row_sma[,min_l],
         y = row_sma[,min_crps],
         col = "blue",
         bg = "blue",
         pch = 21)
  
  
  ## plot for ema ##
  
  plot(x = par_vec, 
       y = row_ema[,-c('min_crps','min_a'),with = FALSE],
       ylim = y_range,
       type = "b",
       col = "blue",
       main = paste0("CRPS for variance estimation by EMA"),
       xlab = "scale parameter",
       ylab = "CRPS"
  )
  
  # highlight minimum and add minimum reference line 
  abline(h = row_ema[,min_crps], lty = "dashed", col = adjustcolor("blue",alpha = .5))
  
  points(x = row_ema[,min_a],
         y = row_ema[,min_crps],
         col = "blue",
         bg = "blue",
         pch = 21)

dev.off()




