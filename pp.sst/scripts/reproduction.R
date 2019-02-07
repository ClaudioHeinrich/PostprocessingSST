

##########################################################################

######################## reproduction script #############################

##########################################################################

# This script reproduces the results of the paper Postprocessing seasonal weather forecasts



### setup ###

rm(list = ls())

setwd('/nr/user/claudio/pkg/paper/PostprocessingSST')

options(max.print = 1e3)

library(devtools)
library(data.table)


install_github('ClaudioHeinrich/PostprocessingSST/pp.sst')

library(pp.sst) # contains all relevant functions as well as an example data set


### set directories for derived data and plots ###

save_dir = "../Derived/"
dir.create(save_dir, showWarnings = FALSE)

plot_dir = "../figures/"
dir.create(plot_dir, showWarnings = FALSE)


############ explore data ###############

DT = copy(example_dt)
DT
save(DT,file = paste0(save_dir,'data.RData'))

# example plot of SST:

plot_diagnostic(DT[year == 2001],var = 'SST_bar',mn = 'SST 06/2001')

# example of residual between observed SST and mean of the raw forecast ensemble:

plot_diagnostic(DT[year == 2001,.(Lon,Lat,SST_bar - Ens_bar)],mn = 'forecast residual 06/2001')

# replace plot_diagnostic by plot_smooth for nicer but slower plots


###########################################################

############## set a couple of parameters #################

ens_size = 9 # number of forecast ensemble

training_years = 1985:2000
validation_years = 2001:2016

mc_cores = 6 # number of cores for parallelization

months = 2 # do not change: only month 2 is included in the example data

####################################################

####### bias correction using moving averages ######



###### bias analysis for simple moving averages ######

past_0 = 5
# The weighting parameter is estimated for each year in the validation period anew.
# Given year y in the validation period, the MSE over the years {min(validation_years) - past_0, ... ,y-1} is computed for a range of weighting parameters.
# The weighting parameter with minimal MSE is then chosen for bias correction in year y. For year y+1 the weighting parameter is estimated anew.


num_years = max(validation_years) - DT[,min(year)]

win_length = 1 : (num_years-1) # range of parameter l considered for simple moving averages

MSE_by_par = function(l) # computes the bias estimate for window length l, return data table with bias estimates and SST estimates
{
  temp = bias_correct_2(dt = DT,
                        method = "sma",
                        par_1 = l,
                        reduced_output = TRUE)
  return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SST_hat])
}

#compute bias correction for all years and all parameters in win_length:

BC = parallel::mclapply(X = win_length, FUN = MSE_by_par,mc.cores = mc_cores)

# restrict to validation years + the past_0 years before:

Bias_est_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                 .(year,month,Lon,Lat,SST_bar)]

for(k in win_length)
{
  Bias_est_dt[,paste0('l',k):= BC[[k]]]
  # compute squared error:
  Bias_est_dt[,paste0('err',k):= (SST_bar - eval(parse(text = paste0('l',k))))^2]
}

#get MSE by year and month, for each year based on all previous years contained in Bias_est_dt

sc_sma = NULL

for(y in validation_years)
{print(y)
  mean_sc_bm = function(m)
  {
    temp = data.table(year = y ,month = m,Bias_est_dt[year < y & month == m,
                                                      lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                      .SDcols = paste0('err',win_length)])
    return(temp)
  }

  MSE_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_bm,mc.cores = mc_cores))
  sc_sma = rbindlist(list(sc_sma,MSE_y))
}

#save:

save(sc_sma, file = paste0(save_dir,"scores.bc.sma.RData"))


###### bias analysis for exponential moving averages ######

# specify the range of scale parameters a considered:

par_vec = seq(0.05,0.4,length.out = 24)

#structure as for simple moving averages:

MSE_by_par = function(a)
{
  temp = bias_correct_2(dt = DT,
                        method = "ema",
                        par_1 = a,
                        reduced_output = TRUE)
  return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SST_hat])
}


BC = parallel::mclapply(X = par_vec, FUN = MSE_by_par,mc.cores = mc_cores)

# restrict to validation years + the past_0 years before

Bias_est_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                 .(year,month,Lon,Lat,SST_bar)]

for(a in par_vec)
{ind = which(par_vec == a )
ra = round(a,4)
Bias_est_dt[,paste0('a',ra):= BC[[ind]]]
Bias_est_dt[,paste0('err',ra):= (SST_bar - eval(parse(text = paste0('a',ra))))^2]
}

#get MSE by year and month, for each year based on all previous years contained in Bias_est_dt

sc_ema = NULL

for(y in validation_years)
{print(y)
  mean_sc_bm = function(m)
  {
    temp = data.table(year = y ,month = m,Bias_est_dt[year < y & month == m,
                                                      lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                      .SDcols = paste0('err',round(par_vec,4))])
    return(temp)
  }

  MSE_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_bm,mc.cores = mc_cores))
  sc_ema = rbindlist(list(sc_ema,MSE_y))
}


save(sc_ema, file = paste0(save_dir,"scores.bc.ema.RData"))

###############################################################

###### plot scores for different ways of bias correction ######

# get means over all months for sma

msc_sma = sc_sma[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_l = NULL

for(i in 1:msc_sma[,.N])
{
  row = msc_sma[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_l = c(min_l,which.min(row))
}

msc_sma[,min_MSE := row_min][,min_l := min_l]


# get means over all months for ema

msc_ema = sc_ema[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_a = NULL

for(i in 1:msc_sma[,.N])
{
  row = msc_ema[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_a = c(min_a,par_vec[which.min(row)])
}

msc_ema[,min_MSE := row_min][,min_a := min_a]


#### plotting ####

# get data for last year, i.e. the MSEs by parameter for the years 1995,...,2016
y = max(validation_years)


row_sma = msc_sma[year == y,-1,with = FALSE]
row_ema = msc_ema[year == y,-1,with = FALSE]

#range for both plots:
y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))



## plot ##

par(mfrow = c(1,2))

plot(x = win_length,
     y = as.vector(row_sma[,-c('min_MSE','min_l'),with = FALSE]),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("SMA"),
     xlab = "window length",
     ylab = "MSE"
)

# highlight minimum and add minimum reference line
abline(h = row_ema[,min_MSE], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_sma[,min_l],
       y = row_sma[,min_MSE],
       col = "blue",
       bg = "blue",
       pch = 21)


## plot for ema ##


plot(x = par_vec,
     y = row_ema[,-c('min_MSE','min_a'),with = FALSE],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("EMA"),
     xlab = "scale parameter",
     ylab = "MSE"
)

# highlight minimum and add minimum reference line
abline(h = row_ema[,min_MSE], lty = "dashed", col = adjustcolor("blue",alpha = .5))

points(x = row_ema[,min_a],
       y = row_ema[,min_MSE],
       col = "blue",
       bg = "blue",
       pch = 21)



###### finding optimal way of bias correction for each year in the validation period ######

# exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

opt_par = data.table(year = validation_years,method = NA_character_,par = NA_real_)


for(y in validation_years)
{
  if(msc_sma[year == y, min_MSE] < 0.95 * msc_ema[year == y, min_MSE]){
    opt_par[year == y, method := 'sma']
    opt_par[year == y,par := msc_sma[year == y,min_l]]

  } else{
    opt_par[year == y, method := 'ema']
    opt_par[year == y,par := msc_ema[year == y,min_a]]
  }
}


#### apply bias correction ####

DT[,Bias_Est := NA_real_][,SST_hat := NA_real_]

for(y in validation_years)
{
  print(y)
  temp = bias_correct_2(dt = DT,
                        method = opt_par[year == y, method],
                        par_1 = opt_par[year == y,par])[year == y,]
  DT[year == y,Bias_Est := temp[,Bias_est]]
  DT[year == y,SST_hat := temp[,SST_hat]]
}

rm(temp)



# For the training years the bias correction considers also the future,
# and uses the parameter and method estimated for the first validation year

method_ty = c(opt_par[year == min(year),method],opt_par[year == min(year),par])

DT = bias_correct_training(dt = DT,
                           method = method_ty,
                           training_years = training_years,
                           save_dir = save_dir)



#######################################################

############## variance estimation ####################

# The variance estimation is largely analogous to the bias correction. The MSE is replaced by CRPS.

###### getting sample variances of ensemble  ######

DT[,var_bar := (SST_hat - SST_bar)^2]

###### running variance estimation for simple moving averages ######

if(mc_cores > 1)
{
  SD_est_by_par = function(k)
  {
    temp = sd_est_2(dt = DT,
                    method = "sma",
                    par_1 = k,
                    reduced_output = TRUE)
    return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SD_hat])
  }


  SD_ls = parallel::mclapply(X = win_length, FUN = SD_est_by_par,mc.cores = mc_cores)
} else {
  SD_ls = list()
  for(k in win_length)
  {
    print(k)

    temp = sd_est_2(dt = DT,
                    method = "sma",
                    par_1 = k,
                    reduced_output = TRUE)

    SD_ls[[k]] = temp[year %between% c(min(validation_years - past_0),max(validation_years)),SD_hat]
  }
}



# restrict to validation years + the past_0 years before

SD_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
           .(year,month,Lon,Lat,SST_bar,SST_hat)]

for(k in win_length)
{
  SD_dt[,paste0('l',k):= SD_ls[[k]]]
  SD_dt[,paste0('crps',k):= crps_na_rm(y = SST_bar,mean = SST_hat,sd = eval(parse(text = paste0('l',k))))]
}

#get mean CRPS by month, for each year averaged over all previous years contained in SD_dt

sc_sd_sma = NULL

for(y in validation_years)
{print(y)
  mean_sc_sd_bm = function(m)
  {
    temp = data.table(year = y ,month = m,SD_dt[year < y & month == m,
                                                lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                .SDcols = paste0('crps',win_length)])
    return(temp)
  }

  CRPS_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_sd_bm,mc.cores = mc_cores))
  sc_sd_sma = rbindlist(list(sc_sd_sma,CRPS_y))
}


save(sc_sd_sma, file = paste0(save_dir,"scores.bc.sd.sma.RData"))


###### variance estimation by exponential moving averages ######

par_vec = seq(0.05,0.4,length.out = 24)

SD_by_par = function(a)
{
  temp = sd_est_2(dt = DT,
                  method = "ema",
                  par_1 = a,
                  saveorgo = FALSE,
                  reduced_output = TRUE)
  return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SD_hat])
}


SD_ls = parallel::mclapply(X = par_vec, FUN = SD_by_par,mc.cores = mc_cores)

# restrict to validation years + the past_0 years before

SD_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
           .(year,month,Lon,Lat,SST_bar,SST_hat)]


for(a in par_vec)
{ind = which(par_vec == a )
print(ind)
ra = round(a,4)
SD_dt[,paste0('a',ra):= SD_ls[[ind]]]
SD_dt[,paste0('crps',ra):= crps_na_rm(y = SST_bar,mean = SST_hat,sd = eval(parse(text = paste0('a',ra))))]
}

#get mean crps by month, averaged over all previous years contained in SD_dt

sc_sd_ema = NULL

for(y in validation_years)
{print(y)
  mean_sc_sd_bm = function(m)
  {
    temp = data.table(year = y ,month = m,SD_dt[year < y & month == m,
                                                lapply(X = .SD,FUN = mean,na.rm = TRUE),
                                                .SDcols = paste0('crps',round(par_vec,4))])
    return(temp)
  }

  CRPS_y = rbindlist(parallel::mclapply(X = months,FUN = mean_sc_sd_bm,mc.cores = mc_cores))
  sc_sd_ema = rbindlist(list(sc_sd_ema,CRPS_y))
}


save(sc_sd_ema, file = paste0(save_dir,"scores.bc.sd.ema.RData"))



###### plotting scores for different ways of variance estimation ######

# get means over all months for sma

msc_sd_sma = sc_sd_sma[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_l = NULL

for(i in 1:msc_sd_sma[,.N])
{
  row = msc_sd_sma[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_l = c(min_l,which.min(row))
}

msc_sd_sma[,min_crps := row_min][,min_l := min_l]


# get means over all months for ema

msc_sd_ema = sc_sd_ema[,lapply(.SD,mean),by = year][,month := NULL]

# get minimum for each row

row_min = NULL
min_a = NULL

for(i in 1:msc_sd_sma[,.N])
{
  row = msc_sd_ema[,-1,with = FALSE][i,]
  row_min = c(row_min,min(row))
  min_a = c(min_a,par_vec[which.min(row)])
}

msc_sd_ema[,min_crps := row_min][,min_a := min_a]


#### plotting ####

# get data for last year
y = max(validation_years)


row_sma = msc_sd_sma[year == y,-1,with = FALSE]
row_ema = msc_sd_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))

####### plotting #######

par(mfrow = c(1,2))

plot(x = win_length,
     y = as.vector(row_sma[,-c('min_crps','min_l'),with = FALSE]),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("SMA"),
     xlab = "window length",
     ylab = "CRPS"
)

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
     main = paste0("EMA"),
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


###### finding optimal way of variance estimation for each year in the validation period ######

# exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

opt_par_sd = data.table(year = validation_years,method = NA_character_,par = NA_real_)


for(y in validation_years)
{
  if(msc_sd_sma[year == y, min_crps] < 0.95 * msc_sd_ema[year == y, min_crps]){
    opt_par_sd[year == y, method := 'sma']
    opt_par_sd[year == y,par := msc_sd_sma[year == y,min_l]]

  } else{
    opt_par_sd[year == y, method := 'ema']
    opt_par_sd[year == y,par := msc_sd_ema[year == y,min_a]]
  }
}


# variance correction year by year

DT[,SD_hat := NA_real_]

for(y in validation_years)
{
  print(y)
  temp = sd_est_2(dt = DT,
                  method = opt_par_sd[year == y, method],
                  par_1 = opt_par_sd[year == y,par],
                  reduced_output = TRUE)[year == y,]
  DT[year == y,SD_hat := temp[,SD_hat]]
}

rm(temp)

# For the training years the variance estimation uses the parameter and method estimated for the first validation year

method_ty = c(opt_par_sd[year == min(year),method],opt_par_sd[year == min(year),par])

temp = sd_est(dt = DT[year %in% training_years],
              method = method_ty[1],
              par_1 = as.numeric(method_ty[2]),
              saveorgo = FALSE)


DT[year %in% training_years, SD_hat:=temp[,SD_hat]]

save(DT,file = paste0(save_dir,'data.RData'))


#######################################################################
#######################################################################
#######################################################################

################## compare moving averages to NGR #####################

# get NGR mean estimates, grouped by month, location and both:

DT = bias_lr_bm(DT,months = months,validation_years = validation_years)
# linear regressions at grid points with no data causes warning:
DT = suppressWarnings(bias_lr_bl(DT,months = months,validation_years = validation_years))
# since we only have one month, NGR_{s,m} = NGR_s
DT[,T_hat_lr_both:=T_hat_lr_loc]


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

# show results:

RMSE_linear_models


###########################################################

########### Compare variance estimation to NGR  ###########


# grouped by month:
DT = var_est_NGR_bm(DT, months = months, validation_years = validation_years)
mean_CRPS_bm = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bm),na.rm = TRUE)]


#grouped by location:
DT = var_est_NGR_bl(DT, months = months, validation_years = validation_years,mc.cores = mc_cores)
mean_CRPS_bl = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bl),na.rm = TRUE)]

#grouped by both (same as by location):
DT[,SD_hat_lr_bb:=SD_hat_lr_bl]
mean_CRPS_bb = DT[,mean(crps_na_rm(SST_bar,SST_hat,SD_hat_lr_bb),na.rm = TRUE)]

# get mean CRPS for moving averages:
CRPS_sma = msc_sd_sma[,mean(min_crps)]
CRPS_ema = msc_sd_ema[,mean(min_crps)]

CRPS_comparison = data.table(mean_CRPS_bm = mean_CRPS_bm,
                             mean_CRPS_bl = mean_CRPS_bl,
                             mean_CRPS_bb = mean_CRPS_bb,
                             mean_CRPS_sma = CRPS_sma,
                             mean_CRPS_ema = CRPS_ema)

CRPS_comparison = round(CRPS_comparison,5)

#print results:

CRPS_comparison

#####################################################
##### conduct permutation tests #####################
#####################################################

# as an example, we conduct a permutation tests comparing bias estimation by EMA
# to the mean model NGR_{s,m}


# number of resamples for permutation test:
N=500

perm_test_dt_MSE = DT[year %in% validation_years & month %in% months,.(year,month,SST_bar,SST_hat,T_hat_lr_both)]

# getting MSEs

perm_test_dt_MSE[,MSE_ma := (SST_bar - SST_hat)^2]
perm_test_dt_MSE[,MSE_lr_bb := (SST_bar - T_hat_lr_both)^2]


### permutation test for MSE_ma ~ MSE_lr_bb ###

pt_MSE = permutation_test_difference(na.omit(perm_test_dt_MSE[,MSE_ma]),na.omit(perm_test_dt_MSE[,MSE_lr_bb]), N = N )

rr = max(abs(1.1*pt_MSE$d_bar),abs(1.1*pt_MSE$D))
rr = c(-rr,rr)

hist(pt_MSE$D, xlim = rr,breaks = 10,
     xlab = '', main = latex2exp::TeX('EMA vs. $NGR_{m,s}$'))

abline(v = pt_MSE$d_bar,col = 'red')


qq = quantile(pt_MSE$D,c(0.05))
abline(v = qq,lty = 2)

# permutation test for spatially averaged MSE:

ptbm = perm_test_dt_MSE[,.('MSE_ma' = mean(MSE_ma,na.rm = TRUE),'MSE_lr_bb' = mean(MSE_lr_bb,na.rm = TRUE)),by = .(year,month)]

pt_MSE_bm = permutation_test_difference(ptbm[,MSE_ma],ptbm[,MSE_lr_bb], N = N  )

rr = max(abs(1.1*pt_MSE_bm$d_bar),abs(1.1*pt_MSE_bm$D))
rr = c(-rr,rr)

hist(pt_MSE_bm$D, xlim = rr,breaks = 20,
     xlab = '', main = latex2exp::TeX('EMA vs. $NGR_{m,s}$, spatially averaged'))

abline(v = pt_MSE_bm$d_bar,col = 'red')

qq = quantile(pt_MSE_bm$D,c(0.05))
abline(v = qq,lty = 2)

############################################################
############################################################
############################################################

####################### PITs ###############################


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
                mn = paste0("PIT mean"))


unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib_1[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0,2*unif_sd),set_white = unif_sd,
                mn = paste0("PIT standard deviation"))


###### for mean estimation by linear regression ######

DT_calib_3 = DT[year %in% validation_years,.(Lon,Lat,year,month,grid_id,SST_bar,T_hat_lr_both,SD_hat,Ens_bar,Ens_sd)]

# PIT for marginally corrected forecast:
DT_calib_3[,"PIT_mc" := dist_fun_tn(SST_bar, mean = T_hat_lr_both, sd = SD_hat)]

DT_calib_3[,"PIT_mc_mean" := mean(PIT_mc), by = grid_id]
DT_calib_3[,"PIT_mc_sd" := sd(PIT_mc), by = grid_id]

########### plot ################

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_mean)],
                rr = c(0,1),
                mn = paste0("PIT mean, linear regression"))

unif_sd = sqrt(1/12) #standard deviation of uniform distribution

plot_diagnostic(DT_calib_3[year == min(year) & Lat %between% Lat_res & month == min(month), .(Lon,Lat,PIT_mc_sd)],
                rr = c(0,2*unif_sd),set_white = unif_sd,
                mn = paste0("PIT sd, linear regression"))


#################################################################
#################################################################
#################################################################

################ set up multivariate methods ####################


##### get weight matrix for tapering and svds of data matrices over the training period #####


Tap_dir = paste0(save_dir,"Tap/")
dir.create(Tap_dir,showWarnings = FALSE)

wm = weight_mat(DT,L = 2500)

num_loc = DT[year == min(year) & month == min(month)][!(is.na(SST_bar) | is.na(SST_hat)) ,.N]


# get full singular value decomposition of the tapered sample covariance matrix.
# This is only used for determining how many principal components should be used.
# for the sake of computational time we only consider the svd for the last year.

for(y in 2016)
{
  print(y)

  svd_by_month = function(m)
  {
    print(paste0("month =",m))

    train_years = DT[month == m][year < y & year > min(year),][,unique(year)]

    data_mat = matrix(DT[month == m][!(is.na(SST_bar) | is.na(SST_hat)) & year %in% train_years,
                                     SST_bar - SST_hat],
                      nrow = num_loc)

    sam_cov_mat = 1/length(train_years) * data_mat %*% t(data_mat)

    tap_cov_mat = wm * sam_cov_mat

    return(svd(tap_cov_mat))
  }
  sin_val_dec = parallel::mclapply(X = months,FUN = svd_by_month,mc.cores = mc_cores)

  save(sin_val_dec,file = paste0(Tap_dir,'svd_y',y,'.RData'))
}




###################################################
###################### PCA  #######################
###################################################

PCA_dir = paste0(save_dir,"PCA/")
dir.create(PCA_dir,showWarnings = FALSE)

# how many PCs should we use?

nPCs = c()

for(m in months)
{
  #plot(sin_val_dec[[m]]$d, main = paste0('month = ',m))
  load(paste0(Tap_dir,'svd_y',2016,'.RData'))

  sum_vec = cumsum(sin_val_dec[[which(months == m)]]$d)
  sum_tot = sum_vec[length(sum_vec)]

  nPCs = c(nPCs,which(sum_vec > 0.9*sum_tot)[1]) #use PCs explaining 90% of the variance in the final year

}


PCA_cov(DT,weight_mat = wm,
        Y = 1985:2016,
        M = months,
        nPCs = nPCs,
        save_years = validation_years,
        save_dir = PCA_dir)


###################################################
################## geostationary ##################
###################################################

GS_dir = paste0(save_dir, "GS/")
dir.create(GS_dir, showWarnings = FALSE)

geostationary_training(dt = DT,
                       training_years = training_years,
                       m = months,
                       save_dir = GS_dir,mc_cores = mc_cores)

########################################
################ ECC  ##################
########################################

ECC_dir = paste0(save_dir, "ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

############################################################################
############################################################################
############################################################################

################## generating multivariate forecasts #######################

# specifications for the desired forecasts:

fc_years = validation_years
fc_months = months
fc_ens_size = 50 # size of the forecast ensemble

mod_vec = c('PCA_mc','PCA_ac','GS','ECC')


###################################################
###################### PCA  #######################
###################################################

# regularization method, multiplicative correction

PCA_fc_mc = forecast_PCA_mult_corr(DT,
                                   Y = fc_years,
                                   M = fc_months,
                                   n = fc_ens_size,
                                   nPCs = nPCs,
                                   cov_dir = PCA_dir)

save(PCA_fc_mc,file = paste0(PCA_dir,"fc_mc.RData"))


# regularization method, additive correction

PCA_fc_ac = forecast_PCA_add_corr(DT,
                                  Y = fc_years,
                                  M = fc_months,
                                  n = fc_ens_size,
                                  nPCs = nPCs,
                                  cov_dir = PCA_dir)

save(PCA_fc_ac,file = paste0(PCA_dir,"fc_ac.RData"))


###################################################
################## geostationary ##################
###################################################

GS_fc = forecast_GS(DT,
                    Y = validation_years,
                    M = months,
                    n = fc_ens_size,
                    var_dir = GS_dir,
                    mc_cores = mc_cores)

save(GS_fc,file = paste0(GS_dir,"fc.RData"))

########################################
################ ECC  ##################
########################################

ECC_fc = forecast_ECC(DT,
                      Y = validation_years,
                      M = months,
                      ens_size = ens_size)

save(ECC_fc,file = paste0(ECC_dir,"fc.RData"))



########################################################################
########################################################################
########################################################################

###################### Validation: variogram scores ####################

# choose power for variogram score:

pp = 0.5

######### PCA #########

vs_PCA_mc = var_sc_par(dt_fc = PCA_fc_mc,  p = pp,
                       years = validation_years, ms = months, n = fc_ens_size,
                       save_dir = NULL,
                       mc_cores = mc_cores)


vs_PCA_ac = var_sc_par(dt_fc = PCA_fc_ac,  p = pp,
                       years = validation_years, ms = months, n = fc_ens_size,
                       save_dir = NULL,
                       mc_cores = mc_cores)



###### GS ######

vs_GS = var_sc_par(dt_fc = GS_fc,  p = pp,
                   years = validation_years, ms = months, n = fc_ens_size,
                   save_dir = NULL,
                   mc_cores = mc_cores)

###### ECC ######

vs_ECC = var_sc_par(dt_fc = ECC_fc,  p = pp,
                    years = validation_years, ms = months, n = ens_size,
                    save_dir = NULL,
                    mc_cores = mc_cores)

###########

# combine and save

vs_dt = data.table(PCA_mc=vs_PCA_mc[,mean(vs)],PCA_ac = vs_PCA_ac[,mean(vs)],GS = vs_GS[,mean(vs)], ECC = vs_ECC[,mean(vs)])

save(vs_dt,vs_PCA_ac,vs_PCA_mc,vs_GS,vs_ECC,
     file = paste0(save_dir,"vs.RData"))


########### permutation tests #####################


mod_vec = c('PCA_mc','PCA_ac','GS','ECC')
n_mod = length(mod_vec)

N = 1000

par( oma=c(0,6,6,0), mfrow=c(n_mod,n_mod), mar=c(3,3,2,2)+0.1 )

for(mod1 in mod_vec)
{

  for(mod2 in mod_vec)
  {
    vs_1 = get(paste0('vs_',mod1))
    vs_2 = get(paste0('vs_',mod2))

    perm_test_dt = merge(vs_1,vs_2,by=c('year','month'))
    setnames(perm_test_dt,c('vs.x','vs.y'), c(mod1,mod2))

    # permutation test for mod1 ~ mod2
    pt_vs = permutation_test_difference(perm_test_dt[,get(mod1)],perm_test_dt[,get(mod2)], N=N)

    x_lim_max = 1.1*max(abs(c(pt_vs$d_bar,pt_vs$D)))

    hist(pt_vs$D, xlim = c(-x_lim_max,x_lim_max),main = paste0(mod1,' vs. ',mod2),breaks = 20)

    qq = quantile(pt_vs$D,c(0.025,0.975))

    abline(v = pt_vs$d_bar,col = 'red')

    abline(v = qq,lty = 2)

    qq = quantile(pt_vs$D,c(0.1,0.9))

  }
}

###################################################################################
###################################################################################
###################################################################################

####################################
####### Case study: route ##########
####################################

# get route:

Bordeaux = c(-0.57,44.8)
Norfolk = c(-76.3,36.9)

p1 = data.table(Lon = Bordeaux[1], Lat = Bordeaux[2], Loc = 'Bordeaux')
p2 = data.table(Lon = Norfolk[1], Lat = Norfolk[2], Loc = "Norfolk")

route_name = "Bordeaux to Norfolk"
file_name = "scores_Bd_to_Nf"


# get grid ids_along this route: fix n and use gcIntermediate to find n coordinates on the route from p1 to p2 as the crow flies:

n = 1000
route = geosphere::gcIntermediate(p1[,.(Lon,Lat)],p2[,.(Lon,Lat)],n = n)

# find the grid_ids in DT closest to the coordinates on the route:

grid_id_dt = unique(DT[,.(Lon,Lat,grid_id)])

point_match = NULL
for(j in 1:dim(route)[1])
{
  a = geosphere::distHaversine(as.vector(route[j,]),as.matrix(grid_id_dt[,.(Lon,Lat)]))
  point_match[j] = which.min(a)
}

dt_route = unique(grid_id_dt[point_match,])
route_ids = dt_route[,grid_id]

###################################

# renaming for convenience

PCA_ac_fc = PCA_fc_ac
PCA_mc_fc = PCA_fc_mc

########################################

mod_vec = c('PCA_mc','PCA_ac','GS','ECC')

mod_vec_all = c('PFC',mod_vec)

# reduce data to route

DT_route = DT[grid_id %in% route_ids & year %in% validation_years,]

for(mod in mod_vec){
  assign(paste0(mod,'_fc_route'),get(paste0(mod,'_fc'))[grid_id %in% route_ids,])
}


########################################

# put considered functionals in a list, the paper shows just the result for the minimum

fun_nlist = list('max','min','mean')
fun_list = lapply(fun_nlist,get)


######################################################
###### Get MSEs for the  maximum along route ########
######################################################

# initialize data table for results:

validation_dt = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt,c("month","year"))

scores_dt = as.data.table(expand.grid(  model = mod_vec_all, fun = fun_nlist,MSE = 0,CRPS = 0))

# get observation

for(name in fun_nlist)
{
  fun = get(name)

  DT_route[,paste0('SST_',name) := fun(SST_bar, na.rm = TRUE), by = .(month,year)]

  temp = DT_route[grid_id == min(grid_id) ,eval(parse(text = paste0('SST_',name)))]
  validation_dt[,paste0('SST_',name) := temp]

  # get value of point forecast

  # the weird column name is for unification with other models
  DT_route[,paste0('PFC_',name,'_fc') := fun(SST_hat,na.rm = TRUE),by = .(month,year)]

  temp = DT_route[grid_id == min(grid_id) ,eval(parse(text = paste0('PFC_',name,'_fc')))]
  validation_dt[,paste0('PFC_',name,'_fc') := temp]

  # get value of multivariate forecasts for the different models:

  for(mod in mod_vec)
  {
    print(c(name,mod))

    # how large is the ensemble?
    es_mod = fc_ens_size
    if( mod == 'ECC')
    {
      es_mod = 9
    }


    fc_temp = get(paste0(mod,'_fc_route'))

    # get the forecasted values:
    for(i in 1:es_mod)
    {
      fc_temp[,paste0(mod,'_',name,'_',i) := lapply(X = .SD,FUN = fun,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]
    }

    # take mean over all forecasts and put this into validation_dt

    fc_temp = fc_temp[grid_id == min(grid_id) ,.SD,.SDcols = paste0(mod,'_',name,'_',1:es_mod)]
    fc_temp[,paste0(mod,'_',name,'_fc') := rowMeans(.SD),.SDcols = paste0(mod,'_',name,'_',1:es_mod)]

    temp = fc_temp[,.SD,.SDcols = paste0(mod,'_',name,'_fc')]

    validation_dt = data.table(validation_dt,temp)

    ### get MSE ###

    obs = validation_dt[,eval(parse(text = paste0('SST_',name)))]
    fcs = validation_dt[,eval(parse(text = paste0(mod,'_',name,'_fc')))]

    validation_dt[,paste0('MSE_',mod,'_',name,'_fc'):=(obs-fcs)^2]

    mse = mean((obs - fcs)^2)

    scores_dt[model == mod & fun == name,MSE:=mse]

    ### get CRPS ###

    data = as.matrix(fc_temp[,.SD,.SDcols = paste0(mod,'_',name,'_',1:es_mod)])

    crps = scoringRules::crps_sample(y = obs,
                                     dat = data)

    validation_dt[,paste0('CRPS_',mod,'_',name,'_fc'):=crps]

    scores_dt[model == mod & fun == name, CRPS := mean(crps)]
  }

  # scores for point forecasts:
  mod = 'PFC'

  obs = validation_dt[,eval(parse(text = paste0('SST_',name)))]
  fcs = validation_dt[,eval(parse(text = paste0(mod,'_',name,'_fc')))]

  mse = mean((obs - fcs)^2)
  crps = mean(abs(obs-fcs))

  scores_dt[model == mod & fun == name,MSE:=mse]
  scores_dt[model == mod & fun == name,CRPS:=crps]
}


save(scores_dt,file = paste0(save_dir,'route_scores.RData'))

### permutation tests ###

funn = 'min'
score = 'MSE'

N = 10000

#pdf(paste0(plot_dir,'perm_test_routescores.pdf'),width = 14,height = 14)

par( oma=c(0,6,6,0), mfrow=c(length(mod_vec),length(mod_vec)), mar=c(3,3,2,2)+0.1 )

for(mod1 in mod_vec)
{

  for(mod2 in mod_vec)
  {
    col_name_1 = paste0(score,'_',mod1,'_',funn,'_fc')
    col_name_2 = paste0(score,'_',mod2,'_',funn,'_fc')

    # permutation test for mod1 ~ mod2
    pt_vs = permutation_test_difference(validation_dt[,get(col_name_1)],validation_dt[,get(col_name_2)], N=N)

    x_lim_max = 1.1*max(abs(c(pt_vs$d_bar,pt_vs$D)))

    hist(pt_vs$D, xlim = c(-x_lim_max,x_lim_max),main = paste0(mod1,' vs. ',mod2),breaks = 20)

    qq_1 = quantile(pt_vs$D,c(0.05))

    abline(v = pt_vs$d_bar,col = 'red')

    abline(v = qq_1,lty = 2)

    qq_2 = quantile(pt_vs$D,c(0.01))

    abline(v = qq_2,lty = 3)
  }
}


### get p-value for all tests ###

library(stats)

p_values_route_scoring = as.data.table(expand.grid(MOD1 = mod_vec,MOD2 = mod_vec,SCORE = c('MSE','CRPS'),FUN = c('max','min')))

p_values_route_scoring[,p_value:=NA_real_]

for(funn in c('max','min') )
{
  for(score in c('CRPS','MSE'))
  {
    print(c(funn,score))

    for(mod1 in mod_vec)
    {

      for(mod2 in mod_vec)
      {
        col_name_1 = paste0(score,'_',mod1,'_',funn,'_fc')
        col_name_2 = paste0(score,'_',mod2,'_',funn,'_fc')

        # permutation test for mod1 ~ mod2
        pt_vs = permutation_test_difference(validation_dt[,get(col_name_1)],validation_dt[,get(col_name_2)], N=N)

        p_val_f = ecdf(pt_vs$D)
        p_val = p_val_f(pt_vs$d_bar)

        p_values_route_scoring[MOD1 == mod1 & MOD2 == mod2 & SCORE == score & FUN == funn, p_value := p_val]

      }
    }

  }
}

save(p_values_route_scoring,file = paste0(save_dir,'p_values.RData'))


#################################
#################################
#################################

 #### multivariate rank histograms ####

##################

brks = 6

### PCA ###

# get forecast:


rks_PCA_mc_route = mv_rank_hist_new(PCA_mc_fc_route, fc_ens_size = fc_ens_size,
                                    mc_cores = mc_cores,
                                    breaks = brks, mn = "PCA_mc rank histograms",
                                    save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_PCA_mc_route")





rks_PCA_ac_route = mv_rank_hist_new(PCA_ac_fc_route, fc_ens_size = fc_ens_size,
                                    mc_cores = mc_cores,
                                    breaks = brks, mn = "PCA_ac rank histograms",
                                    save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_PCA_ac_route")


### geostationary ###

rks_GS_route = mv_rank_hist_new(GS_fc_route, fc_ens_size = fc_ens_size,
                                mc_cores = mc_cores,
                                breaks = brks, mn = "GS rank histograms",
                                save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_GS_route")


### ECC ###

rks_ECC_route = mv_rank_hist_new(ECC_fc_route, fc_ens_size = ens_size,
                                 mc_cores = mc_cores,
                                 breaks = brks, mn = "ECC rank histograms",
                                 save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_ECC_route")





### save entire workspace ###

save.image(file = paste0(save_dir,'results.RData'))

