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

time_s2 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_small" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


####################################################

past_0 = 5 # how many years before validation period do we consider for the RMSEs
num_years = max(validation_years) - 1985


for(m in months)
{
  print( c('month ',m))
  
  DT = get(paste0('DT',m))

  ###### run bias analysis for simple moving averages ######
  
  win_length = 1 : (num_years-1)
  
  MSE_by_par = function(k)
  {
    temp = bias_correct_2(dt = DT, 
                          method = "sma", 
                          par_1 = k,
                          reduced_output = TRUE)
    return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SST_hat])
  }
  
  
  BC = parallel::mclapply(X = win_length, FUN = MSE_by_par, mc.cores = mc_cores)
  
  # restrict to validation years + the past_0 years before 
  
  Bias_est_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                   .(year,month,Lon,Lat,SST_bar)]
  
  for(k in win_length)
  {
    Bias_est_dt[,paste0('l',k):= BC[[k]]]
    Bias_est_dt[,paste0('err',k):= (SST_bar - eval(parse(text = paste0('l',k))))^2]
  }
  
  #get MSE by year and month, for each year based on all previous years contained in Bias_est_dt
  
   mean_sc_by = function(y)
    {
      temp = data.table(year = y,Bias_est_dt[year < y, lapply(X = .SD,FUN = mean,na.rm = TRUE), .SDcols = paste0('err',win_length)])
      return(temp)     
    }
    
  assign(x = paste0('sc_sma_m',m), rbindlist(parallel::mclapply(X = validation_years,FUN = mean_sc_by,mc.cores = mc_cores)))
  
  save(paste0('sc_sma_m',m), file = paste0(save_dir,'scores.bc.sma.m',m,'.RData'))
  
 
  ###### run bias analysis for simple moving averages ######
  
  par_vec = seq(0.05,0.4,length.out = 24) 
  
  MSE_by_par = function(a)
  {
    temp = bias_correct_2(dt = DT, 
                          method = "ema", 
                          par_1 = a,
                          reduced_output = TRUE)
    return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SST_hat])
  }
  
  
  BC = parallel::mclapply(X = par_vec, FUN = MSE_by_par, mc.cores = mc_cores)
  
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
  
  mean_sc_by = function(y)
  {
    temp = data.table(year = y,Bias_est_dt[year < y, lapply(X = .SD,FUN = mean,na.rm = TRUE), .SDcols = paste0('err',round(par_vec,4))])
    return(temp)     
  }
  
  sc_ema_m = rbindlist(parallel::mclapply(X = validation_years,FUN = mean_sc_by,mc.cores = mc_cores))
  
  assign(x = paste0('sc_ema_m',m),value = sc_ema_m)
  
  save(eval(parse(text = paste0('sc_ema_m',m))), file = paste0(save_dir,'scores.bc.ema.m',m,'.RData'))

}
  




###### run bias analysis for exponential moving averages ######

par_vec = seq(0.05,0.4,length.out = 24) 

MSE_by_par = function(a)
{
  temp = bias_correct_2(dt = DT, 
                        method = "ema", 
                        par_1 = a,
                        saveorgo = FALSE,
                        reduced_output = TRUE)
  return(temp[year %between% c(min(validation_years - past_0),max(validation_years)),SST_hat])
}


BC = parallel::mclapply(X = par_vec, FUN = MSE_by_par,mc.cores = mc_cores)

# restrict to validation years + the past_0 years before 

Bias_est_dt = DT[year %between% c(min(validation_years - past_0),max(validation_years)),
                 .(year,month,Lon,Lat,SST_bar)]

for(a in par_vec)
{ind = which(par_vec == a )
  print(ind)
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


###### plotting scores for different ways of bias correction ######

load(paste0(save_dir,"scores.bc.sma.RData"))
load(paste0(save_dir,"scores.bc.ema.RData"))

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

# get data for last year
y = max(validation_years)


row_sma = msc_sma[year == y,-1,with = FALSE]
row_ema = msc_ema[year == y,-1,with = FALSE]

y_range = range(list(row_sma[,-'min_l',with = FALSE],row_ema[,-'min_a',with = FALSE]))  

y_range = sqrt(y_range)

## plot for sma ##

pdf(paste0(plot_dir,"mean_scores_sma.pdf"))
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
dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"mean_scores_ema.pdf"))
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


# bias correction year by year

for(y in validation_years)
{
  print(y) 
  temp = bias_correct_2(dt = DT,
                         method = opt_par[year == y, method],
                         par_1 = opt_par[year == y,par],
                         reduced_output = TRUE)[year == y,]
  DT[year == y,][,Bias_Est := temp[,Bias_Est]][,SST_hat := temp[,SST_hat]]
}

rm(temp)

# For the training years the bias correction considers also the future, 
# and uses the parameter and method estimated for the first validation year

method_ty = c(opt_par[year == min(year),method],opt_par[year == min(year),par])

DT = bias_correct_training(dt = DT,
                           method = method_ty,
                           training_years = training_years,
                           save_dir = save_dir)



time_s2 = proc.time() - time_s2

#### save ####

save.image(file = paste0(save_dir,"setup.RData"))
