
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

time_s3 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

###### getting sample variances of ensemble  ######

DT[,var_bar := (SST_hat - SST_bar)^2]

###### running variance estimation with simple moving averages ######

past_0 = 5 # how many years before validation period do we consider for the 

num_years = max(validation_years) - DT[,min(year)]

win_length = 1 : (num_years-1)

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
{print(k)
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



###### run variance estimation by exponential moving averages ######

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

load(paste0(save_dir,"scores.bc.sd.sma.RData"))
load(paste0(save_dir,"scores.bc.sd.ema.RData"))

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

## plot for sma ##

pdf(paste0(plot_dir,"mean_scores_sd_sma.pdf"))
plot(x = win_length, 
     y = as.vector(row_sma[,-c('min_crps','min_l'),with = FALSE]),
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("CRPS for variance estimation by SMA"),
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
dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"mean_scores_sd_ema.pdf"))
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

for(y in validation_years)
{
  print(y) 
  temp = sd_est_2(dt = DT,
                        method = opt_par_sd[year == y, method],
                        par_1 = opt_par_sd[year == y,par],
                        reduced_output = TRUE)[year == y,]
  DT[year == y,][,SD_hat := temp[,SD_hat]]
}

rm(temp)

# For the training years the variance estimation uses the parameter and method estimated for the first validation year

method_ty = c(opt_par_sd[year == min(year),method],opt_par_sd[year == min(year),par])

temp = sd_est(dt = DT[year %in% training_years],
              method = method_ty[1],
              par_1 = as.numeric(method_ty[2]),
              saveorgo = FALSE)


DT[year %in% training_years, SD_hat:=temp[,SD_hat]]


time_s3 = proc.time() - time_s3

#### save ####

save.image(file = paste0(save_dir,"setup.RData"))


