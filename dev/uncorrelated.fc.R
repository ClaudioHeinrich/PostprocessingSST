

# generate a totally decorrelated prediction and generate the rank histograms for comparison



# Data files: PCA/fc.RData, SE/fc.RData, GS/fc.RData, ECC/fc.Rdata
#
# Requires previous run of 04.master.prep.multivariate.pp.R 
# with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())



setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# specifications for the desired forecasts:

fc_years = validation_years
fc_months = months
fc_ens_size = 500

mod_vec = c('PCA_mc','PCA_ac','GS','ECC')


#################

forecast_uncorrelated = function(DT,Y,M,n)
{
  DT_water = DT[year %in% Y & month %in% M][!is.na(SST_hat)]
  
  
  DT_water_new = data.table()
  for(y in Y)
  {
    for(m in M)
    { print(c(y,m))
      fc_ym = DT_water[year == y & month == m]
      
      nloc = fc_ym[,.N]
      
      mu = fc_ym[,SST_hat]
      sd = fc_ym[,SD_hat]
      
      fc = trc(rnorm(n = n*nloc,mean = mu,sd = sd))
      
      fc_dt = as.data.table(matrix(fc,nrow = nloc))
      setnames(fc_dt,paste0('fc',1:n))
      
      fc_ym = data.table(fc_ym,fc_dt)
      DT_water_new = rbindlist(list(DT_water_new,fc_ym))
    }
  }
  
  DT_new = merge(DT_water_new,DT,by = colnames(DT),all.y = TRUE)
  return(DT_new)
}


#perfectly correlated forecast

forecast_pc = function(DT,Y,M,n)
{
  DT_water = DT[year %in% Y & month %in% M][!is.na(SST_hat)]
  
  
  DT_water_new = data.table()
  for(y in Y)
  {
    for(m in M)
    { print(c(y,m))
      fc_ym = DT_water[year == y & month == m]
      
      nloc = fc_ym[,.N]
      
      mu = fc_ym[,SST_hat]
      sd = fc_ym[,SD_hat]
      
      fc = trc(rnorm(n = n*nloc,mean = mu,sd = sd))
      
      fc_dt = as.data.table(matrix(fc,nrow = nloc))
      setnames(fc_dt,paste0('fc',1:n))
      
      fc_ym = data.table(fc_ym,fc_dt)
      DT_water_new = rbindlist(list(DT_water_new,fc_ym))
    }
  }
  
  DT_new = merge(DT_water_new,DT,by = colnames(DT),all.y = TRUE)
  return(DT_new)
}

###################################################
###################### generate forecast  #######################
###################################################

uc_fc = forecast_uncorrelated(DT[year %in% validation_years], 
                                   Y = validation_years,
                                   M = months,
                                   n = fc_ens_size)


# get multivariate rank histograms

rks_uc = mv_rank_hist_new(uc_fc, fc_ens_size = fc_ens_size, 
                          mc_cores = mc_cores,
                          breaks = brks, mn = "UC rank histograms",
                          save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_uc")

save(uc_fc,rks_uc,file = paste0(save_dir,"fc_uc.RData"))


