
###############################################################################

#############  side script 5.1 - multivariate rank histograms  ################

###############################################################################

# This script compares the multivariate post-processing methods by means of their multivariate rank histograms
# 
# Plots generated: rank_histo_PCA.pdf, rank_histo_SE.pdf, rank_histo_GS.pdf, rank_histo_ECC.pdf,
#   
# Requires previous run of 04.master.multiv.pp.R with the same value of name_abbr as below.

##### setting up ######


rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Full/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

time_s51 = proc.time()

mc_cores = 3

sma_par = data.table(year = validation_years,method = NA_character_,par = NA_real_)


for(y in validation_years)
{
  sma_par[year == y, method := 'sma']
  sma_par[year == y,par := msc_sma[year == y,min_l]]
}


# bias correction year by year

DT[,Bias_Est_sma := 0][,SST_hat_sma := 0]

for(y in validation_years)
{
  print(y) 
  temp = bias_correct_2(dt = DT,
                        method = sma_par[year == y, method],
                        par_1 = sma_par[year == y,par])[year == y,]
  DT[year == y, Bias_Est_sma := temp[,Bias_est]][year == y,SST_hat_sma := temp[,SST_hat]]
}

rm(temp)

save.image(file = paste0(save_dir,"setup.RData"))
