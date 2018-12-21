

##### setting up ######

rm(list = ls())

time_s31 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = 'NAO/lv' 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

mc_cores = 4

###### compare to linear regression models ######

brks = ens_size + 2

### PCA ###

# get forecast:
load(paste0(PCA_dir,"fc_mc.RData"))

rks_PCA_mc = mv_rank_hist_new(PCA_fc_mc, fc_ens_size = fc_ens_size,
                              mc_cores = mc_cores,
                              breaks = brks, mn = "PCA_mc rank histograms",
                              save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_PCA_mc")


save.image(file = paste0(save_dir,"setup.RData"))
