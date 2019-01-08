
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

name_abbr = "NAO_small" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

time_s51 = proc.time()

mc_cores = 5

##################

brks = ens_size + 2

### PCA ###

# get forecast:
load(paste0(PCA_dir,"fc_mc.RData"))

rks_PCA_mc = mv_rank_hist_new(PCA_fc_mc, fc_ens_size = fc_ens_size,
                              mc_cores = mc_cores,
                              breaks = brks, mn = "PCA_mc rank histograms",
                              save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_PCA_mc")



load(paste0(PCA_dir,"fc_ac.RData"))

rks_PCA_ac = mv_rank_hist_new(PCA_fc_ac, fc_ens_size = fc_ens_size, 
                              mc_cores = mc_cores,
                              breaks = brks, mn = "PCA_ac rank histograms",
                              save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_PCA_ac")


### SE ###

# # get forecast:
# load(paste0(SE_dir,"fc.RData"))
# 
# mv_rank_hist(SE_fc, fc_ens_size = fc_ens_size, breaks = brks, mn = "SE rank histograms",
#              save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_SE")

### geostationary ###

# get forecast:
load(paste0(GS_dir,"fc.RData"))

rks_GS = mv_rank_hist_new(GS_fc, fc_ens_size = fc_ens_size, 
                          mc_cores = mc_cores,
                          breaks = brks, mn = "GS rank histograms",
                          save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_GS")


### ECC ###

# get forecast:
load(paste0(ECC_dir,"fc.RData"))


rks_ECC = mv_rank_hist_new(ECC_fc, fc_ens_size = ens_size, 
                           mc_cores = mc_cores,
                           breaks = brks, mn = "ECC rank histograms",
                           save_pdf = TRUE, plot_dir = plot_dir, file_name = "rank_histo_ECC_new")


time_s51 = proc.time() - time_s51

save.image(file = paste0(save_dir,"setup.RData"))
