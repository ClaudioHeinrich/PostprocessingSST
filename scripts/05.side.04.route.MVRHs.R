##### setting up ######


rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

time_s54 = proc.time()

mc_cores = 5

##################

brks = ens_size + 2

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


time_s54 = proc.time() - time_s54

save.image(file = paste0(save_dir,"setup.RData"))

y = 2010

for(m in 1:12){
test = matrixStats::rowRanks(as.matrix(GS_fc_route[year == y &month == m & !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
hist(test[,1],main = c(m,y))
}


test2 = matrixStats::rowRanks(as.matrix(PCA_mc_fc_route[ !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
hist(test2[,1],main = c(m,y))

test3 = matrixStats::rowRanks(as.matrix(PCA_ac_fc_route[ !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
hist(test3[,1],main = c(m,y))

test4 = matrixStats::rowRanks(as.matrix(PCA_mc_fc_route[ !is.na(SST_hat),.SD,.SDcols = c('SST_bar',paste0('fc',1:500))]))
hist(test4[,1],main = c(m,y))


