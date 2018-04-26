
################################################

# this example computes the variogram score for the northern atlantic ocean in January
# for the geostationary model with exponential covariance function.

################################################

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/NAO/")

load(file = paste0(save_dir,"example.snippet.variogram.score.RData"))

### the following is creating and saving the data used in the example and can be skipped
# 
# m = 1 # month
# 
# 
# data_dir = paste0(save_dir,"GeoStat/")
# 
# # get data table
# dt = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_wide_bc.RData"))
# 
# # --- create data table with pair of coordinates: ---
# 
# dt_coor_1 = dt[YM == min(YM), .(Lat,Lon,grid_id)]
# dt_coor_1[,grid_id_ind := match(grid_id,sort(grid_id))]
# setkey(dt_coor_1,Lat,Lon) #ordering needs to be the same as in geostationary_training
# setnames(dt_coor_1,c("Lat1","Lon1","grid_id1","grid_id_ind1"))
# 
# # add dummy key, do outer full merge with a duplicate, and remove key again
# 
# dt_coor_1[,"key" := 1]
# dt_coor_2 = copy(dt_coor_1) 
# setnames(dt_coor_2,c("Lat","Lon2","grid_id2","grid_id_ind2","key"))
# var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
# var_sc_prereq[, "key" := NULL]
# 
# # --- get parameters for the fitted exponential covariance model: ---
# 
# var_file_name = "variogram_exp_m1.RData"
# load(paste0(data_dir,var_file_name))
# 
# psill <- Mod$psill[2]
# range <- Mod$range[2]
# nugget <- Mod$psill[1]
# 
# Sigma <- psill*exp(-Dist/range) # covariance matrix
# 
# save.image(file = paste0(save_dir,"example.snippet.variogram.score.RData"))
# 


##################################################

# the following lines write the variances contained in Sigma at the corresponding location in var_sc_prereq: 
# This takes forever when Sigma is as large as here: 
dim(Sigma)
# Hopefully we can speed this up by using Rcpp...

print("getting variances")

dummy_vec = c()
for(i in 1:var_sc_prereq[,.N])
{if(i %% 10 == 0) print(paste0(i," / ",var_sc_prereq[,.N]))
  dummy_vec = c(dummy_vec, Sigma[var_sc_prereq[,grid_id_ind1][i],var_sc_prereq[,grid_id_ind2][i]])
}
var_sc_prereq[,Var := dummy_vec]
