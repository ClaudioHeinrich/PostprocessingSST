

#############################################################################################

###### master script part 1 - setting up and creating or loading of wide data table #########

#############################################################################################

# This script sets up for a full post-processing analysis 
# for a lon/lat window to be specified below

rm(list = ls())

time_s1 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


# choose an abbreviation for this run

name_abbr = "NAO_small"

# choose the area of the globe to consider
lat_box = c(40,75)
lon_box = c(-30,50)


# # set a couple of parameters

mc_cores = 6

ens_size = 9 # size of forecast ensemble

training_years = 1985:2000
validation_years = 2001:2016 

months = 1:12

##### setting up complete - now move to creation #####

# create directories

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)

### construct load wide data set ###

DT = load_combined_wide()[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]

# tidy up DT:

DT = DT[month %in% months,]
DT[,YM := 12*year + month]
DT = DT[order(year,month,Lon,Lat)]

setcolorder(DT,c("year","month",'Lon','Lat','YM','grid_id','SST_bar','Ens_bar','Ens_sd',paste0('Ens',1:ens_size)))



time_s1 = proc.time() - time_s1
# save everything:

save.image(file = paste0(save_dir,"setup.RData"))


# The following line needs to be run if the full data table has not yet been created, i.e. if new data is available:

# make_combined_wide_dataset(lat_box = lat_box,
#                            lon_box = lon_box,
#                            output_loc = save_dir,
#                            output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))


