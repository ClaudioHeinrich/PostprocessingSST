

#############################################################################################

###### master script part 1 - setting up and creating or loading of wide data table #########

#############################################################################################

# This script sets up for a full post-processing analysis 
# for a lon/lat window to be specified below

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

# choose your favourite area for analysis and give it a name abbreviation

#NAO_2:
 # lat_box = c(40,70)
 # lon_box = c(-60,-30)

#NAO:
# lat_box = c(30,70)
# lon_box = c(-70,-25)

#Pres_Bergen: 
# lat_box = c(50,80)
# lon_box = c(-5,30)

#Europe:
# lat_box = c(30,70)
# lon_box = c(-25,45)



name_abbr = "NAO_3" # for northern atlantic ocean

 lat_box = c(50,85)
 lon_box = c(-20,40)


ens_size = 9 # size of forecast ensemble

training_years = 1985:2011
validation_years = 2012:2016 # all previous years are used for training 
months = 1:12

# create directories

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)

### construct or load wide data set ###

# takes time, avoid if possible: if the data hasn't changed and you're just trying out a new window, just run this:

DT = load_combined_wide()[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]

# prep data further:
DT = DT[month %in% months,]
DT[,YM := 12*year + month]
DT = DT[order(year,month,Lon,Lat)]

setcolorder(DT,c("year","month",'Lon','Lat','YM','grid_id','SST_bar','Ens_bar','Ens_sd',paste0('Ens',1:ens_size)))

#DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"))[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]

# make_combined_wide_dataset(lat_box = lat_box,
#                            lon_box = lon_box,
#                            output_loc = save_dir,
#                            output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))


# save everything:



save.image(file = paste0(save_dir,"setup.RData"))

