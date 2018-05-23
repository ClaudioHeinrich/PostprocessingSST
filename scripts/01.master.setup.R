

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

lat_box = c(40,70)
lon_box = c(-60,-30)

# lat_box = c(-20,0)
# lon_box = c(-110,-75)


name_abbr = "NAO_small" # for northern atlantic ocean

ens_size = 9 # size of forecast ensemble

validation_years = 2001:2010 # all previous years are used for training 


# create directories

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)

### construct or load wide data set ###

# takes time, avoid if possible: if the data hasn't changed and you're just trying a new window, just run this:

DT = load_combined_wide()[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]

# make_combined_wide_dataset(lat_box = lat_box,
#                            lon_box = lon_box,
#                            output_loc = save_dir,
#                            output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))


# save everything:
save.image(file = paste0(save_dir,"setup.RData"))

