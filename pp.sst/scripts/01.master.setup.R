

#############################################################################################

###### master script part 1 - setting up and creating or loading of wide data table #########

#############################################################################################

# This script sets up a run of the post-processing analysis 

# installing the R-package:

rm(list = ls())

install.packages('devtools')
library(devtools)

install_github('ClaudioHeinrich/PostprocessingSST/pp.sst')


#start timer:
time_s1 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(pp.sst)


# choose an abbreviation for this run and give it a description, see the README file for more details. 

name_abbr = "test"

description = 'very small test window. For script testing'

# specify your directories:

data_dir = "~/PostClimDataNoBackup/SFE/Derived/" # where is your prerendered data table stored?

# Directory for derived datasets, this should change when you change name_abbr
save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

# Directory for plots, this should change when you change name_abbr
plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)


# choose the area of the globe to consider:

lat_box = c(55,65)
lon_box = c(-10,5)


 # a couple of parameters:

ens_size = 9 # size of forecast ensemble

training_years = 1985:2000
validation_years = 2001:2016 

months = 1:12

mc_cores = 6 # number of cores for parallelization

### construct or load data set ###

DT = load_combined_wide(data_dir = data_dir)[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]

# tidy up DT:

setkey(x = DT,year,month,Lon,Lat)

DT = DT[month %in% months & year %in% c(training_years,validation_years),]
DT[,YM := 12*year + month]
DT = DT[order(year,month,Lon,Lat)]

setcolorder(DT,c("year","month",'Lon','Lat','YM','grid_id','SST_bar','Ens_bar','Ens_sd'))

#### time, update script counter, save ####

time_s1 = proc.time() - time_s1

script_counter = 1

save.image(file = paste0(save_dir,"setup.RData"))


# The following line needs to be run if the full data table has not yet been created, i.e. if new data is available:

# make_combined_wide_dataset(lat_box = lat_box,
#                            lon_box = lon_box,
#                            output_loc = save_dir,
#                            output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))


