
##################################################################################################

###################  master script part 4 - setting up multivariate models  ######################

##################################################################################################

# This script sets up the different methods of multivariate post-processing, by computing and saving covariance estimates etc.
# As part of this script you should choose whether to consider centered or even standardized variables.
#
# Files generated:
#   
# Data files: PCA/sam_cov_m+_y++.RData, SE/cov_est_SE_m+_y++.RData, 
#             where + labels the months considered and ++ the years in the validation period.
#
# Requires previous run of 03.master.bias.correct 
# with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

time_s4 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO/lv/2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


##### centered variables? #####

# decide whether you want to work with SST, or with SST centered around climatology, or with SST standardized w.r.t. climatology 

SST = ""  # takes 'centered','standardized', or ''

clim_years = training_years


if(SST == "centered")
{
  DT = dt_transform_center(DT,clim_years)
  
  name_abbr = paste0(name_abbr,"/centered" )
  
  save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
  dir.create(save_dir,showWarnings = FALSE)
  
  plot_dir = paste0("./figures/", name_abbr,"/")
  dir.create(plot_dir, showWarnings = FALSE)
  
}
if(SST == "standardized")
{
  DT = dt_transform_stan(DT,clim_years)
  
  name_abbr = paste0(name_abbr,"/standardized")
  
  save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
  dir.create(save_dir,showWarnings = FALSE)
  
  plot_dir = paste0("./figures/", name_abbr,"/")
  dir.create(plot_dir, showWarnings = FALSE)
  
}


###################################################
################## geostationary ##################
###################################################

GS_dir = paste0(save_dir, "GS/")
dir.create(GS_dir, showWarnings = FALSE)

geostationary_training(dt = DT, 
                       training_years = training_years,
                       m = months,
                       save_dir = GS_dir)


# specifications for the desired forecasts:

fc_years = validation_years
fc_months = months
fc_ens_size = 500


###################################################
################## geostationary ##################
###################################################

GS_fc = forecast_GS(DT,
                    Y = validation_years,
                    M = months,
                    n = fc_ens_size,
                    var_dir = GS_dir,
                    mc_cores = mc_cores)

save(GS_fc,file = paste0(GS_dir,"fc.RData"))
rm(GS_fc)

gc()


save.image(file = paste0(save_dir,"setup.RData"))



