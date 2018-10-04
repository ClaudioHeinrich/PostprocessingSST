
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

name_abbr = "Atl" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))


##### centered variables? #####

# decide whether you want to work with SST, or with SST centered around climatology, or with SST standardized w.r.t. climatology 

SST = "standardized"  # takes 'centered','standardized', or ''

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
###################### PCA  #######################
###################################################

PCA_dir = paste0(save_dir,"PCA/")
dir.create(PCA_dir,showWarnings = FALSE)

# get sample covariance matrices

wm = weight_mat(DT) # for weighting down by distance

sam_cov(DT,weight_mat = wm, 
        M = months,
        save_years = validation_years,
        save_dir = PCA_dir)

##################################################
###################### SE  #######################
##################################################

SE_dir = paste0(save_dir,"SE/")
dir.create(SE_dir,showWarnings = FALSE)

# get covariance estimates

cov_est_SE(DT, weight_mat = wm,
           M = months,
           save_years = validation_years,
           save_dir = SE_dir)

###################################################
################## geostationary ##################
###################################################

GS_dir = paste0(save_dir, "GS/")
dir.create(GS_dir, showWarnings = FALSE)

geostationary_training(dt = DT, 
                       training_years = training_years,
                       m = months,
                       save_dir = GS_dir)

########################################
################ ECC  ##################
########################################

ECC_dir = paste0(save_dir, "ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

#####

time_s4 = proc.time() - time_s4

save.image(file = paste0(save_dir,"setup.RData"))
