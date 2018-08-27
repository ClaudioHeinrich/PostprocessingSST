
##########################################################################################

###################  master script part 4 - multivariate modelling  ######################

##########################################################################################

# This script sets up the different methods of multivariate post-processing, by computing and saving covariance estimates etc.
# Generates and saves forecasts for all methods
# 
# Files generated:
#   
# Data files: 
#
# Requires previous run of 03.master.bias.correct 
# with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_3" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# specifications for the desired forecasts:

fc_years = validation_years
fc_months = months
fc_ens_size = 50

##### centered variables? #####

# decide whether you want to work with SST, or with SST centered around climatology, or with SST standardized w.r.t. climatology 

SST = ""  # takes 'centered','standardized', or ''

clim_years = training_years


if(SST == "centered")
{
  DT = dt_transform_center(DT,clim_years)
}
if(SST == "standardized")
{
  DT = dt_transform_stan(DT,clim_years)
}


###################################################
###################### PCA  #######################
###################################################

PCA_dir = paste0(save_dir,"PCA/")
dir.create(PCA_dir,showWarnings = FALSE)

# get sample covariance matrices

wm = weight_mat(DT)

sam_cov(DT,weight_mat = wm, 
        M = months,
        save_years = validation_years,
        save_dir = PCA_dir)

PCA_fc = forecast_PCA(DT, 
                      Y = fc_years,
                      M = fc_months,
                      n = fc_ens_size,
                      cov_dir = PCA_dir)

save(PCA_fc,file = paste0(PCA_dir,"fc.RData"))

rm(PCA_fc)

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

SE_fc = forecast_SE(DT,
                    Y = fc_years,
                    M = fc_months,
                    n = fc_ens_size,
                    cov_dir = SE_dir)

save(SE_fc,file = paste0(SE_dir,"fc.RData"))
rm(SE_fc)

###################################################
################## geostationary ##################
###################################################

GS_dir = paste0(save_dir, "GS/")
dir.create(GS_dir, showWarnings = FALSE)

geostationary_training(dt = DT, 
                       training_years = training_years,
                       m = months,
                       save_dir = GS_dir)

GS_fc = forecast_GS(DT,
                    Y = validation_years,
                    M = months,
                    n = fc_ens_size,
                    var_dir = GS_dir)

save(GS_fc,file = paste0(GS_dir,"fc.RData"))
rm(GS_fc)


########################################
################ ECC  ##################
########################################

ECC_dir = paste0(save_dir, "ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

ECC_fc = forecast_ECC(DT,
                      Y = validation_years,
                      M = months,
                      ens_size = ens_size)

save(ECC_fc,file = paste0(ECC_dir,"fc.RData"))
rm(ECC_fc)


#####

save.image(file = paste0(save_dir,"setup.RData"))
