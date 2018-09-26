
##########################################################################################

###################  master script part 5 - multivariate modelling  ######################

##########################################################################################

# This script generates and saves forecasts for all methods
# 
# Files generated:
#   
# Data files: PCA/fc.RData, SE/fc.RData, GS/fc.RData, ECC/fc.Rdata
#
# Requires previous run of 04.master.prep.multivariate.pp.R 
# with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "Atl" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# specifications for the desired forecasts:

fc_years = validation_years
fc_months = months
fc_ens_size = 100



###################################################
###################### PCA  #######################
###################################################

PCA_fc = forecast_PCA(DT, 
                      Y = fc_years,
                      M = fc_months,
                      n = fc_ens_size,
                      cov_dir = PCA_dir)

save(PCA_fc,file = paste0(PCA_dir,"fc.RData"))

rm(PCA_fc)
gc()

##################################################
###################### SE  #######################
##################################################

SE_fc = forecast_SE(DT,
                    Y = fc_years,
                    M = fc_months,
                    n = fc_ens_size,
                    cov_dir = SE_dir)

save(SE_fc,file = paste0(SE_dir,"fc.RData"))
rm(SE_fc)
gc()

###################################################
################## geostationary ##################
###################################################

GS_fc = forecast_GS(DT,
                    Y = validation_years,
                    M = months,
                    n = fc_ens_size,
                    var_dir = GS_dir)

save(GS_fc,file = paste0(GS_dir,"fc.RData"))
rm(GS_fc)
gc()

########################################
################ ECC  ##################
########################################

ECC_fc = forecast_ECC(DT,
                      Y = validation_years,
                      M = months,
                      ens_size = ens_size)

save(ECC_fc,file = paste0(ECC_dir,"fc.RData"))
rm(ECC_fc)

#####

save.image(file = paste0(save_dir,"setup.RData"))
