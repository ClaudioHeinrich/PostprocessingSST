
###################################################################

#############  side script 4.2 - variogram scores  ################

###################################################################

# This script compares the multivariate post-processing methods based on their variogram score.
# 
# Files generated: 
# data files: vs.RData
# plots: boxplot_vs.pdf
#
# Requires previous run of 04.master.multivar.pp.R with the same value of name_abbr as below.


##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_3" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

##################

# choose power for variogram score:

p = 0.5

### PCA ###

load(paste0(PCA_dir,"fc.RData"))

vs_PCA = var_sc_par(fc_dt = PCA_fc,  p = p, save_dir = NULL)


### SE ###

load(paste0(SE_dir,"fc.RData"))

vs_SE = var_sc_par(fc_dt = SE_fc,  p = p, save_dir = NULL)
  
### GS ###

load(paste0(GS_dir,"fc.RData"))

vs_GS = var_sc_par(fc_dt = GS_fc,  p = p, save_dir = NULL)

### ECC ###

load(paste0(ECC_dir,"fc.RData"))

vs_ECC = var_sc_par(fc_dt = ECC_fc,  p = p, n = ens_size, save_dir = NULL)
  
###########

# combine, boxplot and save

vs_dt = data.table(PCA=vs_PCA[,mean(vs)],SE=vs_SE[,mean(vs)],GS = vs_GS[,mean(vs)], ECC = vs_ECC[,mean(vs)])

pdf(paste0(plot_dir,"boxplot_vs.pdf"))
boxplot( list (vs_PCA[,vs],vs_SE[,vs],vs_GS[,vs], vs_ECC[,vs]),names = c("PCA","SE","GS","ECC"), main = "variogram scores")
dev.off()

save(vs_dt,vs_PCA,vs_SE,vs_GS,vs_ECC,paste0(save_dir,"vs.RData"))
