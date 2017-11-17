rm(list = ls())

##----- Setup -------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
##-------------------

in.dir = "~/PostClimDataNoBackup/SFE/GCFS1"
out.dir = "~PostClimDataNoBackup/SFE/Derived/"
N_ens = 15
verbose = TRUE

make_GCFS1_wide()


