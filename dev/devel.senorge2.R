## Get all the Se norge data together 1991, 2015

rm(list = ls())

##----- Set up ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
##---------------------

##----- Grid Restrictions -----
in.dir = "~/PostClimDataNoBackup/seNorge/"
lon_bound = c(4.5,7.8)
lat_bound = c(57.95, 61.04)
##-----------------------------

  

make_senorge_data()






