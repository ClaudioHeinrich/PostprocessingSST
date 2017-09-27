rm(list = ls())

##---- Setup ----
setwd("~/NR/SFE/")
library(SeasonalForecasting)
options(max.print = 1e3)
##---------------

y_start = 1986
y_stop = 2010
vintage = "mr"
data.dir = "~/PostClimDataNoBackup/"
grid_mapping_loc = "./Data/PostClim/SFE/Derived/"
