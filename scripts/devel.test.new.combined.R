rm(list = ls())

setwd("~/NR/SFE/")
library(SeasonalForecasting)
options(max.print = 1e3)

y_start = 1985
y_stop = 2010
data.dir = "~/PostClimDataNoBackup/"
grid_mapping_loc = "./Data/PostClim/SFE/Derived/"
