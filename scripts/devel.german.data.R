rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)


ncobs = nc_open("~/PostClimDataNoBackup/SFE/GCFS1/Mem01_GCFS1_sst_mm_1981-2015_smon11.nc")
