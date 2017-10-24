rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)

load("./Data/PostClim/SFE/Derived/dt_map.RData")
names(dt_map)= c("Lon_Obs","Lat_Obs","Lon_Ens","Lat_Ens") ## Get rid of this eventually.
write.csv(dt_map, file = "~/grid_mapping.csv", row.names = FALSE)

dt_wide = load_combined_wide()

write.csv(dt_wide,row.names = FALSE, file = "~/wide_data.csv")
