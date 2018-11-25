rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")


X11();plot_smooth(DT_fit_ts[month ==12, .(Lon,Lat,q_500)],exclude_ocean = TRUE, rr = c(-21,21))
X11();plot_smooth(DT_fit_ts[month ==12, .(Lon,Lat,q_9500)],exclude_ocean = TRUE, rr = c(-21,21))

DT_fit_ts[,below_low:=DT_fit_ts[,rowSums(climatology > .SD)/99,.SDcols = paste0("q_",1:99 * 100)]]
plot_smooth(DT_fit_ts[month == 12,.(Lon,Lat,below_low)],exclude_ocean = TRUE)
