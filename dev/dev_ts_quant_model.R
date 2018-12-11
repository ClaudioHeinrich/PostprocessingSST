rm(list = ls())

library(PostProcessing)

load("~/temp/ts_hindcast_cleaned_quantiled.RData")

DT_clean[,ukmo_sd := apply(.SD,1,sd,na.rm=TRUE), .SDcols = paste0("ukmo_anamoly_",1:28)]

