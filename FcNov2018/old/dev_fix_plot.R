rm(list = ls())

library(PostProcessing)
library(data.table)

setwd("~/NR/SFE/")
load("~/PostClimDataNoBackup/SFE/FcNov2018/Forecast_ts.RData")
load("~/PostClimDataNoBackup/SFE/FcNov2018/Forecast_prect.RData")

rank1 = function(X){return(rank(X)[1])}
DT_fit_ts[,r_low := apply(.SD,1,rank1)/100,.SDcols = c("q_low",paste0("q_",1:99))]
dt = DT_fit_ts[,.(Lon,Lat,r_low)]

var = colnames(dt)[3]
mn = var
rr = NULL
theta = 0.5
pixels = 256
col_scheme = "bwr"
set_white = NULL
xlab = ""
ylab = ""
save_pdf = FALSE
save_dir = "./figures/"
file_name = "diag_plot"
stretch_par = NULL
exclude_land = FALSE
exclude_ocean = FALSE


