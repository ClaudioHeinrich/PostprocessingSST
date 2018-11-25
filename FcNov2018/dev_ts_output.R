## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
library(ncdf4)

setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

ncpath <- "./FcNov2018/"
ncname <- "ts"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "ts"  # note: tmp means temperature (not temporary)

lon_all = unique(sort(DT_fit_ts[,Lon]))
lat_all = unique(sort(DT_fit_ts[,Lat]))
                 
londim <- ncdim_def("lon","degrees_east",as.double(lon_all))
latdim <- ncdim_def("lat","degrees_north",as.double(lat_all))
qdim = ncdim_def("q","quantile",as.double(1:99/100))
monthdim = ncdim_def("month","month",as.double(c(11,12,1,2,3,4)))

fillvalue <- 1e32
dlname <- "2m air temperature"
tmp_def <- ncvar_def("ts",
                     "deg_C",
                     list(londim,latdim),
                     fillvalue,
                     dlname,
                     prec="double")
