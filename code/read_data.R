rm(list=ls())

library(ncdf4)
options(max.print = 1e3)

filedir <- "/nr/samba/PostClimDataNoBackup/SFE/NorCPM_Ocean/"
filename <- "ana_19800115_me_20080715_mem01.micom.hm.2009-06.nc"
ncname <- paste(filedir, filename, sep="")

ncin <- nc_open(ncname)
print(ncin)

sst <- ncvar_get(ncin, "sst")
pred.lon <- ncvar_get(ncin, "time")

obsdir <- "/nr/project/stat/SFE/Data/HadiSST2/"
obsname <- "SST_ens_2009_05.nc"
ncnameobs <- paste(obsdir, obsname, sep="")
ncobs <- nc_open(ncnameobs)
print(ncobs)
obs.sst <- ncvar_get(ncobs, "sst")
obs.lon <- ncvar_get(ncobs, "longitude")
obs.lat <- ncvar_get(ncobs, "latitude")

gridname <- "/nr/project/stat/SFE/Data/grid.nc"
ncgrid <- nc_open(gridname)
print(ncgrid)
grid.plon <- ncvar_get(ncgrid, "plon")
grid.plat <- ncvar_get(ncgrid, "plat")
