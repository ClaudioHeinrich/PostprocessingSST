rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")

load("./Fc_201812/Obs_2tm.RData")

DT_ukmo = readr::read_rds("./Fc_201812/ukmo_ts_hindcast_12.rds")
DT_ukmo[,Lon := lon]
DT_ukmo[,Lat := lat]
DT_ukmo[, c("lon","lat") := NULL]
setkey(DT_ukmo, "Lon","Lat", "year","month")
DT_ukmo = merge(DT_ukmo, DT_obs_t2m, by = c("year","month","Lon","Lat"),all.x = TRUE)
A = DT_ukmo[Lat == 58.5 & Lon == 7.5 & month == 1]

A[,ts_median := apply(.SD,1,"median"), .SDcols = paste0("ts",1:28)]
B = A[,apply(.SD,1,"sort"),.SDcols = paste0("ts",1:28)]

plot(B[1,], ylim = range(B), type = "l")
points(B[1,],pch=20)
for(j in 2:dim(B)[1])
{
    lines(B[j,])
    points(B[j,], pch = 20)
}
lines(A[,obs_t2m], col="red")
lines(A[,ts_bar], col="blue", lwd = 2)
abline(h = A[,mean(ts_bar)],lwd = 2, col="blue", lty = 2)
abline(h = median(),lwd = 2, col="blue", lty = 3)
