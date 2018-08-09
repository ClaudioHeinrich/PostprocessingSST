## score maximum, minimum and average temperature along routes

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

training_years = 1985:2000
validation_years = 2001:2010
months = 1:12

PCA_dir = paste0(save_dir,"PCA/")

forecasts = forecast_PCA_newnew(dt = DT, y = validation_years, m = months, n = 1, PCA_depth = 10, saveorgo = FALSE, cov_dir = PCA_dir)

geostat_fc = forecast_geostat(dt = DT, y = validation_years, m = months, n=1, saveorgo = FALSE, data_dir = geostat_dir)

ECC_fc = forecast_ECC (dt = DT, y = validation_years, m = months, saveorgo = FALSE)
                                 

########################################

# pick route:

route_ids = unique(DT[Lat == 50.5, grid_id])
route = "Lat = 50.5"


DT_route = DT[grid_id %in% route_ids,]
fc_route = forecasts[grid_id %in% route_ids,]
fc_gs_route = geostat_fc[grid_id %in% route_ids,]
fc_ecc_route = ECC_fc[grid_id %in% route_ids,]


# initialize data table:

validation_dt = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt,c("month","year"))

# get maximum, minimum, and mean of observation:

fc_route[,SST_max := max(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_max]
validation_dt[,"SST_max" := temp]

fc_route[,SST_min := min(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_min]
validation_dt[,"SST_min" := temp]

fc_route[,SST_mean := mean(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_mean]
validation_dt[,"SST_mean" := temp]


# get maximum, minimum, and mean of PCA forecast:

fc_route[,PCA_max := max(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_max]
validation_dt[,"PCA_max" := temp]

fc_route[,PCA_min := min(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_min]
validation_dt[,"PCA_min" := temp]

fc_route[,PCA_mean := mean(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_mean]
validation_dt[,"PCA_mean" := temp]


# get maximum, minimum, and mean of geostat forecast:

fc_gs_route[,GS_max := max(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_max]
validation_dt[,"GS_max" := temp]

fc_gs_route[,GS_min := min(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_min]
validation_dt[,"GS_min" := temp]

fc_gs_route[,GS_mean := mean(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_mean]
validation_dt[,"GS_mean" := temp]

# get maximum, minimum, and mean of ECC forecast:

fc_ecc_route[,ECC_max := max(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_max]
validation_dt[,"ECC_max" := temp]

fc_ecc_route[,ECC_min := min(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_min]
validation_dt[,"ECC_min" := temp]

fc_ecc_route[,ECC_mean := mean(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_mean]
validation_dt[,"ECC_mean" := temp]



#############################



RMSEs = data.table("route" = route,
                   "PCA_max" = validation_dt[,sqrt(mean((SST_max - PCA_max)^2))],
                   "GS_max" = validation_dt[,sqrt(mean((SST_max - GS_max)^2))],
                   "ECC_max" = validation_dt[,sqrt(mean((SST_max - ECC_max)^2))],
                   "PCA_min" = validation_dt[,sqrt(mean((SST_min - PCA_min)^2))],
                   "GS_min" = validation_dt[,sqrt(mean((SST_min - GS_min)^2))],
                   "ECC_min" = validation_dt[,sqrt(mean((SST_min - ECC_min)^2))],
                   "PCA_mean" = validation_dt[,sqrt(mean((SST_mean - PCA_mean)^2))],
                   "GS_mean" = validation_dt[,sqrt(mean((SST_mean - GS_mean)^2))],
                   "ECC_mean" = validation_dt[,sqrt(mean((SST_mean - ECC_mean)^2))])


################# next route #########################

# pick route:

route_ids = unique(DT[Lat == 60.5, grid_id])
route = "Lat = 60.5"


DT_route = DT[grid_id %in% route_ids,]
fc_route = forecasts[grid_id %in% route_ids,]
fc_gs_route = geostat_fc[grid_id %in% route_ids,]
fc_ecc_route = ECC_fc[grid_id %in% route_ids,]


# initialize data table:

validation_dt = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt,c("month","year"))

# get maximum, minimum, and mean of observation:

fc_route[,SST_max := max(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_max]
validation_dt[,"SST_max" := temp]

fc_route[,SST_min := min(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_min]
validation_dt[,"SST_min" := temp]

fc_route[,SST_mean := mean(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_mean]
validation_dt[,"SST_mean" := temp]


# get maximum, minimum, and mean of PCA forecast:

fc_route[,PCA_max := max(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_max]
validation_dt[,"PCA_max" := temp]

fc_route[,PCA_min := min(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_min]
validation_dt[,"PCA_min" := temp]

fc_route[,PCA_mean := mean(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_mean]
validation_dt[,"PCA_mean" := temp]


# get maximum, minimum, and mean of geostat forecast:

fc_gs_route[,GS_max := max(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_max]
validation_dt[,"GS_max" := temp]

fc_gs_route[,GS_min := min(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_min]
validation_dt[,"GS_min" := temp]

fc_gs_route[,GS_mean := mean(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_mean]
validation_dt[,"GS_mean" := temp]

# get maximum, minimum, and mean of ECC forecast:

fc_ecc_route[,ECC_max := max(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_max]
validation_dt[,"ECC_max" := temp]

fc_ecc_route[,ECC_min := min(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_min]
validation_dt[,"ECC_min" := temp]

fc_ecc_route[,ECC_mean := mean(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_mean]
validation_dt[,"ECC_mean" := temp]



#############################



temp = data.table("route" = route,
                   "PCA_max" = validation_dt[,sqrt(mean((SST_max - PCA_max)^2))],
                   "GS_max" = validation_dt[,sqrt(mean((SST_max - GS_max)^2))],
                   "ECC_max" = validation_dt[,sqrt(mean((SST_max - ECC_max)^2))],
                   "PCA_min" = validation_dt[,sqrt(mean((SST_min - PCA_min)^2))],
                   "GS_min" = validation_dt[,sqrt(mean((SST_min - GS_min)^2))],
                   "ECC_min" = validation_dt[,sqrt(mean((SST_min - ECC_min)^2))],
                   "PCA_mean" = validation_dt[,sqrt(mean((SST_mean - PCA_mean)^2))],
                   "GS_mean" = validation_dt[,sqrt(mean((SST_mean - GS_mean)^2))],
                   "ECC_mean" = validation_dt[,sqrt(mean((SST_mean - ECC_mean)^2))])

RMSEs = rbindlist(list(RMSEs,temp))


################# next route #########################

# pick route:

route_ids = unique(DT[Lon == -55.5, grid_id])
route = "Lon = -55.5"


DT_route = DT[grid_id %in% route_ids,]
fc_route = forecasts[grid_id %in% route_ids,]
fc_gs_route = geostat_fc[grid_id %in% route_ids,]
fc_ecc_route = ECC_fc[grid_id %in% route_ids,]


# initialize data table:

validation_dt = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt,c("month","year"))

# get maximum, minimum, and mean of observation:

fc_route[,SST_max := max(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_max]
validation_dt[,"SST_max" := temp]

fc_route[,SST_min := min(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_min]
validation_dt[,"SST_min" := temp]

fc_route[,SST_mean := mean(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_mean]
validation_dt[,"SST_mean" := temp]


# get maximum, minimum, and mean of PCA forecast:

fc_route[,PCA_max := max(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_max]
validation_dt[,"PCA_max" := temp]

fc_route[,PCA_min := min(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_min]
validation_dt[,"PCA_min" := temp]

fc_route[,PCA_mean := mean(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_mean]
validation_dt[,"PCA_mean" := temp]


# get maximum, minimum, and mean of geostat forecast:

fc_gs_route[,GS_max := max(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_max]
validation_dt[,"GS_max" := temp]

fc_gs_route[,GS_min := min(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_min]
validation_dt[,"GS_min" := temp]

fc_gs_route[,GS_mean := mean(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_mean]
validation_dt[,"GS_mean" := temp]

# get maximum, minimum, and mean of ECC forecast:

fc_ecc_route[,ECC_max := max(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_max]
validation_dt[,"ECC_max" := temp]

fc_ecc_route[,ECC_min := min(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_min]
validation_dt[,"ECC_min" := temp]

fc_ecc_route[,ECC_mean := mean(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_mean]
validation_dt[,"ECC_mean" := temp]



#############################



temp = data.table("route" = route,
                  "PCA_max" = validation_dt[,sqrt(mean((SST_max - PCA_max)^2))],
                  "GS_max" = validation_dt[,sqrt(mean((SST_max - GS_max)^2))],
                  "ECC_max" = validation_dt[,sqrt(mean((SST_max - ECC_max)^2))],
                  "PCA_min" = validation_dt[,sqrt(mean((SST_min - PCA_min)^2))],
                  "GS_min" = validation_dt[,sqrt(mean((SST_min - GS_min)^2))],
                  "ECC_min" = validation_dt[,sqrt(mean((SST_min - ECC_min)^2))],
                  "PCA_mean" = validation_dt[,sqrt(mean((SST_mean - PCA_mean)^2))],
                  "GS_mean" = validation_dt[,sqrt(mean((SST_mean - GS_mean)^2))],
                  "ECC_mean" = validation_dt[,sqrt(mean((SST_mean - ECC_mean)^2))])

RMSEs = rbindlist(list(RMSEs,temp))

################# next route #########################

# pick route:

route_ids = unique(DT[Lon == -45.5, grid_id])
route = "Lon = -45.5"


DT_route = DT[grid_id %in% route_ids,]
fc_route = forecasts[grid_id %in% route_ids,]
fc_gs_route = geostat_fc[grid_id %in% route_ids,]
fc_ecc_route = ECC_fc[grid_id %in% route_ids,]


# initialize data table:

validation_dt = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt,c("month","year"))

# get maximum, minimum, and mean of observation:

fc_route[,SST_max := max(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_max]
validation_dt[,"SST_max" := temp]

fc_route[,SST_min := min(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_min]
validation_dt[,"SST_min" := temp]

fc_route[,SST_mean := mean(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_mean]
validation_dt[,"SST_mean" := temp]


# get maximum, minimum, and mean of PCA forecast:

fc_route[,PCA_max := max(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_max]
validation_dt[,"PCA_max" := temp]

fc_route[,PCA_min := min(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_min]
validation_dt[,"PCA_min" := temp]

fc_route[,PCA_mean := mean(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_mean]
validation_dt[,"PCA_mean" := temp]


# get maximum, minimum, and mean of geostat forecast:

fc_gs_route[,GS_max := max(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_max]
validation_dt[,"GS_max" := temp]

fc_gs_route[,GS_min := min(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_min]
validation_dt[,"GS_min" := temp]

fc_gs_route[,GS_mean := mean(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_mean]
validation_dt[,"GS_mean" := temp]

# get maximum, minimum, and mean of ECC forecast:

fc_ecc_route[,ECC_max := max(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_max]
validation_dt[,"ECC_max" := temp]

fc_ecc_route[,ECC_min := min(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_min]
validation_dt[,"ECC_min" := temp]

fc_ecc_route[,ECC_mean := mean(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_mean]
validation_dt[,"ECC_mean" := temp]



#############################



temp = data.table("route" = route,
                  "PCA_max" = validation_dt[,sqrt(mean((SST_max - PCA_max)^2))],
                  "GS_max" = validation_dt[,sqrt(mean((SST_max - GS_max)^2))],
                  "ECC_max" = validation_dt[,sqrt(mean((SST_max - ECC_max)^2))],
                  "PCA_min" = validation_dt[,sqrt(mean((SST_min - PCA_min)^2))],
                  "GS_min" = validation_dt[,sqrt(mean((SST_min - GS_min)^2))],
                  "ECC_min" = validation_dt[,sqrt(mean((SST_min - ECC_min)^2))],
                  "PCA_mean" = validation_dt[,sqrt(mean((SST_mean - PCA_mean)^2))],
                  "GS_mean" = validation_dt[,sqrt(mean((SST_mean - GS_mean)^2))],
                  "ECC_mean" = validation_dt[,sqrt(mean((SST_mean - ECC_mean)^2))])

RMSEs = rbindlist(list(RMSEs,temp))


################# next route #########################

# pick route:

route_ids = unique(DT[Lon == -35.5, grid_id])
route = "Lon = -35.5"


DT_route = DT[grid_id %in% route_ids,]
fc_route = forecasts[grid_id %in% route_ids,]
fc_gs_route = geostat_fc[grid_id %in% route_ids,]
fc_ecc_route = ECC_fc[grid_id %in% route_ids,]


# initialize data table:

validation_dt = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt,c("month","year"))

# get maximum, minimum, and mean of observation:

fc_route[,SST_max := max(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_max]
validation_dt[,"SST_max" := temp]

fc_route[,SST_min := min(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_min]
validation_dt[,"SST_min" := temp]

fc_route[,SST_mean := mean(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,SST_mean]
validation_dt[,"SST_mean" := temp]


# get maximum, minimum, and mean of PCA forecast:

fc_route[,PCA_max := max(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_max]
validation_dt[,"PCA_max" := temp]

fc_route[,PCA_min := min(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_min]
validation_dt[,"PCA_min" := temp]

fc_route[,PCA_mean := mean(fc1PC10,na.rm = TRUE),by = .(month,year)]
temp = fc_route[grid_id == min(grid_id) ,PCA_mean]
validation_dt[,"PCA_mean" := temp]


# get maximum, minimum, and mean of geostat forecast:

fc_gs_route[,GS_max := max(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_max]
validation_dt[,"GS_max" := temp]

fc_gs_route[,GS_min := min(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_min]
validation_dt[,"GS_min" := temp]

fc_gs_route[,GS_mean := mean(fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_gs_route[grid_id == min(grid_id) ,GS_mean]
validation_dt[,"GS_mean" := temp]

# get maximum, minimum, and mean of ECC forecast:

fc_ecc_route[,ECC_max := max(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_max]
validation_dt[,"ECC_max" := temp]

fc_ecc_route[,ECC_min := min(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_min]
validation_dt[,"ECC_min" := temp]

fc_ecc_route[,ECC_mean := mean(ecc_fc1,na.rm = TRUE),by = .(month,year)]
temp = fc_ecc_route[grid_id == min(grid_id) ,ECC_mean]
validation_dt[,"ECC_mean" := temp]



#############################



temp = data.table("route" = route,
                  "PCA_max" = validation_dt[,sqrt(mean((SST_max - PCA_max)^2))],
                  "GS_max" = validation_dt[,sqrt(mean((SST_max - GS_max)^2))],
                  "ECC_max" = validation_dt[,sqrt(mean((SST_max - ECC_max)^2))],
                  "PCA_min" = validation_dt[,sqrt(mean((SST_min - PCA_min)^2))],
                  "GS_min" = validation_dt[,sqrt(mean((SST_min - GS_min)^2))],
                  "ECC_min" = validation_dt[,sqrt(mean((SST_min - ECC_min)^2))],
                  "PCA_mean" = validation_dt[,sqrt(mean((SST_mean - PCA_mean)^2))],
                  "GS_mean" = validation_dt[,sqrt(mean((SST_mean - GS_mean)^2))],
                  "ECC_mean" = validation_dt[,sqrt(mean((SST_mean - ECC_mean)^2))])

RMSEs = rbindlist(list(RMSEs,temp))


##################################
############# plots ##############
##################################

# max
pdf(paste0(plot_dir,"RMSE_max_along_route.pdf"))

rr = range(RMSEs[,.SD,.SDcols = c("PCA_max","GS_max","ECC_max")])

plot(RMSEs[,PCA_max],type = "b", col = "blue",
     xlab = "route", ylab = "RMSE", main = "RMSE of max. temp. on a route",
     ylim = rr)
lines(RMSEs[,ECC_max],type = "b", col = "darkgreen")
lines(RMSEs[,GS_max],type = "b", col = "darkred")

legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
dev.off()

# min
pdf(paste0(plot_dir,"RMSE_min_along_route.pdf"))

rr = range(RMSEs[,.SD,.SDcols = c("PCA_min","GS_min","ECC_min")])

plot(RMSEs[,PCA_min],type = "b", col = "blue",
     xlab = "route", ylab = "RMSE", main = "RMSE of min. temp. on a route",
     ylim = rr)
lines(RMSEs[,ECC_min],type = "b", col = "darkgreen")
lines(RMSEs[,GS_min],type = "b", col = "darkred")

legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
dev.off()


# mean
pdf(paste0(plot_dir,"RMSE_mean_along_route.pdf"))

rr = range(RMSEs[,.SD,.SDcols = c("PCA_mean","GS_mean","ECC_mean")])

plot(RMSEs[,PCA_mean],type = "b", col = "blue",
     xlab = "route", ylab = "RMSE", main = "RMSE of mean temp. on a route",
     ylim = rr)
lines(RMSEs[,ECC_mean],type = "b", col = "darkgreen")
lines(RMSEs[,GS_mean],type = "b", col = "darkred")

legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
dev.off()

