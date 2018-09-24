########################################################################################################

#############  side script 4.3 - score various functionals that depend on areas of SST  ################

########################################################################################################

# We score max, min, and mean SST over a predefined areal, mostly a route or a square
# 
# Files generated: file name depends on areal, see below
#   
# Requires previous run of 04.master.multivar.pp.R with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_3" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

####################################
####### area specification #########
####################################

# specify the area to score over. Version 1: specify a route:

# Fix two points p1 and p2. These are then connected by a route (as the crow flies, not caring about land in between).
# A grid point in dt is then considered for scoring if it is not on land and is nearest neighbour to a point on the route. 

p1 = data.table(Lon = 5.32, Lat = 60.4, Loc = "Bergen") 
p2 = data.table(Lon = -21.83, Lat = 64.13, Loc = "Reykyavik") 

route_name = "Bergen to Reykyavik"
file_name = "scores_Bergen_to_Reykyavik"

# get grid ids_along this route: fix n and use gcIntermediate to find n coordinates on the route from p1 to p2 as the crow flies:

n = 500
route = geosphere::gcIntermediate(p1[,.(Lon,Lat)],p2[,.(Lon,Lat)],n = n)

# find the grid_ids in DT closest to the coordinates on the route:

grid_id_dt = unique(DT[,.(Lon,Lat,grid_id)])

point_match = NULL
for(j in 1:dim(route)[1])
{
  a = geosphere::distHaversine(as.vector(route[j,]),as.matrix(grid_id_dt[,.(Lon,Lat)]))
  point_match[j] = which.min(a)
}

dt_route = unique(grid_id_dt[point_match,])
route_ids = dt_route[,grid_id]

###################################

# specify the area to score over. Version 2: specify a window:

# option two: Lon/Lat window

Lon_min = -25
Lon_max = -15
Lat_min = 60
Lat_max = 64

grid_id_dt = unique(DT[Lon >= Lon_min & Lon <= Lon_max,][Lat >= Lat_min & Lat <= Lat_max,.(Lon,Lat,grid_id)])
route_ids = grid_id_dt[,grid_id]

route_name = "South of Island"
file_name = "scores_south_of_Island"


######################################

# get forecasts

load(paste0(PCA_dir,"fc.RData"))
load(paste0(SE_dir,"fc.RData"))
load(paste0(GS_dir,"fc.RData"))
load(paste0(ECC_dir,"fc.RData"))

########################################

# reduce data to route

DT_route = DT[grid_id %in% route_ids & year %in% validation_years,]
PCA_fc_route = PCA_fc[grid_id %in% route_ids,]
SE_fc_route = SE_fc[grid_id %in% route_ids,]
GS_fc_route = GS_fc[grid_id %in% route_ids,]
ECC_fc_route = ECC_fc[grid_id %in% route_ids,]



######################################################
###### Get RMSEs for the  maximum along route ########
######################################################

# initialize data table:

validation_dt_max = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt_max,c("month","year"))

# get maximum of observation

DT_route[,SST_max := max(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_max]
validation_dt_max[,"SST_max" := temp]

# get maximum of point forecast

DT_route[,SST_hat_max := max(SST_hat,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_hat_max]
validation_dt_max[,"SST_hat_max" := temp]

# get maximum of PCA forecast:

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_max = unique(PCA_fc_route[,.SD,.SDcols = paste0('PCA_max_',1:fc_ens_size)])

PCA_fc_route[, paste0('PCA_max_',1:fc_ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,PCA_max)

# get maximum of SE forecast:

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_max = unique(SE_fc_route[,.SD,.SDcols = paste0('SE_max_',1:fc_ens_size)])

SE_fc_route[, paste0('SE_max_',1:fc_ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,SE_max)

# get maximum of GS forecast:

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_max = unique(GS_fc_route[,.SD,.SDcols = paste0('GS_max_',1:fc_ens_size)])

GS_fc_route[, paste0('GS_max_',1:fc_ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,GS_max)

# get maximum of ECC forecast:

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_max = unique(ECC_fc_route[,.SD,.SDcols = paste0('ECC_max_',1:ens_size)])

ECC_fc_route[, paste0('ECC_max_',1:ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,ECC_max)

############# putting together ##############

RMSEs_max = data.table("route" = route_name,
                       "RMSE_PCA" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('PCA_max_',1:fc_ens_size)]),
                       "RMSE_SE" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('SE_max_',1:fc_ens_size)]),
                       "RMSE_GS" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('GS_max_',1:fc_ens_size)]),
                       "RMSE_ECC" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('ECC_max_',1:ens_size)]),
                       "RMSE_PFC" = sqrt(validation_dt_max[,mean((SST_max - SST_hat_max)^2)])
                       )



######################################################
###### Get RMSEs for the  minimum along route ########
######################################################


# initialize data table:

validation_dt_min = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt_min,c("month","year"))

# get minimum of observation

DT_route[,SST_min := min(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_min]
validation_dt_min[,"SST_min" := temp]

# get minimum of point forecast

DT_route[,SST_hat_min := min(SST_hat,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_hat_min]
validation_dt_min[,"SST_hat_min" := temp]

# get minimum of PCA forecast:

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_min = unique(PCA_fc_route[,.SD,.SDcols = paste0('PCA_min_',1:fc_ens_size)])

PCA_fc_route[, paste0('PCA_min_',1:fc_ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,PCA_min)

# get minimum of SE forecast:

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_min = unique(SE_fc_route[,.SD,.SDcols = paste0('SE_min_',1:fc_ens_size)])

SE_fc_route[, paste0('SE_min_',1:fc_ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,SE_min)

# get minimum of GS forecast:

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_min = unique(GS_fc_route[,.SD,.SDcols = paste0('GS_min_',1:fc_ens_size)])

GS_fc_route[, paste0('GS_min_',1:fc_ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,GS_min)

# get minimum of ECC forecast:

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_min = unique(ECC_fc_route[,.SD,.SDcols = paste0('ECC_min_',1:ens_size)])

ECC_fc_route[, paste0('ECC_min_',1:ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,ECC_min)


########## putting together ###############

RMSEs_min = data.table("route" = route_name,
                       "RMSE_PCA" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('PCA_min_',1:fc_ens_size)]),
                       "RMSE_SE" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('SE_min_',1:fc_ens_size)]),
                       "RMSE_GS" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('GS_min_',1:fc_ens_size)]),
                       "RMSE_ECC" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('ECC_min_',1:ens_size)]),
                       "RMSE_PFC" = sqrt(validation_dt_min[,mean(SST_min - SST_hat_min)^2])
)




###################################################
###### Get RMSEs for the  mean along route ########
###################################################


# initialize data table:

validation_dt_mean = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt_mean,c("month","year"))

# get mean of observation

DT_route[,SST_mean := mean(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_mean]
validation_dt_mean[,"SST_mean" := temp]

# get mean of point forecast

DT_route[,SST_hat_mean := mean(SST_hat,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_hat_mean]
validation_dt_mean[,"SST_hat_mean" := temp]

# get mean of PCA forecast:

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_mean_',i) := lapply(.SD, mean, na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_mean = unique(PCA_fc_route[,.SD,.SDcols = paste0('PCA_mean_',1:fc_ens_size)])

PCA_fc_route[, paste0('PCA_mean_',1:fc_ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,PCA_mean)

# get mean of SE forecast:

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_mean = unique(SE_fc_route[,.SD,.SDcols = paste0('SE_mean_',1:fc_ens_size)])

SE_fc_route[, paste0('SE_mean_',1:fc_ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,SE_mean)

# get mean of GS forecast:

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_mean = unique(GS_fc_route[,.SD,.SDcols = paste0('GS_mean_',1:fc_ens_size)])

GS_fc_route[, paste0('GS_mean_',1:fc_ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,GS_mean)

# get mean of ECC forecast:

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_mean = unique(ECC_fc_route[,.SD,.SDcols = paste0('ECC_mean_',1:ens_size)])

ECC_fc_route[, paste0('ECC_mean_',1:ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,ECC_mean)

############ putting together ###############

RMSEs_mean = data.table("route" = route_name,
                       "RMSE_PCA" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('PCA_mean_',1:fc_ens_size)]),
                       "RMSE_SE" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('SE_mean_',1:fc_ens_size)]),
                       "RMSE_GS" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('GS_mean_',1:fc_ens_size)]),
                       "RMSE_ECC" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('ECC_mean_',1:ens_size)]),
                       "RMSE_PFC" = sqrt(validation_dt_mean[,mean((SST_mean - SST_hat_mean)^2)])
)


#############################################

### combining all RMSEs ###

RMSEs = rbindlist(list(RMSEs_max[,fun := "max"],RMSEs_min[,fun := "min"],RMSEs_mean[,fun:= 'mean']))

RMSEs = RMSEs[,route := NULL]

setcolorder(RMSEs,'fun')

temp <- round(RMSEs[,-1],3)

RMSEs = data.table(RMSEs[,.(fun)],temp)

##############################################


#########################################
######### Do the same with CRPS #########
#########################################


######################################################
###### Get CRPS for the  maximum along route #########
######################################################

### PFC ###

# for point forecasts, CRPS is just the expected absolute distance:

PFC_CRPS = validation_dt_max[,mean(abs(SST_max - SST_hat_max))]

### PCA ###

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_max = unique(PCA_fc_route[,.SD,.SDcols = paste0('PCA_max_',1:fc_ens_size)])

PCA_mean = rowMeans(PCA_max)

PCA_sd = sqrt(rowMeans((PCA_max - PCA_mean)^2))

PCA_CRPS = crps.na.rm(validation_dt_max[,SST_max], mean = PCA_mean, sd = PCA_sd)


### SE ###

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_max = unique(SE_fc_route[,.SD,.SDcols = paste0('SE_max_',1:fc_ens_size)])

SE_mean = rowMeans(SE_max)

SE_sd = sqrt(rowMeans((SE_max - SE_mean)^2))

SE_CRPS = crps.na.rm(validation_dt_max[,SST_max], mean = SE_mean, sd = SE_sd)

### GS ###

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_max = unique(GS_fc_route[,.SD,.SDcols = paste0('GS_max_',1:fc_ens_size)])

GS_mean = rowMeans(GS_max)

GS_sd = sqrt(rowMeans((GS_max - GS_mean)^2))

GS_CRPS = crps.na.rm(validation_dt_max[,SST_max], mean = GS_mean, sd = GS_sd)

### ECC ###

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_max = unique(ECC_fc_route[,.SD,.SDcols = paste0('ECC_max_',1:ens_size)])

ECC_mean = rowMeans(ECC_max)

ECC_sd = sqrt(rowMeans((ECC_max - ECC_mean)^2))

ECC_CRPS = crps.na.rm(validation_dt_max[,SST_max], mean = ECC_mean, sd = ECC_sd)


######### putting together ###############

CRPS_max = data.table("route" = route_name,
                      "CRPS_PCA" = mean(PCA_CRPS),
                      "CRPS_SE" = mean(SE_CRPS),
                      "CRPS_GS" = mean(GS_CRPS),
                      "CRPS_ECC" = mean(ECC_CRPS),
                      "CRPS_PFC" = PFC_CRPS)




######################################################
###### Get CRPS for the  minimum along route #########
######################################################

### PFC ###

# for point forecasts, CRPS is just the expected absolute distance:

PFC_CRPS = validation_dt_min[,mean(abs(SST_min - SST_hat_min))]

### PCA ###

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_min = unique(PCA_fc_route[,.SD,.SDcols = paste0('PCA_min_',1:fc_ens_size)])

PCA_mean = rowMeans(PCA_min)

PCA_sd = sqrt(rowMeans((PCA_min - PCA_mean)^2))

PCA_CRPS = crps.na.rm(validation_dt_min[,SST_min], mean = PCA_mean, sd = PCA_sd)


### SE ###

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_min = unique(SE_fc_route[,.SD,.SDcols = paste0('SE_min_',1:fc_ens_size)])

SE_mean = rowMeans(SE_min)

SE_sd = sqrt(rowMeans((SE_min - SE_mean)^2))

SE_CRPS = crps.na.rm(validation_dt_min[,SST_min], mean = SE_mean, sd = SE_sd)

### GS ###

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_min = unique(GS_fc_route[,.SD,.SDcols = paste0('GS_min_',1:fc_ens_size)])

GS_mean = rowMeans(GS_min)

GS_sd = sqrt(rowMeans((GS_min - GS_mean)^2))

GS_CRPS = crps.na.rm(validation_dt_min[,SST_min], mean = GS_mean, sd = GS_sd)

### ECC ###

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_min = unique(ECC_fc_route[,.SD,.SDcols = paste0('ECC_min_',1:ens_size)])

ECC_mean = rowMeans(ECC_min)

ECC_sd = sqrt(rowMeans((ECC_min - ECC_mean)^2))

ECC_CRPS = crps.na.rm(validation_dt_min[,SST_min], mean = ECC_mean, sd = ECC_sd)


#  mean CRPS

CRPS_min = data.table("route" = route_name,
                      "CRPS_PCA" = mean(PCA_CRPS),
                      "CRPS_SE" = mean(SE_CRPS),
                      "CRPS_GS" = mean(GS_CRPS),
                      "CRPS_ECC" = mean(ECC_CRPS),
                      "CRPS_PFC" = PFC_CRPS)



######################################################
###### Get CRPS for the  mean along route #########
######################################################

### PFC ###

# for point forecasts, CRPS is just the expected absolute distance:

PFC_CRPS = validation_dt_mean[,mean(abs(SST_mean - SST_hat_mean))]


### PCA ###

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_mean = unique(PCA_fc_route[,.SD,.SDcols = paste0('PCA_mean_',1:fc_ens_size)])

PCA_mn = rowMeans(PCA_mean)

PCA_sd = sqrt(rowMeans((PCA_mean - PCA_mn)^2))

PCA_CRPS = crps.na.rm(validation_dt_mean[,SST_mean], mean = PCA_mn, sd = PCA_sd)


### SE ###

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_mean = unique(SE_fc_route[,.SD,.SDcols = paste0('SE_mean_',1:fc_ens_size)])

SE_mn = rowMeans(SE_mean)

SE_sd = sqrt(rowMeans((SE_mean - SE_mn)^2))

SE_CRPS = crps.na.rm(validation_dt_mean[,SST_mean], mean = SE_mn, sd = SE_sd)

### GS ###

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_mean = unique(GS_fc_route[,.SD,.SDcols = paste0('GS_mean_',1:fc_ens_size)])

GS_mn = rowMeans(GS_mean)

GS_sd = sqrt(rowMeans((GS_mean - GS_mn)^2))

GS_CRPS = crps.na.rm(validation_dt_mean[,SST_mean], mean = GS_mn, sd = GS_sd)

### ECC ###

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_mean = unique(ECC_fc_route[,.SD,.SDcols = paste0('ECC_mean_',1:ens_size)])

ECC_mn = rowMeans(ECC_mean)

ECC_sd = sqrt(rowMeans((ECC_mean - ECC_mn)^2))

ECC_CRPS = crps.na.rm(validation_dt_mean[,SST_mean], mean = ECC_mn, sd = ECC_sd)


########### putting together ##############

CRPS_mean = data.table("route" = route_name,
                      "CRPS_PCA" = mean(PCA_CRPS),
                      "CRPS_SE" = mean(SE_CRPS),
                      "CRPS_GS" = mean(GS_CRPS),
                      "CRPS_ECC" = mean(ECC_CRPS),
                      "CRPS_PFC" = PFC_CRPS)

###########################################

########## combining CRPS #################

CRPS = rbindlist(list(CRPS_max[,fun := "max"],CRPS_min[,fun := "min"],CRPS_mean[,fun:= 'mean']))

CRPS = CRPS[,route := NULL]

setcolorder(CRPS,'fun')

temp <- round(CRPS[,-1],3)

CRPS = data.table(CRPS[,.(fun)],temp)

CRPS

### save everything ####

save(RMSEs,CRPS,file = paste0(save_dir,file_name,".RData"))





# ##################################
# ############# plots ##############
# ##################################
# 
# # max
# pdf(paste0(plot_dir,"RMSE_max_along_route.pdf"))
# 
# rr = range(RMSEs[,.SD,.SDcols = c("PCA_max","GS_max","ECC_max")])
# 
# plot(RMSEs[,PCA_max],type = "b", col = "blue",
#      xlab = "route", ylab = "RMSE", main = "RMSE of max. temp. on a route",
#      ylim = rr)
# lines(RMSEs[,ECC_max],type = "b", col = "darkgreen")
# lines(RMSEs[,GS_max],type = "b", col = "darkred")
# 
# legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
# dev.off()
# 
# # min
# pdf(paste0(plot_dir,"RMSE_min_along_route.pdf"))
# 
# rr = range(RMSEs[,.SD,.SDcols = c("PCA_min","GS_min","ECC_min")])
# 
# plot(RMSEs[,PCA_min],type = "b", col = "blue",
#      xlab = "route", ylab = "RMSE", main = "RMSE of min. temp. on a route",
#      ylim = rr)
# lines(RMSEs[,ECC_min],type = "b", col = "darkgreen")
# lines(RMSEs[,GS_min],type = "b", col = "darkred")
# 
# legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
# dev.off()
# 
# 
# # mean
# pdf(paste0(plot_dir,"RMSE_mean_along_route.pdf"))
# 
# rr = range(RMSEs[,.SD,.SDcols = c("PCA_mean","GS_mean","ECC_mean")])
# 
# plot(RMSEs[,PCA_mean],type = "b", col = "blue",
#      xlab = "route", ylab = "RMSE", main = "RMSE of mean temp. on a route",
#      ylim = rr)
# lines(RMSEs[,ECC_mean],type = "b", col = "darkgreen")
# lines(RMSEs[,GS_mean],type = "b", col = "darkred")
# 
# legend(x = "topright", legend = c("PCA","GS","ECC"),col = c("blue","darkred","darkgreen"),lty = c(1,1,1))
# dev.off()
