## score maximum, minimum and average temperature along routes or collections of locations

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_3" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))

# specify route by getting the corresponding grid ids
# This is done by fixing two points p1 and p2 and connecting them as the crow flies and considering the grid ids lying along this route.
# In particular we assume no land to be between p1 and p2

p1 = data.table(Lon = 5.32, Lat = 60.4, Loc = "Bergen") 
p2 = data.table(Lon = -21.83, Lat = 64.13, Loc = "Reykyavik") 

route_name = "Bergen to Reykyavik"

# get grid ids_along this route: First use gcIntermediate to find n coordinates on the route from p1 to p2 as the bird flies:

n = 500
route = geosphere::gcIntermediate(p1[,.(Lon,Lat)],p2[,.(Lon,Lat)],n = n)

# now, find the grid_ids in DT closest to the coordinates on the route:

grid_id_dt = unique(DT[,.(Lon,Lat,grid_id)])

point_match = NULL
for(j in 1:dim(route)[1])
{
  a = geosphere::distHaversine(as.vector(route[j,]),as.matrix(grid_id_dt[,.(Lon,Lat)]))
  point_match[j] = which.min(a)
}

dt_route = unique(grid_id_dt[point_match,])
route_ids = dt_route[,grid_id]

################

# get forecasts

load(paste0(PCA_dir,"fc.RData"))
load(paste0(SE_dir,"fc.RData"))
load(paste0(GS_dir,"fc.RData"))
load(paste0(ECC_dir,"fc.RData"))

########################################

# pick route:

DT_route = DT[grid_id %in% route_ids & year %in% validation_years,]
PCA_fc_route = PCA_fc[grid_id %in% route_ids,]
SE_fc_route = SE_fc[grid_id %in% route_ids,]
GS_fc_route = GS_fc[grid_id %in% route_ids,]
ECC_fc_route = ECC_fc[grid_id %in% route_ids,]

###### score maximum along route ########

# initialize data table:

validation_dt_max = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt_max,c("month","year"))

# get maximum of observation

DT_route[,SST_max := max(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_max]
validation_dt_max[,"SST_max" := temp]

# get maximum of PCA forecast:

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_max = PCA_fc_route[,.SD,.SDcols = paste0('PCA_max_',1:fc_ens_size)]

PCA_fc_route[, paste0('PCA_max_',1:fc_ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,PCA_max)

# get maximum of SE forecast:

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_max = SE_fc_route[,.SD,.SDcols = paste0('SE_max_',1:fc_ens_size)]

SE_fc_route[, paste0('SE_max_',1:fc_ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,SE_max)

# get maximum of GS forecast:

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_max = GS_fc_route[,.SD,.SDcols = paste0('GS_max_',1:fc_ens_size)]

GS_fc_route[, paste0('GS_max_',1:fc_ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,GS_max)

# get maximum of ECC forecast:

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_max_',i) := max(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_max = ECC_fc_route[,.SD,.SDcols = paste0('ECC_max_',1:ens_size)]

ECC_fc_route[, paste0('ECC_max_',1:ens_size) := NULL]

validation_dt_max = data.table(validation_dt_max,ECC_max)

###########################



RMSEs_max = data.table("route" = route_name,
                       "RMSE_PCA" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('PCA_max_',1:fc_ens_size)]),
                       "RMSE_SE" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('SE_max_',1:fc_ens_size)]),
                       "RMSE_GS" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('GS_max_',1:fc_ens_size)]),
                       "RMSE_ECC" = sqrt(validation_dt_max[,mean((SST_max - .SD)^2),.SDcols = paste0('ECC_max_',1:ens_size)])
                       )


############# RMSEs for minimum temperature along route #####################


# initialize data table:

validation_dt_min = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt_min,c("month","year"))

# get minimum of observation

DT_route[,SST_min := min(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_min]
validation_dt_min[,"SST_min" := temp]

# get minimum of PCA forecast:

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_min = PCA_fc_route[,.SD,.SDcols = paste0('PCA_min_',1:fc_ens_size)]

PCA_fc_route[, paste0('PCA_min_',1:fc_ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,PCA_min)

# get minimum of SE forecast:

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_min = SE_fc_route[,.SD,.SDcols = paste0('SE_min_',1:fc_ens_size)]

SE_fc_route[, paste0('SE_min_',1:fc_ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,SE_min)

# get minimum of GS forecast:

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_min = GS_fc_route[,.SD,.SDcols = paste0('GS_min_',1:fc_ens_size)]

GS_fc_route[, paste0('GS_min_',1:fc_ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,GS_min)

# get minimum of ECC forecast:

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_min_',i) := min(.SD,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_min = ECC_fc_route[,.SD,.SDcols = paste0('ECC_min_',1:ens_size)]

ECC_fc_route[, paste0('ECC_min_',1:ens_size) := NULL]

validation_dt_min = data.table(validation_dt_min,ECC_min)

###########################



RMSEs_min = data.table("route" = route_name,
                       "RMSE_PCA" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('PCA_min_',1:fc_ens_size)]),
                       "RMSE_SE" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('SE_min_',1:fc_ens_size)]),
                       "RMSE_GS" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('GS_min_',1:fc_ens_size)]),
                       "RMSE_ECC" = sqrt(validation_dt_min[,mean((SST_min - .SD)^2),.SDcols = paste0('ECC_min_',1:ens_size)])
)


############# scoring the mean along the route #########################


# initialize data table:

validation_dt_mean = as.data.table(expand.grid(months,validation_years))
setnames(validation_dt_mean,c("month","year"))

# get mean of observation

DT_route[,SST_mean := mean(SST_bar,na.rm = TRUE),by = .(month,year)]
temp = DT_route[grid_id == min(grid_id) ,SST_mean]
validation_dt_mean[,"SST_mean" := temp]

# get meanimum of PCA forecast:

for(i in 1:fc_ens_size)
{
  PCA_fc_route[,paste0('PCA_mean_',i) := lapply(.SD, mean, na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

PCA_mean = PCA_fc_route[,.SD,.SDcols = paste0('PCA_mean_',1:fc_ens_size)]

PCA_fc_route[, paste0('PCA_mean_',1:fc_ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,PCA_mean)

# get mean of SE forecast:

for(i in 1:fc_ens_size)
{
  SE_fc_route[,paste0('SE_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

SE_mean = SE_fc_route[,.SD,.SDcols = paste0('SE_mean_',1:fc_ens_size)]

SE_fc_route[, paste0('SE_mean_',1:fc_ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,SE_mean)

# get meanimum of GS forecast:

for(i in 1:fc_ens_size)
{
  GS_fc_route[,paste0('GS_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

GS_mean = GS_fc_route[,.SD,.SDcols = paste0('GS_mean_',1:fc_ens_size)]

GS_fc_route[, paste0('GS_mean_',1:fc_ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,GS_mean)

# get meanimum of ECC forecast:

for(i in 1:ens_size)
{
  ECC_fc_route[,paste0('ECC_mean_',i) := lapply(.SD,mean,na.rm = TRUE),by = .(month,year),.SDcols = paste0('fc',i)]  
}

ECC_mean = ECC_fc_route[,.SD,.SDcols = paste0('ECC_mean_',1:ens_size)]

ECC_fc_route[, paste0('ECC_mean_',1:ens_size) := NULL]

validation_dt_mean = data.table(validation_dt_mean,ECC_mean)

###########################



RMSEs_mean = data.table("route" = route_name,
                       "RMSE_PCA" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('PCA_mean_',1:fc_ens_size)]),
                       "RMSE_SE" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('SE_mean_',1:fc_ens_size)]),
                       "RMSE_GS" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('GS_mean_',1:fc_ens_size)]),
                       "RMSE_ECC" = sqrt(validation_dt_mean[,mean((SST_mean - .SD)^2),.SDcols = paste0('ECC_mean_',1:ens_size)])
)


# pugging together, bringing in shape
RMSEs = rbindlist(list(RMSEs_max[,fun := "max"],RMSEs_min[,fun := "min"],RMSEs_mean[,fun:= 'mean']))

RMSEs = RMSEs[,route := NULL]

setcolorder(RMSEs,'fun')

temp <- round(RMSEs[,-1],3)

RMSEs = data.table(RMSEs[,.(fun)],temp)


#############################



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
