

############# Stuff for Kristoffer #########################

#### setting up, you don't need to look at this ####

rm(list = ls())

library(data.table)

# 
# setwd("~/NR/SFE")
# options(max.print = 1e3)
# 
# library(PostProcessing)
# 
# 
# 
# Kris_dir = paste0("./Kristoffer/")
# dir.create(Kris_dir, showWarnings = FALSE)
# 
# # loading and modifying data:
# 
# data_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/NAO/")
# DT = load_combined_wide(data_dir = data_dir, output_name = paste0("dt_combine_wide_bc_var.RData"))
# 
# 
#   DT = DT[,  "Res1" := trc(Ens1 + Bias_Est) - SST_bar]
# 
# 
# mon=1
# 
# DT = DT[month == mon, .SD,.SDcols = c("Lon","Lat","year","Res1")]
# setnames(DT,c("Lon","Lat","year","X"))
# 
# # save everything:
# save.image(file = paste0(Kris_dir,"setup.RData"))

#######################

setwd(".") # set your wd to the path of the downloaded folder

load(file = "./setup.RData")

# look at the data:

DT

# We have a data variable X at spatial locations given by Lons/Lats over a range of years.
# We want to use these observations to estimate a spatial variogram gamma such that for a given year 
# E[(X_i-X_j)^2] \approx gamma(d_ij), where d_ij is the spatial distance of the locations associated with X_i and X_j


### 1.: bring your data table into STFDF format: ###

# get spatial coordinates

sp <- sp::SpatialPoints(cbind(x=DT[year == min(year), Lon],
                              y=DT[year == min(year), Lat]), 
                        proj4string = sp::CRS("+proj=longlat +datum=WGS84"))

# get time stamps for your data:  
# Convert year Y into date, namely 1st of January

time_convert = function(Y){
  as.Date(paste0(Y,"-01-01"))
}

time = as.POSIXct( time_convert(unique(DT[,year])), tz = "GMT")

# get your data

setkey(DT,year,Lat,Lon) # for creating STFDFs the data should be ordered such that the 'spatial index is moving fastest'
data = DT[,.(X)] 

# create stfdf:  
  
stfdf = spacetime::STFDF(sp, time, data)
  
  
#### Calculate the empirical semi-variogram ####
#### -------------------------------------- ####
  
  nintv = 100 # set, how many different distance bins we consider

  ## calculate the distance matrix [km], full symmetric matrix
  Dist <- sp::spDists(sp, longlat = TRUE)
  
  ## set the intervals
  up_Dist <- Dist[upper.tri(Dist, diag = FALSE)] 
  sort_up <- sort(up_Dist)
  
  bound_id = seq(1,length(up_Dist),length.out = nintv + 1)
  boundaries <- sort_up[bound_id]
  
  plot(boundaries)
  
  #compute empirical variogram
  empVgm <- gstat::variogramST( X ~ 1, stfdf, tlags=0, boundaries = boundaries,assumeRegular = TRUE, na.omit = TRUE) 
  
  plot(empVgm$gamma)
  
  #### Fit an the Exponential semi-variogram function
  #### -------------------------------------------------- 
  
  
  ## set the cutoff - in our case, usually at least the last 10% of the variogram fit look quite bad
    cutoff_ind  = ceiling(0.9*nintv)
    cutoff = empVgm$dist[cutoff_ind]
  
    ## prepare empirical variograms for fitting
  spEmpVgm <- empVgm[empVgm$dist<=cutoff,1:3]  
  class(spEmpVgm) <- c("gstatVariogram", "data.frame")
  spEmpVgm$dir.hor <- 0
  spEmpVgm$dir.ver <- 0
  
  
  ## Exponential semi-variogram function with nugget, see the 1st argument. 
  ## Fixing "psill", fit "nugget" and "range", the 2nd arg.
  ## Using fit.method = 7, the 3rd arg.
  Mod <- gstat::fit.variogram(spEmpVgm, gstat::vgm("Exp"), fit.sills = c(T,F), fit.method = 7)
  
  psill <- Mod$psill[2]
  range <- Mod$range[2]
  nugget <- Mod$psill[1]
  
  # getting covariance matrix
  Sigma <- psill*exp(-Dist/range)
  sills <- diag(Sigma) + nugget
  diag(Sigma) <- sills
  

