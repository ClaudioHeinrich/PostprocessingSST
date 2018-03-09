
#######################################################################
##  Estimate a monthly Exponential semi-variogram function.
#######################################################################

rm(list = ls())

library(SeasonalForecasting)


setwd("~/NR/SFE")
options(max.print = 500)


#------- set up ------

m = 7
training_years = 1985:2000


#fitting area for the variogram and area for plotting:

Lon_min = -60
Lon_max = 15
Lat_min = 30
Lat_max = 70

#loading data

DT_complete = load_combined_wide(bias = TRUE)

# --- modify data for fitting variogram with spacetime package ---

library(sp) # SpatialPoints(), spDists(), SpatialPointsDataFrame(), spplot
library(spacetime)  # "STFDF" object
library(gstat)  # variogramST(), fit.variogram()

DT = DT_complete[month == m & year %in% training_years & 
                            Lon >= Lon_min & Lon <= Lon_max & Lat >= Lat_min & Lat <= Lat_max,
                            .(Lon,Lat,year,month,YM, Res = SST_hat - SST_bar)]


sp <- SpatialPoints(cbind(x=DT[YM == min(YM), Lon],y=DT[YM == min(YM), Lat]), proj4string = CRS("+proj=longlat +datum=WGS84"))

# Convert YM into date

time_convert = function(YM){
  M = YM %% 12
  M[M == 0] = 12
  Y = (YM - M)/ 12
  M[M < 10]  = paste0(0,M[M < 10])
  as.Date(paste0(Y,"-",M,"-15"))
}

time = as.POSIXct( time_convert(unique(DT[,YM])), tz = "GMT")

data = DT[,.(Res)]

stfdf = STFDF(sp, time, data)


#### Calculate the empirical semi-variogram
#### ----------------------------------------------

## calculate the distance matrix [km], full symmetric matrix
Dist <- spDists(sp, longlat = T)

## set the intervals
up.Dist <- Dist[upper.tri(Dist, diag = F)] 
sort.up <- sort(up.Dist)

nintv <- 150
bound.id <- seq(1, length(up.Dist), length(up.Dist)/nintv)
boundaries <- sort.up[bound.id]
boundaries <- c(boundaries, max(up.Dist)) 

png("./figures/boundaries.png")
plot(boundaries, main = nintv)
dev.off()


empVgm <- variogramST(Res~1, stfdf, tlags=0, boundaries = boundaries,na.omit = TRUE) #"gstat"
show(empVgm[1:3,])


dim(empVgm)
## 140  7
head(empVgm, 3)
##         np     dist      gamma   id timelag spacelag  avgDist
## 2 16434600 2.593084 0.01181281 lag0       0 1.802588 2.593084
## 3 15930300 4.567351 0.01421445 lag0       0 4.351693 4.567351
## 4 16504200 5.896617 0.01516245 lag0       0 5.749505 5.896617

#### ----------------------------------------------



#### Fit to the Exponential semi-variogram function
#### -------------------------------------------------- 

## set the cutoff, here I use no cutoff
cutoff <- max(empVgm$dist) # 4300

## prepare empirical variograms for fitting
spEmpVgm <- empVgm[empVgm$dist<=cutoff,]  
spEmpVgm <- spEmpVgm[spEmpVgm$timelag==0,]
sSpEmpVgm <- spEmpVgm[spEmpVgm$np!=0,] 
spEmpVgm <- sSpEmpVgm[,1:3] 
class(spEmpVgm) <- c("gstatVariogram", "data.frame")
spEmpVgm$dir.hor <- 0
spEmpVgm$dir.ver <- 0


## Exponential semi-variogram function with nugget, see the 1st argument. 
## Fixing "psill", fit "nugget" and "range", the 2nd arg.
## Using fit.method = 7, the 3rd arg.
Mod <- fit.variogram(spEmpVgm, vgm("Exp"), fit.sills = c(T,F), fit.method = 7)
psill <- Mod$psill[2]
range <- Mod$range[2]
nugget <- Mod$psill[1]
c(psill, range, nugget)
## 9.650506e-02 1.011401e+02 8.779565e-03
str(Mod)

#### -------------------------------------------------



#### Plot the fitted variogram on top of the empirical variograms
#### --------------------------------------------------------------

range(spEmpVgm$gamma)
## 0.01181281 0.11760182
range(spEmpVgm$dist)
## 2.593084 98.518692

ycuts <- seq(0, range(spEmpVgm$gamma)[2], 0.07)
xcuts <- seq(0, range(spEmpVgm$dist)[2] , 150)

png("./figures/Vgm_uncensored.png")
plot(spEmpVgm$dist,spEmpVgm$gamma, ylim=c(0,range(spEmpVgm$gamma)[2]), xlim = c(0,range(spEmpVgm$dist)[2]),
     ylab = "", xlab="", axes = FALSE, 
     main = paste0("Empirical and fitted semi-variograms for month",m) )
axis(1, at = xcuts)
axis(2, at = ycuts)
lines(variogramLine(vgm(psill,"Exp",range, nugget), maxdist=max(spEmpVgm$dist)),
      col="blue", lwd=2)
dev.off()
#### ------------------------------------------------------------



#### Simulation and plotting
#### ----------------------------------------------



Sigma <- psill*exp(-Dist/range)
sills <- diag(Sigma) + nugget
diag(Sigma) <- sills

ns <- length(sp)

# plot and save 5 example residuals with exponential covariance

plotting_DT = DT_complete[month == m & year == min(training_years) & 
                                                Lon >= Lon_min & Lon <= Lon_max & Lat >= Lat_min & Lat <= Lat_max,
                                              .(Lon,Lat,year,month,YM, Res = SST_hat - SST_bar)]

for(Ex in 1:5){
  no <- mvrnorm(n=1, mu=rep(0,ns), Sigma=Sigma)
  
  plotting_DT[,noise := no]
  plotting_DT[is.na(Res), noise := NA]
  
  
  file.name =paste0( "Exp_cov_mon",m,"_ex",Ex)
  plot_diagnostic(plotting_DT[,.(Lon,Lat,noise)],
                  mn = paste0("exponential covariance res for June"),
                  save.pdf = TRUE,file.name = file.name)
  print(paste0("Example ",Ex," complete"))
}

# plot and save example residuals observed in different (oos)-years

ex_years = 2001:2006

DT_oos = DT_complete[month == m & year %in% ex_years & 
                                           Lon >= Lon_min & Lon <= Lon_max & Lat >= Lat_min & Lat <= Lat_max,
                                         .(Lon,Lat,year,month,YM, Res = SST_hat - SST_bar)]

for(Ex in ex_years){
  
  file.name =paste0( "Obs_res_mon",m,"_y",Ex)
  plot_diagnostic(DT_oos[year == Ex,.(Lon,Lat,Res)],
                  mn = paste0("Observed res for ",m," / ",Ex),
                  save.pdf = TRUE,file.name = file.name)
  print(paste0("Example ",Ex," complete"))
}

# plot and save example residuals generated by PCA

PCs = 15

setup_PCA(dt = DT_complete, m = m, y = ex_years, max_PCA_depth = 50)

for(Ex in ex_years){
  
  Noise = forecast_PCA(y = Ex, m = m, PCA_depth = PCs,saveorgo = FALSE)[Lon >= Lon_min & 
                                                                          Lon <= Lon_max & 
                                                                          Lat >= Lat_min & 
                                                                          Lat <= Lat_max,.(Lon,Lat,noise)]
  
  
  file.name =paste0( "PCA_res_mon",m,"_y",Ex)
  plot_diagnostic(Noise,
                  mn = paste0("PCA res for ",m," / ",Ex),
                  save.pdf = TRUE,file.name = file.name)
  print(paste0("Example ",Ex," complete"))
}

 
  
#### Plot simulated residuals
#### ------------------------------------------------------- 
library(RColorBrewer)
my.palette <- rev(brewer.pal(n = 11, name = "RdBu"))
cont.palette <- colorRampPalette(my.palette)(1000000)

range(wStar) 
# -0.232786  1.1516399

asSp <- SpatialPointsDataFrame(coords = sp, data = data.frame(wStar))
out <- "wStar_June.png"
png(out, width = 800, height = 500)
spplot(asSp, col.regions = cont.palette, col="transparent", colorkey=T,
       main="wStar in June",as.table=T, cuts=seq(-1.152,1.152,0.01))
dev.off()

#### --------------------------------------------------------
