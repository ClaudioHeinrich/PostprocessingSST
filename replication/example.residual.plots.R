rm(list = ls())

library(SeasonalForecasting)
library(data.table)
library(maps)
library(fields)
setwd("~/NR/SFE")
options(max.print = 1e3)


#---- get example residual plots for the years and months specified in Yvec and Mvec ----

dt = load_combined_wide(bias = TRUE)
dt[,"res":= SST_hat - SST_bar]

Yvec = 2004:2008
Mvec = 9

#--- get maximum range for uniform plotting scale ----

rr = c()
rr_sd = c()

for(Y in Yvec){
  for(M in Mvec){
    rr = c(rr,range(dt[year == Y & month == M,res],na.rm = TRUE))
    rr_sd = c(rr_sd,range(dt[year == Y & month == M,Ens_sd],na.rm = TRUE))
  }
}
  

rr = range(rr)
rr_sd = range(rr_sd)


for(Y in Yvec){
  for(M in Mvec){
    
    plot_system(Y = Y, M=M, type = "res", print_figs = TRUE, plot_title = paste0("Residual 0",M," / ",Y))
  }
}




for(Y in Yvec){
  for(M in Mvec){
    plot_diagnostic(dt[year == Y & month == M,.(Lon,Lat,Ens_sd)], mn = paste0("Ensemble spread 0",M," / ",Y),rr = rr_sd, save.pdf = TRUE, file.name = paste0("Ens_spread_",M,"_",Y) )
  }
}

