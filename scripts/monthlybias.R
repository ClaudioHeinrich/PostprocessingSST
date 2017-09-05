
## plots the averaged monthly bias for different vintages on the globe and the 2 hemispheres

rm(list = ls())

##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"

options(max.print = 1e3)

##------ Set up -------

dt_complete = list()
k = 1
for(vin in c("mr","Jan","Apr","Jul","Oct"))
{   print(vin)
    temp <- load_combined_vin(vintage = vin) 
    dt_complete[[k]] = rbindlist(temp)
    dt_complete[[k]][,vintage := vin]
    k = k + 1
}

dt_complete <- rbindlist(dt_complete)

# compute biases

dt_complete[,bias := SST_bar - SST_hat_grid]

##---------------------


# get monthly bias averaged over the globe and southern and northern hemisphere

monbias <- dt_complete[,.(mon_bias = mean(bias, na.rm = TRUE)),by = .(month,vintage)]
monbias_sh <- dt_complete[ Lat<0 ,.(mon_bias = mean(na.omit(bias))),by = .(month,vintage)]
monbias_nh <- dt_complete[ Lat>0 ,.(mon_bias = mean(na.omit(bias))),by = .(month,vintage)]


##---------plotting--------

titles <- c("most recent","January","April","July","October")
save.dir <- "./figures/vintagebias/"

k = 0
for(vin in c("mr","Jan","Apr","Jul","Oct"))
  {
  # global
  jpeg(file = paste0(save.dir,"glob",vin,".jpeg"))
  plot(1:12,monbias[,mon_bias[(k*12+1):(k*12+12)]], xlab = "month", ylab = "mean bias",main = paste0("Mean global bias, ",titles[k+1]," vintage"),type = "b")
  dev.off()
  
  # southern hemisphere
  jpeg(file = paste0(save.dir,"sh",vin,".jpeg"))
  plot(1:12,monbias_sh[,mon_bias[(k*12+1):(k*12+12)]], xlab = "month", ylab = "mean bias",main = paste0("Mean SH bias, ",titles[k+1]," vintage"),type = "b")
  dev.off()
  
  #  northern hemisphere
  jpeg(file = paste0(save.dir,"nh",vin,".jpeg"))
  plot(1:12,monbias_nh[,mon_bias[(k*12+1):(k*12+12)]], xlab = "month", ylab = "mean bias",main = paste0("Mean NH bias, ",titles[k+1]," vintage"),type = "b")
  dev.off()
  
    k = k+1
  }

