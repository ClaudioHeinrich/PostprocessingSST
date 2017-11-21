rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)







month = 1:12
y = 1985:2010


setup_PCA(m = month, y = y, oriented = TRUE)

vec = 1:4
for(k in vec){

calib = forecast_PCA(y = y, m = month, PCA_depth = k, output_opts = "mar_sd" , saveorgo = FALSE)

calib[,"variance" := forecast^2]

calib[,"PIT" := pnorm(SST_bar, mean = SST_hat, sd = sqrt(variance))]


for(moment in 1:2){
  
  if(moment == 1) calib[,"moment":= mean(PIT,na.rm = TRUE),by = grid_id]
  if(moment == 2) calib[,"moment":= var(PIT,na.rm = TRUE),by = grid_id]


  save.dir="./Data/PostClim/SFE/Derived/PCA"
  save(calib, file = paste0(save.dir,"/cal",k,"_mom",moment,".RData"))
}

plot_system(type = "cal", moment = 1, depth = k, plot_title = paste0("PIT mean, ",k, " prin.com."))
plot_system(type = "cal", moment = 2, depth = k, plot_title = paste0("PIT variance, ",k, " prin.com."))


}