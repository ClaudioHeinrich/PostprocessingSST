rm(list = ls())

library(PostProcessing)

setwd("~/NR/SFE/")

load("./Data/MultiEnsemble/DT.RData")
DT = DT[!is.na(SST_bar)]

setkey(DT,grid_id,year,month)

DT[,global_mean:=mean(SST_bar),grid_id]


mm = DT[,median(global_mean)]

gg = 58178

DT[,mean_y:=mean(SST_bar),.(grid_id,year)]
DT[,mean_m:=mean(SST_bar),.(grid_id,month)]
DT[,mean_e_1:=mean(NPCM_Ens1,na.rm=TRUE),.(grid_id,month)]
DT_i = DT[grid_id == gg]

plot(DT_i[,SST_bar], type="l")
points(DT_i[,SST_bar], pch = 20)
lines(DT_i[,mean_y], col="blue")
lines(DT_i[,NPCM_Ens1], col="red")
points(DT_i[,NPCM_Ens1], col="red", pch=20)

lines(DT_i[,GCFS_Ens1], col="green")
points(DT_i[,GCFS_Ens1], col="green", pch=20)

lines(DT_i[,mean_m], col="purple")
points(DT_i[,mean_m], col="purple", pch=20)


X11();plot(DT_i[,SST_bar - mean_m],type="l",ylim = c(-3,4))
lines(DT_i[,NPCM_Ens1 - mean_e_1], col="red")

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

DT[,climatology := ff_mean(SST_bar), .(grid_id, month)]
DT[,obs_anamoly:=SST_bar - climatology]
DT[,obs_anamoly_1 := shift(obs_anamoly,1,NA,"lag"),.(grid_id,month)] 
        
DT[year > 1995,.("RMSE" = sqrt( mean( (climatology - SST_bar)^2, na.rm=TRUE)),
      "RMSE_1" = sqrt( mean( (climatology + .27 * obs_anamoly_1 + .208 - SST_bar)^2, na.rm=TRUE))),
   month]
