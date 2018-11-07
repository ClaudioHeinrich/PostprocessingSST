rm(list = ls())

library(PostProcessing)
library(data.table)
setwd("~/NR/SFE/")

load("./Data/MultiEnsemble/DT.RData")
DT = DT[!is.na(SST_bar)]

setkey(DT,grid_id,year,month)

ff_mean = function(x){
    x_s = cumsum(x)
    x_s = x_s - x
    return(c(0,x_s[-1]/1:(length(x_s) - 1)))
}

DT[,climatology := ff_mean(SST_bar), .(grid_id, month)]
DT[,obs_anamoly := SST_bar - climatology]
DT[,obs_anamoly_1 := shift(obs_anamoly,1,NA,"lag"), .(grid_id, month)] 
DT[,NPCM_mean:= rowMeans(.SD), .SDcols = paste0("NPCM_Ens",1:9)]
DT[,NPCM_climatology := ff_mean(NPCM_mean), .(grid_id, month)]
DT[,NCPM_anamoly := NPCM_mean - NPCM_climatology]
DT[,GCFS_mean:= rowMeans(.SD), .SDcols = paste0("GCFS_Ens",1:15)]
DT[,GCFS_climatology := ff_mean(GCFS_mean), .(grid_id, month)]
DT[,GCFS_anamoly := GCFS_mean - GCFS_climatology]
Years_test = 1996:2017
coef = NULL
for(i in 1:length(Years_test)){
    y = Years_test[i]
    DT_train = DT[year < y]
    mod0 = lm(obs_anamoly ~  climatology, data = DT_train)
    mod1 = lm(obs_anamoly ~  obs_anamoly_1, data = DT_train)
    mod2 = lm(obs_anamoly ~  climatology + obs_anamoly_1, data = DT_train)
    mod3 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + NCPM_anamoly , data = DT_train)
    mod4 = lm(obs_anamoly ~  NCPM_anamoly , data = DT_train)
    mod5 = lm(obs_anamoly ~  NCPM_anamoly + obs_anamoly_1 , data = DT_train)
    mod6 = lm(obs_anamoly ~  climatology + NCPM_anamoly , data = DT_train)
    mod7 = lm(obs_anamoly ~  climatology + obs_anamoly_1 + NCPM_anamoly + GCFS_anamoly , data = DT_train)
    coef = rbind(coef, mod7$coef)
    DT[year == y, pred0 := climatology + predict(mod0,newdata = DT[year == y])]
    DT[year == y, pred1 := climatology + predict(mod1,newdata = DT[year == y])]
    DT[year == y, pred2 := climatology + predict(mod2,newdata = DT[year == y])]
    DT[year == y, pred3 := climatology + predict(mod3,newdata = DT[year == y])]
    DT[year == y, pred4 := climatology + predict(mod4,newdata = DT[year == y])]
    DT[year == y, pred5 := climatology + predict(mod5,newdata = DT[year == y])]
    DT[year == y, pred6 := climatology + predict(mod6,newdata = DT[year == y])]
    DT[year == y, pred7 := climatology + predict(mod7,newdata = DT[year == y])]
}

gg = function(a,b){return( mean( (a - b)^2))}

Score = DT[year > 1995 & !is.na(pred7),lapply(.SD,gg,SST_bar),.SDcols = c("climatology", paste0("pred",0:7))]

plot(Years_test, coef[,4], ylim  = c(.2,.6), xlab = "Year", ylab = "Coef", pch = 19)
lines(Years_test, coef[,4], lty = 1)
lines(Years_test, coef[,5], col="red")
points(Years_test, coef[,5], pch = 19, col="red")


