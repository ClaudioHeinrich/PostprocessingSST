## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(data.table)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "/home/alex/NR/SFE/Presentations/20181120_TLT/fig_alex/"
print_figs = FALSE

load("./FcNov2018/ts_hindcast_slimmed.RData")

ens_names = c("norcpm_ts_bar", "ecmwf_ts_bar","mf_ts_bar","ukmo_ts_bar")
DT_final[,equal_weighted := rowMeans(.SD), .SDcol = ens_names]
DT_final[,persistence := obs_anamoly_1]
for(y in 2006:2017){

    print(y)
    DT_train = DT_final[year < y]
    mod = lm(obs_anamoly ~  climatology + obs_anamoly_1 + norcpm_anamoly + ecmwf_anamoly + mf_anamoly + ukmo_anamoly, data = DT_train)
    DT_final[year == y,
             trained_model := climatology + predict(mod,newdata = DT_final[year == y])]
}

gg = function(a,b){return( sqrt(mean( (a - b)^2)))}

Lat_nordic = c(55,80)
Lon_nordic = c(5,30)
Score = DT_final[between(Lat,Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]) & year > 2005 & !is.na(trained_model),
                 lapply(.SD,gg,obs_erai_ts),
                 month,
                 .SDcols = c("climatology","equal_weighted","trained_model")]


if(print_figs){pdf(paste0(path_out,"/Model_to_Basic.pdf"))}else{X11()}
rr = range(Score[,-1])
plot(1:12,Score[[2]], type  = "l", xlab = "Month", ylab = "RMSE")
points(1:12,Score[[2]], pch = 20)
for(j in 3:4){
    lines(1:12,Score[[j]], col=j - 1)
    points(1:12,Score[[j]], pch = 20, col = j - 1)
}
legend("topleft", lty = 1, pch = 20, legend = c("Climatology","Equal Weighted", "Trained"), col=1:3)
if(print_figs)dev.off()
