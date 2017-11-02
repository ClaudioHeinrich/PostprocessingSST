rm(list = ls())



##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/SFE/PCACov/"
options(max.print = 1e3)


#---- for bias correction with exponential moving averages----

exp_mov_av_bc = function (vec,ratio){
  a = c(0,vec[1])   
  for(k in 3:length(vec)){
    temp=ratio * vec[k-1] + (1-ratio)*a[length(a)]
    a = c(a,temp)
  }
  return(a)
}




bias_correct = function(dt = NULL,
                        Y = 1985:2010,
                        M = 1:12,
                        method = "ema",  # "gwa" for 'gliding window average', i.e. simple moving average, 
                                         # "ema" for 'exponential moving average'
                        par_1 = .2,  # if method = gwa, par_1 is the length of gliding window
                                    # if method = ema, par_1 is the ratio of the exp. mov. av.
                        global_mean_scores = FALSE,
                        reduced_output = FALSE
                        ){
  
  if(is.null(dt)) dt = load_combined_wide()
  
  dt = dt[year %in% Y & month %in% M,]
  dt_new = copy(dt)
  
  if(method == "gwa"){
    dt_new = dt_new[!is.na(Ens_bar) & !is.na(SST_bar) ,
                    "Bias_Est" := par_1/(par_1-1) * SMA(SST_bar - Ens_bar,n = par_1) - (SST_bar - Ens_bar)/(par_1-1),
                    by = .(grid_id, month)]
  
    dt_new = dt_new[year < min(Y) + par_1, 
                    Bias_Est := (cumsum(SST_bar-Ens_bar) - (SST_bar-Ens_bar)) / (year - min(year)+1),
                    by = .(grid_id, month)]
  }
  if (method == "ema"){
    dt_new = dt_new[!is.na(Ens_bar) & !is.na(SST_bar) ,
                    "Bias_Est" := exp_mov_av_bc(SST_bar - Ens_bar, ratio = par_1),
                    by = .(grid_id, month)]
  }
  
  
  dt_new[,"SST_hat" := Ens_bar + Bias_Est]
  
  if(global_mean_scores){
    glob_mean_sc = dt_new[,.( "RMSE" = 
                                sqrt(mean( (SST_bar - SST_hat)^2, 
                                           na.rm=TRUE)),
                              "MAE" = mean(abs(SST_bar - SST_hat),
                                           na.rm=TRUE)),
                            keyby = YM]
  return(glob_mean_sc)
    }else if(reduced_output)
    {result = dt_new[,.(grid_id,year,month,YM,SST_hat)]
    return(result)
  } else return(dt_new)
}




#---- get global mean scores for a range of parameters and save them ---

test_bias_correct = function(method = "gwa", # also accepts ema
                             num.years = 26,
                             saveorgo = TRUE,
                             save.dir = "~/PostClimDataNoBackup/SFE/Derived/"
                             ){
  
  dt = load_combined_wide()
  
  
  sc = list()
  
  if(method == "gwa"){
    for (k in 2:num.years){
      print(paste0("gwa length = ",k))
      temp = bias_correct(dt = dt, par_1 = k, global_mean_scores = TRUE)
      sc[[k]] = temp[,"win_length" := k]
    }
    sc = rbindlist(sc)
    if(saveorgo) save(sc, file = paste0(save.dir,"glob.scores.bc.gwa.RData"))
  }
  
  if(method == "ema"){
    ratio.vec = seq(1/num.years,.8,length.out = 20)
    ind = 1
    for (k in ratio.vec){
      print(paste0("ratio = ",k))
      temp = bias_correct(dt = dt, method = "ema",par_1 = k, global_mean_scores = TRUE)
      sc[[ind]] = temp[,"ratio" := k]
      ind=ind + 1
    }
    sc = rbindlist(sc)
    if(saveorgo) save(sc, file = paste0(save.dir,"glob.scores.bc.ema.RData"))
  }
  
}


#--------plot mean scores for gwa --------

load("~/PostClimDataNoBackup/SFE/Derived/glob.scores.bc.gwa.RData")

mean_scores = sc[,.(RMSE_bar = mean(RMSE), MAE_bar = mean(MAE)), by = win_length]

yrange = c(mean_scores[,min(range(RMSE_bar),range(MAE_bar))],
           mean_scores[,max(range(RMSE_bar),range(MAE_bar))])
  
setwd("~/NR/SFE/")
plot.dir = "./figures/"
pdf(paste0(plot.dir,"mean_scores_gwa.pdf"))
plot(x = mean_scores[,win_length],
     y = mean_scores[,RMSE_bar],
     ylim = yrange,
     type = "b",
     col = "blue",
     main = "mean scores for SMA by length of window",
     xlab = "window length",
     ylab = "mean score"
     )
lines( mean_scores[,win_length], 
       mean_scores[,MAE_bar],
       type = "b",
       col = "darkgreen")

#---- add minima -----
abline(h = mean_scores[,min(RMSE_bar)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
abline(h = mean_scores[,min(MAE_bar)], lty = "dashed", col = adjustcolor("darkgreen",alpha = .5))

min_loc_RMSE = mean_scores[,which.min(RMSE_bar)]
min_loc_MAE = mean_scores[,which.min(MAE_bar)]

points(x = mean_scores[,win_length][min_loc_RMSE],
       y = mean_scores[,RMSE_bar][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)
points(x = mean_scores[,win_length][min_loc_MAE],
       y = mean_scores[,MAE_bar][min_loc_MAE],
       col = "darkgreen",
       bg = "darkgreen",
       pch = 21)
         
legend("topright",legend = c("RMSE","MAE"),lty = c(1,1),col = c("blue","darkgreen"))
dev.off()


#--------plot mean scores for ema --------

load("~/PostClimDataNoBackup/SFE/Derived/glob.scores.bc.ema.RData")

mean_scores = sc[,.(RMSE_bar = mean(RMSE), MAE_bar = mean(MAE)), by = ratio]

yrange = c(mean_scores[,min(range(RMSE_bar),range(MAE_bar))],
           mean_scores[,max(range(RMSE_bar),range(MAE_bar))])

setwd("~/NR/SFE/")
plot.dir = "./figures/"
pdf(paste0(plot.dir,"mean_scores_ema.pdf"))
plot(x = mean_scores[,ratio],
     y = mean_scores[,RMSE_bar],
     ylim = yrange,
     type = "b",
     col = "blue",
     main = "Mean scores for EMAs by ratio",
     xlab = "ratio",
     ylab = "mean score"
)
lines( mean_scores[,ratio], 
       mean_scores[,MAE_bar],
       type = "b",
       col = "darkgreen")

#---- add minima -----
abline(h = mean_scores[,min(RMSE_bar)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
abline(h = mean_scores[,min(MAE_bar)], lty = "dashed", col = adjustcolor("darkgreen",alpha = .5))

min_loc_RMSE = mean_scores[,which.min(RMSE_bar)]
min_loc_MAE = mean_scores[,which.min(MAE_bar)]

points(x = mean_scores[,ratio][min_loc_RMSE],
       y = mean_scores[,RMSE_bar][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)
points(x = mean_scores[,ratio][min_loc_MAE],
       y = mean_scores[,MAE_bar][min_loc_MAE],
       col = "darkgreen",
       bg = "darkgreen",
       pch = 21)

legend("right",legend = c("RMSE","MAE"),lty = c(1,1),col = c("blue","darkgreen"))
dev.off()


