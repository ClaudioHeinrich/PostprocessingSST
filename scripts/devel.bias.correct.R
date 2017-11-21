rm(list = ls())



##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)





#---- get global mean scores for a range of parameters and save them ---

test_bias_correct = function(model = "NorESM",
                             method = "gwa", # also accepts ema
                             saveorgo = TRUE,
                             save.dir = "~/PostClimDataNoBackup/SFE/Derived/"
                             ){
  
  
  dt = load_combined_wide(model = model)
  
  if(model == "senorge") setnames(dt,"senorge_grid_id","grid_id")
  
  file.out = paste0(save.dir,"glob.scores.bc.",method,".",model,".RData")
  
  num.years = dt[,range(year)][2] - dt[,range(year)][1] + 1
  
  sc = list()
  
  if(method == "gwa"){
    for (k in 2:num.years){
      print(paste0("gwa length = ",k))
      temp = bias_correct(dt = dt, model = model, par_1 = k, global_mean_scores = TRUE)
      sc[[k]] = temp[,"win_length" := k]
    }
    sc = rbindlist(sc)
    if(saveorgo) save(sc, file = file.out)
  }
  
  if(method == "ema"){
    ratio.vec = seq(1/num.years,.8,length.out = 20)
    ind = 1
    for (k in ratio.vec){
      print(paste0("ratio = ",k))
      temp = bias_correct(dt = dt, method = "ema", model = model, par_1 = k, global_mean_scores = TRUE)
      sc[[ind]] = temp[,"ratio" := k]
      ind=ind + 1
    }
    sc = rbindlist(sc)
    if(saveorgo) save(sc, file = file.out )
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


