
rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)


#---- get global mean scores for a range of parameters and save them ---

test_bias_correct = function(method = "gwa", # also accepts ema
                             saveorgo = TRUE,
                             save.dir = "~/PostClimDataNoBackup/SFE/Derived/"
){

  dt = load_combined_wide()
  
  file.out = paste0(save.dir,"glob.scores.bc.",method,".","NorESM",".RData")
  
  num.years = dt[,range(year)][2] - dt[,range(year)][1] + 1
  
  if(method == "gwa"){
    sc = list()
    dummy_function = function(k){
      temp = bias_correct(dt = dt, 
                          model = model, 
                          optimal = FALSE, 
                          method = "gwa", 
                          par_1 = k, 
                          global_mean_scores = TRUE)
      sc[[k]] = temp[,"win_length" := k]
    }
    sc = mclapply(X = 2:num.years, FUN = dummy_function, mc.cores = 8)
    sc = rbindlist(sc)
    if(saveorgo) save(sc, file = file.out)
  }
  
  if(method == "ema"){
    
    ratio.vec = seq(1/num.years,.8,length.out = 25)
    sc = list()
    dummy_function = function(k){
      print(paste0("ratio = ",ratio.vec[k]))
      temp = bias_correct(dt = dt, optimal = FALSE, method = "ema",par_1 = ratio.vec[k], global_mean_scores = TRUE)
      sc[[k]] = temp[,"ratio" := ratio.vec[k]]
    }
    
    sc = mclapply(X = 1:length(ratio.vec), FUN = dummy_function,mc.cores = 8)
    sc = rbindlist(sc)
    if(saveorgo) save(sc, file = file.out )
  }
  
}


#----- these lines take their time -----

#test_bias_correct(method = "gwa")
#test_bias_correct(method = "ema")


#--------- plot RMSE for bias correction with SMA and EMA 

load("~/PostClimDataNoBackup/SFE/Derived/glob.scores.bc.gwa.RData")
sc_sma = sc[YM>2001*12, .( MRMSE = mean(RMSE)), by = win_length]
load("~/PostClimDataNoBackup/SFE/Derived/glob.scores.bc.ema.RData")
sc_ema = sc[YM>2001*12, .( MRMSE = mean(RMSE)), by = ratio]

plot.dir = "./figures/"

# ensure that they are plotted vs the same range
y_range = range(c(sc_sma[,MRMSE],sc_ema[,MRMSE]))  

pdf(paste0(plot.dir,"mean_scores_sma.pdf"))
plot(x = sc_sma[,win_length],
     y = sc_sma[,MRMSE],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = "RMSE for bias correction by SMA",
     xlab = "window length",
     ylab = "mean RMSE"
)

#---- add minimum -----
abline(h = sc_sma[,min(MRMSE)], lty = "dashed", col = adjustcolor("blue",alpha = .5))

min_loc_RMSE = sc_sma[,which.min(MRMSE)]

points(x = sc_sma[,win_length][min_loc_RMSE],
       y = sc_sma[,MRMSE][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

#plot for ema

pdf(paste0(plot.dir,"mean_scores_ema.pdf"))
plot(x = sc_ema[,ratio],
     y = sc_ema[,MRMSE],
     ylim = y_range,
     type = "b",
     col = "darkgreen",
     main = "RMSE for bias correction by EMA",
     xlab = "ratio",
     ylab = "mean RMSE"
)

#---- add minimum -----
abline(h = sc_ema[,min(MRMSE)], lty = "dashed", col = adjustcolor("darkgreen",alpha = .5))

min_loc_RMSE = sc_ema[,which.min(MRMSE)]

points(x = sc_ema[,ratio][min_loc_RMSE],
       y = sc_ema[,MRMSE][min_loc_RMSE],
       col = "darkgreen",
       bg = "darkgreen",
       pch = 21)

dev.off()





load(file = "~/PostClimDataNoBackup/SFE/Derived/dt_combine_wide_bias.RData")

num.years = 26

dt[,"oos_clim" := mean(SST_bar,na.rm = TRUE) - SST_bar/num.years, by = .(month,grid_id)]
dt[,"oos_clim_past" := (cumsum(SST_bar) - SST_bar)/(year - min(year) ),
           by = .(month,grid_id)]

save(dt, file = "~/PostClimDataNoBackup/SFE/Derived/dt_combine_wide_bias.RData")



glob_mean_sc = dt[,.( month,"RMSE_raw" = sqrt(mean( (SST_bar - Ens_bar)^2, na.rm=TRUE)),
                              "RMSE_bc" = sqrt(mean( (SST_bar - SST_hat)^2, na.rm=TRUE)),
                              "MAE" = mean(abs(SST_bar - SST_hat), na.rm=TRUE),
                              "RMSE_clim" = sqrt(mean( (SST_bar - oos_clim)^2, na.rm=TRUE)),
                              "MAE_clim" = mean(abs(SST_bar - oos_clim), na.rm=TRUE),
                              "RMSE_clim_past" = sqrt(mean( (SST_bar - oos_clim_past)^2, na.rm=TRUE)),
                              "MAE_clim_past" = mean(abs(SST_bar - oos_clim_past), na.rm=TRUE)),
                          keyby = YM]

reduced_mean_sc = glob_mean_sc[,.("RMSE_raw" = mean(RMSE_raw), 
                                  "RMSE_bc" = mean(RMSE_bc), 
                                  "RMSE_clim" = mean(RMSE_clim)),by = month]


rr = range(reduced_mean_sc[,.(RMSE_raw,RMSE_bc,RMSE_clim)])

plot.dir = "./figures/"
pdf(paste0(plot.dir,"RMSE_by_month_3.pdf"))

plot(x = reduced_mean_sc[,month],
     y = reduced_mean_sc[,RMSE_raw],
     col = "blue",
     type = "b",
     ylim = rr,
     main = "RMSE by month",
     ylab = "RMSE",
     xlab = "month")
lines(x = reduced_mean_sc[,month],
      y = reduced_mean_sc[,RMSE_bc],
      col = "darkred",
      type = "b")
lines(x = reduced_mean_sc[,month],
      y = reduced_mean_sc[,RMSE_clim],
      col = "darkgreen",
      type = "b")
legend("topleft", 
       legend = c("raw forecast","climatology","bias corrected forecast"), 
       col = c("blue","darkgreen","darkred"),
       lty = c(1,1,1))
dev.off()


load("~/PostClimDataNoBackup/SFE/Derived/glob.scores.bc.ema.RData")

mean_scores = sc[,.(RMSE_bar = mean(RMSE,na.rm = TRUE), MAE_bar = mean(MAE,na.rm = TRUE)), by = ratio]


#---- plot ---

plot.dir = "./figures/"
pdf(paste0(plot.dir,"mean_scores_ema_clim.pdf"))


yrange = c(mean_scores[,range(RMSE_bar)], mean_scores[,range(MAE_bar)], range(reduced_mean_sc))
yrange = range(yrange)

plot(x = mean_scores[,ratio],
     y = mean_scores[,RMSE_bar],
     ylim = yrange,
     type = "b",
     col = "blue",
     main = "mean scores for SMA by length of window",
     xlab = "window length",
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

# ---- add climatology scores ---------

abline(h = reduced_mean_sc[,RMSE_clim], lty = "dashed", col = "darkred")
abline(h = reduced_mean_sc[,MAE_clim], lty = "dashed", col = "darkred")
abline(h = reduced_mean_sc[,RMSE_clim_past], lty = "dashed", col = "yellow4")
abline(h = reduced_mean_sc[,MAE_clim_past], lty = "dashed", col = "yellow4")


legend("topright",legend = c("RMSE","MAE","full clim.", "past clim."),lty = c(1,1),col = c("blue","darkgreen","darkred","yellow4"))
dev.off()




