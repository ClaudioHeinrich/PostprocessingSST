rm(list = ls())



##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/SFE/PCACov/"
options(max.print = 1e3)



bias_correct = function(dt = NULL,
                        dt.load = TRUE,
                        Y = 1985:2010,
                        M = 1:12,
                        method = "gwa",  # for 'gliding window average'
                        par_1 = 5,  # length of gliding window
                        global_mean_scores = FALSE,
                        reduced_output = FALSE
                        ){
  
  if(dt.load) dt = load_combined_wide()
  
  dt = dt[year %in% Y & month %in% M,]
  dt_new = copy(dt)
  
  
  dt_new = dt_new[!is.na(Ens_bar) & !is.na(SST_bar) ,
                  "Bias_Est" := SMA(SST_bar - Ens_bar,n = par_1),
                  by = .(grid_id, month)]
  
  dt_new = dt_new[year < min(Y) + par_1, Bias_Est := 
                  (cumsum(SST_bar-Ens_bar) - (SST_bar-Ens_bar)) / (year - min(year)+1),
              by = .(grid_id, month)]
  
  
  
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


test_bias_correct = function(saveorgo = TRUE,
                             save.dir = "~/PostClimDataNoBackup/SFE/Derived/"
                             ){
  
  dt = load_combined_wide()
  
  
  sc = list()
  for (k in 1:26){
    print(paste0("gwa length = ",k))
    temp = bias_correct(dt = dt,
                                    dt.load = FALSE,
                                    par_1 = k,
                                    global_mean_scores = TRUE)
    sc[[k]] = temp[,"win_length" := k]
  }
  rbindlist(sc)
  
  if(saveorgo) save(sc, file = paste0(save.dir,"glob.scores.bc.gwa.RData"))
}








