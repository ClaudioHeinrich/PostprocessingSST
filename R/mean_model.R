#' exponential and simple moving averages for bias correction, deals with missing values by skipping them
#' @export
exp_mov_av_bc = function (vec,ratio){
  if(all(is.na(vec))){as.numeric(rep(NA,times = length(vec)))
    }else{
      na_loc = is.na(vec)
      vec_1 = vec[!is.na(vec)]
      a = c(0,vec_1[1]) 
      for(k in 3:length(vec_1)){
        temp=ratio * vec_1[k-1] + (1-ratio)*a[length(a)]
        a = c(a,temp)
      }
      vec[!is.na(vec)] = a
    return(vec)
    }
}

#' @export
sim_mov_av_bc = function (vec,win_length){
  
  if(all(is.na(vec))){return(as.numeric(rep(NA,times = length(vec))))
  } else {
    a = c(0) 
    for(k in 2:win_length){
      a = c(a,mean(vec[1:(k-1)], na.rm = TRUE))
    }
    if(win_length < length(vec)){
      for(k in (win_length + 1):length(vec)){
        a = c(a,mean(vec[(k-win_length):(k-1)],na.rm = TRUE))
        }
      }
    return(a)
  }
}

#' @export
bias_correct = function(dt = NULL,
                        method = "gwa", # "gwa" for 'gliding window average', i.e. simple moving average, 
                                        # "ema" for 'exponential moving average'
                        par_1 = 16,   # if method = gwa, par_1 is the length of gliding window
                                      # if method = ema, par_1 is the ratio of the exp. mov. av.
                        global_mean_scores = FALSE,
                        reduced_output = FALSE,
                        saveorgo = FALSE,
                        data.dir = "~/PostClimDataNoBackup/"
){
  
  
  if(is.null(dt)) dt = load_combined_wide() 
  
  
  Y = dt[,min(year)]:dt[,max(year)]
    M = 1:12
  
  dt_new = copy(dt)
  
  if(method == "gwa"){
    dt_new = dt_new[,"Bias_Est" := sim_mov_av_bc(SST_bar - Ens_bar, win_length = par_1),
                    by = .(grid_id, month)]
  }
  if (method == "ema"){
    dt_new = dt_new[,"Bias_Est" := exp_mov_av_bc(SST_bar - Ens_bar, ratio = par_1),
                    by = .(grid_id, month)]
  }
  
  
  dt_new[,"SST_hat" := Ens_bar + Bias_Est]
  
  if(saveorgo){dt = dt_new
    save(dt, file = paste0(data.dir,"SFE/Derived/dt_combine_wide_bias.RData"))
  }
  
  
  
  if(global_mean_scores){
    glob_mean_sc = dt_new[,.( "RMSE" = sqrt(mean( (SST_bar - SST_hat)^2, na.rm=TRUE)),
                             "MAE" = mean(abs(SST_bar - SST_hat),na.rm=TRUE)),
                          keyby = YM]
    return(glob_mean_sc)
  }else if(reduced_output)
  {result = dt_new[,.(grid_id,year,month,YM,SST_hat)]
    return(result)
  } else return(dt_new)
  
}

