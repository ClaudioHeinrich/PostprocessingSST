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




#' vec and years are vectors of the same length

#' @export
#' 
sim_mov_av = function( l,vec, years ){
  
  all_years = min(years):max(years)
  
  sv = rep(0,length(vec))
  
  for (i in 2:length(vec)){
    year_ind = which(all_years == years[i])
    weight_vec = (year_ind - which(all_years %in% years[1:i-1])) <= l
    if(TRUE %in% weight_vec) weight_vec = weight_vec/sum(weight_vec)
    sv[i] = sum(weight_vec*vec[1:i-1]^2)
  }
  return(sv)
}


crps.na.nrm = function(y, mean, sd){
  
  na_loc = which( is.na(y) | is.na(mean) | is.na(sd) | sd == 0)
  
  if(identical(na_loc,integer(0))){
    x = crps(y, family = "normal", mean = mean, sd = sd)
  }else{
    x = rep(0,length(y))
    x[-na_loc] = crps(y[-na_loc], family = "normal", mean = mean[-na_loc], sd = sd[-na_loc])
    x[na_loc] = NA
  }
  return(x)
}


global_scores = function (DT, eval_years = 2001:2010, averaging = TRUE){
  glob_sc = DT[year %in% eval_years, 
               .("CRPS" = crps.na.nrm(SST_bar,SST_hat,sd = Ens_sd),SST_bar,SST_hat),
               by = .(YM, year, month, grid_id,Lon,Lat)]
  if(averaging){
    glob_sc = glob_sc[,.("mean_CRPS" = mean(CRPS,na.rm = TRUE),
                         "RMSE" = sqrt(mean( (SST_bar - SST_hat)^2, na.rm=TRUE))),
                      by = .(YM,year,month)]
  }
  return(glob_sc)
}

#' @export
bias_correct = function(dt = NULL,
                        method = "sma", # "sma" for 'simple moving average',
                                        # "ema" for 'exponential moving average'
                        par_1 = 16,   # if method == sma, par_1 is the length of window for the sma
                                      # if method == ema, par_1 is the ratio of the exp. mov. av.
                        scores = FALSE,
                        eval_years = 2001:2010,
                        saveorgo = TRUE,
                        data.dir = "~/PostClimDataNoBackup/"){
  if(is.null(dt)) dt = load_combined_wide(bias = FALSE) 
  
  Y = dt[,unique(year)]
  M = dt[,unique(month)]
  
  if(method == "sma"){
    dt = dt[,"Bias_Est" := sim_mov_av(l = par_1, 
                                              vec = SST_bar - Ens_bar, 
                                              years = year),
                    by = .(grid_id, month)]
  }
  if (method == "ema"){
    dt_new = dt_new[,"Bias_Est" := exp_mov_av_bc(SST_bar - Ens_bar, ratio = par_1),
                    by = .(grid_id, month)]
  }
  
  dt[,"SST_hat" := Ens_bar + Bias_Est]
  
  if(saveorgo){
    save(dt, file = paste0(data.dir,"SFE/Derived/dt_combine_wide_bias.RData"))
  }
  
  if(scores){
    mean_sc = global_scores(dt)
    return(glob_sc)
  } else return(dt_new)
  
}

