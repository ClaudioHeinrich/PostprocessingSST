

#' computes simple moving averages, is resistant to missing values  
#'
#' @param l length of the averaging window.
#' @param vec,years vectors of the same length, vec[i] contains the value corresponding to year years[i]
#' @skip Integer. If skip = n > 0 the moving average skips the most recent n years. Useful if realizations in the most recent past are missing.
#'                                    
#' @return a vector of the same length as the input vectors. At location j it contains the average of the values contained in vec that fall into the period of the last l years (excluding the present). First entry is 0.
#'
#' @author Claudio Heinrich
#' @examples sim_mov_av(5, rnorm(10), c(1990,1993,1995:2002))
#' 
#' 
#' @export

sim_mov_av = function(l,vec, years, skip = 0){
  
  all_years = min(years):max(years)
  
  sma = rep(0,length(vec))
  
  for (i in (2 + skip) :length(vec)){
    year_ind = which(all_years == years[i])
    weight_vec = ((year_ind - which(all_years %in% years[1:(i-1)])) <= l ) & ((year_ind - which(all_years %in% years[1:(i-1)])) > skip)
    if(TRUE %in% weight_vec) weight_vec = weight_vec/sum(weight_vec)
    sma[i] = sum(weight_vec*vec[1:(i-1)])
  }
  return(sma)
}

#' computes exponentially weighted moving averages, is resistant to missing values  
#'
#' @param a weight parameter.
#' @param vec,years vectors of the same length, vec[i] contains the value corresponding to year years[i]
#'                  
#' @return a vector of the same length as the input vectors. At location j it contains the average of the past entries of vec, weighted by exp(-a*d) where d is the distance to the current year. First entry is 0.
#' @skip Integer. If skip = n > 0 the moving average skips the most recent n years. Useful if realizations in the most recent past are missing.
#'
#' @author Claudio Heinrich
#' @examples exp_mov_av(.1, rnorm(10), c(1990,1993,1995:2002))
#' 
#' 
#' @export


exp_mov_av = function( a,vec, years, skip = 0 ){
  
  all_years = min(years):max(years)
  
  exp_weights = rep(0,length(all_years))
  for(i in 1:length(all_years)){
    exp_weights[i] = exp(-a*(i-1))
  }
  
  ema = rep(0,length(vec))
  
  for (i in (skip + 2):length(vec)){
    weight_vec = rev(exp_weights[which(all_years %in% years[1:(i-1-skip)])])
    weight_vec = weight_vec/sum(weight_vec)
    ema[i] = sum(weight_vec*vec[1:(i-1-skip)])
  }
  return(ema)
}



#' crps score for normal distribution, resistant to missing values and 0 standard deviation  
#'
#' @param y vector of observations.
#' @param mean,sd mean and sd of the forecasting distribution
#'                   
#' @return vector of the same length as y containing crps scores.
#'
#' @author Claudio Heinrich
#' @examples crps.na.rm(c(NA,rnorm(10)), 1,1)
#' 
#' @importFrom scoringRules crps
#' 
#' @export


crps.na.rm = function(y, mean, sd){
  
  na_loc = which( is.na(y) | is.na(mean) | is.na(sd) | sd == 0)
  
  if(identical(na_loc,integer(0))){
    x = scoringRules::crps(y, family = "normal", mean = mean, sd = sd)
  }else{
    x = rep(0,length(y))
    x[-na_loc] = scoringRules::crps(y[-na_loc], family = "normal", mean = mean[-na_loc], sd = sd[-na_loc])
    x[na_loc] = NA
  }
  return(x)
}

#' Computes scores for bias correction and variance estimation
#'
#' @param DT The data table
#' @param ens_size Size of the ensemble.
#' @param eval_years The years to compute the scores for
#' @param var Logical. If TRUE the CRPS is computed, else the MSE of the bias corrected forecast is computed
#'                   
#' @return A one-row data table containing the mean score.
#'
#' @author Claudio Heinrich
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' global_mean_scores(DT)}
#' 
#' 
#' @export


global_mean_scores_EnsBMA = function (DT, ens_size = 9, eval_years = 2001:2010, var = TRUE)
{
  if(var){
    # marginal distribution is a mixture distribution of normal distributions centered at the ensemble members + Bias_Est, all with variance SD_hat^2.
    
    na_ids = DT[year %in% eval_years,which(is.na(Ens_bar) | is.na(Bias_Est) | is.na(SD_hat))]
    DT_temp = DT[year %in% eval_years,][-na_ids,]
    
    obs = DT_temp[,SST_bar]
    means = as.matrix(DT_temp[,trc(.SD + Bias_Est),.SDcols = paste0("Ens",1:ens_size)])
    sds = rep(DT_temp[,SD_hat], ens_size)
    sds = matrix(sds, ncol = ens_size)
    # get CRPS for each location for Gaussian mixture distribution
    mar_CRPS = scoringRules::crps_mixnorm(obs,means,sds)
    glob_mean_sc = data.table("CRPS" = mean(mar_CRPS))
    
  } else {
    glob_mean_sc = DT[year %in% eval_years, 
                      .("MSE" = mean( (SST_bar - SST_hat)^2, na.rm=TRUE))]
  }
  
  return(glob_mean_sc)
}


global_mean_scores = function (DT, eval_years = 2001:2010, var = TRUE){
  
  if(var){
    glob_mean_sc = DT[year %in% eval_years, 
                    .("CRPS" = mean (crps.na.rm(SST_bar, mean = SST_hat,sd = SD_hat),na.rm = TRUE))]
  } else glob_mean_sc = DT[year %in% eval_years, 
                           .("MSE" = mean( (SST_bar - SST_hat)^2, na.rm=TRUE))]
  
  return(glob_mean_sc)
}

#' sets up the estimation of marginal variance by computing the ensemble variance by location, month and year.
#'
#' @param dt The data table containing the estimated bias.
#' @param ens_size integer. Size of the ensemble.
#' @param saveorgo logical. Do we save the data table?
#' @param save_dir,file_name character strings. Directory and file name for saving.
#' 
#' @return the data table dt with a new column labelled 'var_bar'
#' 
#' @author Claudio Heinrich
#' 
#' @examples 
#' \dontrun{load_combined_wide()
#'          ens_sd_est(dt)}
#'          
#' @export

ens_sd_est = function(dt, 
                      ens_size = 9,
                      saveorgo = TRUE,
                      mean_est = "bcf",
                      save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                      file_name = "dt_combine_wide_bias_sd.RData")
{# get the average variation of the bias corrected truncated ensemble members around the true value by year, month and grid_id:
  if(mean_est == "bcf"){
    get_variances = function(i)
    {
      var_setup_dt = as.data.table(dt[, (.SD + Bias_Est - SST_bar)^2/ens_size,.SDcols = paste0("Ens",i)])
      #var_setup_dt = as.data.table(var_setup_dt)
      return(var_setup_dt)
    }
  } else if (mean_est == "sm"){
    get_variances = function(i)
    {
    var_setup_dt = as.data.table(dt[, (.SD - Ens_bar)^2/(ens_size-1),.SDcols = paste0("Ens",i)])
    #var_setup_dt = as.data.table(var_setup_dt)
    return(var_setup_dt)
    }
  }
  
  var_list = parallel::mclapply(1:ens_size,get_variances,mc.cores = ens_size)
  var_list = as.data.table(var_list)
  setnames(var_list,paste0("temp",1:ens_size))
  
  dt[,"var_bar" := rowSums(var_list)]
  
return(dt)
}



#' Applies bias correction with a specified method to the data and saves or returns scores.
#'
#' @param dt The data table.
#' @param method Method of bias correction. Takes "sma" for simple moving average and "ema" for exponential moving average. 
#' @param par_1 Numeric. If method == "sma", par_1 is the (integer) length of the moving average window, if method == "ema", par_1 is the scale parameter for the exponential downscaling, typically in (0,1).
#' @param scores Logical. If true, the MSE is returned.
#' @param eval_years Numerical vector. The years for evaluating the score.
#' @param saveorgo Logical. If TRUE, the data table with corrected SST_hat and new column Bias_Est is saved.
#' @param save_dir,file_name Directory and name for the saved file.
#' @param skip Integer. Passed on to sim_mov_av or exp_mov_av.
#'                   
#'                   
#' @return The data table with corrected SST_hat and new column Bias_Est.
#'
#' @author Claudio Heinrich
#' @examples \dontrun{bias_correct(saveorgo = FALSE)}
#' 
#' 
#' @export

bias_correct = function(dt = NULL,
                        method = "sma", # "sma" for 'simple moving average',
                                        # "ema" for 'exponential moving average'
                        par_1 = 16,   # if method == sma, par_1 is the length of window for the sma
                                      # if method == ema, par_1 is the ratio of the exp. mov. av.
                        scores = FALSE,
                        eval_years = 2001:2010,
                        saveorgo = TRUE,
                        save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                        file_name = "dt_combine_wide_bias.RData",
                        skip = 0){
  
  if(is.null(dt)) dt = load_combined_wide(bias = FALSE) 
  
  if(method == "sma"){
    dt = dt[,"Bias_Est" := sim_mov_av(l = par_1, 
                                      vec = SST_bar - Ens_bar, 
                                      years = year,
                                      skip = skip),
                    by = .(grid_id, month)]
  }
  if (method == "ema"){
    dt = dt[,"Bias_Est" := exp_mov_av(a = par_1,
                                      vec = SST_bar - Ens_bar,
                                      years = year,
                                      skip = skip),
                    by = .(grid_id, month)]
  }
  
  dt[,"SST_hat" := trc(Ens_bar + Bias_Est)]
  
  if(saveorgo){
    save(dt, file = paste0(save_dir,file_name))
  }
  
  if(scores){
    mean_sc = global_mean_scores(dt, eval_years = eval_years, var = FALSE)
    return(mean_sc)
  } else return(dt)
  
}



#' Estimates standard deviation with a specified method to the data and saves or returns scores.
#'
#' @param dt The data table.
#' @param method Method of variance estimation. Takes "sma" for simple moving average and "ema" for exponential moving average. 
#' @param par_1 Numeric. If method == "sma", par_1 is the (integer) length of the moving average window, if method == "ema", par_1 is the scale parameter for the exponential downscaling, typically in (0,1).
#' @param scores Logical. If true, the CRPS is returned.
#' @param eval_years Numerical vector. The years for evaluating the score.
#' @param saveorgo Logical. If TRUE, the data table with new column SD_hat is saved.
#' @param save_dir,file_name Directory and name for the saved file.
#'                   
#' @return The data table dt containing a new column SD_hat.
#'
#' @author Claudio Heinrich
#' @examples \dontrun{sd_est(saveorgo = FALSE)}
#' 
#' 
#' @export

sd_est = function(dt = NULL,
                  method = "sma", # "sma" for 'simple moving average',
                        # "ema" for 'exponential moving average'
                  par_1 = 16,   # if method == sma, par_1 is the length of window for the sma
                        # if method == ema, par_1 is the ratio of the exp. mov. av.
                  scores = FALSE,
                  ens_mean = FALSE,
                  ens_size = 9,
                  eval_years = 2001:2010,
                  saveorgo = TRUE,
                  save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                  file_name = "dt_combine_wide_bc_var.RData"){
  
  if(is.null(dt)) 
    {
    dt = load_combined_wide(bias = TRUE)
    }
  
  
  if(method == "sma"){
    dt[year != min(year),"SD_hat" := sqrt(sim_mov_av(l = par_1, 
                                                     vec = var_bar, 
                                                     years = year)),
       by = .(grid_id, month)]
  }
  if (method == "ema"){
    dt = dt[year != min(year),"SD_hat" := sqrt(exp_mov_av(a = par_1,
                              vec = var_bar, 
                              years = year)),
            by = .(grid_id, month)]
  }
  
  if(saveorgo){
    save(dt, file = paste0(save_dir,file_name))
  }
  
  if(scores){
    mean_sc = global_mean_scores(dt, eval_years = eval_years, var = TRUE)
    return(mean_sc)
  } else return(dt)
  
}

