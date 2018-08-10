

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

sim_mov_av = function(l,vec, years, skip = 0, twosided = FALSE){
  
  #take out NA values:
  na_loc = which(is.na(vec))
  if(!identical(na_loc,integer(0)))
  {
    vec = vec[-na_loc]
    years = years[-na_loc]
  }
  
  
  all_years = min(years):max(years)
  
  sma = rep(0,length(vec))
  
  if(!twosided)
  {
    for (i in 1:length(vec)){
      y=years[i]
      weight_vec = rep(0,length(vec))
      rel_loc = (y > years) & ((y - years) <= l)
      skip_loc =  y-years <= skip
      weight_vec = rel_loc & !skip_loc
      if(TRUE %in% weight_vec) weight_vec = weight_vec/sum(weight_vec)
      sma[i] = sum(weight_vec*vec)
    }
    return(sma)  
  } else if(twosided)
  {
    for (i in 1:length(vec)){
      y=years[i]
      weight_vec = rep(0,length(vec))
      rel_loc = abs(years - y) <= l 
      skip_loc = abs(years - y) <= skip
      weight_vec = rel_loc & !skip_loc
      if(TRUE %in% weight_vec) weight_vec = weight_vec/sum(weight_vec)
      sma[i] = sum(weight_vec*vec)
    }
    return(sma)
  }
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


exp_mov_av = function( a,vec, years, skip = 0,twosided = FALSE){
  
  
  #take out NA values:
  na_loc = which(is.na(vec))
  if(!identical(na_loc,integer(0)))
  {
    vec = vec[-na_loc]
    years = years[-na_loc]
  }
  
  lv = length(vec)
  
  all_years = min(years):max(years)
  
  exp_weights = exp(-a * (0:(length(all_years)-1)))
  
  ema = rep(0,lv)
  
  if(!twosided)
  {
    for (i in (skip + 2):lv){
      weight_vec = rev(exp_weights[which(all_years %in% years[1:(i-1-skip)])])
      weight_vec = weight_vec/sum(weight_vec)
      ema[i] = sum(weight_vec*vec[1:(i-1-skip)])
    }
  } else if(twosided)
  {
      for (i in 1:lv){
        year_dist = abs(years[1 : lv] - years[i])
        weights = rep(0,lv)
        weights[!(year_dist <= skip)] = exp_weights[year_dist]
        weights = weights/sum(weights)
        ema[i] = sum(weights*vec)
    }
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
#' @param mean_est
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
                      mean_est = "sv",
                      save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                      file_name = "dt_combine_wide_bias_sd.RData")
{# get the average variation of the bias corrected truncated ensemble members around the true value by year, month and grid_id:
  if(mean_est == "bcf"){
    get_variances = function(i)
    {
      var_setup_dt = as.data.table(dt[, (trc(.SD + Bias_Est) - SST_bar)^2/ens_size,.SDcols = paste0("Ens",i)])
      #var_setup_dt = as.data.table(var_setup_dt)
      return(var_setup_dt)
    }
    
    var_list = parallel::mclapply(1:ens_size,get_variances,mc.cores = ens_size)
    var_list = as.data.table(var_list)
    setnames(var_list,paste0("temp",1:ens_size))
  
    dt[,"var_bar" := rowSums(var_list)]
      
  } else if (mean_est == "sm"){
    get_variances = function(i)
    {
    var_setup_dt = as.data.table(dt[, (.SD - Ens_bar)^2/(ens_size-1),.SDcols = paste0("Ens",i)])
    #var_setup_dt = as.data.table(var_setup_dt)
    return(var_setup_dt)
    }
    
    var_list = parallel::mclapply(1:ens_size,get_variances,mc.cores = ens_size)
    var_list = as.data.table(var_list)
    setnames(var_list,paste0("temp",1:ens_size))
  
    dt[,"var_bar" := rowSums(var_list)]
  } else if (mean_est == "sv") 
  {
    dt[,"var_bar" := (SST_bar-SST_hat)^2]
  }
  
  
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

bias_correct = function(dt, method, par_1,
                        scores = FALSE,
                        eval_years = 2001:2010,
                        saveorgo = TRUE,
                        save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                        file_name = "dt_combine_wide_bias.RData",
                        skip = 0){
  
  
  if(method == "sma"){
       dt = dt[,"Bias_Est" := sim_mov_av(l = par_1, 
                                      vec = SST_bar - Ens_bar, 
                                      years = year,
                                      skip = skip),
                    by = .(grid_id, month)]
  }
  if (method == "ema"){
       dt[,"Bias_Est" := exp_mov_av(a = par_1,
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



#' computes the RMSE when instead of bias correction the linear model SST_hat = a + b Ens_bar is used (as in standard NGR).
#' 
#' @param DT the data table.
#' @param months The considered months.
#' @param training_years,validation_years Training years and validation years.
#' 
#' @return data table containing the RMSEs for the model above, where the coefficients are estimated in three different ways: grouped by month, location and by both.
#'
#' @examples \dontrun{ DT = load_combined_wide()
#'                     bias_lr(DT = DT)}
#'   
#' @author Claudio Heinrich
#' 
#' @export

bias_lr = function(DT,
                   months = 1:12,
                   validation_years = 2001:2010)
{
  
  for(y in validation_years)
  {
    print(y)
    # grouped by month
    
    fits_by_month = lme4::lmList(formula = SST_bar ~ 1 + Ens_bar | month,
                                 data=DT[year < y,])
    months = as.character(DT[,month])
    DT[,c("a","b"):= coef(fits_by_month)[months,]]
    DT[year == y, T_hat_lr_m := a + b * Ens_bar]
    
    # grouped by location
    
    fits_by_loc = lme4::lmList(formula = SST_bar ~ 1 + Ens_bar | grid_id,
                         data=DT[year < y,])
    
    grid_ids = as.character(DT[,grid_id])
    DT[,c("a","b"):= coef(fits_by_loc)[grid_ids,]]
    DT[year == y,T_hat_lr_loc := a + b * Ens_bar]
    
    # by both
    
    DT[,m_gid := interaction(month,grid_id)]
    fits_by_both = lme4::lmList(formula = SST_bar ~ 1 + Ens_bar | m_gid,
                          data=DT[year < y,])
    m_gids = as.character(DT[,m_gid])
    DT[,c("a","b"):= coef(fits_by_both)[m_gids,]]
    DT[year == y,T_hat_lr_both := a + b * Ens_bar]
    
  }
  
  
  RMSE_lr_m = sqrt(DT[year %in% validation_years,mean((T_hat_lr_m-SST_bar)^2, na.rm = TRUE)])
  RMSE_lr_loc = sqrt(DT[year %in% validation_years,mean((T_hat_lr_loc-SST_bar)^2, na.rm = TRUE)])
  RMSE_lr_both = sqrt(DT[year %in% validation_years,mean((T_hat_lr_both-SST_bar)^2, na.rm = TRUE)])
  
  new_col_names = c("a","b","T_hat_lr_m","T_hat_lr_loc","T_hat_lr_both","m_gid")
  DT[,(new_col_names) := NULL]
  
  RMSE_linear_models = data.table(RMSE_lr_m = RMSE_lr_m, RMSE_lr_loc = RMSE_lr_loc, RMSE_lr_both = RMSE_lr_both)
  
  return(RMSE_linear_models)
  
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

bias_correct_training = function(dt = NULL,
                                method,
                                training_years = 1985:2000,
                                eval_years = 2001:2010,
                                saveorgo = TRUE,
                                save_dir = "~/PostClimDataNoBackup/SFE/Derived/",
                                file_name = "dt_combine_wide_bias.RData",
                                skip = 0){
  
  if(is.null(dt)) dt = load_combined_wide(bias = TRUE) 
  
  if(method[1] == "sma")
    {
      dt[year %in% training_years,"Bias_Est" :=   sim_mov_av(l = as.double(method[2]),
                                                           vec = SST_bar - Ens_bar,
                                                           years = year,
                                                           skip = skip,
                                                           twosided = TRUE),
         by = .(grid_id, month)]
  }
  if (method[1] == "ema"){
    dt[year %in% training_years,"Bias_Est" :=   exp_mov_av(a = as.double(method[2]),
                                                           vec = SST_bar - Ens_bar,
                                                           years = year,
                                                           skip = skip,
                                                           twosided = TRUE),
       by = .(grid_id, month)]
    
  }
  
  dt[,"SST_hat" := trc(Ens_bar + Bias_Est)]
  
  if(saveorgo){
    save(dt, file = paste0(save_dir,file_name))
  } 
    
  return(dt)
}

# Variance as in NGR

var_est_NGR = function(dt,
                       months = 1:12,
                       validation_years = 2001:2010)
{
  na_loc = which(dt[,is.na(SST_bar) | is.na(Ens_bar)])
  dt_new = dt[-na_loc,]
  
  #grouped by month:
  
  var_est_by_month = as.data.table(expand.grid(validation_years,months))
  setnames(var_est_by_month,c("year","month"))
  
  
  print("minimize CRPS for data grouped by month:")
  for(y in validation_years)
  {
    print(y)
    # grouped by month
    for(m in months)
    {
      temp = dt_new[year < y & month == m,]
      CRPS_by_month = function(cd)
      {
        return(mean(crps.na.rm(temp[,SST_bar], temp[,SST_hat], cd[1]^2 + cd[2]^2 * temp[,Ens_sd]), na.rm = TRUE))
      } 
      opt_par = optim(par = c(0,1),fn = CRPS_by_month)
      var_est_by_month[year == y & month == m, "c":= opt_par$par[1]]
      var_est_by_month[year == y & month == m, "d":= opt_par$par[2]]
      }
  }
  
  dt = merge(dt,var_est_by_month,by = c("year","month"),all.x = TRUE)
  dt[,SD_hat_lr_bm := c^2 + d^2*Ens_sd]
  dt[,c("c","d"):=NULL]
  
  
  #grouped by location:
  
  print("minimize CRPS for data grouped by location:")
  
  vebl_parallel = function(y)
  { return_DT = data.table(grid_id = dt[,unique(grid_id)], year = y)
    for(gid in dt_new[,unique(grid_id)])
    {
      temp = dt_new[year < y & grid_id == gid,]
      CRPS_by_gid = function(cd)
      {
        return(mean(crps.na.rm(temp[,SST_bar], temp[,SST_hat], cd[1]^2 + cd[2]^2 * temp[,Ens_sd]), na.rm = TRUE))
      } 
      opt_par = optim(par = c(0,1),fn = CRPS_by_gid)
      return_DT[grid_id == gid, "c" := opt_par$par[1]]
      return_DT[grid_id == gid, "d" := opt_par$par[2]]
    }
  return(return_DT)
  }
  
  var_est_by_loc = parallel::mclapply(validation_years,vebl_parallel,mc.cores =  min(12,length(validation_years)))
  var_est_by_loc = rbindlist(var_est_by_loc)
  
  dt = merge(dt,var_est_by_loc,by = c("year","grid_id"),all.x = TRUE)
  dt[,SD_hat_lr_bl := c^2 + d^2*Ens_sd]
  dt[,c("c","d"):=NULL]
  
  #grouped by both:
  
  print("minimize CRPS for data grouped by both:")
  
  vebb_parallel = function(y)
  { return_DT = as.data.table(expand.grid(dt_new[,unique(grid_id)],months))
    setnames(return_DT,c("grid_id","month"))
    return_DT[,year := y]
    
    for(m in months)
    {
      temp = dt_new[year < y & month == m,]
      
      for(gid in dt_new[,unique(grid_id)])
      {
        temp_2 = temp[grid_id == gid,]
        CRPS_by_both = function(cd)
        {
          return(mean(crps.na.rm(temp_2[,SST_bar], temp_2[,SST_hat], cd[1]^2 + cd[2]^2 * temp_2[,Ens_sd]), na.rm = TRUE))
        } 
        opt_par = optim(par = c(0,1),fn = CRPS_by_both)
        return_DT[month == m & grid_id == gid, "c":= opt_par$par[1]]
        return_DT[month == m & grid_id == gid, "d":= opt_par$par[2]]
      }  
    }
    return(return_DT)
  }
  
  var_est_by_both = parallel::mclapply(validation_years,vebb_parallel,mc.cores =  min(12,length(validation_years)))
  var_est_by_both = rbindlist(var_est_by_both)
  
  dt = merge(dt,var_est_by_both,by = c("year","grid_id","month"),all.x = TRUE)
  dt[,SD_hat_lr_bb := c^2 + d^2*Ens_sd]
  dt[,c("c","d"):=NULL]

  return(dt)
  
}




