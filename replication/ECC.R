
############################################################################
############### Creates a postprocessed forecast by ECC ####################
############################################################################


#' ECC forecasts 
#' 
#' @description Generates forecasts post-processed by ensemble copula coupling (ECC)
#' 
#' @param dt the data table containing bias and variance estimates.
#' @param y,m vectors contining the years and months for the forecast, default are all years and months in dt.
#' @param ens_size Integer. Size of the NWP ensemble.
#' @param method takes "R" or "Q". Method of ECC, see Schefzik et al. 2013.
#' @param truncate logical, whether the forecasted temperature is truncated at freezing temperature.
#' @param saveorgo Logical, whether we save or not. 
#' @param save_dir,file_name The directory to save in and the name of the saved file.
#' 
#' @return data table containing \code{ens_size} columns that are marginally calibrated, labelled "fc"1:ens_size, and the calibrated forecast by ECC, labelled "ecc_fc"1:ens_size
#' 
#' @examples \dontrun{ DT = load_combined_wide(data_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO", output_name = "dt_combine_NAO_wide_bc_var.RData") \cr
#'                     forecast_ECC(dt = DT,save_dir = "./Data/PostClim/SFE/Derived/NAO/")}
#' 
#' @author Claudio Heinrich        
#' 
#' @importFrom matrixStats rowRanks
#' 
#' @export



forecast_ECC = function(dt,
                        y = NULL,
                        m = NULL,
                        ens_size = 9,
                        method = "Q",
                        saveorgo = TRUE,
                        save_dir = "./Data/PostClim/SFE/Derived/ECC/",
                        file_name = "ECC_fc.RData")
{
  # prepare data table:
  
  if(!is.null(y))
  {
    dt = dt[year %in% y,]
  }
  
  
  if(!is.null(m))
  {
    dt = dt[month %in% m,]
  }  
  
  na_rows = dt[ ,which(is.na(Bias_Est) |  is.na(Ens_bar) | is.na(SD_hat) )]
  
  
  
  if(!identical(na_rows,integer(0)))
  {
    dt_temp = dt[-na_rows,]
  } else {
    dt_temp = dt
  }
  
  # do univariate postprocessing:  
  
  length_norm = dt_temp[,.N]
  
  if(method == "R")
  {
    for(i in 1:ens_size)
    {
      dt_temp[ ,paste0("fc",i) := rnorm(n = length_norm, 
                                        mean = dt_temp[,SST_hat], 
                                        sd = dt_temp[,SD_hat])]
    }
  }
  
  if(method == "Q")
  {
    for(i in 1:ens_size)
    {
      dt_temp[ ,paste0("fc",i) := qnorm(i/(ens_size + 1), 
                                        mean = dt_temp[,SST_hat], 
                                        sd = dt_temp[,SD_hat])]
    }
  }
  
  
  # get rank order of the ensemble and reorder the post-processed ensemble to mach the rank order statistic of the ensemble
  
  rks_ens = dt_temp[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("Ens",1:9)]
  
  rks_fc = dt_temp[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("fc",1:9)]
  
  fcs = as.matrix(dt_temp[,.SD,.SDcols = paste0("fc",1:ens_size)])
  
  num_row = nrow(fcs)
  
  ecc_fcs = matrix(NA,nrow = num_row,ncol = ens_size)
  for(i in 1:num_row){
    if( i %% 100000 == 0)
    {
      print(paste0(i," / ",num_row))
    }
    if(!is.na(fcs[i,1]))
    { 
      ecc_fcs[i,] = fcs[i, match(rks_ens[i,],rks_fc[i,])]
    }
  }
  
  ecc_fcs = data.table(ecc_fcs)
  setnames(ecc_fcs,paste0("ecc_fc",1:ens_size))
  
  dt_temp = data.table(dt_temp,ecc_fcs)
  
  # put back missing locations

  if(!identical(na_rows,integer(0)))
  {
    key_dt = key(dt)
    dt = rbindlist(list(dt_temp,dt[na_rows,]),fill = TRUE)
    setkeyv(dt,key_dt)
  } else {
    dt = dt_temp
  }
  
  
  if(saveorgo)
  {
    save(dt, file = paste0(save_dir,file_name))
  }
  
  return(dt)
}





