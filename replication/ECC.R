
############################################################################
############### Creates a postprocessed forecast by ECC ####################
############################################################################

rm(list = ls())

library(PostProcessing)

setwd("~/NR/SFE")
options(max.print = 500)




#' ECC forecasts 
#' 
#' @description Generates forecasts post-processed by ensemble copula coupling (ECC)
#' 
#' @param dt the data table containing bias and variance estimates, if NULL, it is loaded.
#' @param ens_size Integer. Size of the NWP ensemble.
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



forecast_ECC = function(dt = NULL,
                        ens_size = 9,
                        saveorgo = TRUE,
                        save_dir = "./Data/PostClim/SFE/Derived/ECC/",
                        file_name = "ECC_fc.RData"
                        ){
  
if(is.null(dt))
{
  load_combined_wide(var = TRUE)
}
  
# generate noise and marginally post-process each ensemble member
  
na_cols = DT[ ,which(is.na(Bias_Est) |  is.na(Ens_bar) | is.na(var_bar) )]

length_norm = DT[ -na_cols,.N]

for(i in 1:ens_size)
  {
  norm_rv = rnorm(n = length_norm, mean = DT[-na_cols,eval(parse(text = paste0("Ens",i))) + Bias_Est], sd = DT[-na_cols,SD_hat])
  DT[ -na_cols,paste0("fc",i) := norm_rv]
  }

#get rank order of the ensemble and reorder the post-processed ensemble to mach the rank order statistic of the ensemble


rks_ens = DT[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("Ens",1:9)]

rks_fc = DT[,matrixStats::rowRanks(as.matrix(.SD)),.SDcols = paste0("fc",1:9)]

fcs = as.matrix(DT[,.SD,.SDcols = paste0("fc",1:ens_size)])

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

DT = data.table(DT,ecc_fcs)

if(saveorgo)
{
  save(DT, file = paste0(save_dir,file_name))
}

return(DT)

}




for(i in 1:ens_size){
  name_1 = paste0("ecc_fc",i)
  plot_diagnostic(DT[year == 2010 & month == 2,.SD,.SDcols = c("Lon","Lat",name_1)],mn = name_1)
  name_2 = paste0("fc",i)
  plot_diagnostic(DT[year == 2010 & month == 2,.SD,.SDcols = c("Lon","Lat",name_2)],mn = name_2)
}


