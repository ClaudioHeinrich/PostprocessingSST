
#############################################################################################

###### master script part 1 - setting up and creating or loading of wide data table #########

#############################################################################################

# This script sets up for a full post-processing analysis 
# for a lon/lat window to be specified below

rm(list = ls())

time_s1 = proc.time()

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


# choose an abbreviation for this run

name_abbr = "NAO/lv"

# choose the area under consideration
lat_box = c(5,75)
lon_box = c(-100,50)


# # set a couple of parameters

ens_size = 9 # size of forecast ensemble

training_years = 1985:2000
validation_years = 2001:2016 

months = 1:12

mc_cores = 3

##### setting up complete - now move to creation #####

# create directories

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)

### construct load wide data set ###

DT = load_combined_wide()[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]

# tidy up DT:

setkey(x = DT,year,month,Lon,Lat)

DT = DT[month %in% months,]
DT[,YM := 12*year + month]
DT = DT[order(year,month,Lon,Lat)]

DT[,paste0('Ens',1:ens_size):=NULL]


setcolorder(DT,c("year","month",'Lon','Lat','YM','grid_id','SST_bar','Ens_bar','Ens_sd'))

# split DT by months

# for(m in months)
# {
#   assign(x = paste0('DT',m),value = DT[month == m])
# }
# 
# rm(DT)
# gc()

time_s1 = proc.time() - time_s1
# save everything:

save.image(file = paste0(save_dir,"setup.RData"))


# The following line needs to be run if the full data table has not yet been created, i.e. if new data is available:

# make_combined_wide_dataset(lat_box = lat_box,
#                            lon_box = lon_box,
#                            output_loc = save_dir,
#                            output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))





###### run bias analysis for simple moving averages ######

num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1

sc_sma = list()
dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "sma", 
                      par_1 = k,
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE)
  sc_sma[[k]] = temp[,"win_length" := k]
}

sc_sma = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = mc_cores)
sc_sma = rbindlist(sc_sma)

save(sc_sma, file = paste0(save_dir,"scores.bc.sma.RData"))


###### run bias analysis for exponential moving averages ######

par_vec = seq(0.05,0.4,length.out = 24) 

sc_ema = list()

dummy_function = function(k){
  temp = bias_correct(dt = DT, 
                      method = "ema",
                      par_1 = par_vec[k], 
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE)
  sc_ema[[k]] = temp[,"a" := par_vec[k]]
}

sc_ema = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = mc_cores)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save_dir,"scores.bc.ema.RData"))


###### plotting scores for different ways of bias correction ######

load(paste0(save_dir,"scores.bc.sma.RData"))
load(paste0(save_dir,"scores.bc.ema.RData"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma[,sqrt(MSE)],sc_ema[,sqrt(MSE)]))  

## plot for sma ##

pdf(paste0(plot_dir,"mean_scores_sma.pdf"))
plot(x = sc_sma[,win_length],
     y = sc_sma[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for bias correction by SMA"),
     xlab = "window length",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = sc_sma[,min(sqrt(MSE))], lty = "dashed", col = adjustcolor("blue",alpha = .5))
min_loc_RMSE = sc_sma[,which.min(MSE)]
points(x = sc_sma[,win_length][min_loc_RMSE],
       y = sc_sma[,sqrt(MSE)][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"mean_scores_ema.pdf"))
plot(x = sc_ema[,a],
     y = sc_ema[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for bias correction by EMA"),
     xlab = "weight parameter a",
     ylab = "RMSE"
)

# highlight minimum and add minimum reference line 
abline(h = sc_ema[,min(sqrt(MSE))], lty = "dashed", col = adjustcolor("blue",alpha = .5))
min_loc_RMSE = sc_ema[,which.min(MSE)]
points(x = sc_ema[,a][min_loc_RMSE],
       y = sc_ema[,sqrt(MSE)][min_loc_RMSE],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()


###### finding and applying optimal way of bias correcion, exponential moving averages are preferred, as the parameter estimation for them typically is more stable ######

if(sc_sma[,min(MSE)] < 0.9 * sc_ema[,min(MSE)]){
  opt_par = c("sma",sc_sma[,win_length][,which.min(MSE)])
} else{
  opt_par = c("ema",sc_ema[,a][sc_ema[,which.min(MSE)]])
}

DT = bias_correct(dt = DT,
                  method = opt_par[1],
                  par_1 = as.double(opt_par[2]),
                  save_dir = save_dir)

# For the training years the bias correction considers also the future

DT = bias_correct_training(dt = DT,
                           method = opt_par,
                           training_years = training_years,
                           save_dir = save_dir)




#### save ####

save.image(file = paste0(save_dir,"setup.RData"))


###### getting sample variances of ensemble  ######

DT = ens_sd_est(dt = DT,
                ens_size = ens_size,
                save_dir = save_dir,
                mean_est = "sv",
                file_name = paste0("dt_combine_wide_bc_var.RData"))


###### getting scores for simple moving averages ######

num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1

sc_sma_var = list()
dummy_function = function(k){
  temp = sd_est(dt = DT[year > min(year) ,], 
                method = "sma", 
                par_1 = k,
                scores = TRUE,
                eval_years = validation_years,
                saveorgo = FALSE)
  sc_sma_var[[k]] = temp[,"win_length" := k]
}

sc_sma_var = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = mc_cores)
sc_sma_var = rbindlist(sc_sma_var)

save(sc_sma_var, file = paste0(save_dir,"scores.bc.sd.sma.Rdata"))


###### getting scores for exponential moving averages ######

par_vec = seq(0.01,0.4,length.out = 24)

sc_ema_var = list()

dummy_function = function(k){
  temp = sd_est(dt = DT[year > min(year) ,], 
                method = "ema",
                par_1 = par_vec[k], 
                scores = TRUE,
                eval_years = validation_years,
                saveorgo = FALSE)
  sc_ema_var[[k]] = temp[,"a" := par_vec[k]]
}

sc_ema_var = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = 8)
sc_ema_var = rbindlist(sc_ema_var)

save(sc_ema_var, file = paste0(save_dir,"scores.bc.sd.ema.Rdata"))


#### plotting CRPS for different ways of variance estimation ####

load(paste0(save_dir,"/scores.bc.sd.sma.Rdata"))
load(paste0(save_dir,"/scores.bc.sd.ema.Rdata"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma_var[,CRPS],sc_ema_var[,CRPS]))  

## plot for sma ##

pdf(paste0(plot_dir,"/mean_scores_sd_sma.pdf"))
plot(x = sc_sma_var[,win_length],
     y = sc_sma_var[,CRPS],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("CRPS for ",name_abbr," SD estimation by SMA"),
     xlab = "window length",
     ylab = "variance score"
)

# highlight minimum and add minimum reference line 
abline(h = sc_sma_var[,min(CRPS)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
min_loc_CRPS = sc_sma_var[,which.min(CRPS)]
points(x = sc_sma_var[,win_length][min_loc_CRPS],
       y = sc_sma_var[,CRPS][min_loc_CRPS],
       col = "blue",
       bg = "blue",
       pch = 21)

dev.off()

## plot for ema ##

pdf(paste0(plot_dir,"/mean_scores_sd_ema.pdf"))
plot(x = sc_ema_var[,a],
     y = sc_ema_var[,CRPS],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("CRPS for ",name_abbr," SD estimation by EMA"),
     xlab = "weight parameter a",
     ylab = "variance score"
)

# highlight minimum and add minimum reference line 
abline(h = sc_ema_var[,min(CRPS)], lty = "dashed", col = adjustcolor("blue",alpha = .5))
min_loc_CRPS = sc_ema_var[,which.min(CRPS)]
points(x = sc_ema_var[,a][min_loc_CRPS],
       y = sc_ema_var[,CRPS][min_loc_CRPS],
       col = "blue",
       bg = "blue",
       pch = 21)

dev.off()

### finding optimal way of variance modelling ###

# if the optimal scores for simple moving averages and exponential moving acerages are essentially the same (less than 1 % difference),
# we pick exponential moving averages, since they are more stable.

if(sc_sma_var[,min(CRPS)] < sc_ema_var[,min(CRPS)]){
  opt_par_var = c("sma",sc_sma_var[,win_length][sc_sma_var[,which.min(CRPS)]])
} else{
  opt_par_var = c("ema",sc_ema_var[,a][sc_ema_var[,which.min(CRPS)]])
}

### variance estimation ###

DT = sd_est(dt = DT,
            method = opt_par_var[1],
            par_1 = as.double(opt_par_var[2]),
            save_dir = save_dir,
            file_name = paste0("dt_combine_wide_bc_var.RData"))



#####

save.image(file = paste0(save_dir,"setup.RData"))


##### centered variables? #####

# decide whether you want to work with SST, or with SST centered around climatology, or with SST standardized w.r.t. climatology 

SST = ""  # takes 'centered','standardized', or ''

clim_years = training_years


if(SST == "centered")
{
  DT = dt_transform_center(DT,clim_years)
  
  name_abbr = paste0(name_abbr,"/centered" )
  
  save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
  dir.create(save_dir,showWarnings = FALSE)
  
  plot_dir = paste0("./figures/", name_abbr,"/")
  dir.create(plot_dir, showWarnings = FALSE)
  
}
if(SST == "standardized")
{
  DT = dt_transform_stan(DT,clim_years)
  
  name_abbr = paste0(name_abbr,"/standardized")
  
  save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
  dir.create(save_dir,showWarnings = FALSE)
  
  plot_dir = paste0("./figures/", name_abbr,"/")
  dir.create(plot_dir, showWarnings = FALSE)
  
}


##### get weight matrix for tapering and svds of data matrices over the training period #####


Tap_dir = paste0(save_dir,"Tap/")
dir.create(Tap_dir,showWarnings = FALSE)

wm = weight_mat(DT,L = 2500) 

num_loc = DT[year == min(year) & month == min(month)][!(is.na(SST_bar) | is.na(SST_hat)) ,.N]


for(y in validation_years)
{
  
  print(y)
  
  svd_by_month = function(m)
  {
    print(paste0("month =",m))  
    
    train_years = DT[month == m][year < y & year > min(year),][,unique(year)]
    
    data_mat = matrix(DT[month == m][!(is.na(SST_bar) | is.na(SST_hat)) & year %in% train_years,
                                     SST_bar - SST_hat],
                      nrow = num_loc)
    
    sam_cov_mat = 1/length(train_years) * data_mat %*% t(data_mat) 
    
    tap_cov_mat = wm * sam_cov_mat
    
    return(svd(tap_cov_mat))
  }
  sin_val_dec = parallel::mclapply(X = months,FUN = svd_by_month,mc.cores = mc_cores)
  
  save(sin_val_dec,file = paste0(Tap_dir,'svd_y',y,'.RData'))
}




###################################################
###################### PCA  #######################
###################################################

PCA_dir = paste0(save_dir,"PCA/")
dir.create(PCA_dir,showWarnings = FALSE)


# how many PCs should we use?

nPCs = c()

for(m in months)
{
  #plot(sin_val_dec[[m]]$d, main = paste0('month = ',m))
  
  sum_vec = cumsum(sin_val_dec[[m]]$d)
  sum_tot = sum_vec[length(sum_vec)]
  
  nPCs = c(nPCs,which(sum_vec > 0.9*sum_tot)[1])
  
}


PCA_cov(DT,weight_mat = wm, 
        Y = 1985:2016,
        M = months,
        nPCs = nPCs,
        save_years = validation_years,
        save_dir = PCA_dir)

##################################################
###################### SE  #######################
##################################################
#
# SE_dir = paste0(save_dir,"SE/")
# dir.create(SE_dir,showWarnings = FALSE)
# 
# # get covariance estimates
# 
# cov_est_SE(DT, weight_mat = wm,
#            M = months,
#            save_years = validation_years,
#            save_dir = SE_dir)
#
###################################################
################## geostationary ##################
###################################################

GS_dir = paste0(save_dir, "GS/")
dir.create(GS_dir, showWarnings = FALSE)

geostationary_training(dt = DT, 
                       training_years = training_years,
                       m = months,
                       save_dir = GS_dir,mc_cores)

########################################
################ ECC  ##################
########################################

ECC_dir = paste0(save_dir, "ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

#####

time_s4 = proc.time() - time_s4

save.image(file = paste0(save_dir,"setup.RData"))



