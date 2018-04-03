
#############################################################################################

############## master script for the analysis of the northern atlantic ocean ################

#############################################################################################


rm(list = ls())

library(PostProcessing)
library(data.table)
library(parallel)

setwd("~/NR/SFE")
options(max.print = 1e3)

# box for analysis, name abbreviation for saved files and directories

lat_box = c(30,70)
lon_box = c(-60,15)

name_abbr = "NAO" 

save.dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr)
dir.create(save.dir, showWarnings = FALSE)

plot.dir = paste0("./figures/", name_abbr)
dir.create(plot.dir, showWarnings = FALSE)

### construct and load wide data set ###

# this one takes time, avoid if possible:
# make_combined_wide_dataset(lat_box = lat_box, lon_box = lon_box, output_loc = save.dir, output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))
                           
DT = load_combined_wide(data.dir = save.dir, output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))



########################################
####### univariate modelling ###########
########################################


############### bias ###################


### run bias analysis ###
 
validation_years = 2001:2010

# for simple moving averages

num.years = DT[,range(year)][2] - DT[,range(year)][1] + 1

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

sc_sma = mclapply(X = 1:(num.years-1), FUN = dummy_function, mc.cores = 8)
sc_sma = rbindlist(sc_sma)

save(sc_sma, file = paste0(save.dir,"/scores.bc.sma.Rdata"))
    
# for exponential moving averages

par.vec = seq(0.05,0.4,length.out = 24)
      
sc_ema = list()

dummy_function = function(k){
            temp = bias_correct(dt = DT, 
                                method = "ema",
                                par_1 = par.vec[k], 
                                scores = TRUE,
                                eval_years = validation_years,
                                saveorgo = FALSE)
      sc_ema[[k]] = temp[,"a" := par.vec[k]]
      }

sc_ema = mclapply(X = 1:length(par.vec), FUN = dummy_function,mc.cores = 8)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save.dir,"/scores.bc.ema.Rdata"))


#### plotting scores for different ways of bias correction ####

load(paste0(save.dir,"/scores.bc.sma.Rdata"))
load(paste0(save.dir,"/scores.bc.ema.Rdata"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma[,sqrt(MSE)],sc_ema[,sqrt(MSE)]))  

## plot for sma ##

pdf(paste0(plot.dir,"/mean_scores_sma.pdf"))
plot(x = sc_sma[,win_length],
     y = sc_sma[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for ",name_abbr," bias correction by SMA"),
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

pdf(paste0(plot.dir,"/mean_scores_ema.pdf"))
plot(x = sc_ema[,a],
     y = sc_ema[,sqrt(MSE)],
     ylim = y_range,
     type = "b",
     col = "blue",
     main = paste0("RMSE for ",name_abbr," bias correction by EMA"),
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

### finding optimal way of bias correcion ###

if(sc_sma[,min(MSE)] < sc_ema[,min(MSE)]){
  print(paste0("optimal bias correction uses simple moving averages with window length of ",sc_sma[,which.min(MSE)], " years, and achieves a RMSE of ",round(sc_sma[,sqrt(min(MSE))],3),
               ". Best correction with exponential moving averages achieves a RMSE of ",round(sc_ema[,sqrt(min(MSE))],3),"."))
  opt_par = c("sma",sc_sma[,which.min(MSE)])
  } else{
  print(paste0("optimal bias correction uses exponential moving averages with parameter a = ",round(sc_ema[,a][sc_sma[,which.min(MSE)]],3),
              ", and achieves a RMSE of ",round(sc_ema[,sqrt(min(MSE))],3),". Best correction with simple moving averages achieves a RMSE of ",round(sc_sma[,sqrt(min(MSE))],3),"."))
  opt_par = c("ema",sc_ema[,a][sc_sma[,which.min(MSE)]])
  }

### bias correction ###

bias_correct(dt = DT,
             method = opt_par[1],
             par_1 = as.double(opt_par[2]),
             save.dir = paste0(save.dir,"/"),
             file.name = paste0("dt_combine_",name_abbr,"_wide_bc.RData")
             )

DT = load_combined_wide(data.dir = save.dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"))



############### variance ###################
    
                    
### variance analysis ###

validation_years = 2001:2010


# using simple moving averages

num.years = DT[,range(year)][2] - DT[,range(year)][1] + 1

sc_sma_var = list()
dummy_function = function(k){
  temp = sd_est(dt = DT[year > min(year),], # for the first year the forecast is not bias corrected and (obs - fc)^2 is not a good approximation for the variance
                      method = "sma", 
                      par_1 = k,
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE)
  sc_sma_var[[k]] = temp[,"win_length" := k]
}

sc_sma_var = mclapply(X = 1:(num.years-1), FUN = dummy_function, mc.cores = 8)
sc_sma_var = rbindlist(sc_sma_var)

save(sc_sma_var, file = paste0(save.dir,"/scores.bc.sd.sma.Rdata"))


# for exponential moving averages

par.vec = seq(0.01,0.4,length.out = 24)

sc_ema_var = list()

dummy_function = function(k){
  temp = sd_est(dt = DT[year > min(year),], 
                      method = "ema",
                      par_1 = par.vec[k], 
                      scores = TRUE,
                      eval_years = validation_years,
                      saveorgo = FALSE)
  sc_ema_var[[k]] = temp[,"a" := par.vec[k]]
}

sc_ema_var = mclapply(X = 1:length(par.vec), FUN = dummy_function,mc.cores = 8)
sc_ema_var = rbindlist(sc_ema_var)

save(sc_ema_var, file = paste0(save.dir,"/scores.bc.sd.ema.Rdata"))


#### plotting scores for different ways of bias correction ####

load(paste0(save.dir,"/scores.bc.sd.sma.Rdata"))
load(paste0(save.dir,"/scores.bc.sd.ema.Rdata"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma_var[,CRPS],sc_ema_var[,CRPS]))  

## plot for sma ##

pdf(paste0(plot.dir,"/mean_scores_sd_sma.pdf"))
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

pdf(paste0(plot.dir,"/mean_scores_sd_ema.pdf"))
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

if(sc_sma_var[,min(CRPS)] < sc_ema_var[,min(CRPS)]){
  print(paste0("optimal variance estimation uses sample variance with window length of ",sc_sma_var[,which.min(CRPS)], " years, and achieves a CRPS of ",round(sc_sma_var[,min(CRPS)],3),
               ". Best estimation with exponentially weighted sample variance achieves a RMSE of ",round(sc_ema_var[,min(CRPS)],3),"."))
  opt_par = c("sma",sc_sma_var[,which.min(CRPS)])
} else{
  print(paste0("optimal variance estimation uses exponentially weighted sample variance with parameter a = ",round(sc_ema_var[,a][sc_ema_var[,which.min(CRPS)]],3),
               ", and achieves a CRPS of ",round(sc_ema_var[,min(CRPS)],3),". Best estimation with simple moving averages achieves a CRPS of ",round(sc_sma_var[,min(CRPS)],3),"."))
  opt_par = c("ema",sc_ema_var[,a][sc_sma_var[,which.min(CRPS)]])
}

### variance estimation ###

sd_est(dt = DT,
             method = opt_par[1],
             par_1 = as.double(opt_par[2]),
             save.dir = paste0(save.dir,"/"),
             file.name = paste0("dt_combine_",name_abbr,"_wide_bc_sd.RData")
)

DT = load_combined_wide(data.dir = save.dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc_sd.RData"))



#################################################
########## multivariate modelling ###############
#################################################


###################  PCA  #######################

training_years = 1985:2000
eval_years = 2001:2010

ens.size = 9

##### compute and save principal components #####

cov.dir = paste0(save.dir,"/PCACov")
dir.create(cov.dir, showWarnings = FALSE)

for_res_cov(Y = training_years,dt = DT, save.dir = cov.dir,ens.size = ens.size)

############## set up PCA #######################

setup_PCA(dt = DT, y = eval_years, cov.dir = cov.dir)

#### a brief example ####

ex_depth = 15
ex_month = 1
ex_year = 2001

result = forecast_PCA(y=ex_year,m=ex_month,PCA_depth = ex_depth, saveorgo = FALSE)

plot_diagnostic(result[,.(Lon,Lat,noise)], mn = paste0("Noise for ",ex_month,"/",ex_year," generated by ",d," PCs"))


########### compute variogram scores ############

dvec = c(1:50) # the number of principal components the variogram score is computed for
p = 0.5 # power for variogram score

months = 1:12
year = eval_years

## --------------------

var_sc = list() 

score_by_month = function(month){
  
  print(paste("starting month",month))
  
  # for computing pth moments we require the PCA data matrix, and we save it in form of a datatable:
  
  setup_PCA(dt = DT, m=month,cov.dir = cov.dir)
  PCA <- eval(parse(text = paste0("PCA",month)))
  PCA_DT = fc[year == min(year)][month == min(month),.(Lon,Lat,grid_id)]
  
  for(d in  dvec){
    PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
  }
  
  
  
  # we construct a data table that contains pairs of coordinates within the region specified 
  # by lat_box and lon_box defined in the beginning
  
  dummy_dt = fc[year == min(year) & month == min(month),  .(Lon,Lat,grid_id)]
  setnames(dummy_dt,c("Lon1","Lat1","grid_id1"))
  
  # add dummy key, do outer full merge with a duplicate, and remove key again:
  
  dummy_dt[,"key" := 1]
  dummy_2 = copy(dummy_dt) 
  setnames(dummy_2,c("Lon2","Lat2","grid_id2","key"))
  
  var_sc_prereq = merge(dummy_dt,dummy_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
  ##--- For each pair of coordinates (i,j), compute the variance of X_i-X_j 
  ##--- for the predictive distribution with d principal components
  
  var_d = function(d,coor1,coor2){
    i = match(coor1,PCA_DT[,grid_id])
    j = match(coor2,PCA_DT[,grid_id])
    charvec = paste0("PC",1:d)
    sq_diff_sum=0
    for(sum_ind in 1:d){
      sq_diff_sum = sq_diff_sum +(PCA_DT[,eval(parse(text = charvec[sum_ind]))][i] - 
                                    PCA_DT[,eval(parse(text = charvec[sum_ind]))][j]  )^2
    }
    return(sq_diff_sum)
  }
  
  
  for(d in dvec){
    var_sc_prereq[,(paste0("Var",d)) := var_d(d,grid_id1,grid_id2)]
    print( paste0("computing variances for ",d," principal components complete"))
  }
  
  
  # --- now get the differences of the residuals required for the variogram score
  
  obs_diff = function(y,coor1,coor2){
    i = match(coor1,DT[,grid_id])
    j = match(coor2,DT[,grid_id])
    
    abs_diff_obs = (abs(DT[year == y,SST_bar-SST_hat][i] - 
                          DT[y == year,SST_bar-SST_hat][j]  ))
    return(abs_diff_obs)
  }
  
  
  for(y in eval_years){
    var_sc_prereq[,(paste0("ObsDif",y)) := obs_diff(y,grid_id1,grid_id2)]
  }
  
  
  variogram_score = function(var_sc_prereq,
                             p ,
                             dvec,
                             eval_years){
    
    var_sc = as.data.table(expand.grid(dvec,eval_years))
    setnames(var_sc,c("PCs","year"))
    setkey(var_sc,PCs)
    
    # get the p-th moment of a standard normal distribution
    p_mom_sn = 2^(p/2)*gamma((p+1)/2)/sqrt(pi) 
    
    sc = c()
    for (d in dvec){
      print(paste0(d," principal components"))
      for(y in eval_years){
        print(paste0("year ",y))
        sc = c(sc,sum(((var_sc_prereq[,get(paste0("Var",d))])^(p/2) * p_mom_sn - 
                         (var_sc_prereq[,get(paste0("ObsDif",y))])^(p/2))^2))
      }
    }
    
    var_sc[,"v_sc" := sc]
    return(var_sc)
  }
  
  
  
  var_sc[[month]] = variogram_score(var_sc_prereq = var_sc_prereq,
                                    p = p,
                                    dvec = dvec,
                                    eval_years = eval_years)[,"month" := month]
  
  
  #print(paste0("month ",month," complete"))
  
}

var_sc = mclapply(X = months, FUN = score_by_month, mc.cores = 12)

var_sc = rbindlist(var_sc)

save(var_sc, file = paste0(save.dir,"/variogram_scores.RData"))


# --- plot variograms ---

load(file = paste0(save.dir,"/variogram_scores.RData"))

mean_sc = var_sc[,mean(v_sc), by = PCs]

pdf(paste0(plot.dir,"/mean_variogram_scores.pdf"))
plot(x = mean_sc[[1]],
     y = mean_sc[[2]],
     type = "b",
     col = "blue",
     main = paste0("mean variogram scores for ",name_abbr),
     xlab = "number of principal components",
     ylab = "mean score"
)
#---- add minima: -----
abline(h = min(mean_sc[[2]]), lty = "dashed", col = adjustcolor("blue",alpha = .5))

min_loc = mean_sc[,which.min(V1)]

points(x = mean_sc[[1]][min_loc],
       y = mean_sc[[2]][min_loc],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

############################################################




