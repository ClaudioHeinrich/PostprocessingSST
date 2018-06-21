

#############################################################################################

###### BIG test script for a bunch of models #########

#############################################################################################

# This script tests variants of the PCA method on three different 

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

# choose your favourite area for analysis and give it a name abbreviation

#NAO_2:

 lat_box = c(35,75)
 lon_box = c(-70,40)

name_abbr = "NAO_large" # for northern atlantic ocean

ens_size = 9 # size of forecast ensemble

validation_years = 2001:2010 # all previous years are used for training 
months = 1:12

# create directories

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr,"/")
dir.create(plot_dir, showWarnings = FALSE)

### construct or load wide data set ###

DT = load_combined_wide()[Lon >= lon_box[1] & Lon <= lon_box[2] & Lat >= lat_box[1] & Lat <= lat_box[2]]


training_year_index = !(DT[,unique(year)] %in% validation_years) 
training_years = DT[,unique(year)][training_year_index]



##########################################################
###### run bias analysis for simple moving averages ######
##########################################################


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

sc_sma = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = 8)
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

sc_ema = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = 8)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save_dir,"scores.bc.ema.RData"))

###### plotting scores for different ways of bias correction ######

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





###### finding and applying optimal way of bias correcion ######

if(sc_sma[,min(MSE)] < sc_ema[,min(MSE)]){
  opt_par = c("sma",sc_sma[,which.min(MSE)])
} else{
  opt_par = c("ema",sc_ema[,a][sc_sma[,which.min(MSE)]])
}

DT = bias_correct(dt = DT,
                  method = opt_par[1],
                  par_1 = as.double(opt_par[2]),
                  save_dir = save_dir)


DT = bias_correct_training(dt = DT,
                           method = opt_par,
                           training_years = training_years,
                           save_dir = save_dir)

###### compare to linear regression models ######

RMSE_linear_models = bias_lr(DT,validation_years = validation_years)

RMSE_linear_models



#######################################################################################

variance_methods = c("sv","bcf")


for(i in 1:length(variance_methods))
{
  
  print(paste0("variance estimation method: ",variance_methods[i]))
  
  DT = ens_sd_est(dt = DT,
                  ens_size = ens_size,
                  save_dir = save_dir,
                  mean_est = variance_methods[i],
                  file_name = paste0("dt_combine_wide_bc_var.RData"))
  
  
  ###### getting scores for simple moving averages ######
  
  num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1
  
  sc_sma_var = list()
  dummy_function = function(k){
    temp = sd_est(dt = DT, 
                  method = "sma", 
                  par_1 = k,
                  scores = TRUE,
                  eval_years = validation_years,
                  saveorgo = FALSE)
    sc_sma_var[[k]] = temp[,"win_length" := k]
  }
  
  sc_sma_var = parallel::mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = 8)
  sc_sma_var = rbindlist(sc_sma_var)
  
    ###### getting scores for exponential moving averages ######
  
  par_vec = seq(0.01,0.4,length.out = 24)
  
  sc_ema_var = list()
  
  dummy_function = function(k){
    temp = sd_est(dt = DT[year > min(year),], 
                  method = "ema",
                  par_1 = par_vec[k], 
                  scores = TRUE,
                  eval_years = validation_years,
                  saveorgo = FALSE)
    sc_ema_var[[k]] = temp[,"a" := par_vec[k]]
  }
  
  sc_ema_var = parallel::mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = 8)
  sc_ema_var = rbindlist(sc_ema_var)
  
  
  
  ### finding optimal way of variance modelling ###
  
  if(sc_sma_var[,min(CRPS)] < sc_ema_var[,min(CRPS)]){
    opt_par = c("sma",sc_sma_var[,which.min(CRPS)])
  } else{
    opt_par = c("ema",sc_ema_var[,a][sc_sma_var[,which.min(CRPS)]])
  }
  
  DT = sd_est(dt = DT,
              method = opt_par[1],
              par_1 = as.double(opt_par[2]),
              save_dir =save_dir,
              file_name = paste0("dt_combine_wide_bc_var.RData")
  )
  
  setnames(DT,"SD_hat",paste0("SD_hat_",variance_methods[i]))
  
  ############# saving #####################
  
  save(sc_sma_var,sc_ema_var,opt_par, file = paste0(save_dir,"scores.bc.sd.",variance_methods[i],".Rdata"))

}



###### Compare to SD estimation in NGR fashion: #####

DT = var_est_NGR(DT, months = months, validation_years = validation_years)

sc_sma = data.table(win_length = 1:(num_years-1))
sc_ema = data.table(a = par_vec)

for(i in 1:length(variance_methods))
{
  var_met = variance_methods[i]
  
  load(file = paste0(save_dir,"scores.bc.sd.",var_met,".Rdata"))
  sc_sma = merge(sc_sma,sc_sma_var)
  setnames(sc_sma,"CRPS",paste0("sma_CRPS_",var_met))
  
  sc_ema = merge(sc_ema,sc_ema_var)
  setnames(sc_ema,"CRPS",paste0("ema_CRPS_",var_met))
  
}



CRPS_comparison = data.table(mean_CRPS_bm = DT[,mean(crps.na.rm(SST_bar,SST_hat,SD_hat_lr_bm),na.rm = TRUE)],
                             mean_CRPS_bl = DT[,mean(crps.na.rm(SST_bar,SST_hat,SD_hat_lr_bl),na.rm = TRUE)],
                             mean_CRPS_bb = DT[,mean(crps.na.rm(SST_bar,SST_hat,SD_hat_lr_bb),na.rm = TRUE)],
                             sc_ema[,lapply(.SD,min),.SDcols = 2:ncol(sc_ema)],
                             sc_sma[,lapply(.SD,min),.SDcols = 2:ncol(sc_ema)])

round(CRPS_comparison,3)


##### CRPS plots for different versions of marginal variance estimation: #####

col_vec = c("blue","darkgreen","darkred","black","yellow2")
leg_vec = variance_methods

y_range = range(sc_sma[,2:( 1+length(variance_methods))], sc_ema[,2:( 1+length(variance_methods))])  


# sma plot

pdf(paste0(plot_dir,"/mean_scores_sd_sma.pdf"))
  plot(sc_sma[,.SD,.SDcols = c("win_length",paste0("CRPS_",variance_methods[1]))],
       ylim = y_range,
       type = "l",
       col = col_vec[1],
       main = paste0("CRPS for ",name_abbr," SD estimation by SMA"),
       xlab = "window length",
       ylab = "variance score"
  )
  
  # highlight minimum 
  min_loc_CRPS = which.min(sc_sma[,eval(parse(text = paste0("CRPS_",variance_methods[1])))])
  points(x = sc_sma[,win_length][min_loc_CRPS],
         y = sc_sma[,2][min_loc_CRPS],
         col = col_vec[1],
         bg = col_vec[1],
         pch = 21)
  
  for(i in 2:length(variance_methods))
  {
    lines(sc_sma[,.SD,.SDcols = c("win_length",paste0("CRPS_",variance_methods[i]))],
          type = "l",col = col_vec[i]) 
    
    min_loc_CRPS = which.min(sc_sma[,eval(parse(text = paste0("CRPS_",variance_methods[i])))])
    points(x = sc_sma[,win_length][min_loc_CRPS],
           y = sc_sma[,1+i, with = FALSE][min_loc_CRPS],
           col = col_vec[i],
           bg = col_vec[i],
           pch = 21)
  }
  
  # highlight minimum and add minimum reference line 
  abline(h = min(sc_sma[,-1, with = FALSE]), lty = "dashed", col = adjustcolor("black",alpha = .5))
  
  legend("topright",legend = leg_vec,lty = 1,col = col_vec[1:length(variance_methods)])
  
  dev.off()
  

# ema plot
  
  
  pdf(paste0(plot_dir,"/mean_scores_sd_ema.pdf"))
  plot(sc_ema[,.SD,.SDcols = c("a",paste0("CRPS_",variance_methods[1]))],
       ylim = y_range,
       type = "l",
       col = col_vec[1],
       main = paste0("CRPS for ",name_abbr," SD estimation by EMA"),
       xlab = "a",
       ylab = "CRPS"
  )
  
  # highlight minimum 
  min_loc_CRPS = which.min(sc_ema[,eval(parse(text = paste0("CRPS_",variance_methods[1])))])
  points(x = sc_ema[,a][min_loc_CRPS],
         y = sc_ema[,2][min_loc_CRPS],
         col = col_vec[1],
         bg = col_vec[1],
         pch = 21)
  
  for(i in 2:length(variance_methods))
  {
    lines(sc_ema[,.SD,.SDcols = c("a",paste0("CRPS_",variance_methods[i]))],
          type = "l",col = col_vec[i]) 
    
    min_loc_CRPS = which.min(sc_ema[,eval(parse(text = paste0("CRPS_",variance_methods[i])))])
    points(x = sc_ema[,a][min_loc_CRPS],
           y = sc_ema[,1+i, with = FALSE][min_loc_CRPS],
           col = col_vec[i],
           bg = col_vec[i],
           pch = 21)
  }
  
  # highlight minimum and add minimum reference line 
  abline(h = min(sc_ema[,-1, with = FALSE]), lty = "dashed", col = adjustcolor("black",alpha = .5))
  
  legend("topright",legend = leg_vec,lty = 1,col = col_vec[1:length(variance_methods)])
  
  dev.off()
  
  
#########################################################
######## univariate rank histograms #####################
#########################################################

methods_sdnames = c(paste0("SD_hat_",variance_methods),paste0("SD_hat_lr_b",c("m","l","b")))
method_names = c(paste0("WMA w.r.t. ",variance_methods),paste0("NGR est. aggr. by ",c("month","location","both")))

na_loc = which(DT[,is.na(SST_bar) | is.na(Ens_bar)])  

n=5000
nbreaks = 20

for(i in 1: length(methods_sdnames))
{
  m_sdn = methods_sdnames[i]
  print(method_names[i])
  
  means = DT[-na_loc,][year %in% validation_years & month %in% months,SST_hat]
  sds = DT[-na_loc,][year %in% validation_years & month %in% months,eval(parse(text = m_sdn))]
  obs = DT[-na_loc,][year %in% validation_years & month %in% months,SST_bar]
  
  fcs = rnorm(n*length(means),mean = means,sd = sds)

  rh_mat = matrix(c(obs,fcs),ncol = n+1)
  
  rank_mat = matrixStats::rowRanks(rh_mat)
  
  pdf(file=paste0(plot_dir,"rkh_",m_sdn,".pdf"))
  hist(rank_mat[,1],
       breaks=seq(0, n+1, length.out = nbreaks),
       main = method_names[i],
       xlab = "", ylab = "", axes=FALSE, col="gray80", border="gray60")
  abline(a=dim(rank_mat)[1]/(nbreaks-1), b=0, lty=2, col="gray30")
  dev.off()
}

DT[,SD_hat := SD_hat_sv]

save.image(file = paste0(save_dir,"univariate.model.comparison.RData"))




##############################################################

############## Multivariate Methods ##########################

##############################################################

name_abbr = "NAO_2" # for northern atlantic ocean

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"univariate.model.comparison.RData"))


##########################################
##### testing different PCA versions #####
##########################################
  
versions = c("aggr_by_season","sum_of_squares","wrt_ens_mean")
version_titles = c("aggr. by season","w.r.t. ensemble members","w.r.t. ensemble mean")
version_abbrs = c("abs","ssq","scm")
  
for(ver_ind in 1:length(versions))
{
  version = versions[ver_ind]
  version_title = version_titles[ver_ind]
  version_abbr = version_abbrs[ver_ind]
 
  print(paste0("version = ",version))
   
  PCA_dir = paste0(save_dir,"PCA/",version_abbr,"/")
  dir.create(PCA_dir, showWarnings = FALSE)
  
  training_years = DT[!(year %in% validation_years),unique(year)]
  
  for_res_cov(Y = training_years,
              dt = DT, 
              save_dir = PCA_dir,
              ens_size = ens_size,
              version = version)
  
  
  
  ### range of PCs to test ###
  
  
  if(version == "aggr_by_season")
  {
    PCs = c(1:10,15,20)  
  }
  if(version == "sum_of_squares")
  {
    PCs = c(1:10,20,50)  
  }
  if(version == "wrt_ens_mean")
  {
    PCs = c(1:10)  
  }
  
  #### variogram scores computation: ####
  
  
  for(m in months)
  { 
    # setup: get principal components and marginal variances for the given month m:
    
    print(paste0("m = ",m))
    load(file = paste0(PCA_dir,"CovRes_mon",m,".RData"))
    PCA = irlba::irlba(res_cov, nv = max(PCs),fastpath = FALSE)
    
    land_ids = which(DT[year == min(year) & month == min(month),is.na(Ens_bar) | is.na(SST_bar)])
    
    PCA_DT = DT[year == min(year) & month == min(month),][-land_ids,.(Lon,Lat,grid_id)]
    
    for(d in  min(PCs):max(PCs))
    {
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
    } 
    
    # also get marginal variances
    variances = list()
    d = min(PCs) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = vec^2
    
    if(length(PCs)>1)
    {
      for(d in (min(PCs)+1):max(PCs)){
        vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
        variances[[d]] = variances[[d-1]] + vec^2
      }
    }
    names(variances) = paste0("var",min(PCs):max(PCs))
    PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))    
    
    # now move to getting the variogram scores:
    
    # without marginal correction:
    
    dummy_fct = function(y)
    {
      var_sc_PCA_old(m, y, DT, PCA = PCA, PCA_DT = PCA_DT, dvec = PCs, ens_size = ens_size,
                     weighted = FALSE,
                     file_name = paste0("var_sc_by_PC", weight_name_addition),
                     marginal_correction = FALSE, 
                     cov_dir = PCA_dir, save_dir = PCA_dir)
    }
    parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
    
    # with marginal correction:
    
    dummy_fct = function(y)
    {
      print(c(paste0("y = ",y),paste0("m = ",m)))
      var_sc_PCA_old(m, y, DT, PCA = PCA, PCA_DT = PCA_DT, dvec = PCs, ens_size = ens_size,
                     weighted = FALSE,
                     file_name = paste0("var_sc_by_PC", weight_name_addition),
                     cov_dir = PCA_dir, save_dir = PCA_dir)
    }
    parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
  }
  
  
  
  ###### combine: ######
  
  scores_nmc = list()
  k=0
  for(m in months){
    for(y in validation_years)
    {k=k+1
    load(file = paste0(PCA_dir,"var_sc_by_PC",weight_name_addition,"_nmc_m",m,"_y",y,".RData"))
    scores_nmc[[k]] = scores
    }
  }
  scores_nmc = rbindlist(scores_nmc)
  
  scores_mc = list()
  k=0
  for(m in months){
    for(y in validation_years)
    {k=k+1
    load(file = paste0(PCA_dir,"var_sc_by_PC",weight_name_addition,"_m",m,"_y",y,".RData"))
    scores_mc[[k]] = scores
    }
  }
  scores_mc = rbindlist(scores_mc)
  
  mean_sc = scores_mc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
  mean_sc = unique(mean_sc)
  
  
  mean_sc_nmc = scores_nmc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
  mean_sc_nmc = unique(mean_sc_nmc)
  
  # save
  
  save(scores_mc,scores_nmc,mean_sc,mean_sc_nmc,file = paste0(PCA_dir,"variogram_scores.RData"))
  
}

##### geostationary #####

geostat_dir = paste0(save_dir, "GeoStat/")
dir.create(geostat_dir, showWarnings = FALSE)

geostationary_training(dt = DT, save_dir = geostat_dir,m = months)


for(m in months)
{ 
  # setup: get principal components and marginal variances for the given month m:
  
  print(paste0("m = ",m))
  
  load(paste0(geostat_dir,"variogram_exp_m",m,".RData")) 
  
  
  
  # variogram scores without marginal correction:
  
  dummy_fct = function(y)
  {
    var_sc_geoStat_new( DT, m, y, Mod = Mod, Dist = Dist,
                        weighted = weight_a_minute, weight_fct = weight_fct, 
                        file_name = paste0("var_sc",weight_name_addition),
                        save_dir = geostat_dir, data_dir = geostat_dir,
                        mar_var_cor = FALSE)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
  
  # with marginal correction:
  
  dummy_fct = function(y)
  {
    var_sc_geoStat_new( DT,m, y, Mod = Mod, Dist = Dist,
                        weighted = weight_a_minute, weight_fct = weight_fct, 
                        file_name = paste0("var_sc",weight_name_addition),
                        save_dir = geostat_dir, data_dir = geostat_dir,
                        mar_var_cor = TRUE)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
}


####### combine #########

scores_geostat_nmc = list()
k=0
for(m in months){
  for(y in validation_years)
  {k=k+1
  load(file = paste0(geostat_dir,"var_sc",weight_name_addition,"_nmc_m",m,"_y",y,".RData"))
  scores_geostat_nmc[[k]] = scores
  }
}
scores_geostat_nmc = rbindlist(scores_geostat_nmc)

scores_geostat_mc = list()
k=0
for(m in months){
  for(y in validation_years)
  {k=k+1
  load(file = paste0(geostat_dir,"var_sc",weight_name_addition,"_m",m,"_y",y,".RData"))
  scores_geostat_mc[[k]] = scores
  }
}
scores_geostat_mc = rbindlist(scores_geostat_mc)


mean_geostat_sc = scores_geostat_mc[,mean_sc := mean(sc)][,.(mean_sc)]
mean_geostat_sc = unique(mean_geostat_sc)


mean_geostat_sc_nmc = scores_geostat_nmc[,mean_sc := mean(sc)][,.(mean_sc)]
mean_geostat_sc_nmc = unique(mean_geostat_sc_nmc)

# save

save(scores_geostat_mc,scores_geostat_nmc,mean_geostat_sc,mean_geostat_sc_nmc,file = paste0(geostat_dir,"variogram_scores",weight_name_addition,".RData"))



### plotting and comparison ###

for(ver_ind in 1:length(versions))
{
  version = versions[ver_ind]
  version_title = version_titles[ver_ind]
  version_abbr = version_abbrs[ver_ind]

  PCA_dir = paste0(save_dir,"PCA/",version_abbr,"/")
  
  load(file = paste0(PCA_dir,"variogram_scores.RData"))
  
  assign(paste0("mean_sc_",version_abbr),mean_sc)
  assign(paste0("mean_sc_nmc_",version_abbr),mean_sc_nmc)
}


##################################
# plot with marginal correction: #
##################################

pdf(paste0(plot_dir,"/mean_variogram_scores.pdf"))
rr = range(c(mean_sc_abs[[2]],mean_sc_ssq[[2]],mean_sc_scm[[2]],mean_geostat_sc))

plot(x = mean_sc_ssq[[1]],
     y = mean_sc_ssq[[2]],
     ylim = rr,
     type = "b",
     col = "blue",
     main = paste0("mean variogram scores"),
     xlab = "number of principal components",
     ylab = "mean score")

#---- add minima: -----


lines(x = mean_sc_abs[[1]],
      y = mean_sc_abs[[2]],
      type = "b",
      col = "darkred")

lines(x = mean_sc_scm[[1]],
      y = mean_sc_scm[[2]],
      type = "b",
      col = "darkgreen")


# --- add geostat value: ----

abline(h = mean_geostat_sc[,mean_sc], lty = "dashed", col = adjustcolor("black"))

#abline(h = sc_ECC[,mean(sc)], lty = "dashed", col = adjustcolor("pink"))


legend("topright",legend = c("w.r.t. ens. mem.","aggr. by season","w.r.t ens. mean","geostat"),
       col = c("blue","darkgreen","darkred","black"),lty = c(1,1,1,2))
# legend("topright",legend = c("PCA, m.c.v.","PCA","geostat, m.c.v.","geostat","ECC"),
#        col = c("blue","darkgreen","darkred","black","pink"),lty = c(1,1,2,2,2))
dev.off()



#####################################
# plot without marginal correction: #
#####################################

pdf(paste0(plot_dir,"/mean_variogram_scores_nmc.pdf"))
rr = range(c(mean_sc_nmc_abs[[2]],mean_sc_nmc_ssq[[2]],mean_sc_nmc_scm[[2]],mean_geostat_sc_nmc,mean_geostat_sc))

plot(x = mean_sc_nmc_ssq[[1]],
     y = mean_sc_nmc_ssq[[2]],
     ylim = rr,
     type = "b",
     col = "blue",
     main = paste0("mean variogram scores without mar. cor."),
     xlab = "number of principal components",
     ylab = "mean score")

#---- add minima: -----


lines(x = mean_sc_nmc_abs[[1]],
      y = mean_sc_nmc_abs[[2]],
      type = "b",
      col = "darkred")

lines(x = mean_sc_nmc_scm[[1]],
      y = mean_sc_nmc_scm[[2]],
      type = "b",
      col = "darkgreen")


# --- add geostat value: ----

abline(h = mean_geostat_sc_nmc[,mean_sc], lty = "dashed", col = adjustcolor("black"))
abline(h = mean_geostat_sc[,mean_sc], lty = "dashed", col = adjustcolor("black"))

#abline(h = sc_ECC[,mean(sc)], lty = "dashed", col = adjustcolor("pink"))


legend("topright",legend = c("w.r.t. ens. mem.","aggr. by season","w.r.t ens. mean","geostat"),
       col = c("blue","darkgreen","darkred","black"),lty = c(1,1,1,2))
# legend("topright",legend = c("PCA, m.c.v.","PCA","geostat, m.c.v.","geostat","ECC"),
#        col = c("blue","darkgreen","darkred","black","pink"),lty = c(1,1,2,2,2))
dev.off()




