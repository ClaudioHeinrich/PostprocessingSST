
#############################################################################################

############## master script for the analysis of the northern atlantic ocean ################

#############################################################################################


rm(list = ls())

library(PostProcessing)
library(data.table)
library(parallel)
library(fExtremes)
library(fields)
library(vegan)

setwd("~/NR/SFE")
options(max.print = 1e3)

# box for analysis, name abbreviation for saved files and directories

lat_box = c(30,70)
lon_box = c(-60,15)

name_abbr = "NAO" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr)
dir.create(save_dir, showWarnings = FALSE)

plot_dir = paste0("./figures/", name_abbr)
dir.create(plot_dir, showWarnings = FALSE)

### construct and load wide data set ###

# this one takes time, avoid if possible:
# make_combined_wide_dataset(lat_box = lat_box, lon_box = lon_box, output_loc = save_dir, output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))
                           
DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_",name_abbr,"_wide.RData"))



########################################
####### univariate modelling ###########
########################################


############### bias ###################


### run bias analysis ###
 
validation_years = 2001:2010

# for simple moving averages

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

sc_sma = mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = 8)
sc_sma = rbindlist(sc_sma)

save(sc_sma, file = paste0(save_dir,"/scores.bc.sma.Rdata"))
    
# for exponential moving averages

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

sc_ema = mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = 8)
sc_ema = rbindlist(sc_ema)

save(sc_ema, file = paste0(save_dir,"/scores.bc.ema.Rdata"))


#### plotting scores for different ways of bias correction ####

load(paste0(save_dir,"/scores.bc.sma.Rdata"))
load(paste0(save_dir,"/scores.bc.ema.Rdata"))

# ensure that they are plotted on the the same range
y_range = range(c(sc_sma[,sqrt(MSE)],sc_ema[,sqrt(MSE)]))  

## plot for sma ##

pdf(paste0(plot_dir,"/mean_scores_sma.pdf"))
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

pdf(paste0(plot_dir,"/mean_scores_ema.pdf"))
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
             save_dir = paste0(save_dir,"/"),
             file_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData")
             )

DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc.RData"))



############### variance ###################
    
                    
### variance analysis  - incomplete ###

validation_years = 2001:2010


# using simple moving averages

num_years = DT[,range(year)][2] - DT[,range(year)][1] + 1

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

sc_sma_var = mclapply(X = 1:(num_years-1), FUN = dummy_function, mc.cores = 8)
sc_sma_var = rbindlist(sc_sma_var)

save(sc_sma_var, file = paste0(save_dir,"/scores.bc.sd.sma.Rdata"))


# for exponential moving averages

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

sc_ema_var = mclapply(X = 1:length(par_vec), FUN = dummy_function,mc.cores = 8)
sc_ema_var = rbindlist(sc_ema_var)

save(sc_ema_var, file = paste0(save_dir,"/scores.bc.sd.ema.Rdata"))


#### plotting scores for different ways of bias correction ####

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
             save_dir = paste0(save_dir,"/"),
             file_name = paste0("dt_combine_",name_abbr,"_wide_bc_sd.RData")
)

DT = load_combined_wide(data_dir = save_dir, output_name = paste0("dt_combine_",name_abbr,"_wide_bc_sd.RData"))



#################################################
########## multivariate modelling ###############
#################################################


###################  PCA  #######################

training_years = 1985:2000
eval_years = 2001:2010

ens_size = 9

##### compute and save principal components #####

cov_dir = paste0(save_dir,"/PCACov")
dir.create(cov_dir, showWarnings = FALSE)

for_res_cov(Y = training_years,dt = DT, save_dir = cov_dir,ens_size = ens_size)

############## set up PCA #######################

setup_PCA(dt = DT, y = eval_years, cov_dir = cov_dir)

#### a brief example ####

ex_depth = 15
ex_month = 1
ex_year = 2001

result = forecast_PCA(y=ex_year,m=ex_month,PCA_depth = ex_depth, saveorgo = FALSE)

plot_diagnostic(result[,.(Lon,Lat,noise)], mn = paste0("Noise for ",ex_month,"/",ex_year," generated by ",ex_depth," PCs"))


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
  
  setup_PCA(dt = DT, m=month,cov_dir = cov_dir)
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

save(var_sc, file = paste0(save_dir,"/variogram_scores.RData"))


# --- plot variograms ---

load(file = paste0(save_dir,"/variogram_scores.RData"))

mean_sc = var_sc[,mean(v_sc), by = PCs]

pdf(paste0(plot_dir,"/mean_variogram_scores.pdf"))
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

opt_num_PCs = mean_sc[,which.min(V1)]

points(x = mean_sc[[1]][opt_num_PCs],
       y = mean_sc[[2]][opt_num_PCs],
       col = "blue",
       bg = "blue",
       pch = 21)
dev.off()

############################################################


###### multivariate rank histograms ######

## Minimum spanning tree ranks 
mst.rank <- function (x) {
  l.mst <- NULL
  for(f in 1:(dim(x)[2])) {
    euc.dist <- rdist(t(x[,-f]))
    l.mst <- c(l.mst,sum(spantree(euc.dist)$dist))
  }
  x.rank <- rank(l.mst,ties="random")
  return(x.rank)
}

## Multivariate ranks 
mv.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- numeric(d[2])
  for(i in 1:d[2]) {
    x.prerank[i] <- sum(apply(x<=x[,i],2,all))
  }
  x.rank <- rank(x.prerank,ties="random")
  return(x.rank)
}

## Average ranks
avg.rank <- function(x)  {
  x.ranks <- apply(x,1,rank)
  x.preranks <- apply(x.ranks,1,mean)
  x.rank <- rank(x.preranks,ties="random")
  return(x.rank)
}

## Band depth ranks
bd.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- array(NA,dim=d)
  for(i in 1:d[1]) {
    tmp.ranks <- rank(x[i,])
    x.prerank[i,] <- (d[2] - tmp.ranks) * (tmp.ranks - 1)
  }
  x.rank <- apply(x.prerank,2,mean) + d[2] - 1
  x.rank <- rank(x.rank,ties="random")
  return(x.rank)
} 


## rank histograms for data tables

# The data table should have the key variable (most commonly YM) as first column and the ranks of the observations as second

rhist.dt <- function(B, ens_size, breaks = seq(0, ens_size + 1, length.out = min(ens_size + 1,20)), hist_xlab="", hist_ylab="", hist_ylim=NULL)
{
  hist(as.vector(B[[2]]),breaks = breaks, main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray80",border="gray60",ylim=hist_ylim)
  abline(a=length(B[[1]])/length(breaks), b=0, lty=2, col="gray30")
}


##################################################
#### compute rank histograms for PCA forecast ####
##################################################




validation_years = 2001:2010  # validation years
validation_months = 1:12       # validation months
MC_sample_size = 100   # how often we generate PCA noise
ens_size = 9   #size of forecast ensemble

PCvec = c(5,opt_num_PCs,50)      # number of considered principal components

for(PCs in PCvec){

  print(paste0("computing RHs for ",PCs," principal components."))
  no_dt = list()
  for(i in 1:MC_sample_size){
    print(paste0("generating Monte Carlo sample ",i,"/",MC_sample_size))
    no_dt[[i]] = forecast_PCA(m = validation_months, y = validation_years, PCA_depth = PCs, saveorgo = FALSE)[,noise]
  }
  no_dt = as.data.table(no_dt)
  
  setnames(no_dt,paste0("no",1:MC_sample_size))

  DT_pca = no_dt[,c("year","month","YM","SST_bar","Bias_Est",paste0("Ens",1:ens_size)) := DT[year %in% validation_years & month %in% validation_months,c("year","month","YM","SST_bar","Bias_Est",paste0("Ens",1:ens_size)),with = FALSE]]


  #choose random ensemble members (REM) and generate forecast as REM + bias + noise

  ens_mem = sample.int(ens_size,MC_sample_size,replace = TRUE)
  for(i in 1:MC_sample_size){
    dummy_dt = DT_pca[,.SD,.SDcols = c(paste0("no",i),paste0("Ens",ens_mem[i]),"Bias_Est")]
    forecast = dummy_dt[[1]] + dummy_dt[[2]] + dummy_dt[[3]]
    DT_pca = DT_pca[,paste0("fc",i) := forecast]
  }
  
  ym = unique(DT_pca[,YM])
  
  ranks.matrix = matrix(ym,nrow = length(ym),ncol = 1+3*(MC_sample_size +1)) 
  #ncol: 1 col for YM, 3 methods of ranking, for each we get ranks for observation and each Monte Carlo sample
  
  drm = dim(ranks.matrix)
  
  YM_ind = 0
  
  for(yearmonth in ym){
    print(paste0("YM = ",yearmonth,"/",ym[length(ym)]))
    YM_ind = YM_ind + 1
    fc_obs_mat = na.omit(DT_pca[YM == yearmonth,.SD,.SDcols = c("SST_bar",paste0("fc",1:MC_sample_size))])
    
    # get ranks
    ranks.matrix[YM_ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
    rm(fc_obs_mat)
  }
  
  names_vec = c("YM","mst.rk.obs",paste0("mst.r.",1:MC_sample_size),
                "av.rk.obs",paste0("av.r.",1:MC_sample_size),
                "bd.rk.obs",paste0("bd.rk",1:MC_sample_size))
  
  ranks = data.table(ranks.matrix)
  
  setnames(ranks, names_vec)
  
  # --- save ---
  
  save(ranks,file = paste0(save_dir,"/ranks_pca_em_",PCs,"pcs.Rdata"))
  
  # ---- plotting ----
  pdf(file=paste0(plot_dir,"/rks_pca_em_",PCs,"pcs.pdf"),width=8,height=2,points=12)
  par(mfrow=c(1,3),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  rhist.dt(ranks[,.(YM,mst.rk.obs)], ens_size = MC_sample_size,  hist_xlab = "minimum spanning tree")
  rhist.dt(ranks[,.(YM,av.rk.obs)], ens_size = MC_sample_size, hist_xlab = "average")
  rhist.dt(ranks[,.(YM,bd.rk.obs)], ens_size = MC_sample_size, hist_xlab = "band depth")
  
  title(paste0("RHs for forecast by ",PCs," pcs"),outer = TRUE)
  dev.off()

}

###################


###### Example plots of forecasted SST and anomalies w.r.t climatology #######

ex_depth = opt_num_PCs
ex_months = 4:9
ex_month_names = c("April","May","June","July","August","September")
ex_year = 2010

clim_years = 1985:2009 # the years to compute the climatology from

MC_sample_size = 10   # number of plots with independently generated noise

ens_size = 9   #size of forecast ensemble

PCs = opt_num_PCs      # number of considered principal components


#compute climatology
climatology = DT[year %in% clim_years, clim := mean(SST_bar),by = .(grid_id,month)][year == min(year) ,.(Lon,Lat,month,clim)]



for(m in ex_months){
  print(paste0("month = ",m))
    
    #generate noise:
    no_dt = list()
    for(i in 1:MC_sample_size){
      no_dt[[i]] = forecast_PCA(m = m, y = ex_year, PCA_depth = PCs, saveorgo = FALSE)[,.(Lon,Lat,noise), keyby = .(Lon,Lat)][,noise]
    }
    no_dt = as.data.table(no_dt)
    setnames(no_dt,paste0("no",1:MC_sample_size))
    
    DT_pca_plot = no_dt[,c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens_size)) := 
                          DT[year == ex_year & month == m, c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens_size)),
                             with = FALSE]]
    DT_pca_plot[,clim := climatology[month == m, clim]]
    
    # choose random ensemble members (REM) and generate forecast as REM + bias + noise
    
    ens_mem = sample.int(ens_size,MC_sample_size,replace = TRUE)
    for(i in 1:MC_sample_size){
      dummy_dt = DT_pca_plot[,.SD,.SDcols = c(paste0("no",i),paste0("Ens",ens_mem[i]),"Bias_Est")]
      forecast = dummy_dt[[1]] + dummy_dt[[2]] + dummy_dt[[3]]
      DT_pca_plot = DT_pca_plot[,paste0("fc",i) := forecast]
    }
    
    #forecast plots:
    
    rr_sst = range(na.omit(DT_pca_plot[,.SD,.SDcols = c(paste0("fc",1:MC_sample_size))]))
    
    for(i in 1:MC_sample_size){
      plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i))))],
                      rr = rr_sst,
                      mn = paste0("SST forecast for ",ex_month_names[which(ex_months == m)]),
                      save_pdf = TRUE, 
                      save_dir = paste0(plot_dir,"/"),
                      file_name = paste0("m",m,"_fc",i),
                      stretch_par = .8)
    }
    
    #anomaly plot 
    rr_clim = range(na.omit(DT_pca_plot[,.SD - clim,.SDcols = c(paste0("fc",1:MC_sample_size))]))
    
    for(i in 1:MC_sample_size){
      plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i)))-clim)],
                      rr = rr_clim,
                      mn = paste0("Anomaly forecast for ",ex_month_names[which(ex_months == m)]),
                      save_pdf = TRUE, 
                      save_dir = paste0(plot_dir,"/"),
                      file_name = paste0("m",m,"_afc",i),
                      stretch_par = .8)
    }
    
  }



