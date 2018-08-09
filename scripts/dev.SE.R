

########### shrinking of eigenvalues ###############

for_res_cov_SE = function(dt = NULL,
                       Y = 1985:2000,
                       M = 1:12,
                       save_dir = "~/PostClimDataNoBackup/SFE/SECov/",
                       ens_size = 9,
                       version = "wrt_ens_mean")
{  
  
  if(is.null(dt))
  {
    print("loading data")
    dt = load_combined_wide(var = TRUE)
  }
  
  land_ids <- which(dt[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt = dt[-land_ids,]
  }
  
  
    for(mon in M)
    {
      print(paste0("month =",mon))  
      
      dt_SE = copy(dt[month == mon & year %in% Y,])
      
      dt_SE[,paste0("Res",1:ens_size):= .SD-SST_bar,.SDcols = paste0("Ens",1:ens_size)]
      
      
      # if(version == "sum_of_squares")
      # {
      #   dt_PCA[,"center" := - Bias_Est]
      #   decor_factor = 1 / sqrt( length(Y) * ens_size )
      #   
      #   sqrt_cov_mat = as.matrix(na.omit(dt_PCA[,.SD - center,.SDcols = paste0("Res",1:9)]))
      #   res_cov = decor_factor * matrix(sqrt_cov_mat,ncol = length(Y) * ens_size)
      # }
      if(version == "wrt_ens_mean")
      {
        dt_SE[,"Res":= SST_bar - trc(Ens_bar + Bias_Est)]  
        cor_factor = 1 / sqrt( length(Y))
        
        sqrt_cov_mat = as.matrix(na.omit(dt_SE[,Res]))
        
        sqrt_cov_mat = cor_factor * matrix(sqrt_cov_mat,ncol = length(Y) )
      
        sc_mat = sqrt_cov_mat %*% t(sqrt_cov_mat)  
      
        p = dim(sc_mat)[1]
        n = length(Y)
        
        m_n = sum(diag(sc_mat))/p
        d_n_sq = sum(diag((sc_mat - diag(rep(m_n,p))) %*% (sc_mat - diag(rep(m_n,p)))))/p
        
        b_n_sq_bar = 0
        for(k in 1:n)
        { x_k_sq = (sqrt_cov_mat[,k]/cor_factor)^2
          mat_temp = t(diag(x_k_sq) - sc_mat) %*% (diag(x_k_sq) - sc_mat)
          b_n_sq_bar = b_n_sq_bar + (1/n^2)*sum(diag(mat_temp))/p
        }
        
        b_n_sq = min(d_n_sq,b_n_sq_bar)
        a_n_sq = d_n_sq - b_n_sq
          
        Sigma = diag(rep(b_n_sq * m_n / d_n_sq,p)) + a_n_sq/d_n_sq * sc_mat
      
        save(Sigma, file = paste0(save_dir,"CovRes_mon",mon,".RData"))
      
      rm(dt_SE)
    }
  }
  
}



var_sc_SE_new = function(m, y, dt, weight_mat,
                          Sigma, 
                          ens_size = 9,
                          save_dir = "~/PostClimDataNoBackup/SFE/Derived/SE/",
                          file_name = "var_sc",
                          weighted = FALSE,
                          weight_fct = NULL)
{
  # setting up:
  
  dt = dt[year %in% y & month %in% m,]
  
  land_ids <- which(dt[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt = dt[-land_ids,]
  }
  
  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt:
  
  dt_coor_1 = dt[year == min(year) & month ==  min(month),.(Lon,Lat,grid_id)]
  setnames(dt_coor_1,c("Lon1","Lat1","grid_id1"))
  # add dummy key, do outer full merge with a duplicate, and remove key again:
  dt_coor_1[,"key" := 1]
  dt_coor_2 = copy(dt_coor_1) 
  setnames(dt_coor_2,c("Lon2","Lat2","grid_id2","key"))
  var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
  # get indices for pairs of locations
  
  n_loc_pair = dt_coor_1[,.N]
  id_1 = as.vector(mapply(rep,1:n_loc,times = n_loc_pair))
  id_2 = rep(1:n_loc_pair, times = n_loc_pair)
  
  
  # find variance and covariance locations in Sigma:
  n_loc = dt[,.N]
  var_id_1 = id_1 + (id_1-1)*n_loc
  var_id_2 = id_2 + (id_2-1)*n_loc
  cov_id = id_1 + (id_2 - 1)*n_loc
  
  
  #get weights for variogram score
  
  if(weighted)
  {
    if(is.null(weight_fct))
    {
      weight_fct = function(x)
      { y = rep(1, length(x))
      y[abs(x)>1] = (1/x[abs(x)>1]^2)
      return(y)
      }
    }
    
    var_sc_prereq[,"weights" := weight_fct(weight_dt[,mean_SST][id_1] - weight_dt[,mean_SST][id_2])]
    #normalizing:
    var_sc_prereq[,"weights" := weights/sum(weights)]
  } else {
    var_sc_prereq[,"weights" := 1/.N]
  }
  
  
  
  ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
  ##--- for the predictive distribution with d principal components
  
  print("getting variances of differences:")
  
  diff_var = list()
  
  
  Sigma_wt = weight_mat * Sigma
  
  var_sc_prereq[,var:= Sigma_wt[var_id_1] + Sigma_wt[var_id_2] - 2*Sigma_wt[cov_id]]
  
  #complement var_sc_prereq by the squared differences of the mean vectors and squared differences of observation:
  
  var_sc_prereq[,sq_m_diff := (dt[,SST_hat][id_1]-dt[,SST_hat][id_2])^2]
  
  var_sc_prereq[,sq_obs_diff := (dt[,SST_bar][id_1]-dt[,SST_bar][id_2])^2]
  
  
  # get variogram scores
  print("done - compute variogram scores:")
  
  
  #var_sc_by_PC = function(d){
  
  return_data =  data.table(year = y)
    
  dt_temp = var_sc_prereq[,.SD,.SDcols = c("var","sq_m_diff","sq_obs_diff","weights")]
  setnames(dt_temp,c("dvar","sq_mean_diff", "sq_obs_diff","weights"))
  return_data[,sc := variogram_score_nrm_p2(dt_temp)]
    
    
  scores = return_data[,month := m]  
    
  save(scores, file = paste0(save_dir,file_name,"_m",m,"_y",y,".RData"))
}


### unfinished work. Implementing the SE forecast seems bad, since its marginal variance is off? Maybe do both 
forecast_SE = function(dt, y, m, n = 1,
                       saveorgo = FALSE,
                       save_dir = "./Data/PostClim/SFE/Derived/SE/",
                       file_name = "forecastSE",
                       cov_dir)
{ 
  #find land grid ids:

  land_ids <- dt[year %in% y & month %in% m, which(is.na(Ens_bar) | is.na(SST_bar))]

  SD_cols = c("Lon","Lat","grid_id","month","year","YM",
              "SST_hat","SST_bar",paste0("Ens",1:ens_size),"Ens_bar","Bias_Est","var_bar","SD_hat")
  fc = dt[year %in% y & month %in% m ,.SD,.SDcols = SD_cols]

  if(!identical(land_ids,integer(0)))
  {
    fc_water <-  dt[year %in% y & month %in% m ,][-land_ids,.SD,.SDcols = SD_cols]
  } else {
    fc_water = fc
  }
  
  # get weight matrix 
  sp <- sp::SpatialPoints(cbind(x=fc_water[YM == min(YM), Lon],
                                y=fc_water[YM == min(YM), Lat]), 
                          proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  Dist <- sp::spDists(sp, longlat = TRUE)
  
  weights = GneitingWeightFct(Dist, L=L)
  weight_mat = matrix(weights, nrow = dim(Dist)[1])
  
  

  for(mon in m)
    {
      load(file = paste0(cov_dir,"CovRes_mon",mon,".RData"))
      
      Sigma_hat = weight_mat * Sigma
      
      for(y_0 in y)
        {
          mcf = fc_water[year == y_0 & month == mon, SD_hat]/sqrt(Sigma[1])
          Sigma_hat = diag(mcf) %*% Sigma %*% diag(mcf)
          no <- matrix(MASS::mvrnorm(n=n, mu=rep(0,length(sp)), Sigma=Sigma_hat),nrow = n)
          for (i in 1:n)
            {
            fc_water[year == y_0 & month == mon, paste0("no",i):= no[i,]]
            fc_water[year == y_0 & month == mon,paste0("fc",i):= trc(Ens_bar + Bias_Est) + .SD, 
                     .SDcols = paste0("no",i)]
            }
        }
    }
  
  # add land:
  
  fc = merge(fc,fc_water, by = colnames(fc), all.x = TRUE)
  setkeyv(fc,key(DT))
  
  if(saveorgo){
    save(fc,file = paste0(save_dir,file_name))
  }


  return(fc)
}




###############################################
########### SE variogram score ################
###############################################


# to skip this section, run instead:

# SE_dir = paste0(save_dir, "SE/")
# load(file = paste0(SE_dir,"variogram_scores.RData"))


##### setting up ######

SE_dir = paste0(save_dir,"SE/")
dir.create(SE_dir, showWarnings = FALSE)

training_years = DT[!(year %in% validation_years),unique(year)]

for_res_cov_SE(Y = training_years,
               dt = DT, 
               save_dir = SE_dir,
               ens_size = ens_size)



L = 2500

NA_rows = DT[,which(is.na(SST_bar) | is.na(SST_hat) | is.na(SD_hat))]

DT_NA_free = DT[-NA_rows,]

sp <- sp::SpatialPoints(cbind(x=DT_NA_free[YM == min(YM), Lon],
                              y=DT_NA_free[YM == min(YM), Lat]), 
                        proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
Dist <- sp::spDists(sp, longlat = TRUE)

weights = weight_fct(Dist)
weight_mat = matrix(weights, nrow = dim(Dist)[1])



#### variogram score computation: ####


for(m in months)
{ 
  # setup: get principal components and marginal variances for the given month m:
  
  print(paste0("m = ",m))
  load(file = paste0(SE_dir,"CovRes_mon",m,".RData"))
  
  land_ids = which(DT[year == min(year) & month == min(month),is.na(Ens_bar) | is.na(SST_bar)])
  
  
  dummy_fct = function(y)
  {
    var_sc_SE_new(m, y, DT,weight_mat = weight_mat, Sigma = Sigma, ens_size = ens_size,
                   weighted = FALSE,
                   save_dir = SE_dir)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
  
}



###### combine: ######

scores_SE = list()
k=0
for(m in months){
  for(y in validation_years)
  {k=k+1
  load(file = paste0(SE_dir,"var_sc_m",m,"_y",y,".RData"))
  scores_SE[[k]] = scores
  }
}
scores_SE = rbindlist(scores_SE)

mean_sc_SE = scores_SE[,mean_sc := mean(sc)]


# save

save(scores_SE_mc,scores_SE_nmc,mean_sc,mean_sc_nmc,file = paste0(PCA_dir,"variogram_scores",weight_name_addition,".RData"))

##########################################
