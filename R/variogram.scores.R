
######################################################
###### Functions for computing variogram scores ######
######################################################


#' Computing variogram scores
#' 
#' Computes the variogram for a normal forecasting distribution
#'
#' @param dt a data table containing the columns .(sq_mean_diff, dvar, sq_obs_diff), containing the squared differences of forecast means, the variances for differences of the forecast distribution,
#' as well as the square of pairwise differences of observations.
#' @param p power for the variogram score.
#' 
#' @return a double var_sc
#' 
#' @author Claudio Heinrich
#' 
#' @export


variogram_score_nrm_p2 = function(dt,weights = TRUE)
  {
  
  # get the p-th moment of the forecasting distribution, which is non-central Gaussian
  if(!weights)
  {
    var_sc = sum(dt[,(sq_mean_diff + dvar - sq_obs_diff)^2])/dt[,.N]  
  }
  if(weights)
  {
    var_sc = sum(dt[,weights * (sq_mean_diff + dvar - sq_obs_diff)^2])
  }
  
  return(var_sc)
}


get_moment_nrm = function(dt,p = 0.5)
{
  fc_mean = dt[,fc_mean]
  fc_var = dt[,fc_var]
  
  # get the p-th abs. moment of the forecasting distribution, which is non-central Gaussian
  f = 2^(p/2) * gamma((p+1)/2)/sqrt(pi)
  p_mom := f * fc_var^(p/2) * hypergeo::genhypergeo(-p/2,1/2,-1/2*(fc_mean/sqrt(fc_var))^2)
  return(p_mom)
}

#' Computing variogram scores
#'
#' @param dt a data table containing the columns .(year, grid_id1,grid_id2,fc_var, obs_diff) 
#' @param p power for the variogram score.
#' @param eval_years integer vector containing the years you want to use for validation, need to be contained in dt[,year]
#' 
#' @return a double var_sc
#' 
#' @author Claudio Heinrich
#' 
#' @export


variogram_score_old = function(dt,
                           p = 0.5,
                           eval_years=2001:2010){
  
  
  # get the p-th moment of a standard normal distribution
  p_mom_sn = 2^(p/2)*gamma((p+1)/2)/sqrt(pi) 
  
  var_sc = sum(((dt[,fc_var])^(p/2) * p_mom_sn - 
                  (dt[,obs_diff])^p)^2)
  
  return(var_sc)
}


#####################################################
####################### PCA: ########################
#####################################################


#' Computing variogram scores for PCA post-processing 
#'
#' @description  This function computes the variogram scores for
#' the PCA post-processing method for a given month and year for a range of considered numbers of PCs.
#'
#' @param m,y The month and year.
#' @param dt The data table.
#' @param p The power for the variogram score.
#' @param PCA,PCA_DT Optional. These objects depend only on the month, not the year. Therefore, to speed things up, they can be computed prior to running the function for a range of years.
#' @param dvec Integer vector. Contains the numbers of principal components we want to test for.
#' @param marginal_correction Logical. If TRUE, PCA with corrected marginal variance is considered (see paper).
#' @param cov_dir String. The covariance directory for the PCA.
#' @param ens_size The number of ensemble members.
#' @param save_dir The directory to save in.
#' @param file_name The file is saved as \code{file_name(_nmc)_m#_y##.RData}, where # is the month and ## the year and _nmc is added if no marginal correction is chosen, i.e. if marginal_correction == F.
#' 
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' var_sc_PCA(m = 1,y = 2000, dt = DT)}
#'
#' @author Claudio Heinrich
#' 
#' @export


var_sc_PCA = function(m, y, dt,
                      p = 0.5,
                      PCA = NULL, PCA_DT = NULL,
                      dvec = c(1:10,12,14,16,18,20,25,30,50,70),
                      marginal_correction = TRUE,
                      cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                      ens_size = 9,
                      save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                      file_name = "var_sc_by_PC",
                      weighted = TRUE,
                      training_years = 1985:2000,
                      weight_fct = NULL)
{
  
  if(weighted)
  {
    # the weights are a function of the SST difference averaged over training_years
    weight_dt = dt[year %in% training_years & month == m,][,mean(SST_bar), by = grid_id]
    setnames(weight_dt,c("grid_id","mean_SST"))
  }
  
  
  # for computing pth moments we require the PCA data matrix, and we save it in form of a datatable:
  
  dt = dt[year %in% y & month %in% m,]
  
  land_ids <- which(dt[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt = dt[-land_ids,]
    
    if(weighted)
    {
      weight_dt = weight_dt[-land_ids,]
    }
  }
  
  #get covariance matrix
  
  if(is.null(PCA))
  {
    load(file = paste0(cov_dir,"CovRes_mon",m,".RData"))
    
    PCA = irlba::irlba(res_cov, nv = max(dvec))
  }
  
  
  if(is.null(PCA_DT))
  {
    PCA_DT = dt[,.(Lon,Lat,grid_id)]
    
    for(d in  min(dvec):max(dvec)){
      if(d ==0) PCA_DT [,paste0("PC",d) := 0]
      if(d>0)
      {
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
      }
    } 
    
    # also get marginal variances
    variances = list()
    d = min(dvec) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = vec^2
    
    if(length(dvec)>1)
    {
      for(d in (min(dvec)+1) : max(dvec)){
        vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
        variances[[d]] = variances[[d-1]] + vec^2
      }  
    }
    
    names(variances) = paste0("var",min(dvec):max(dvec))
    PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))
    PCA_DT[,var0 := 0]
  }
  
  
  # marginal SD correction factor:
  SD = dt[,SD_hat]
  
  
  
  PCA_DT[,paste0("marSDcf",dvec) := SD/sqrt(.SD),.SDcols = paste0("var",dvec)]
  
  # at critical indices we just use the SD of the PCA:
  crit_ind  = which(SD < 0.001 | PCA_DT[,var1] < 1e-20)
  PCA_DT[crit_ind, paste0("marSDcf",dvec) := 1]
  
  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt:
  
  dt_coor_1 = PCA_DT[,.(Lon,Lat,grid_id)]
  setnames(dt_coor_1,c("Lon1","Lat1","grid_id1"))
  # add dummy key, do outer full merge with a duplicate, and remove key again:
  dt_coor_1[,"key" := 1]
  dt_coor_2 = copy(dt_coor_1) 
  setnames(dt_coor_2,c("Lon2","Lat2","grid_id2","key"))
  var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
  
  ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
  ##--- for the predictive distribution with d principal components, first without marginal variance correction
  
  n_loc = dt_coor_1[,.N]
  id_1 = as.vector(mapply(rep,1:n_loc,times = n_loc))
  id_2 = rep(1:n_loc, times = n_loc)
  
  if(weighted)
  {
    if(is.null(weight_fct))
    {
      weight_fct = function(x)
      { y = rep(1, length(x))
      y[abs(x)>1] = (1/x[abs(x)>1]^4)
      return(y)
      }
    }
    
    var_sc_prereq[,"weights" := weight_fct(weight_dt[,mean_SST][id_1] - weight_dt[,mean_SST][id_2])]
    #normalizing:
    var_sc_prereq[,"weights" := weights/sqrt(sum(weights^2))]
  } else {
    var_sc_prereq[,"weights" := 1/sqrt(.N)]
  }
  
  
  print("getting variances of differences:")
  
  diff_var = list()
  
  if(marginal_correction)
  {
    list_page = 1
    for(d in dvec)
    {
      print(paste0("d = ",d))
      mar_cor = PCA_DT[,eval(parse(text = paste0("marSDcf",d)))]
      temp = mar_cor * PCA_DT[,.SD,.SDcols = paste0("PC",1:d)]
      diff_var[[list_page]] =  rowSums((temp[id_1]-temp[id_2])^2)
      list_page = list_page + 1
    }  
    names(diff_var) = paste0("dvar",dvec)
    
  } else if(!marginal_correction)
  {
    d = min(dvec) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    diff_var[[1]] = (vec[id_1]-vec[id_2])^2
    
    list_page = 2
    for(d in (min(dvec)+1):max(dvec)){
      vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
      diff_var[[list_page]] = diff_var[[list_page-1]] + (vec[id_1]-vec[id_2])^2
      list_page = list_page +1
    }
    names(diff_var) = paste0("dvar",min(dvec):max(dvec))
  }
  
  diff_var = rbindlist(list(diff_var))
  
  var_sc_prereq = data.table(var_sc_prereq,diff_var)
  
  #complement var_sc_prereq by the suqared differences of the mean vectors and squared differences of observation:
  
  sq_m_diff = 0
  for(k in  1:ens_size)
  {
    vec = dt[,eval(parse(text = paste0("EMF",k)))]
    sq_m_diff = sq_m_diff +(vec[id_1]-vec[id_2])^2
  }
  
  var_sc_prereq[,sq_m_diff := sq_m_diff/ens_size]
  
  var_sc_prereq[,"sq_obs_diff" := (dt[,SST_bar][id_1]-dt[,SST_bar][id_2])^2]
  
  
  # get variogram scores
  print("done - compute variogram scores:")
  scores = list()
  
  #var_sc_by_PC = function(d){
  list_page = 1
  
  for(d in dvec){
    print(paste0("d = ",d))
    return_data =  data.table(year = y, d = d)
    
    dt_temp = var_sc_prereq[,.SD,.SDcols = c("sq_m_diff",paste0("dvar",d),"sq_obs_diff","weights")]
    setnames(dt_temp,c("sq_mean_diff", "dvar", "sq_obs_diff","weights"))
    return_data[,sc := variogram_score_nrm_p2(dt_temp)]
    
    
    scores[[list_page]] = return_data[,month := m]  
    list_page = list_page +1
  }
  
  scores = rbindlist(scores)
  
  if(!marginal_correction)
  {
    file_name = paste0(file_name,"_nmc")
  }
  
  save(scores, file = paste0(save_dir,file_name,"_m",m,"_y",y,".RData"))
}


#' Computing variogram scores for PCA post-processing 
#'
#' @description  This function computes the variogram scores for
#' the PCA post-processing method for a given month and year for a range of considered numbers of PCs.
#'
#' @param m,y The month and year.
#' @param dt The data table.
#' @param p The power for the variogram score.
#' @param PCA,PCA_DT Optional. These objects depend only on the month, not the year. Therefore, to speed things up, they can be computed prior to running the function for a range of years.
#' @param dvec Integer vector. Contains the numbers of principal components we want to test for.
#' @param marginal_correction Logical. If TRUE, PCA with corrected marginal variance is considered (see paper).
#' @param cov_dir String. The covariance directory for the PCA.
#' @param ens_size The number of ensemble members.
#' @param save_dir The directory to save in.
#' @param file_name The file is saved as \code{file_name(_nmc)_m#_y##.RData}, where # is the month and ## the year and _nmc is added if no marginal correction is chosen, i.e. if marginal_correction == F.
#' 
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' var_sc_PCA(m = 1,y = 2000, dt = DT)}
#'
#' @author Claudio Heinrich
#' 
#' @export


var_sc_PCA_old = function(m, y, dt,
                          p = 0.5,
                          PCA = NULL, PCA_DT = NULL,
                          dvec = c(1:10,12,14,16,18,20,25,30,50,70),
                          marginal_correction = TRUE,
                          cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                          ens_size = 9,
                          save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                          file_name = "var_sc_by_PC")
{
  # for computing pth moments we require the PCA data matrix, and we save it in form of a datatable:
  
  dt = dt[year %in% y & month %in% m,]
  
  land_ids <- which(dt[, is.na(Ens_bar) | is.na(SST_bar)])
  if(!identical(land_ids,integer(0)))
  {
    dt = dt[-land_ids,]
  }
  
  #get covariance matrix
  
  if(is.null(PCA))
  {
    load(file = paste0(cov_dir,"CovRes_mon",m,".RData"))
    
    PCA = irlba::irlba(res_cov, nv = max(dvec))
  }
  
  
  if(is.null(PCA_DT))
  {
    PCA_DT = dt[,.(Lon,Lat,grid_id)]
    
    for(d in  min(dvec):max(dvec)){
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
    } 
  
    # also get marginal variances
    variances = list()
    d = min(dvec) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = vec^2
    
    if(length(dvec)>1)
    {
      for(d in (min(dvec)+1) : max(dvec)){
        vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
        variances[[d]] = variances[[d-1]] + vec^2
      }  
    }
    
    names(variances) = paste0("var",min(dvec):max(dvec))
    PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))  
  }
  
  
  # marginal SD correction factor:
  SD = dt[,SD_hat]
  
  crit_ind  = which(SD < 0.001 | PCA_DT[,var1] < 1e-20)
  
  PCA_DT[,paste0("marSDcf",dvec) := SD/sqrt(.SD),.SDcols = paste0("var",dvec)]
  
  # at critical indices we just use the SD of the PCA:
  PCA_DT[crit_ind, paste0("marSDcf",dvec) := 1]

  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt:
  
  dt_coor_1 = PCA_DT[,.(Lon,Lat,grid_id)]
  setnames(dt_coor_1,c("Lon1","Lat1","grid_id1"))
  # add dummy key, do outer full merge with a duplicate, and remove key again:
  dt_coor_1[,"key" := 1]
  dt_coor_2 = copy(dt_coor_1) 
  setnames(dt_coor_2,c("Lon2","Lat2","grid_id2","key"))
  var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
  
  ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
  ##--- for the predictive distribution with d principal components, first without marginal variance correction
  
  n_loc = dt_coor_1[,.N]
  id_1 = as.vector(mapply(rep,1:n_loc,times = n_loc))
  id_2 = rep(1:n_loc, times = n_loc)
  
  print("getting variances of differences:")
  
  diff_var = list()
  
  if(marginal_correction)
  {
    list_page = 1
    for(d in dvec)
    {
      print(paste0("d = ",d))
      mar_cor = PCA_DT[,eval(parse(text = paste0("marSDcf",d)))]
      temp = mar_cor * PCA_DT[,.SD,.SDcols = paste0("PC",1:d)]
      diff_var[[list_page]] =  rowSums((temp[id_1]-temp[id_2])^2)
      list_page = list_page + 1
    }  
    names(diff_var) = paste0("dvar",dvec)
    
  } else if(!marginal_correction)
  {
    d = min(dvec) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    diff_var[[1]] = (vec[id_1]-vec[id_2])^2
    
    list_page = 2
    for(d in (min(dvec)+1):max(dvec)){
      vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
      diff_var[[list_page]] = diff_var[[list_page-1]] + (vec[id_1]-vec[id_2])^2
      list_page = list_page +1
    }
    names(diff_var) = paste0("dvar",min(dvec):max(dvec))
  }
  
  diff_var = rbindlist(list(diff_var))
  
  var_sc_prereq = data.table(var_sc_prereq,diff_var)
  
  #complement var_sc_prereq by difference of observation:
  
  k  = sample(ens_size,1)
  dt[,obs_res := trc(.SD+Bias_Est)-SST_bar,.SDcols = paste0("Ens",k)]
  vec = dt[,obs_res]
  var_sc_prereq[,(paste0("ObsDiff_y",y)) := abs(vec[id_1]-vec[id_2])]
  
  
  # get variogram scores
  print("done - compute variogram scores:")
  scores = list()
  
  #var_sc_by_PC = function(d){
  list_page = 1
  
  for(d in dvec){
    print(paste0("d = ",d))
    return_data =  data.table(year = y, d = d)
    
    dt_temp = var_sc_prereq[,.SD,.SDcols = c(paste0("dvar",d),paste0("ObsDiff_y",y))]
    setnames(dt_temp,c("fc_var","obs_diff"))
    return_data[, sc := variogram_score(dt_temp, p=p, eval_years = eval_years)]
    scores[[list_page]] = return_data[,month := m]
    list_page = list_page +1
  }
  
  scores = rbindlist(scores)
  
  if(!marginal_correction)
  {
    file_name = paste0(file_name,"_nmc")
  }
  
  save(scores, file = paste0(save_dir,file_name,"_m",m,"_y",y,".RData"))
}




########################################################
################## Geostationary: ######################
########################################################

#' Setup for computing variogram scores 
#'
#' @description the variogram score for a multivariate normal forecast distribution X requires the variances of all differences of 
#' univariate projections X_i-X_j, as well as the differences of the corresponding observations obs_i-obs_j. This function computes 
#' these values for the (exponential) geostationary post-processing method.
#'
#' @param dt optional. The wide data table.
#' @param months the months considered.
#' @param eval_years integer vector contining the years you want to use for validation.
#' @param savedir the directory to save in.
#' @param file_name the name of the saved file.
#' @param data_dir where the data of the fitted geostationary model is stored.
#' @param var_file_names character vector. Contains the names of the files containing the fitted variograms.
#' @param finite_time logical. If finite_time == TRUE, the computation is sped up by randomly sampling \code{sample_size} (non-land) locations and only considering pairs of these locations.
#' @param sample_size Integer. See \code{finite_time}.
#' 
#' @return A data table var_sc_prereq
#' 
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' geostat_dir = paste0(save_dir,"/GeoStat")
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' setup_var_sc_geoStat(dt = DT,data_dir = geostat_dir,save_dir = geostat_dir)
#' }
#'
#' @author Claudio Heinrich
#' 
#' @export
   
setup_var_sc_geoStat = function(dt = NULL,
                          months = 1:12,
                          eval_years = 2001:2010,
                          save_dir = "~/PostClimDataNoBackup/SFE/Derived/GeoStat/",
                          file_name = "diff_var_geoStat_m",
                          data_dir = "./Data/PostClim/SFE/Derived/GeoStat/",
                          var_file_names = paste0("variogram_exp_m",months,".RData"),
                          mar_var_cor = TRUE)
{
  if (is.null(dt))
  {
  load_combined_wide(var = TRUE)
  }
  
  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt that are not land_id
  land_grid_id = which(dt[,is.na(SST_bar) | is.na(Ens_bar)])
  if(!identical(land_grid_id,integer(0)))
  {
    dt = dt[-land_grid_id,]  
  }
  
  
  dt_coor_1 = dt[YM == min(YM), .(Lat,Lon,grid_id)]
  dt_coor_1[,grid_id_ind := match(grid_id,sort(grid_id))]
  setkey(dt_coor_1,Lat,Lon) #ordering needs to be the same as in geostationary_training
  setnames(dt_coor_1,c("Lat1","Lon1","grid_id1","grid_id_ind1"))
  
  # add dummy key, do outer full merge with a duplicate, and remove key again
  
  dt_coor_1[,"key" := 1]
  dt_coor_2 = copy(dt_coor_1) 
  setnames(dt_coor_2,c("Lat","Lon2","grid_id2","grid_id_ind2","key"))
  var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
    # parallelize rest
  setup_by_month = function(m)
  {
    ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
    load(paste0(data_dir,var_file_names[which(m == months)]))
    
    psill <- Mod$psill[2]
    range <- Mod$range[2]
    nugget <- Mod$psill[1]
    
    # full covariance matrix is used if finite time == FALSE, else the distance matrix Dist_vs has been computed above
    Sigma <- psill*exp(-Dist/range)
    sills <- diag(Sigma) + nugget
    diag(Sigma) <- sills
    mar_var = Sigma[1]
    
    print("getting variances")
    
    id_1 = var_sc_prereq[,grid_id_ind1]
    id_2 = var_sc_prereq[,grid_id_ind2]
    ff = id_1 + (id_2 - 1)*dim(Sigma)[1]
    var_sc_prereq[,covar:=Sigma[ff]] 
    
    
    print("getting variances of forecast distribution done, moving to obs differences")
    
    # attach differences of observed residuals
    
    dt = dt[month == m ,]
    
    for(y in eval_years)
    {
      print(paste0("year ",y))
      k  = sample(ens_size,1)
      dt_temp = dt[year == y,]
      dt_temp[,obs_res := trc(.SD+Bias_Est)-SST_bar,.SDcols = paste0("Ens",k)]
      if(mar_var_cor == TRUE)
      {
        mar_var_corr_fac = dt_temp[,SD_hat/sqrt(mar_var)]
        c_1 = mar_var_corr_fac[id_1]
        c_2 = mar_var_corr_fac[id_2]
        # Use Var(X_i-X_j) = Var(X_i)+Var(X_j) - 2Covar(X_i,X_j):
        var_sc_prereq[,(paste0("covar_y",y)) := (c_1^2 +c_2^2)*mar_var - 2*c_1*c_2*covar] 
      } 
      if(!mar_var_cor)
      {
        var_sc_prereq[,(paste0("covar_y",y)) := 2*(mar_var - covar)]
      }
      vec = dt_temp[,obs_res]
      var_sc_prereq[,(paste0("ObsDiff_y",y)) := abs(vec[id_1]-vec[id_2])]
      
    }
    save(var_sc_prereq,file = paste0(save_dir,file_name,m,".RData"))
    
  }
    # run the parallelized prepfunction by month:
    parallel::mclapply(months, setup_by_month,mc.cores = 12) # setup_by_month saves a file for each run  
}



#' Computing variogram scores for geostationary post-processing 
#' 
#' @param p The power used in the variogram score.
#' @param months,eval_years Integer vectors containing the months and years considered.
#' @param saveorgo Logical. Whether or not the results are saved.
#' @param save_dir The directory to save in and to load the prepared data tables from (see setup_var_sc_geoStat)
#' @param file_name Name of the saved file.
#' 
#' 
#' @examples \dontrun{
#' geostat_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/GeoStat/"
#' var_sc_geoStat(save_dir = geostat_dir)
#' }
#'
#' @author Claudio Heinrich
#' 
#' @export

  
var_sc_geoStat = function(p = 0.5,
                      months = 1:12,
                      eval_years = 2001:2010,
                      saveorgo = TRUE,
                      save_dir = "~/PostClimDataNoBackup/SFE/Derived/GeoStat/",
                      data_name = "diff_var_geoStat_m",
                      file_name = "var_sc.RData")
{ sc = as.data.table(expand.grid(months,eval_years))
  setnames(sc,c("month","year"))
  
  for(m in months)
  { print(paste0("month = ",m))
    load(file = paste0(save_dir,data_name,m,".RData"))
      for(y in eval_years)
      {dt_temp = var_sc_prereq[,.SD,.SDcols = c(paste0("covar_y",y),paste0("ObsDiff_y",y))]
      setnames(dt_temp,c("fc_var","obs_diff"))
      sc[month == m &year == y, sc := variogram_score(dt_temp, p=p, eval_years = eval_years)]
      }
  }
  setkey(sc,year,month)
  
  if(saveorgo) 
  {
    save(sc, file = paste0(save_dir,file_name))
  }else{
    return(sc)
  }
}   

########################################################
####################### ECC: ###########################
########################################################

#' Variogram scores for ECC rely on estimate of the pth moment for the forecast distribution and not on assumed normality. Thus, we need a different function to compute variograms.
#'
#' @param dt a data table containing the columns .(year, grid_id1,grid_id2,fc_var, obs_diff) 
#' @param eval_years integer vector containing the years you want to use for validation, need to be contained in dt[,year]
#' 
#' @return double var_sc
#' 
#' @author Claudio Heinrich
#' 
#' @export


variogram_score_ECC = function(dt,
                           eval_years=2001:2010)
{
  var_sc = sum((dt[,fc_var]- dt[,obs_diff])^2)
  return(var_sc)
}

#' Setup for computing variogram scores for ECC
#'
#' @description the variogram score requires the pth moment of all differences of 
#' univariate projections of the forecast distribution X_i-X_j, as well as the differences of the corresponding observations obs_i-obs_j. 
#' This function computes the moment estimate for X_i-X_j for ECC post-processed forecasts and the differences in the observation.
#'
#' @param dt optional. The wide data table (generated by forecast_ECC).
#' @param p The power for the variogram score.
#' @param months the months considered.
#' @param eval_years integer vector contining the years you want to use for validation.
#' @param ens_size Size of the forecast ensemble.
#' @param savedir the directory to save in.
#' @param file_name the name of the saved files.
#' @param data_dir Only required if dt == NULL. Where we load the ECC forecast from.
#' 
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' ECC_dir = paste0(save_dir,"/ECC")
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' setup_var_sc_ECC(data_dir = ECC_dir,save_dir = ECC_dir)
#' }
#'
#' @author Claudio Heinrich
#' 
#' @export

setup_var_sc_ECC = function(dt = NULL,
                            p = 0.5,
                            months = 1:12,
                            eval_years = 2001:2010,
                            ens_size = 9,
                            save_dir = "~/PostClimDataNoBackup/SFE/Derived/ECC/",
                            file_name = "diff_var_ECC_m",
                            data_dir = "./Data/PostClim/SFE/Derived/ECC/")
{
  if (is.null(dt))
  {
    load(file = paste0(data_dir,"ECC_fc.RData"))
  }
  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt that are not land_id
  land_grid_id = which(dt[,is.na(SST_bar) | is.na(Ens_bar)])
  dt = dt[-land_grid_id,]
  
  dt_coor_1 = dt[YM == min(YM), .(Lat,Lon,grid_id)]
  dt_coor_1[,grid_id_ind := match(grid_id,sort(grid_id))]
  setnames(dt_coor_1,c("Lat1","Lon1","grid_id1","grid_id_ind1"))
  
  # add dummy key, do outer full merge with a duplicate, and remove key again
  
  dt_coor_1[,"key" := 1]
  dt_coor_2 = copy(dt_coor_1) 
  setnames(dt_coor_2,c("Lat2","Lon2","grid_id2","grid_id_ind2","key"))
  var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
  var_sc_prereq[, "key" := NULL]
  
  # parallelize rest
  setup_by_month = function(m)
  {
    print(paste0("month =  ",m))
    
    dt_temp_1 = dt[month == m,]
    print("getting variances")
    
    id_1 = var_sc_prereq[,grid_id_ind1]
    id_2 = var_sc_prereq[,grid_id_ind2]
    
    for(y in eval_years[2:length(eval_years)])
    {
      print(paste0("year = ",y))
      dt_temp_2 = dt_temp_1[year == y,.SD,
                            .SDcols = c(paste0("ecc_fc",1:ens_size),
                                        paste0("Ens",1:ens_size),
                                        "Bias_Est","SST_bar")]
      #get observation:
      k  = sample(ens_size,1)
      dt_temp_2[,obs_res := trc(.SD+Bias_Est)-SST_bar,.SDcols = paste0("Ens",k)]
      vec = dt_temp_2[,obs_res]
      var_sc_prereq[,(paste0("ObsDiff_y",y)) := abs(vec[id_1]-vec[id_2])^p]
      
      # get moment estimate :
      
      mom_vec = 0
      for(i in (1:ens_size)[-k])
      {
        vec = dt_temp_3[,trc(eval(parse(text = paste0("ecc_fc",i)))-SST_bar)]
        mom_vec = mom_vec + abs(vec[id_1]-vec[id_2])^p
      }
      mom_vec = mom_vec/(ens_size-1)
      var_sc_prereq[,(paste0("mom_est_y",y)) := mom_vec]
      
      
      
    }
    
    save(var_sc_prereq,file = paste0(save_dir,file_name,m,".RData"))
    
  }
  # run the parallelized prepfunction by month:
  parallel::mclapply(months, setup_by_month,mc.cores = 12,mc.silent = FALSE) # setup_by_month saves a file for each run  
}



#' Computing variogram scores for ECC post-processing 
#' 
#' @param months,eval_years integer vectors containing the months and years considered.
#' @param saveorgo logical. Whether or not the results are saved.
#' @param save_dir the directory to save in and to load the prepared data tables from (see setup_var_sc_ECC)
#' @param file_name name of the saved file.
#' 
#' 
#' @examples \dontrun{
#' ECC_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/ECC"
#' var_sc_geoStat(save_dir = geostat_dir)
#' }
#'
#' @author Claudio Heinrich
#' 
#' @export


var_sc_ECC = function(months = 1:12,
                      eval_years = 2001:2010,
                      saveorgo = TRUE,
                      save_dir = "~/PostClimDataNoBackup/SFE/Derived/ECC/",
                      file_name = "var_sc.RData")
{ sc = as.data.table(expand.grid(months,eval_years[2:length(eval_years)]))
  setnames(sc,c("month","year"))
  
  for(m in months)
  { print(paste0("month = ",m))
    load(file = paste0(save_dir,"diff_var_ECC_m",m,".RData"))
    for(y in eval_years[2:length(eval_years)])
    {dt_temp = var_sc_prereq[,.SD,.SDcols = c(paste0("mom_est_y",y),paste0("ObsDiff_y",y))]
    setnames(dt_temp,c("fc_var","obs_diff"))
    sc[month == m &year == y, sc := variogram_score_ECC(dt_temp,eval_years = eval_years)]
    }
  }
  setkey(sc,year,month)
  
  if(saveorgo) 
  {
    save(sc, file = paste0(save_dir,file_name))
  }else{
    return(sc)
  }
}   

