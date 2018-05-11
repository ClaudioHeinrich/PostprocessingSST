
######################################################
###### Functions for computing variogram scores ######
######################################################


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


variogram_score = function(dt,
                           p = 0.5,
                           eval_years=2001:2010){
  
  
  # get the p-th moment of a standard normal distribution
  p_mom_sn = 2^(p/2)*gamma((p+1)/2)/sqrt(pi) 
  
  var_sc = sum(((dt[,fc_var])^(p/2) * p_mom_sn - 
                  (dt[,obs_diff])^(p/2))^2)
  
  return(var_sc)
}


#####################################################
####################### PCA: ########################
#####################################################

#' Setup for computing variogram scores 
#'
#' @description the variogram score for a multivariate normal forecast distribution X requires the variances of all differences of 
#' univariate projections X_i-X_j, as well as the differences of the corresponding observations obs_i-obs_j. This function computes these values for 
#' the PCA post-processing method for a range of considered numbers of PCs.
#'
#' @param m the months considered.
#' @param dt optional. The wide data table.
#' @param dvec integer vector. Contains the numbers of principal components we want to test for.
#' @param cov_dir string. The covariance directory for the PCA.
#' @param ens_size The number of ensemble members.
#' @param eval_years integer vector contining the years you want to use for validation.
#' @param saveorgo Logical. Should the results be saved?
#' @param savedir the directory to save in.
#' @param finite_time logical. If finite_time == TRUE, the computation is sped up by randomly sampling \code{sample_size} (non-land) locations and only considering pairs of these locations.
#' @param sample_size Integer. See \code{finite_time}.
#' 
#' @return A data table var_sc_prereq
#' 
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")
#' setup_var_sc_PCA(m = 1, dt = DT)}
#'
#' @author Claudio Heinrich
#' 
#' @export


setup_var_sc_PCA = function(m,
                              dt = NULL,
                              dvec = c(1:50),
                              cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                              ens_size = 9,
                              eval_years = 2001:2010,
                              saveorgo = TRUE,
                              save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                              finite_time = FALSE,
                              sample_size = 500){
  
  print(paste("starting month",m))
  
  # for computing pth moments we require the PCA data matrix, and we save it in form of a datatable:
  
  setup_PCA(dt = dt, y = 2010, m=m, cov_dir = cov_dir,max_PCA_depth = max(dvec)) # the year is irrelevant
  PCA <- eval(parse(text = paste0("PCA",m)))
  PCA_DT = fc[year == min(year)][month == min(month),.(Lon,Lat,grid_id)]
  
  for(d in  dvec){
    PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
  }
  
  if (finite_time)
  {
    # if dt contains too many locations computing variogram scores takes a long time. We compute a variogram for randomly sampled sample_size locations that are not land
    na_grid_ids = dt[YM == min(YM)][is.na(SST_bar) | is.na(Ens_bar),grid_id]
    
    if(dt[YM == min(YM)][!(grid_id %in% na_grid_ids),.N] <= sample_size)
    {
      dt = dt[!(grid_id %in% na_grid_ids),]
      fc = fc[!(grid_id %in% na_grid_ids),]
      PCA_DT = PCA_DT[!(grid_id %in% na_grid_ids),]
    } else {
      grid_id_sample = sample(dt[YM == min(YM)][!(grid_id %in% na_grid_ids),grid_id],sample_size)
      
      dt = dt[grid_id %in% grid_id_sample,]
      fc = fc[grid_id %in% grid_id_sample,]
      PCA_DT = PCA_DT[grid_id %in% grid_id_sample,]
    }
    
  }
  
  print("getting variances:")
  
   # build data table that contains pairs of coordinates for all coordinates contained in dt
  
   dt_coor_1 = fc[year == min(year)][month == min(month),.(Lon,Lat,grid_id)]
   setnames(dt_coor_1,c("Lon1","Lat1","grid_id1"))
   
   # add dummy key, do outer full merge with a duplicate, and remove key again
   
   dt_coor_1[,"key" := 1]
   dt_coor_2 = copy(dt_coor_1) 
   setnames(dt_coor_2,c("Lon2","Lat2","grid_id2","key"))
   var_sc_prereq = merge(dt_coor_1,dt_coor_2, by = "key",allow.cartesian = TRUE)
   var_sc_prereq[, "key" := NULL]
   
   
   ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
   ##--- for the predictive distribution with d principal components
    
   d = min(dvec) 
   
   
   n_loc = dt_coor_1[,.N]
   id_1 = mapply(rep,1:n_loc,times = n_loc)
   id_2 = rep(1:n_loc, times = n_loc)
   
   variances = list()
   
   d = min(dvec) 
   
   vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
   variances[[d]] = (vec[id_1]-vec[id_2])^2
   
   for(d in dvec[2:length(dvec)]){
     if(d %% 10 == 0)
     {
       print(paste0("d = ",d))
     }
     vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
     variances[[d]] = variances[[d-1]] + (vec[id_1]-vec[id_2])^2
   }
   names(variances) = paste0("var",dvec)
   
   variances = list(variances)
   variances = rbindlist(variances)
   
   var_sc_prereq = data.table(var_sc_prereq,variances)
   
  
  
    print("getting variances of forecast distribution done, moving to differences of observed residuals and marginal SDs")
    
    # attach differences of observed residuals
    
    dt = dt[month == m & !(grid_id %in% land_grid_id[,unique(grid_id)]),]
    
    for(y in eval_years)
      {
      print(paste0("year ",y))
      k  = sample(ens_size,1)
      dt[year == y,obs_res := trc(.SD+Bias_Est)-SST_bar,.SDcols = paste0("Ens",k)]
      vec = dt[year == y,obs_res]
      var_sc_prereq[,(paste0("ObsDiff_y",y)) := abs(vec[id_1]-vec[id_2])]
      vec2 = dt[year == y,SD_hat]
      var_sc_prereq[,(paste0("MarSD_y",y)) := vec2[id_1]*vec2[id_2]]
      }
    
    # we need the marginal variances for PCA with different number of principal components:
    
    d = min(dvec)
    PCA_DT[,paste0("var",d) := eval(parse(text = paste0("PC",d)))^2]
    for(d in (min(dvec)+1):max(dvec))
    {
      PCA_DT[,paste0("var",d) := eval(parse(text = paste0("var",d-1))) + eval(parse(text = paste0("PC",d)))^2]
    }
    PCA_DT[,paste0("PC",dvec) := NULL]
    
    
    if(saveorgo)
    {
      print("saving...")
      save(var_sc_prereq,PCA_DT,file = paste0(save_dir,"diff_var_by_PC_m",m,".RData"))
    } else
    {
    return(var_sc_prereq)
    }
    print(paste0("month ",m," done."))
}



#' Computing variogram scores for PCA post-processing 
#' 
#' @param dvec integer vector. Contains the numbers of principal components we want to test for.
#' @param p double. The power used in the variogram score.
#' @param months,eval_years integer vectors containing the months and years considered.
#' @param save_dir the directory to save in and to load the prepared data tables from (see setup_var_sc_PCA)
#' @param file_name name of the saved file.
#' 
#' 
#' @examples \dontrun{
#' save_dir = "~/PostClimDataNoBackup/SFE/Derived/NAO/"
#' var_sc_PCA(months = 1, save_dir = save_dir)}
#'
#' @author Claudio Heinrich
#' 
#' @export

   
var_sc_PCA = function(dvec = 1:50,
                      p = 0.5,
                     months = 1:12,
                     eval_years = 2001:2010,
                     save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                     file_name = "var_sc_by_PC.RData",
                     mar_var_corr = TRUE)
{
   sc = list()
   
   for(m in months)
   { print(paste0("month = ",m))
     load(file = paste0(save_dir,"diff_var_by_PC_m",m,".RData"))
     
     #setup for marginal variance correction factor
     n_loc = PCA_DT[,.N]
     id_1 = mapply(rep,1:n_loc,times = n_loc)
     id_2 = rep(1:n_loc, times = n_loc)
     
     var_sc_by_PC = function(d){
       return_data =  data.table(year = eval_years, d = d)
       
       if(mar_var_corr)
       {
       # marginal variance correction factor:
       sdvec = sqrt(PCA_DT[,eval(parse(text = paste0("var",d)))])
       mvcf = sdvec[id_1] * sdvec[id_2]
       }
       
       for(y in eval_years)
       {dt_temp = var_sc_prereq[,.SD,.SDcols = c(paste0("var",d),paste0("ObsDiff_y",y))]
        
        # marginal variance correction:
       if(mar_var_corr)
       {
       MarSD = var_sc_prereq[,.SD,.SDcols = paste0("MarSD_y",y)]
       dt_temp[,1] = dt_temp[,1]*MarSD/mvcf
       }
       
         setnames(dt_temp,c("fc_var","obs_diff"))
         return_data[year == y, sc := variogram_score(dt_temp, p=p, eval_years = eval_years)]
      }
     return_data[,month := m]
     return(return_data)
     }
     sc_temp = parallel::mclapply(dvec,var_sc_by_PC,mc.cores = 12)
     sc = c(sc,sc_temp)
   }
   sc = rbindlist(sc)
   setkey(sc,d,year,month)
   
   save(sc, file = paste0(save_dir,file_name))
   
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
                          finite_time = FALSE,
                          sample_size = 500)
{
  if (is.null(dt))
  {
  load_combined_wide(var = TRUE)
  }
  
  
  if (finite_time)
  {
    # if dt contains more than 500 locations (i.e. more than 250000 pairs of locations) 
    # it just takes a long time and we compute a variogram for randomly sampled 500 locations that are not land
    na_grid_ids = dt[YM == min(YM)][is.na(SST_bar) | is.na(Ens_bar),grid_id]
    
    if(dt[YM == min(YM)][!(grid_id %in% na_grid_ids),.N] <= sample_size)
    {
      dt = dt[!(grid_id %in% na_grid_ids),]
    } else {
      grid_id_sample = sample(dt[YM == min(YM)][!(grid_id %in% na_grid_ids),grid_id],sample_size)
      dt = dt[grid_id %in% grid_id_sample,]
    }
    
    # get distance matrix:
    
    sp <- SpatialPoints(cbind(x=dt[YM == min(YM), Lon], y=dt[YM == min(YM), Lat]), proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    Dist_vs <- sp::spDists(sp, longlat = TRUE)
  }
  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt that are not land_id
  land_grid_id = which(dt[,is.na(SST_bar) | is.na(Ens_bar)])
  dt = dt[-land_grid_id,]
  
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
    if(!finite_time)
    {
      Dist_vs = Dist
    }
    
    Sigma <- psill*exp(-Dist_vs/range)
    
    print("getting variances")
    
    id_1 = var_sc_prereq[,grid_id_ind1]
    id_2 = var_sc_prereq[,grid_id_ind2]
    ff = id_1 + (id_2 - 1)*dim(Sigma)[1]
    var_sc_prereq[,Var:=Sigma[ff]]
    
    
    print("getting variances of forecast distribution done, moving to obs differences")
    
    # attach differences of observed residuals
    
    dt = dt[month == m ,]
    
    for(y in eval_years)
    {
      print(paste0("year ",y))
      k  = sample(ens_size,1)
      dt_temp = dt[year == y,]
      dt_temp[,obs_res := trc(.SD+Bias_Est)-SST_bar,.SDcols = paste0("Ens",k)]
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
#' @param p double. The power used in the variogram score.
#' @param months,eval_years integer vectors containing the months and years considered.
#' @param saveorgo logical. Whether or not the results are saved.
#' @param save_dir the directory to save in and to load the prepared data tables from (see setup_var_sc_geoStat)
#' @param file_name name of the saved file.
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
                      file_name = "var_sc.RData")
{ sc = as.data.table(expand.grid(months,eval_years))
  setnames(sc,c("month","year"))
  
  for(m in months)
  { print(paste0("month = ",m))
    load(file = paste0(save_dir,"diff_var_geoStat_m",m,".RData"))
      for(y in eval_years)
      {dt_temp = var_sc_prereq[,.SD,.SDcols = c("Var",paste0("ObsDiff_y",y))]
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

#' Variogram scores for ECC rely on estimate of the pth moment for the forecast distribution and not on assumed normality. Thus, we need a differentfct to compute variograms.
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
                           eval_years=2001:2010){
  
  var_sc = sum((dt[,fc_var]- dt[,obs_diff])^2)
  return(var_sc)
}

#' Setup for computing variogram scores for ECC
#'
#' @description the variogram score requires the pth moment of all differences of 
#' univariate projections of the forecast distribution X_i-X_j, as well as the differences of the corresponding observations obs_i-obs_j. 
#' This function computes the moment estimate for X_i-X_j for ECC post-processed forecasts and the differences in the observation.
#'
#' @param dt optional. The wide data table.
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
    
    for(y in eval_years)
    {print(paste0("year = ",y))
      dt_temp_2 = dt_temp_1[year == y,.SD,
                   .SDcols = c(paste0("ecc_fc",1:ens_size),
                               paste0("Ens",1:ens_size),
                               "Bias_Est","SST_bar")]
      mom_vec = 0
      # get moment estimate:
      for(i in 1:ens_size)
      {
        vec = dt_temp_2[,eval(parse(text = paste0("ecc_fc",i)))]
        mom_vec = mom_vec + abs(vec[id_1]-vec[id_2])^p
      }
      mom_vec = mom_vec/ens_size
      var_sc_prereq[,(paste0("mom_est_y",y)) := mom_vec]
      
      #get observation:
      k  = sample(ens_size,1)
      dt_temp_2[,obs_res := trc(.SD+Bias_Est)-SST_bar,.SDcols = paste0("Ens",k)]
      vec = dt_temp_2[,obs_res]
      var_sc_prereq[,(paste0("ObsDiff_y",y)) := abs(vec[id_1]-vec[id_2])^p]
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
{ sc = as.data.table(expand.grid(months,eval_years))
setnames(sc,c("month","year"))

for(m in months)
{ print(paste0("month = ",m))
  load(file = paste0(save_dir,"diff_var_ECC_m",m,".RData"))
  for(y in eval_years)
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

