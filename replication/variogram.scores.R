
rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)



# --- set parameters

months = 1:12
year = 1985:2010

# specify window for scoring

lon.min = -59.5
lon.max = 14.5
lat.min = 30.5
lat.max = 69.5

 # the number of principal components the variogram score is computed for

eval_years = 2001:2010 # the years for which the score is computed

p = 0.5 # power for variogram score

# var_sc_pc()
# the function takes a data table with yearmonth obs and 



#########################

dvec = 1:50

for(m in 1:12)
{
  setup_var_sc_PCA(m,DT,dvec = dvec,eval_years = eval_years,cov_dir = cov_dir,save_dir = save_dir)
}
var_sc_PCA(dvec = dvec,months = 1:12,eval_years = eval_years,save_dir = save_dir)


#########################


setup_var_sc_PCA = function(m,
                              dt = NULL,
                              dvec = c(1:50),
                              cov_dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                              ens_size = 9,
                              eval_years = 2001:2010,
                              saveorgo = TRUE,
                              save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                              finite_time = TRUE,
                              sample_size = 500){
  
  print(paste("starting month",m))

# for computing pth moments we require the PCA data matrix, and we save it in form of a datatable:

setup_PCA(dt = dt, y = 2010, m=m, cov_dir = cov_dir) # the year is irrelevant
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
  
  dummy_fct = function(d){
    return( var_d(d,var_sc_prereq[,grid_id1],var_sc_prereq[,grid_id2]))
    }
  
  dummy_list = mclapply(dvec,dummy_fct,mc.cores = 10)
  names(dummy_list) = paste0("Var",dvec)
  dummy_list = rbindlist(list(dummy_list))
  
  var_sc_prereq = data.table(var_sc_prereq,dummy_list)
  
  print("getting variances of forecast distribution done, moving to obs differences")
  
  # attach differences of observed residuals
  
  dt = dt[month == m,]
  
  
  obs_diff = function(y,coor1,coor2,k){
    i = match(coor1,dt[year == y,grid_id])
    j = match(coor2,dt[year == y,grid_id])
    
    abs_diff_obs = (abs(dt[year == y,.SD+Bias_Est-SST_bar,.SDcols = paste0("Ens",k)][i] - 
                        dt[year == y,.SD+Bias_Est-SST_bar,.SDcols = paste0("Ens",k)][j]  ))
    return(abs_diff_obs)
  }
  
 
  for(y in eval_years)
  {
    print(paste0("year ",y))
    k  = sample(ens_size,1)
    var_sc_prereq[,(paste0("ObsDiff_y",y)) := obs_diff(y,grid_id1,grid_id2,k)]
  }
  
  
  if(saveorgo)
  {
    save(var_sc_prereq,file = paste0(save_dir,"diff_var_by_PC_m",m,".RData"))
  } else
  {
  return(var_sc_prereq)
  }
}



  
  # --- now get the differences of residuals required for the variogram score
   
   
  #dt should be data table containing the columns .(year, grid_id1,grid_id2,fc_var, obs_diff) where year %in% eval_years and 

   variogram_score = function(dt,
                              p = 0.5,
                              eval_years=2001:2010){
     
     
     # get the p-th moment of a standard normal distribution
     p_mom_sn = 2^(p/2)*gamma((p+1)/2)/sqrt(pi) 
     
     var_sc = sum(((dt[,fc_var])^(p/2) * p_mom_sn - 
                      (dt[,obs_diff])^(p/2))^2)
    
     return(var_sc)
   }
     

   
var_sc_PCA = function(dvec = 1:50,
                      p = 0.5,
                     months = 1:12,
                     eval_years = 2001:2010,
                     saveorgo = TRUE,
                     save_dir = "~/PostClimDataNoBackup/SFE/Derived/PCA/",
                     file_name = "var_sc_by_PC.RData")
{
   sc = list()
   
   for(m in months)
   { print(paste0("month = ",m))
     load(file = paste0(save_dir,"diff_var_by_PC_m",m,".RData"))
     var_sc_by_PC = function(d){
       return_data =  data.table(year = eval_years, d = d)
       for(y in eval_years)
       {dt_temp = var_sc_prereq[,.SD,.SDcols = c(paste0("Var",d),paste0("ObsDiff_y",y))]
         setnames(dt_temp,c("fc_var","obs_diff"))
         return_data[year == y, sc := variogram_score(dt_temp, p=p)]
         return_data[year == y, sc := variogram_score(dt_temp, p=p)]
       }
     return_data[,month := m]
     return(return_data)
     }
     sc_temp = mclapply(dvec,var_sc_by_PC,mc.cores = 12)
     sc = c(sc,sc_temp)
   }
   sc = rbindlist(sc)
   setkey(sc,d,year,month)
   
   if(saveorgo) 
   {
     save(sc, file = paste0(save_dir,file_name))
   }
}
   
   
   
setup_var_sc_geoStat = function(dt = NULL,
                          months = 1:12,
                          eval_years = 2001:2010,
                          save_dir = "~/PostClimDataNoBackup/SFE/Derived/GeoStat/",
                          file_name = "diff_var_geoStat_m",
                          data_dir = "./Data/PostClim/SFE/Derived/GeoStat/",
                          var_file_names = paste0("variogram_exp_m",months,".RData"),
                          finite_time = TRUE,
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
  
  
  # build data table that contains pairs of coordinates for all coordinates contained in dt
  
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
    
    dummy_vec = c()
    for(i in 1:var_sc_prereq[,.N])
    {if(i %% 1000 == 0) print(paste0(i," / ",var_sc_prereq[,.N]))
      dummy_vec = c(dummy_vec, Sigma[var_sc_prereq[,grid_id_ind1][i],var_sc_prereq[,grid_id_ind2][i]])
    }
    var_sc_prereq[,Var := dummy_vec]
    
    
    print("getting variances of forecast distribution done, moving to obs differences")
    
    # attach differences of observed residuals
    
    dt = dt[month == m,]
    
    
    obs_diff = function(y,coor1,coor2,k){
      i = match(coor1,dt[year == y,grid_id])
      j = match(coor2,dt[year == y,grid_id])
      
      abs_diff_obs = (abs(dt[year == y,.SD+Bias_Est-SST_bar,.SDcols = paste0("Ens",k)][i] - 
                            dt[year == y,.SD+Bias_Est-SST_bar,.SDcols = paste0("Ens",k)][j]  ))
      return(abs_diff_obs)
    }
    
    
    for(y in eval_years)
    {
      print(paste0("year ",y))
      k  = sample(ens_size,1)
      var_sc_prereq[,(paste0("ObsDiff_y",y)) := obs_diff(y,grid_id1,grid_id2,k)]
    }
    
    
    
      save(var_sc_prereq,file = paste0(save_dir,file_name,m,".RData"))
    
  }
  
    
  # run the parallelized prepfunction by month:
    mclapply(months, setup_by_month,mc.cores = 12) # setup_by_month saves a file for each run  
  
}
  
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
      sc[month == m &year == y, sc := variogram_score(dt_temp, p=p)]
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




# --- plots ---


load(file = paste0(save.dir,"variogram_scores.RData"))

plot.dir = "./figures/"

mean_sc = var_sc[,mean(v_sc), by = PCs]
#mean_sc_all_years = var_sc_all_years[,mean(v_sc), by = PCs]


yrange = range(c(range(mean_sc[[2]])))
#yrange = range(c(range(mean_sc[[2]]),range(mean_sc_all_years[[2]])))

pdf(paste0(plot.dir,"mean_variogram_scores.pdf"))
plot(x = mean_sc[[1]],
     y = mean_sc[[2]],
     ylim = yrange,
     type = "b",
     col = "blue",
     main = "mean variogram score",
     xlab = "number of principal components",
     ylab = "mean score"
)
# lines( mean_sc_all_years[[1]], 
#        mean_sc_all_years[[2]],
#        type = "b",
#        col = "darkgreen")

#---- add minima -----
abline(h = min(mean_sc[[2]]), lty = "dashed", col = adjustcolor("blue",alpha = .5))
#abline(h = min(mean_sc_all_years[[2]]), lty = "dashed", col = adjustcolor("darkgreen",alpha = .5))

min_loc = mean_sc[,which.min(V1)]
#min_loc_all_years = mean_sc_all_years[,which.min(V1)]

points(x = mean_sc[[1]][min_loc],
       y = mean_sc[[2]][min_loc],
       col = "blue",
       bg = "blue",
       pch = 21)
# points(x = mean_sc_all_years[[1]][min_loc_all_years],
#        y = mean_sc_all_years[[2]][min_loc_all_years],
#        col = "darkgreen",
#        bg = "darkgreen",
#        pch = 21)

#legend("topright",legend = c("1985 - 2000","all years"),lty = c(1,1),col = c("blue","darkgreen"))
dev.off()


