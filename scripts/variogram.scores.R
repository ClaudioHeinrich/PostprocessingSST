
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

dvec = c(1:35,40,45,50) # the number of principal components the variogram score is computed for

p = 0.5 # power for variogram score

eval_years = 2001:2010 # the years for which the score is computed

## --------------------

var_sc = list() 

score_by_month = function(month){
  
  print(paste("starting month",month))

# for computing pth moments we require the PCA data matrix, and we save it in form of a datatable:

setup_PCA(y=year, m=month, oriented = TRUE)
PCA <- eval(parse(text = paste0("PCA",month)))
PCA_DT = fc[year == min(year)][month == min(month),.(Lon,Lat,grid_id)]

for(d in  dvec){
  PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
}

PCA_DT_small = PCA_DT[lon.min <= Lon][Lon <= lon.max][lat.min <= Lat][Lat <= lat.max,]
  
 # build data table that contains pairs of coordinates within the region specified by lon.min/max and lat.min/max

 dummy_dt = fc[year == min(year)][month == min(month)][lon.min <= Lon][Lon <= lon.max & 
                      lat.min <= Lat][Lat <= lat.max,  .(Lon,Lat,grid_id)]
 setnames(dummy_dt,c("Lon1","Lat1","grid_id1"))
 
 # add dummy key, do outer full merge with a duplicate, and remove key again
 
 dummy_dt[,"key" := 1]
 dummy_2 = copy(dummy_dt) 
 setnames(dummy_2,c("Lon2","Lat2","grid_id2","key"))
 
 var_sc_prereq = merge(dummy_dt,dummy_2, by = "key",allow.cartesian = TRUE)
 var_sc_prereq[, "key" := NULL]
 
 
 ##--- For each pair of coordinates (i,j) compute the variance of X_i-X_j 
 ##--- for the predictive distribution with d principal components
  
  var_d = function(d,coor1,coor2){
    i = match(coor1,PCA_DT_small[,grid_id])
    j = match(coor2,PCA_DT_small[,grid_id])
    charvec = paste0("PC",1:d)
    sq_diff_sum=0
    for(sum_ind in 1:d){
      sq_diff_sum = sq_diff_sum +(PCA_DT_small[,eval(parse(text = charvec[sum_ind]))][i] - 
                                  PCA_DT_small[,eval(parse(text = charvec[sum_ind]))][j]  )^2
    }
    return(sq_diff_sum)
  }
  
  
  for(d in dvec){
    var_sc_prereq[,(paste0("Var",d)) := var_d(d,grid_id1,grid_id2)]
    print( paste0("computing variances for ",d," principal components complete"))
  }
  
  
  # --- now get the differences of residuals required for the variogram score
   
   obs_fc_data = load_combined_wide(bias = TRUE)
   
   obs_fc_data_small = obs_fc_data[lon.min <= Lon &
                                     Lon <= lon.max &
                                     lat.min <= Lat &
                                     Lat <= lat.max &
                                     month == month]

   obs_diff = function(y,coor1,coor2){
     i = match(coor1,obs_fc_data_small[,grid_id])
     j = match(coor2,obs_fc_data_small[,grid_id])
     
     abs_diff_obs = (abs(obs_fc_data_small[year == y,SST_bar-SST_hat][i] - 
                                     obs_fc_data_small[y == year,SST_bar-SST_hat][j]  ))
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

save.dir="./Data/PostClim/SFE/Derived/PCA/"
save(var_sc, file = paste0(save.dir,"variogram_scores.RData"))

# --- plots ---

plot.dir = "./figures/"

mean_sc = var_sc[,mean(v_sc), by = PCs]
mean_sc_all_years = var_sc_all_years[,mean(v_sc), by = PCs]

yrange = range(c(range(mean_sc[[2]]),range(mean_sc_all_years[[2]])))

pdf(paste0(plot.dir,"mean_variogram_scores.pdf"))
plot(x = mean_sc[[1]],
     y = mean_sc[[2]],
     ylim = yrange,
     type = "b",
     col = "blue",
     main = "mean variogram scores for p=0.5",
     xlab = "Principal Components",
     ylab = "mean score"
)
lines( mean_sc_all_years[[1]], 
       mean_sc_all_years[[2]],
       type = "b",
       col = "darkgreen")

#---- add minima -----
abline(h = min(mean_sc[[2]]), lty = "dashed", col = adjustcolor("blue",alpha = .5))
abline(h = min(mean_sc_all_years[[2]]), lty = "dashed", col = adjustcolor("darkgreen",alpha = .5))

min_loc = mean_sc[,which.min(V1)]
min_loc_all_years = mean_sc_all_years[,which.min(V1)]

points(x = mean_sc[[1]][min_loc],
       y = mean_sc[[2]][min_loc],
       col = "blue",
       bg = "blue",
       pch = 21)
points(x = mean_sc_all_years[[1]][min_loc_all_years],
       y = mean_sc_all_years[[2]][min_loc_all_years],
       col = "darkgreen",
       bg = "darkgreen",
       pch = 21)

legend("topright",legend = c("1985 - 2000","all years"),lty = c(1,1),col = c("blue","darkgreen"))
dev.off()


