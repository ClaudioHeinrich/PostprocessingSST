#######################################################################################

###################  master script part 4 - variogram scores  #########################

#######################################################################################

# This script compares several multivariate post-processing models based on their variogram scores.
#
# 
# Files generated:
# 
# Directories: PCACov, GeoStat, ECC
# Data files (* in 1:12): PCA/CovRes_mon*.RData, PCA/diff_var_by_PC_m*.RData, PCA/var_sc_by_PC.RData, PCA/var_sc_by_PC_no_marg_corr.RData
#                         GeoStat/diff_var_geoStat_m*.RData, GeoStat/variogram_exp_m*.RData, GeoStat/var_sc.RData,
#                         ECC/ECC_fc.RData, ECC/diff_var_ECC_m*.RData, ECC/var_sc.RData
#
# Plots: mean_variogram_scores.pdf
#
# Requires previous run of 03.master.var.est.R
# with the same value of name_abbr as below.


##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_small" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))
DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")



###############################
########### PCA ###############
###############################


##### setting up ######

PCA_dir = paste0(save_dir,"PCA_new/")
dir.create(PCA_dir, showWarnings = FALSE)

training_years = DT[!(year %in% validation_years),unique(year)]

for_res_cov_new(Y = training_years,
            dt = DT, 
            save_dir = PCA_dir,
            ens_size = ens_size)

PCs = 1:50 # range of PCs to test

months = 1:12


#### variogram score computation: ####

# the variogram score computation takes long and can be skipped (go to the 'combine'-section) if ran previously


for(m in months)
  { 
   # setup: get principal components and marginal variances for the given month m:
  
    print(paste0("m = ",m))
    load(file = paste0(PCA_dir,"CovRes_mon",m,".RData"))
    PCA = irlba::irlba(res_cov, nv = max(PCs))
    
    land_ids = which(DT[year == min(year) & month == min(month),is.na(Ens_bar) | is.na(SST_bar)])
    
    PCA_DT = DT[year == min(year) & month == min(month),][-land_ids,.(Lon,Lat,grid_id)]
    
    for(d in  PCs){
      PCA_DT [,paste0("PC",d) := PCA$d[d]*PCA$u[,d]]
    } 
    
    # also get marginal variances
    variances = list()
    d = min(PCs) 
    vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
    variances[[d]] = vec^2
    
    for(d in PCs[2:length(PCs)]){
      vec = PCA_DT[,eval(parse(text = paste0("PC",d)))]
      variances[[d]] = variances[[d-1]] + vec^2
    }
    names(variances) = paste0("var",PCs)
    PCA_DT = data.table(PCA_DT,rbindlist(list(variances)))  
  
  # now move to getting the variogram scores:
    
  # without marginal correction:
    
  dummy_fct = function(y)
  {
    var_sc_PCA_new_no_marg_corr(m,y,DT,PCA = PCA,PCA_DT = PCA_DT,dvec = PCs,cov_dir = PCA_dir,ens_size = ens_size,save_dir = PCA_dir)
  }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))

  # with marginal correction:
  
  dummy_fct = function(y)
    {
    print(c(paste0("y = ",y),paste0("m = ",m)))
    var_sc_PCA_new(m,y,DT,PCA = PCA,PCA_DT = PCA_DT,dvec = PCs,cov_dir = PCA_dir,ens_size = ens_size,save_dir = PCA_dir)
    }
  parallel::mclapply(validation_years,dummy_fct,mc.cores = length(validation_years))
}



###### combine: ######

scores_nmc = list()
k=0
for(m in months){
  for(y in validation_years)
  {k=k+1
  load(file = paste0(PCA_dir,"var_sc_by_PC_no_marg_corr_m",m,"_y",y,".RData"))
  scores_nmc[[k]] = scores
  }
}
scores_nmc = rbindlist(scores_nmc)

scores_mc = list()
k=0
for(m in months){
  for(y in validation_years)
  {k=k+1
  load(file = paste0(PCA_dir,"var_sc_by_PC_m",m,"_y",y,".RData"))
  scores_mc[[k]] = scores
  }
}
scores_mc = rbindlist(scores_mc)


#############################################

mean_sc = scores_mc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
mean_sc = unique(mean_sc)


mean_sc_nmc = scores_nmc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
mean_sc_nmc = unique(mean_sc_nmc)

print(paste0("Minimal variogram score is achieved for ",mean_sc[mean_sc == min(mean_sc),d][1],
             " principal components with a score of ",mean_sc[,min(mean_sc)],
             ". Without marginal correction it is achieved for ",mean_sc_nmc[mean_sc == min(mean_sc),d][1],
             " principal components with a score of ",mean_sc_nmc[,min(mean_sc)],"."))


#########################################
########### Geostationary ###############
#########################################

# to skip this section, run instead:

# geostat_dir = paste0(save_dir, "GeoStat/")
# 
# load(file = paste0(geostat_dir,"var_sc_nmc.RData"))
# sc_geostat_nmc = sc
# 
# load(file = paste0(geostat_dir,"var_sc.RData"))
# sc_geostat = sc


#########################################

geostat_dir = paste0(save_dir, "GeoStat/")
dir.create(geostat_dir, showWarnings = FALSE)

geostationary_training(dt = DT, save_dir = geostat_dir)

# without marginal variance correction:

setup_var_sc_geoStat(dt = DT,
                     eval_years = validation_years,
                     data_dir = geostat_dir,
                     save_dir = geostat_dir,
                     file_name = "diff_var_geoStat_nmc_m",
                     mar_var_cor = FALSE)

var_sc_geoStat(eval_years = validation_years,
               save_dir = geostat_dir,
               data_name = "diff_var_geoStat_nmc_m",
               file_name = "var_sc_nmc.RData")

# with marginal variance correction:

setup_var_sc_geoStat(dt = DT,
                     eval_years = validation_years,
                     data_dir = geostat_dir,
                     save_dir = geostat_dir)

var_sc_geoStat(eval_years = validation_years,
               save_dir = geostat_dir)



load(file = paste0(geostat_dir,"var_sc_nmc.RData"))
sc_geostat_nmc = sc

load(file = paste0(geostat_dir,"var_sc.RData"))
sc_geostat = sc


#########################################
################ ECC ####################
#########################################

# to skip this section, run instead:

# ECC_dir = paste0(save_dir,"ECC/")
# load(file = paste0(ECC_dir,"var_sc.RData"))
# sc_ECC = sc

#########################################

ECC_dir = paste0(save_dir,"ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

forecast_ECC(dt = DT,save_dir = ECC_dir)

setup_var_sc_ECC(eval_years = validation_years,data_dir = ECC_dir,save_dir = ECC_dir)

var_sc_ECC(eval_years = validation_years,save_dir = ECC_dir)

load(file = paste0(ECC_dir,"var_sc.RData"))
sc_ECC = sc


###########################################
########### plotting scores: ##############
###########################################

# we leave out ECC, because its score is just too bad:
print(paste0("variogram score for ECC is ",sc_ECC[,mean(sc)]))


pdf(paste0(plot_dir,"/mean_variogram_scores.pdf"))
  rr = range(c(mean_sc[[2]],mean_sc_nmc[[2]],sc_geostat[,mean(sc)],sc_geostat_nmc[,mean(sc)]))
  
  #### with marginal correction ####
  
  plot(x = mean_sc[[1]],
       y = mean_sc[[2]],
       ylim = rr,
       type = "b",
       col = "blue",
       main = paste0("mean variogram scores for ",name_abbr),
       xlab = "number of principal components",
       ylab = "mean score")
  
  #---- add minima: -----
  abline(h = min(mean_sc[[2]]), lty = "dashed", col = adjustcolor("blue",alpha = .5))
  
  opt_num_PCs = mean_sc[,which.min(mean_sc)]
  
  points(x = mean_sc[[1]][opt_num_PCs],
         y = mean_sc[[2]][opt_num_PCs],
         col = "blue",
         bg = "blue",
         pch = 21)
  
  #### without marginal correction ####
  
  lines(x = mean_sc_nmc[[1]],
       y = mean_sc_nmc[[2]],
       type = "b",
       col = "darkgreen")
  
  #---- add minima: -----
  abline(h = min(mean_sc_nmc[[2]]), lty = "dashed", col = adjustcolor("darkgreen",alpha = .5))
  
  opt_num_PCs = mean_sc_nmc[,which.min(mean_sc)]
  
  points(x = mean_sc_nmc[[1]][opt_num_PCs],
         y = mean_sc_nmc[[2]][opt_num_PCs],
         col = "darkgreen",
         bg = "darkgreen",
         pch = 21)
  
  # --- add geostat value: ----
  
  abline(h = sc_geostat[,mean(sc)], lty = "dashed", col = adjustcolor("darkred"))
  abline(h = sc_geostat_nmc[,mean(sc)], lty = "dashed", col = adjustcolor("black"))
  abline(h = sc_ECC[,mean(sc)], lty = "dashed")
  
  legend("topright",legend = c("PCA, mar.cor.var.","just PCA","geostat,m.c.v.","geostat"),col = c("blue","darkgreen","darkred"),lty = c(1,1,2,2))
dev.off()

