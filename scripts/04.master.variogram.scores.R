#######################################################################################

###################  master script part 4 - variogram scores  #########################

#######################################################################################

# This script compares several multivariate post-processing models based on their variogram scores
#
# 
# Files generated:
#   
# Data files: scores.bc.sd.sma.Rdata, scores.bc.sd.ema.Rdata
# Plots: mean_scores_sd_sma.pdf, mean_scores_sd_ema.pdf 
#
# Requires previous run of 03.master.var.est.R
# with the same value of name_abbr as below.


##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))
DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")



###############################
########### PCA ###############
###############################


##### setting up ######

training_years = DT[!(year %in% validation_years),unique(year)]

cov_dir = paste0(save_dir,"/PCACov/")
dir.create(cov_dir, showWarnings = FALSE)

for_res_cov(Y = training_years,dt = DT, save_dir = cov_dir,ens_size = ens_size)

PCs = 1:50 # range of PCs to test

# the variogram score computation takes long and can be skipped by running the following lines

# load(file = paste0(save_dir,"var_sc_by_PC_no_marg_corr.RData"))
# sc_nmc = sc
# load(file = paste0(save_dir,"var_sc_by_PC.RData"))

#### variogram score computation: ####

dummy_fct = function(m)
{
   setup_var_sc_PCA(m,DT,dvec = PCs,eval_years = validation_years,cov_dir = cov_dir,save_dir = save_dir)
}

parallel::mclapply(1:12,dummy_fct,mc.cores = 12)

var_sc_PCA(dvec = PCs,
           months = 1:12,
           eval_years = validation_years,
           save_dir = save_dir)

# for comparison do the same without marginal correction:

var_sc_PCA(dvec = PCs,
           months = 1:12,
           eval_years = validation_years,
           save_dir = save_dir,
           file_name = "var_sc_by_PC_no_marg_corr.RData",
           mar_var_corr = FALSE)


load(file = paste0(save_dir,"var_sc_by_PC_no_marg_corr.RData"))
sc_nmc = sc
load(file = paste0(save_dir,"var_sc_by_PC.RData"))

#############################################

mean_sc = sc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
mean_sc = unique(mean_sc)


mean_sc_nmc = sc_nmc[,mean_sc := mean(sc), by = d][,.(d,mean_sc)]
mean_sc_nmc = unique(mean_sc_nmc)

print(paste0("Minimal variogram score is achieved for ",mean_sc[mean_sc == min(mean_sc),d][1],
             " principal components with a score of ",mean_sc[,min(mean_sc)],
             ". Without marginal correction it is achieved for ",mean_sc_nmc[mean_sc == min(mean_sc),d][1],
             " principal components with a score of ",mean_sc_nmc[,min(mean_sc)],"."))


#########################################
########### Geostationary ###############
#########################################

# to skip this section, run instead:

# load(file = paste0(geostat_dir,"var_sc.RData"))

#########################################

geostat_dir = paste0(save_dir, "GeoStat/")
dir.create(geostat_dir, showWarnings = FALSE)

geostationary_training(dt = DT, save_dir = geostat_dir)

setup_var_sc_geoStat(dt = DT,
                     eval_years = validation_years,
                     data_dir = geostat_dir,
                     save_dir = geostat_dir)

var_sc_geoStat(eval_years = validation_years,
               save_dir = geostat_dir)

load(file = paste0(geostat_dir,"var_sc.RData"))


#########################################
################ ECC ####################
#########################################

ECC_dir = paste0(save_dir,"ECC/")
dir.create(ECC_dir, showWarnings = FALSE)

forecast_ECC(dt = DT,save_dir = ECC_dir)

setup_var_sc_ECC(eval_years = validation_years,data_dir = ECC_dir)


###########################################
########### plotting scores: ##############
###########################################

pdf(paste0(plot_dir,"/mean_variogram_scores.pdf"))

rr = range(c(mean_sc[[2]],mean_sc_nmc[[2]]))

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

abline(h = sc[,mean(sc)], lty = "dashed", col = adjustcolor("darkred"))

legend("topright",legend = c("PCA, mar.cor.var.","just PCA","geostat"),col = c("blue","darkgreen","darkred"),lty = c(1,1,2))

dev.off()

