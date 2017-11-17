## Now try to align ensemble
rm(list = ls())

##----- Set up ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
library(rgdal)
library(parallel)
##---------------------

##----- Get the SeNorge Data -----
load("~/PostClimDataNoBackup/SFE/Derived/dt_senorge_1957_1990.RData")
load("~/PostClimDataNoBackup/SFE/Derived/dt_senorge_1991_2015.RData")
load("~/PostClimDataNoBackup/SFE/Derived/GCFS1_wide.RData")
load("~/PostClimDataNoBackup/SFE/Derived/senorge2_gfsc1_map.RData")
##--------------------------------

##---- Let's get this going ----
dt_senorge = rbind(dt_senorge_early[year >= 1991], dt_senorge)
##------------------------------

YM_all = unique(dt_senorge[,.(year,month)])

dt_senorge_all = list()
helper = function(ym)
{

  print(ym)
  ##------ Extract and organize observations ----
  y = YM_all[ym,year]
  m = YM_all[ym,month]
  dt_senorge_j = dt_senorge[ (year == y) & (month == m),]
  setkey(dt_senorge_j, senorge_grid_id, lon, lat)
  setkey(dt_senorge_grid, senorge_grid_id, lon, lat)
  dt_senorge_j = merge(dt_senorge_j, dt_senorge_grid)
  ##----------------------------------------------

  ##----- Get Ensemble Information ----
  dt_ens_j = dt_ens[ (year == y) & (month == m)]
  setkey(dt_ens_j, lon, lat, GCFS1_id)
  setkey(dt_ens_grid, lon, lat, GCFS1_id)
  dt_ens_j = merge(dt_ens_j, dt_ens_grid)
  setkey(dt_ens_j, GCFS1_id)
  ##------------------------------------

  ##----- Now Loop Through ---------
  N_ens = 15
  for(k in 1:N_ens)
  {
    E = matrix(NA, dt_senorge_j[,.N],4)
    D = matrix(NA, dt_senorge_j[,.N],4)
    for(j in 1:4)
    {
      E[,j] = dt_ens_j[ .(dt_senorge_j[, get(paste0("ens_id_",j))]),
                       get(paste0("Ens_",j))]
      D[,j] = 1/dt_senorge_j[,get(paste0("dist_ens_",j))]
    }
    D_tot = rowSums(D)
    Ens = rowSums(E * D) / D_tot
    dt_senorge_j[,paste0("Ens_",k):= Ens]
  }
  ##-----------------------------------

  ##----- Clean up and submit ----
  for(j in 1:4)
  {
    dt_senorge_j[,paste0("ens_id_",j) := NULL]
    dt_senorge_j[,paste0("dist_ens_",j) := NULL]
  }
  dt_senorge_j[,"Ens_bar":=rowMeans(dt_senorge_j[,.SD,.SDcols = paste0("Ens_",1:N_ens)])]
  return(dt_senorge_j)
  ##-------------------------------

}

dt_senorge_all = mclapply(1:YM_all[,.N], "helper", mc.cores = 30, mc.silent = FALSE)
dt_senorge = rbindlist(dt_senorge_all)

save(dt_senorge, file = "~/PostClimDataNoBackup/SFE/Derived/senorge2_gfsc1_combined.RData")
