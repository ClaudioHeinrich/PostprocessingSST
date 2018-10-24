rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)


name_abbr_1 = "Aut_2018_small" 
save_dir_1= paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr_1,"/")

name_abbr_2 = "Aut_2018_small/GCFS1" 
save_dir_2 = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr_2,"/")


### load data ###

load(file = paste0(save_dir_1,"setup.RData"))

nms_keep = c("Lon","Lat","year","month","grid_id","YM",paste0("Ens",1:9),"SST_bar")
nms_new = c("Lon","Lat","year","month","grid_id","YM",paste0("NPCM_Ens",1:9),"SST_bar")

DT_1 = DT[,.SD,.SDcols = nms_keep]
names(DT_1) = nms_new
    
load(file = paste0(save_dir_2,"setup.RData"))

nms_keep = c("grid_id","YM",paste0("Ens",1:15))
nms_new = c("grid_id","YM",paste0("GCFS_Ens",1:15))

DT_2 = DT[,.SD,.SDcols=nms_keep]
names(DT_2) = nms_new

DT = merge(DT_2,DT_1, by=c("grid_id","YM"))
DT = DT[year != 2018]
save(DT, file = "./Data/MultiEnsemble/DT.RData")

