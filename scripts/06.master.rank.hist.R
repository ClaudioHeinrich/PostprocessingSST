
# rm(list = ls())
# 
# setwd("~/NR/SFE")
# options(max.print = 1e3)
# 
# library(PostProcessing)
# library(data.table)

name_abbr = "NAO_2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

#load(file = paste0(save_dir,"setup.RData"))
DT = load_combined_wide(data_dir = save_dir, output_name = "dt_combine_wide_bc_var.RData")


##################################################
#### compute rank histograms for PCA forecast ####
##################################################

PCA_dir = paste0(save_dir,"PCA_new/")


eval_years = 2001:2010  # validation years
months = 1:12       # validation months
MC_sample_size = 100   # how often we generate PCA noise
ens_size = 9   #size of forecast ensemble

PCvec = c(5,20,50,70)      # number of considered principal components

DT_PCA = forecast_PCA_new(dt = DT, y = eval_years, m = months, n = MC_sample_size, PCA_depth = PCvec,
                          ens_member = FALSE, saveorgo = FALSE, cov_dir = PCA_dir)

land_ids = DT_PCA[,which(is.na(Ens_bar) | is.na(SST_bar))]

DT_PCA = DT_PCA[-land_ids,]

for(PCs in PCvec){
  
  print(paste0("Consider ",PCs," principal components:"))
  
  
  fc_dt = DT_PCA[,.SD,.SDcols = paste0("fc",1:MC_sample_size,"PC",PCs)]
  
  setnames(fc_dt,paste0("fc",1:MC_sample_size))
  
  DT_RH = data.table( DT_PCA[,c("year","month","YM","SST_bar","Bias_Est",paste0("Ens",1:ens_size)), with = FALSE], fc_dt)
                      
   
  ym = unique(DT_RH[,YM])
  
  ranks.matrix = matrix(ym,nrow = length(ym),ncol = 1+3*(MC_sample_size +1)) 
  #ncol: 1 col for YM, 3 methods of ranking, for each we get ranks for observation and each Monte Carlo sample
  
  drm = dim(ranks.matrix)
  
  YM_ind = 0
  
  for(yearmonth in ym){
    
      print(paste0("YM = ",yearmonth-min(ym)+1,"/",ym[length(ym)]-min(ym)+1))  
    
    YM_ind = YM_ind + 1
    fc_obs_mat = DT_RH[YM == yearmonth,.SD,.SDcols = c("SST_bar",paste0("fc",1:MC_sample_size))]
    
    # get ranks
    ranks.matrix[YM_ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
    }
  
  names_vec = c("YM","mst.rk.obs",paste0("mst.r.",1:MC_sample_size),
                "av.rk.obs",paste0("av.r.",1:MC_sample_size),
                "bd.rk.obs",paste0("bd.rk",1:MC_sample_size))
  
  ranks = data.table(ranks.matrix)
  
  setnames(ranks, names_vec)
  
  # --- save ---
  
  save(ranks,file = paste0(PCA_dir,"/ranks_pca_em_",PCs,"pcs.Rdata"))
  
}

for(PCs in PCvec)
{
  load(file = paste0(PCA_dir,"/ranks_pca_em_",PCs,"pcs.Rdata"))
  # ---- plotting ----
  pdf(file=paste0(plot_dir,"/rks_pca_",PCs,"pcs.pdf"),width=8,height=2,points=12)
  par(mfrow=c(1,3),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  rhist_dt(ranks[,.(YM,mst.rk.obs)], ens_size = MC_sample_size, breaks = 15,  hist_xlab = "minimum spanning tree")
  rhist_dt(ranks[,.(YM,av.rk.obs)], ens_size = MC_sample_size, breaks = 15, hist_xlab = "average")
  rhist_dt(ranks[,.(YM,bd.rk.obs)], ens_size = MC_sample_size, breaks = 15, hist_xlab = "band depth")
  
  title(paste0("RHs for forecast by ",PCs," pcs"),outer = TRUE)
  dev.off()
  
}

###################
