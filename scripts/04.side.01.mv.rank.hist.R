
rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_3" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))



##################################################################
##### function for plotting rank histograms for PCA forecast #####
##################################################################


mv_rank_hist = function(dt,
                        fc_ens_size,
                        mn = "",
                        breaks = 15,
                        save_pdf = FALSE, plot_dir = "", file_name = "")
{
  ym = unique(dt[,YM])
    
  ranks_matrix = matrix(ym,nrow = length(ym),ncol = 1+3*(fc_ens_size +1)) 
    #ncol: 1 col for YM, 3 methods of ranking, for each we get ranks for observation and each Monte Carlo sample
    
  drm = dim(ranks_matrix)
    
  YM_ind = 0
    
  for(yearmonth in ym)
    {
    print(paste0("YM = ",yearmonth-min(ym)+1,"/",ym[length(ym)]-min(ym)+1))  
    
    YM_ind = YM_ind + 1
    fc_obs_mat = dt[YM == yearmonth,.SD,.SDcols = c("SST_bar",paste0("fc",1:fc_ens_size))]
      
    # get ranks
    ranks_matrix[YM_ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
    }
    
  names_vec = c("YM","mst.rk.obs",paste0("mst.r.",1:fc_ens_size),
                  "av.rk.obs",paste0("av.r.",1:fc_ens_size),
                  "bd.rk.obs",paste0("bd.rk",1:fc_ens_size))
    
  ranks = data.table(ranks_matrix)
    
  setnames(ranks, names_vec)
  
  if(save_pdf)
  {
    pdf(file=paste0(plot_dir,file_name,".pdf"),width=8,height=2,points=12)
  }
    
  par(mfrow=c(1,3),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  
  rhist_dt(ranks[,.(YM,mst.rk.obs)], ens_size = fc_ens_size, breaks = breaks,  hist_xlab = "minimum spanning tree")
  rhist_dt(ranks[,.(YM,av.rk.obs)], ens_size = fc_ens_size, breaks = breaks, hist_xlab = "average")
  rhist_dt(ranks[,.(YM,bd.rk.obs)], ens_size = fc_ens_size, breaks = breaks, hist_xlab = "band depth")
  
  title(mn,outer = TRUE)
  
  if(save_pdf)
  {
    dev.off()
  }
  
}


##################

### PCA ###


# get forecast:
load(paste0(PCA_dir,"fc.RData"))

PCA_plot_dir = paste0(plot_dir,"PCA/")
dir.create(PCA_plot_dir,showWarnings = FALSE)

mv_rank_hist(PCA_fc, fc_ens_size = fc_ens_size, mn = "PCA rank histograms",
             save_pdf = TRUE, plot_dir = PCA_plot_dir, file_name = "rank_histo_PCA")


### SE ###


# get forecast:
load(paste0(SE_dir,"fc.RData"))

SE_plot_dir = paste0(plot_dir,"SE/")
dir.create(SE_plot_dir,showWarnings = FALSE)

mv_rank_hist(SE_fc, fc_ens_size = fc_ens_size, mn = "SE rank histograms",
             save_pdf = TRUE, plot_dir = SE_plot_dir, file_name = "rank_histo_SE")


### geostationary ###

# get forecast:
load(paste0(GS_dir,"fc.RData"))

GS_plot_dir = paste0(plot_dir,"GS/")
dir.create(GS_plot_dir,showWarnings = FALSE)

mv_rank_hist(GS_fc, fc_ens_size = fc_ens_size, mn = "GS rank histograms",
             save_pdf = TRUE, plot_dir = GS_plot_dir, file_name = "rank_histo_GS")


### ECC ###

# get forecast:
load(paste0(ECC_dir,"fc.RData"))

ECC_plot_dir = paste0(plot_dir,"ECC/")
dir.create(ECC_plot_dir,showWarnings = FALSE)

mv_rank_hist(GS_fc, fc_ens_size = ens_size, breaks = 11, mn = "ECC rank histograms",
             save_pdf = TRUE, plot_dir = ECC_plot_dir, file_name = "rank_histo_ECC")
