

#### This Script is currently stripped for scraps and NOT operational ###


### example plots ###

ex_years = 2002:2004
ex_month = 3

ex_PCs = 25

DT_PCA = forecast_PCA_new(dt = DT,y=ex_years, m = ex_month, n= 1, PCA_depth = ex_PCs,saveorgo = FALSE,cov_dir = cov_dir)

for(y in ex_years){
  plot_diagnostic(DT_PCA[year == y,.(Lon,Lat,SST_bar)],
                  mn = paste0("SST 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("SST",y))

  plot_diagnostic(DT_PCA[year == y,.(Lon,Lat,Ens1)],
                  mn = paste0("Ens1 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("raw_fc",y))
  
  plot_diagnostic(DT_PCA[year == y,.(Lon,Lat,Ens1+Bias_Est)],
                  mn = paste0("Ens1 bc 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("bc_fc",y))
  plot_diagnostic(DT_PCA[year == y,.(Lon,Lat,Ens1 + Bias_Est + rnorm(DT_PCA[year == y,.N],mean = 0,sd = DT_PCA[year == y,SD_hat]))],
                  mn = paste0("Ens1 marg. cal. 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("marcal_fc",y))
  plot_diagnostic(DT_PCA[year == y,.(Lon,Lat,fc1PC25)],
                  mn = paste0("PCA fc 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("pca_fc",y))
  
  
}

# geostationary:
DT_geostat = forecast_geostat(dt = DT,n=1,y = ex_years,m = ex_month,saveorgo = FALSE, data_dir = geostat_dir)

for(y in ex_years){
  plot_diagnostic(DT_geostat[year == y,.(Lon,Lat,fc1)],
                  mn = paste0("SST 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("fc_geostat",y))
  
 
  
}

# ECC

DT_ECC = forecast_ECC(dt = DT[year %in% ex_years & month == ex_month,],saveorgo = FALSE)

for(y in ex_years){
  plot_diagnostic(DT_ECC[year == y,.(Lon,Lat,ecc_fc1)],
                  mn = paste0("ECC forecast 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("fc_ecc",y))
  
  
  
}


############################################################


###### multivariate rank histograms ######

## Minimum spanning tree ranks 
mst.rank <- function (x) {
  l.mst <- NULL
  for(f in 1:(dim(x)[2])) {
    euc.dist <- rdist(t(x[,-f]))
    l.mst <- c(l.mst,sum(spantree(euc.dist)$dist))
  }
  x.rank <- rank(l.mst,ties="random")
  return(x.rank)
}

## Multivariate ranks 
mv.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- numeric(d[2])
  for(i in 1:d[2]) {
    x.prerank[i] <- sum(apply(x<=x[,i],2,all))
  }
  x.rank <- rank(x.prerank,ties="random")
  return(x.rank)
}

## Average ranks
avg.rank <- function(x)  {
  x.ranks <- apply(x,1,rank)
  x.preranks <- apply(x.ranks,1,mean)
  x.rank <- rank(x.preranks,ties="random")
  return(x.rank)
}

## Band depth ranks
bd.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- array(NA,dim=d)
  for(i in 1:d[1]) {
    tmp.ranks <- rank(x[i,])
    x.prerank[i,] <- (d[2] - tmp.ranks) * (tmp.ranks - 1)
  }
  x.rank <- apply(x.prerank,2,mean) + d[2] - 1
  x.rank <- rank(x.rank,ties="random")
  return(x.rank)
} 


## rank histograms for data tables

# The data table should have the key variable (most commonly YM) as first column and the ranks of the observations as second

rhist.dt <- function(B, ens_size, breaks = seq(0, ens_size + 1, length.out = min(ens_size + 1,20)), hist_xlab="", hist_ylab="", hist_ylim=NULL)
{
  hist(as.vector(B[[2]]),breaks = breaks, main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray80",border="gray60",ylim=hist_ylim)
  abline(a=length(B[[1]])/length(breaks), b=0, lty=2, col="gray30")
}


##################################################
#### compute rank histograms for PCA forecast ####
##################################################




validation_years = 2001:2010  # validation years
validation_months = 1:12       # validation months
MC_sample_size = 100   # how often we generate PCA noise
ens_size = 9   #size of forecast ensemble

PCvec = c(5,opt_num_PCs,50)      # number of considered principal components

for(PCs in PCvec){

  print(paste0("computing RHs for ",PCs," principal components."))
  no_dt = list()
  for(i in 1:MC_sample_size){
    print(paste0("generating Monte Carlo sample ",i,"/",MC_sample_size))
    no_dt[[i]] = forecast_PCA(m = validation_months, y = validation_years, PCA_depth = PCs, saveorgo = FALSE)[,noise]
  }
  no_dt = as.data.table(no_dt)
  
  setnames(no_dt,paste0("no",1:MC_sample_size))

  DT_pca = no_dt[,c("year","month","YM","SST_bar","Bias_Est",paste0("Ens",1:ens_size)) := DT[year %in% validation_years & month %in% validation_months,c("year","month","YM","SST_bar","Bias_Est",paste0("Ens",1:ens_size)),with = FALSE]]


  #choose random ensemble members (REM) and generate forecast as REM + bias + noise

  ens_mem = sample.int(ens_size,MC_sample_size,replace = TRUE)
  for(i in 1:MC_sample_size){
    dummy_dt = DT_pca[,.SD,.SDcols = c(paste0("no",i),paste0("Ens",ens_mem[i]),"Bias_Est")]
    forecast = dummy_dt[[1]] + dummy_dt[[2]] + dummy_dt[[3]]
    DT_pca = DT_pca[,paste0("fc",i) := forecast]
  }
  
  ym = unique(DT_pca[,YM])
  
  ranks.matrix = matrix(ym,nrow = length(ym),ncol = 1+3*(MC_sample_size +1)) 
  #ncol: 1 col for YM, 3 methods of ranking, for each we get ranks for observation and each Monte Carlo sample
  
  drm = dim(ranks.matrix)
  
  YM_ind = 0
  
  for(yearmonth in ym){
    print(paste0("YM = ",yearmonth,"/",ym[length(ym)]))
    YM_ind = YM_ind + 1
    fc_obs_mat = na.omit(DT_pca[YM == yearmonth,.SD,.SDcols = c("SST_bar",paste0("fc",1:MC_sample_size))])
    
    # get ranks
    ranks.matrix[YM_ind,2:drm[2]] = c(mst.rank(as.matrix(fc_obs_mat)),
                                      avg.rank(as.matrix(fc_obs_mat)),
                                      bd.rank(as.matrix(fc_obs_mat)))
    rm(fc_obs_mat)
  }
  
  names_vec = c("YM","mst.rk.obs",paste0("mst.r.",1:MC_sample_size),
                "av.rk.obs",paste0("av.r.",1:MC_sample_size),
                "bd.rk.obs",paste0("bd.rk",1:MC_sample_size))
  
  ranks = data.table(ranks.matrix)
  
  setnames(ranks, names_vec)
  
  # --- save ---
  
  save(ranks,file = paste0(save_dir,"/ranks_pca_em_",PCs,"pcs.Rdata"))
  
  # ---- plotting ----
  pdf(file=paste0(plot_dir,"/rks_pca_em_",PCs,"pcs.pdf"),width=8,height=2,points=12)
  par(mfrow=c(1,3),mex=0.5,oma = c(0,0,2.5,0),mar=c(2.5,2.5,2.5,2.5)+0.1,mgp=c(0.5,0,0))
  rhist.dt(ranks[,.(YM,mst.rk.obs)], ens_size = MC_sample_size,  hist_xlab = "minimum spanning tree")
  rhist.dt(ranks[,.(YM,av.rk.obs)], ens_size = MC_sample_size, hist_xlab = "average")
  rhist.dt(ranks[,.(YM,bd.rk.obs)], ens_size = MC_sample_size, hist_xlab = "band depth")
  
  title(paste0("RHs for forecast by ",PCs," pcs"),outer = TRUE)
  dev.off()

}

###################


###### Example plots of forecasted SST and anomalies w.r.t climatology #######

ex_depth = opt_num_PCs
ex_months = 4:9
ex_month_names = c("April","May","June","July","August","September")
ex_year = 2010

clim_years = 1985:2009 # the years to compute the climatology from

MC_sample_size = 10   # number of plots with independently generated noise

ens_size = 9   #size of forecast ensemble

PCs = opt_num_PCs      # number of considered principal components


#compute climatology
climatology = DT[year %in% clim_years, clim := mean(SST_bar),by = .(grid_id,month)][year == min(year) ,.(Lon,Lat,month,clim)]



for(m in ex_months){
  print(paste0("month = ",m))
    
    #generate noise:
    no_dt = list()
    for(i in 1:MC_sample_size){
      no_dt[[i]] = forecast_PCA(m = m, y = ex_year, PCA_depth = PCs, saveorgo = FALSE)[,.(Lon,Lat,noise), keyby = .(Lon,Lat)][,noise]
    }
    no_dt = as.data.table(no_dt)
    setnames(no_dt,paste0("no",1:MC_sample_size))
    
    DT_pca_plot = no_dt[,c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens_size)) := 
                          DT[year == ex_year & month == m, c("year","month","YM","Lat","Lon","SST_bar","Bias_Est",paste0("Ens",1:ens_size)),
                             with = FALSE]]
    DT_pca_plot[,clim := climatology[month == m, clim]]
    
    # choose random ensemble members (REM) and generate forecast as REM + bias + noise
    
    ens_mem = sample.int(ens_size,MC_sample_size,replace = TRUE)
    for(i in 1:MC_sample_size){
      dummy_dt = DT_pca_plot[,.SD,.SDcols = c(paste0("no",i),paste0("Ens",ens_mem[i]),"Bias_Est")]
      forecast = dummy_dt[[1]] + dummy_dt[[2]] + dummy_dt[[3]]
      DT_pca_plot = DT_pca_plot[,paste0("fc",i) := forecast]
    }
    
    #forecast plots:
    
    rr_sst = range(na.omit(DT_pca_plot[,.SD,.SDcols = c(paste0("fc",1:MC_sample_size))]))
    
    for(i in 1:MC_sample_size){
      plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i))))],
                      rr = rr_sst,
                      mn = paste0("SST forecast for ",ex_month_names[which(ex_months == m)]),
                      save_pdf = TRUE, 
                      save_dir = paste0(plot_dir,"/"),
                      file_name = paste0("m",m,"_fc",i),
                      stretch_par = .8)
    }
    
    #anomaly plot 
    rr_clim = range(na.omit(DT_pca_plot[,.SD - clim,.SDcols = c(paste0("fc",1:MC_sample_size))]))
    
    for(i in 1:MC_sample_size){
      plot_diagnostic(DT_pca_plot[,.(Lon,Lat,eval(parse(text = paste0("fc",i)))-clim)],
                      rr = rr_clim,
                      mn = paste0("Anomaly forecast for ",ex_month_names[which(ex_months == m)]),
                      save_pdf = TRUE, 
                      save_dir = paste0(plot_dir,"/"),
                      file_name = paste0("m",m,"_afc",i),
                      stretch_par = .8)
    }
    
  }



