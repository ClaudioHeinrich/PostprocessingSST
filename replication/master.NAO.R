

#### This Script is currently stripped for scraps and NOT operational ###


### example plots ###

ex_years = 2002:2004
ex_month = 3

ex_PCs = 25

DT_PCA = forecast_PCA_new(dt = DT,y = ex_years, m = ex_month, n= 1, PCA_depth = ex_PCs,saveorgo = FALSE, cov_dir = PCA_dir)



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
  
  plot_diagnostic(DT_PCA[year == y,.(Lon,Lat,no1PC25)],
                  mn = paste0("PCA noise 0",ex_month," / ",y),
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("pca_no",y))
  
  
}


# show principle components:

PCs = 6

prin_comp_dt = get_PCs(dt = DT, y = ex_years, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)

#without marginal correction

for(y in ex_years){
  rr_PC_raw = range(prin_comp_dt[year == y,.SD,.SDcols = c(paste0("PC",1:PCs))],na.rm = TRUE)
  rr_PC_raw = c(-max(abs(rr_PC_raw)), max(abs(rr_PC_raw)))
  
  for(d in 1:PCs)
  {
    plot_diagnostic(prin_comp_dt[year == y,.SD,.SDcols = c("Lon","Lat",paste0("PC",d))],
                    rr = rr_PC_raw,
                    mn = paste0("PC ",d,", 0",ex_month," / ",y),
                    save_pdf = TRUE,
                    save_dir = paste0(plot_dir,"/"),
                    file_name = paste0("PC",d,"_y",y,"_raw"))
    
  }
}

#with marginal correction:

for(y in ex_years){
  rr_PC_raw = range(prin_comp_dt[year == y,.SD,.SDcols = c(paste0("PC_marcor_",1:PCs))],na.rm = TRUE)
  rr_PC_raw = c(-max(abs(rr_PC_raw)), max(abs(rr_PC_raw)))
  
  for(d in 1:PCs)
  {
    plot_diagnostic(prin_comp_dt[year == y,.SD,.SDcols = c("Lon","Lat",paste0("PC_marcor_",d))],
                    rr = rr_PC_raw,
                    mn = paste0("PC ",d,", 0",ex_month," / ",y),
                    save_pdf = TRUE,
                    save_dir = paste0(plot_dir,"/"),
                    file_name = paste0("PC",d,"_y",y,"_mc"))
    
  }
}

# plotting marginal SDs

ex_years_sd = 2002:2007
prin_comp_dt = get_PCs(dt = DT, y = ex_years_sd, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)

for(y in ex_years_sd){
  rr = range(prin_comp_dt[,SD_hat],na.rm = TRUE)
  rr = c(0,max(rr))
  plot_diagnostic(prin_comp_dt[year == y,.(Lon,Lat,SD_hat)],
                  mn = paste0("SD 0",ex_month," / ",y),
                  rr = rr,
                  col_scheme = "wr",
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("SD_hat",y))
  
}

# plotting difference of SDs

ex_years_sd = 2002:2007
prin_comp_dt = get_PCs(dt = DT, y = ex_years_sd, m = ex_month, PCA_depth = PCs, cov_dir = PCA_dir)

y = min(ex_years_sd)

rr = range(prin_comp_dt[,SD_hat],na.rm = TRUE)
rr = c(0,max(rr))
plot_diagnostic(prin_comp_dt[year == y,.(Lon,Lat,SD_hat)],
                mn = paste0("SD 0",ex_month," / ",y),
                rr = rr,
                col_scheme = "wr",
                save_pdf = TRUE,
                save_dir = paste0(plot_dir,"/"),
                file_name = paste0("SD_hat_diff",y))


for(y in ex_years_sd[2:length(ex_years_sd)]){
  rr = range(prin_comp_dt[,SD_hat],na.rm = TRUE)
  rr = c(0,max(rr))
  plot_diagnostic(prin_comp_dt[year == y,.(Lon,Lat,SD_hat)],
                  mn = paste0("SD 0",ex_month," / ",y),
                  rr = rr,
                  col_scheme = "wr",
                  save_pdf = TRUE,
                  save_dir = paste0(plot_dir,"/"),
                  file_name = paste0("SD_hat",y))
  
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



