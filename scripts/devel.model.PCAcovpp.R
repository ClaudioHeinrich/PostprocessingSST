rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)


load(file = paste0("~/PostClimDataNoBackup/SFE/Derived/dt_reduced_PCA.RData"))


setup_PCA = function(m=7,
                     y = 1999,
                     max_PCA_depth = 200,
                     cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                     data.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                     obs_unc = FALSE){
  
  if(! exists("dt_reduced")) { 
      load(file = paste0(data.dir,"dt_reduced_PCA.RData"))
      dt_reduced <<- dt_reduced
      } 
  print("loading complete")
  
  
  #find land grid ids:
  
  land_grid_id <<- dt_reduced[variable == "SST1" & 
                                year == min(year) & 
                                month == min(m) & 
                                is.na(Bias_Est),
                              .(Lon,Lat,grid_id,month,year)]
  print("finding land complete")
  
  
  fc <<- na.omit( dt_reduced[variable == "SST1" 
                            & year %in% y 
                            & month %in% m,
                            .(Lon,Lat,grid_id,month,year,SST_hat)])
  
  print("creating fc complete")
  

  #get covariance matrices
  
  for(mon in m){
    if(obs_unc == TRUE){
    load(file = paste0(cov.dir,"CovOU_mon",mon,".RData"))  
    load(file = paste0(cov.dir,"Cov_FU_mon",mon,".RData"))
    }else load(file = paste0(cov.dir,"CovRes_mon",mon,".RData"))
    
    if(obs_unc == TRUE){assign(paste0("A",mon),
           matrix(c(cov_for_unc,cov_obs_unc), nrow = dim(cov_obs_unc)[1]),
           envir = globalenv()
           )
    }
    if(!obs_unc){assign(paste0("A",mon),
                        matrix(res_cov, 
                               nrow = dim(res_cov)[1]),
                        envir = globalenv()  )
    }
    
    
    assign(paste0("PCA",mon),
           irlba(eval(parse(text = paste0 ("A",mon))), nv = max_PCA_depth),
           envir = globalenv())
    print(paste0("Month ",mon," complete"))
  
  }
}


setup_PCA()


forecast_PCA = function(y = 1999, 
                        m = 7, 
                        PCA_depth = 5,  #accepts 0, then the mean of the bias corrected 
                                        #ensemble is returned. If saveorgo = TRUE, PCA_depth can be
                                        #a vector
                        save.dir="./Data/PostClim/SFE/Derived/PCA",
                        cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                        data.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                        output_opts = "Forecast", # also takes "mar_sd", then the approximate marginal 
                                                  # variance is returned, or "PC" where the dth Eigenvector
                                                  # is returned (d=PCA_depth)
                        saveorgo = TRUE,
                        truncate = TRUE,
                        max_PCA_depth = 100) {
  
  # Check whether setup_PCA has been run
  
  if(! exists("land_grid_id")) setup_PCA(m = m, max_PCA_depth = max_PCA_depth)
  print("setup complete")
  
  
  
  
  #----- generate noise ------------
  
  for (d in PCA_depth){
    print(d)
    no <- c()
    
    for(mon in m){
      A <- eval(parse(text = paste0("A",mon)))
      PCA <- eval(parse(text = paste0("PCA",mon)))
      eigen_vectors <- PCA$u[,1:d]
      sing_values <- diag(x = PCA$d[1:d], nrow = length(PCA$d[1:d]))
        
      for(year in y){
        
          if(output_opts == "forecast") a <- eigen_vectors  %*% sing_values %*% rnorm(d)
          if(output_opts == "mar_sd") a <- sqrt(( eigen_vectors %*% sing_values)^2  %*% rep(1,d))
          if(output_opts == "PC") {
            if(d == 1) a <-  PCA$d[d]*eigen_vectors else a <-  PCA$d[d]*eigen_vectors[,d]
          }
          
          no <- c(no, a)
        
          if(output_opts == "forecast")  {
            fc = fc[,"noise":= no]
            fc=fc[,"forecast":=SST_hat+noise]
          }
          
          if(output_opts == "mar_sd" | output_opts == "PC")   {
            fc = fc[year == min(year), "forecast" := no]
            #because the plotting functions plot the forecast
            }
          
          #-------- add land --------------
          
          fc_land <- list()
          fc_land[[1]] = fc
          fc_land[[2]] = land_grid_id
          fc_land = rbindlist(fc_land, fill = TRUE)
          fc_land = fc_land[order(fc_land[["grid_id"]])]
          
          #-------- truncate negative temperatures
          
          if(truncate & output_opts %in% c("forecast","PC")) fc_land[forecast < -1.62, forecast := -1.62]
          if(truncate & output_opts == "mar_sd") fc_land[SST_hat < -1.61, forecast := 0]
          
          #-------- save -------
          
          if(saveorgo & output_opts == "forecast"){ 
            save(fc_land, file = paste0(save.dir,"/fc_",d,"pc_",y,"_",m,".RData"))
          }
          if(saveorgo & output_opts == "mar_sd"){ 
            save(fc_land, file = paste0(save.dir,"/PCA_mar_sd",d,".RData"))
          }
          
          if(saveorgo & output_opts == "PC"){ 
            save(fc_land, file = paste0(save.dir,"/PCA_",d,"PC.RData"))
          }
          
        }
        
      }
    }
  
  
  return(fc_land)
}


vec1 = c(1:5,10,15,25,50,100,200)
forecast_PCA(y=1999, m=7, PCA_depth = vec1, max_PCA_depth = 100, output_opts = "mar_sd" )
for(d in vec1){
  plot_EVs(M=7,depth = d)
}


vec2 = c(1:25)
forecast_PCA(y=1999, m=7, PCA_depth = vec2, max_PCA_depth = 100, output_opts = "PC" )
for(d in vec2){
  plot_EVs(M=7,depth = d,type = "PC")
}
  
