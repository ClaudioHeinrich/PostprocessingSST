rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)


load(file = paste0("~/PostClimDataNoBackup/SFE/Derived/dt_reduced_PCA.RData"))


setup_PCA = function(m=7,
                     max_PCA_depth = 500,
                     cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                     data.dir = "~/PostClimDataNoBackup/SFE/Derived/"){
  
  if(! exists("dt_reduced")) { load(file = paste0(data.dir,"dt_reduced_PCA.RData"))
      dt_reduced <<- dt_reduced
      } print("loading complete")
  
  
  #find land grid ids:
  
  land_grid_id <<- dt_reduced[variable == "SST1" & 
                                year == min(year) & 
                                month == min(m) & 
                                is.na(Bias_Est),
                              .(Lon,Lat,grid_id,month,year)]
  print("finding land complete")

  #get covariance matrices
  
  for(mon in m){
    load(file = paste0(cov.dir,"CovOU_mon",mon,".RData"))  
    load(file = paste0(cov.dir,"Cov_FU_mon",mon,".RData"))  
    
    assign(paste0("A",mon),
           matrix(c(cov_for_unc,cov_obs_unc), nrow = dim(cov_obs_unc)[1]),
           envir = globalenv()
           )
    
    assign(paste0("PCA",mon),
           irlba(crossprod(eval(parse(text = paste0 ("A",mon)))), nv = max_PCA_depth),
           envir = globalenv())
    print(paste0("Month ",mon," complete"))
  
  }
}


forecast_PCA = function(y = 1999, 
                        m = 7, 
                        PCA_depth = 5,  #accepts 0, then the mean of the bias corrected 
                                        #ensemble is returned. If saveorgo = TRUE, PCA_depth can be
                                        #a vector
                        save.dir="./Data/PostClim/SFE/Derived/PCA",
                        cov.dir = "~/PostClimDataNoBackup/SFE/PCACov/",
                        data.dir = "~/PostClimDataNoBackup/SFE/Derived/",
                        showEVs = FALSE, #If TRUE, the eigenspaces are plotted
                        saveorgo = TRUE,
                        truncate = TRUE,
                        max_PCA_depth = 500) {
  
  # Check whether setup_PCA has been run
  
  if(! exists("land_grid_id")) setup_PCA(m = m, max_PCA_depth = max_PCA_depth)
  print("setup complete")
  
  
  fc <- na.omit( dt_reduced[variable == "SST1" 
                            & year %in% y 
                            & month %in% m,
                            .(Lon,Lat,grid_id,month,year,SST_hat)])
  
  
  
  print("creating fc completed")
  
  
  #----- generate noise ------------
  
  for (d in PCA_depth){
    print(d)
    no <- c()
    
    for(mon in m){
      A <- eval(parse(text = paste0("A",mon)))
      PCA <- eval(parse(text = paste0("PCA",mon)))
        
      for(year in y){
        
          if(!showEVs) a <- A %*% PCA$v  %*% c(rnorm(d),rep(0,max_PCA_depth-d))
          if(showEVs) a <- A %*% PCA$v  %*% c(rep(1,d),rep(0,max_PCA_depth-d))
          no <- c(no, a)
        
          if(!showEVs)  {
            fc = fc[,"noise":= no]
            fc=fc[,"forecast":=SST_hat+noise]
          }
          
          if(showEVs)   {
            fc = fc[year == min(year),"forecast":= abs(no)]
            #the plot function plots the forecast, 
            #so the forecast should contain the abs value of the Eigenspaces
          }
          
          #-------- add land --------------
          
          fc_land <- list()
          fc_land[[1]] = fc
          fc_land[[2]] = land_grid_id
          fc_land = rbindlist(fc_land, fill = TRUE)
          fc_land = fc_land[order(fc_land[["grid_id"]])]
          
          #-------- truncate negative temperatures
          
          if(truncate & !showEVs) fc_land[forecast < -1.62, forecast := -1.62]
          if(truncate & showEVs) fc_land[SST_hat < -1.61, forecast := 0]
          
          #-------- save -------
          
          if(saveorgo & !showEVs){ 
            save(fc_land, file = paste0(save.dir,"/fc_",d,"pc_",y,"_",m,".RData"))
          }
          if(saveorgo & showEVs){ 
            save(fc_land, file = paste0(save.dir,"/PCA_",d,"evs.RData"))
          }
          
        }
        
      }
    }
  
  
  return(fc_land)
}



vec1 = c(1:5,10,15,25,50,100,200)
for(obs in c(1,2,"mean")){
print(paste0("Obs = ",obs))
  for(d in vec1){
  print(d)
  forecast_PCA(dt_reduced = dt_reduced,y=1999, m=7, PCA_depth = d)
  plot_residuals(Y=1999,M=7,depth = d,obs_num = obs)
  }}

vec1 = c(1:5,10,15,25,50,100,200,500)
forecast_PCA(y=1999, m=7, PCA_depth = vec1, showEVs = TRUE)
for(d in vec1){
  plot_EVs(M=7,depth = d)
}

  
