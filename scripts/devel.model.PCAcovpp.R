## Plots examples for PCA
rm(list = ls())

library(SeasonalForecasting)
setwd("~/NR/SFE")
options(max.print = 1e3)



setup_PCA(y=1985:2010,oriented = TRUE)





#---- plot uncertainty things ----

mon = 7
load("~/PostClimDataNoBackup/SFE/Derived/range_sd_res.RData")
plot_range = c(0,rr[month == mon,max_sd_res])
vec1 = c(1:5)
forecast_PCA(m=mon, PCA_depth = vec1,  output_opts = "mar_sd" )
for(d in vec1){
  plot_system(M=mon, 
              type = "mar_sd",
              depth = d, 
              rr = plot_range, 
              plot_title = paste0("ssd for res, month7, k = ",d))
}


#---- plot PCs ----

mon = 7
vec2 = c(1:25)
# get range:
PCA = eval(parse(text = paste0("PCA",mon)))
plot_range = range(PCA$u %*% diag(PCA$d))
forecast_PCA(y=1999, m=mon, PCA_depth = vec2, max_PCA_depth = 100, output_opts = "PC" )
for(d in vec2){
  plot_system(M=mon,depth = d,type = "PC", rr = plot_range)
}
  
#---- plot summed PCs

mon = 7
vec2 = c(1:15)
forecast_PCA(y=1999, m=mon, PCA_depth = vec2, max_PCA_depth = 100, output_opts = "PCsum" )

# get range:
PCA = eval(parse(text = paste0("PCA",mon)))
range_temp = c()
for(d in vec2){
  if (d == 1) range_temp = c(range_temp,range(PCA$u[,1] * PCA$d[1] ))
  if (d > 1) range_temp = c(range_temp,range(PCA$u[,1:d] %*% diag(PCA$d[1:d]) %*% rep(1,d)))
}
range = c(min(range_temp),max(range_temp))

for(d in vec2){
    plot_system(M=mon,depth = d,type = "PCsum", rr = range, plot_title = paste0("scaled sum of first ",d," eigenvectors"))
}

