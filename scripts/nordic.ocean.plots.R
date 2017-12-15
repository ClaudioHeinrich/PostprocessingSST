

data = load_combined_wide(bias = TRUE)

data[,"clim" := mean(SST_bar,na.rm = TRUE), by = .(grid_id,month)]


#--- set lon lat window ---
lon_r = c(-20,40)
lat_r = c(50,80)

lon_min = floor(lon_r[1]+0.5)-0.5
lon_max = floor(lon_r[2]+0.5)-0.5
Lons = seq(lon_min,lon_max,by = 1)  

lat_min = floor(lat_r[1]+0.5)-0.5
lat_max = floor(lat_r[2]+0.5)-0.5
Lats = seq(lat_min,lat_max,by = 1)  



#---show that climatology is questionable---

mean_temp = data[Lon %in% Lons][Lat %in% Lats, .(year, mean(SST_bar,na.rm = TRUE)), by = year]
 
pdf("./figures/Nordic Ocean/climatetrend.pdf")
plot(mean_temp[,year],mean_temp[,V2], type = 'b', main = "mean SST in the Nordic Sea", xlab = "year", ylab = "mean SST")
lines(x = c(1985,2015), y = c(5.45,6.1), col = "blue", lty = "dashed")
abline(h = 5.778, col = "darkgreen", lty = "dashed")
legend("bottomright", legend = c("trend","est. climatology"), col = c("blue","darkgreen"),lty = "dashed")
dev.off()



#--- plot forecast ---

r= data[Lon %in% Lons][Lat %in% Lats][YM %in% 24216:24219,range(SST_hat,na.rm = TRUE)]

plot_diagnostic(dt = data[year == 2017 & month == 12, 
                          .(Lon,Lat,SST_hat)],
                lons = lon_r,
                lats = lat_r,
                save.pdf = TRUE,
                rr=r,
                file.name = "bc_fc_y2017_m12",
                mn = "SST forecast 12 / 2017")

plot_diagnostic(dt = data[year == 2018 & month == 1, 
                          .(Lon,Lat,SST_hat)],
                lons = lon_r,
                lats = lat_r,
                save.pdf = TRUE,
                rr = r,
                file.name = "bc_fc_y2018_m01",
                mn = "SST forecast 01 / 2018")
plot_diagnostic(dt = data[year == 2018 & month == 2, 
                          .(Lon,Lat,SST_hat)],
                lons = c(-20,40),
                lats = c(50,80),
                save.pdf = TRUE,
                rr = r,
                file.name = "bc_fc_y2018_m02",
                mn = "SST forecast 02 / 2018")

plot_diagnostic(dt = data[year == 2018 & month == 3, 
                          .(Lon,Lat,SST_hat)],
                lons = c(-20,40),
                lats = c(50,80),
                save.pdf = TRUE,
                file.name = "bc_fc_y2018_m03",
                mn = "SST forecast 03 / 2018")

#----- plot clim. anomalies------

r= data[Lon %in% Lons][Lat %in% Lats][YM %in% 24216:24219,range(SST_hat-clim,na.rm = TRUE)]

# looks cold, even though its warm, make the scale 'fair'

r[1] = -1.69


plot_diagnostic(dt = data[year == 2017 & month == 12, 
                          .(Lon,Lat,SST_hat-clim)],
                lons = lon_r,
                lats = lat_r,
                save.pdf = TRUE,
                rr=r,
                file.name = "ano_y2017_m12",
                mn = "SST anomaly forecast 12 / 2017",
                set.white = 0)

plot_diagnostic(dt = data[year == 2018 & month == 1, 
                          .(Lon,Lat,SST_hat-clim)],
                lons = lon_r,
                lats = lat_r,
                save.pdf = TRUE,
                rr=r,
                file.name = "ano_y2018_m01",
                mn = "SST anomaly forecast 01 / 2018",
                set.white = 0)

plot_diagnostic(dt = data[year == 2018 & month == 2, 
                          .(Lon,Lat,SST_hat-clim)],
                lons = lon_r,
                lats = lat_r,
                save.pdf = TRUE,
                rr=r,
                file.name = "ano_y2018_m02",
                mn = "SST anomaly forecast 02 / 2018",
                set.white = 0)

plot_diagnostic(dt = data[year == 2018 & month == 3, 
                          .(Lon,Lat,SST_hat-clim)],
                lons = lon_r,
                lats = lat_r,
                save.pdf = TRUE,
                rr=r,
                file.name = "ano_y2018_m03",
                mn = "SST anomaly forecast 03 / 2018",
                set.white = 0)



#---- plot marginal sd ----

a = load("~/PostClimDataNoBackup/SFE/Derived/dtcombine_mr_wide_sd.RData")

months = c(1:3,12)

rr = dt_sd[month %in% months & Lon>= lon_r[1] & Lon <= lon_r[2] & Lat>= lat_r[1] & Lat <= lat_r[2], range(res_sd_by_loc, na.rm = TRUE)] 

global_uncertainty_plot(dt_sd = dt_sd, M = months,
                        file_out = "./figures/sd_res_nordic_sea",
                        lons = lon_r,
                        lats = lat_r,
                        rr=rr)



#---- plot summed PCs


vec2 = c(20)
forecast_PCA(y=1999, m=months, PCA_depth = vec2, max_PCA_depth = 100, output_opts = "PCsum" )

# get range:
PCA = eval(parse(text = paste0("PCA",mon)))
range_temp = c()
for(d in vec2){
  if (d == 1) range_temp = c(range_temp,range(PCA$u[,1] * PCA$d[1] ))
  if (d > 1) range_temp = c(range_temp,range(PCA$u[,1:d] %*% diag(PCA$d[1:d]) %*% rep(1,d)))
}
range = c(min(range_temp),max(range_temp))

for(mon in months){
for(d in vec2){
  plot_system(M=mon,
              depth = d,
              type = "PCsum", 
              rr = range, 
              lons = lon_r,
              lats = lat_r,
              plot_title = paste0("residual error covariance"),
              file_dir = "./figures/Nordic_sea_")
}
}


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





