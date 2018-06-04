rm(list = ls())

library(PostProcessing)
library(data.table)
library(fields)
setwd("~/NR/SFE/")

our_plot = function (Lon,Lat,z,
                     mn = "",
                     save_name = NULL,
                     lons = NULL,
                     lats = NULL,
                     rr = NULL, 
                     set_white = NULL,
                     col_scheme = "bwr",
                     stretch_par = NULL) 
{
    Lons = unique(Lon)
    Lats = unique(Lat)
    n_lon = length(Lons)
    n_lat = length(Lats)
    if(is.null(rr)) rr = range(z, na.rm = TRUE)
    A = matrix(z, n_lon, n_lat, byrow=TRUE)
    brk = seq(rr[1],rr[2], length = 500)
    brk.ind = round(seq(1, length(brk), length = 10))
    brk.lab = round(brk[brk.ind], 2)
    brk.at = brk[brk.ind]
    col_scheme = "bwr"
    set_white = 0.0
    zero.ind = min(which(brk > set_white))/length(brk)
    color <- designer.colors(n = length(brk) - 1, col = c("darkblue", 
                "white", "darkred"), x = c(0, zero.ind, 1))
    stretch_par = n_lat/n_lon
    if (!is.null(save_name)){
        pdf(save_name, width = 7, height = stretch_par * 7)
    }else{
        dev.new()
    }
    image.plot(Lons, Lats, A, xlab = "Longitude", ylab = "Latitude", 
               zlim = rr, xlim = range(Lons), ylim = range(Lats), breaks = brk, 
               col = color, main = mn, cex.main = 1.8, cex.lab = 1.4, 
               cex.axis = 1, axis.args = list(cex.axis = 1, at = brk.at, 
                                              label = brk.lab))
    map("world", add = TRUE)

    if (!is.null(save_name)){ 
        dev.off()
    }
}

DT = load_combined_wide()

##for(y in 2001:2010){
y = 2009
    DT_sub = DT[year < (y + 1) & month == 7 & between(Lon,-5,30) & between(Lat, 50, 80)]
    setkey(DT_sub,Lon,Lat,year)
    DT_sub[, climatology:= (cumsum(SST_bar) - SST_bar)/(year - min(year)), grid_id]
    DT_sub[, climatology:= (cumsum(SST_bar) - SST_bar)/(year - min(year)), grid_id]
    DT_sub[, observed_anomaly := SST_bar - climatology, grid_id]
    DT_sub[, forecast_anomaly := (Ens_bar - climatology)]
    DT_sub[, error := observed_anomaly - forecast_anomaly]
    DT_sub = DT_sub[!(year == min(year))]
    DT_sub[, bias_corrected := (Ens_bar - climatology) + (cumsum(error) - error) / (year - min(year)), grid_id]
    
    Lon  = DT_sub[,Lon]
    Lat = DT_sub[,Lat]
    rr = range(DT_sub[year == y, .(observed_anomaly, forecast_anomaly, bias_corrected)], na.rm=TRUE)
    rr = c(-5,5)
    our_plot(Lon, Lat,DT_sub[year == 2001,observed_anomaly],mn = paste0("Observed Anomaly 7/",y),rr = rr, save_name = "~/NR/SFE/figures/July_Observed.pdf")
    our_plot(Lon, Lat,DT_sub[year == 2001,forecast_anomaly],mn = paste0("Forecasted Anomaly 7/",y),rr = rr, save_name = "~/NR/SFE/figures/July_Forecast.pdf")
    our_plot(Lon, Lat,DT_sub[year == 2001,bias_corrected],mn = paste0("Bias Corrected 7/",y),rr = rr, save_name = "~/NR/SFE/figures/July_Bias_Corrected.pdf")

## 1. Raw Ensemble in Anomalie (x - y_bar)
## 2. Anomalie Forecast (y - y_bar)
## 3. Bias Corrected Fore (x - y_bar) + bar( (x - y ))
