

aa = NULL
for(j in 1:99){
    aa[j] = mod_q[[1]][[j]]$coefficients[2]
}

DT  = DT_forecast[month == 12]
DT[, pred := predict(mod_q[[1]][[50]],newdata = DT)]
plot_smooth(DT[,.(Lon, Lat, pred)], exclude_ocean = TRUE)

rr = c(-1,5.1,5)
pdf("~/Assessments.pdf")
plot_smooth(DT[,.(Lon, Lat, q_50 - obs_ts_climatology)], exclude_ocean = TRUE, mn = "Median - Mean temp", rr = c(-1.5,1.5))
plot_smooth(DT[,.(Lon, Lat, pred)], exclude_ocean = TRUE, mn = "50% Forecasted Anomaly", rr = c(-1.5,1.5))
plot_smooth(DT[,.(Lon, Lat, q_hat_50 - q_50)], exclude_ocean = TRUE, mn = "Difference")
plot(1:99/1e2, aa, xlab = "Quantile", ylab= "Regression coefficient on ECMWF", pch = 20)
lines(1:99/1e2, aa, xlab = "Quantile", ylab = "Slope on ECMWF anomaly")
dev.off()
