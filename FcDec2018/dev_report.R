rm(list = ls())

library(PostProcessing)
library(data.table)

setwd("~/NR/SFE/")
load("~/PostClimDataNoBackup/SFE/FcNov2018/Forecast_ts.RData")
load("~/PostClimDataNoBackup/SFE/FcNov2018/Forecast_prect.RData")
print_figs = TRUE

latex_table = function(X, nms, rnd = 2){
    for(i in 1:dim(X)[1]){
        cat(nms[i])
        for(j in 1:dim(X)[2]){
            cat(" & ")
            cat(round(X[i,j], rnd))
        }
        cat("\\\\")
        cat("\n")
    }
}

A = t(Score_ts)
nms = rownames(A)
latex_table(A, nms, 2)

A = t(Score_prect[month != 4])
nms = rownames(A)
latex_table(A, nms, 2)


months = c(12,1,2,3,4)

rank1 = function(X){return(rank(X)[1])}
DT_fit_ts[,p_low := apply(.SD,1,"rank1")/100,.SDcols = c("q_low", paste0("q_",1:99))]
DT_fit_ts[,p_median := apply(.SD,1,"rank1")/100,.SDcols = c("climatology", paste0("q_",1:99))]
DT_fit_ts[,p_high := apply(.SD,1,"rank1")/100,.SDcols = c("q_high", paste0("q_",1:99))]

for(m in months){
    if(print_figs){pdf(paste0("./figures/FcNov2018_ts_", m,".pdf"))}else{X11()}
    r_m = max(abs(DT_fit_ts[month == m, model_mean]), na.rm=TRUE)
    plot_smooth(DT_fit_ts[month == m, .(Lon,Lat,model_mean)],
                exclude_ocean = TRUE,
                rr = c(-r_m,r_m),
                mn = paste0("Temp Month ",m))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_ts_anamoly_", m,".pdf"))}else{X11()}
    r_m = max(abs(DT_fit_ts[month == m, model_mean - climatology]), na.rm=TRUE)
    plot_smooth(DT_fit_ts[month == m, .(Lon,Lat,model_mean - climatology)],
                exclude_ocean = TRUE,
                rr = c(-r_m,r_m),
                mn = paste0("Temp Anamoly Month ",m))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_ts_plow_", m,".pdf"))}else{X11()}
    plot_smooth(DT_fit_ts[month == m, .(Lon,Lat,p_low)],
                exclude_ocean = TRUE,
                mn = paste0("Probability Colder than 25 percentile, Month ",m),
                rr = c(0,0.5))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_ts_pmed_", m,".pdf"))}else{X11()}
    plot_smooth(DT_fit_ts[month == m, .(Lon,Lat,p_median)],
                exclude_ocean = TRUE,
                mn = paste0("Probability Colder than 50 percentile, Month ",m),
                rr = c(.25,.75))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_ts_phigh_", m,".pdf"))}else{X11()}
    plot_smooth(DT_fit_ts[month == m, .(Lon,Lat,p_high)],
                exclude_ocean = TRUE,
                mn = paste0("Probability Colder than 75 percentile, Month ",m),
                rr = c(0.5,1))
    if(print_figs)dev.off()

}


DT_fit_prect[,p_low := apply(.SD,1,"rank1")/100,.SDcols = c("q_low", paste0("q_",1:99))]
DT_fit_prect[,p_median := apply(.SD,1,"rank1")/100,.SDcols = c("climatology", paste0("q_",1:99))]
DT_fit_prect[,p_high := apply(.SD,1,"rank1")/100,.SDcols = c("q_high", paste0("q_",1:99))]

for(m in months){
    if(print_figs){pdf(paste0("./figures/FcNov2018_prect_", m,".pdf"))}else{X11()}
    r_m = max(abs(DT_fit_prect[month == m, model_mean]), na.rm=TRUE)
    plot_smooth(DT_fit_prect[month == m, .(Lon,Lat,model_mean)],
                exclude_ocean = TRUE,
                rr = c(0,r_m),
                mn = paste0("Precip Month ",m))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_prect_anamoly_", m,".pdf"))}else{X11()}
    r_m = max(abs(DT_fit_prect[month == m, model_mean - climatology]), na.rm=TRUE)
    plot_smooth(DT_fit_prect[month == m, .(Lon,Lat,model_mean - climatology)],
                exclude_ocean = TRUE,
                rr = c(-r_m,r_m),
                mn = paste0("Precip Anamoly Month ",m))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_prect_plow_", m,".pdf"))}else{X11()}
    plot_smooth(DT_fit_prect[month == m, .(Lon,Lat,p_low)],
                exclude_ocean = TRUE,
                mn = paste0("Probability Drier than 25 percentile, Month ",m),
                rr = c(0,0.5))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_prect_pmed_", m,".pdf"))}else{X11()}
    plot_smooth(DT_fit_prect[month == m, .(Lon,Lat,p_median)],
                exclude_ocean = TRUE,
                mn = paste0("Probability Drier than 50 percentile, Month ",m),
                rr = c(.25,.75))
    if(print_figs)dev.off()

    if(print_figs){pdf(paste0("./figures/FcNov2018_prect_phigh_", m,".pdf"))}else{X11()}
    plot_smooth(DT_fit_prect[month == m, .(Lon,Lat,p_high)],
                exclude_ocean = TRUE,
                mn = paste0("Probability Drier than 75 percentile, Month ",m),
                rr = c(0.5,1))
    if(print_figs)dev.off()

}


