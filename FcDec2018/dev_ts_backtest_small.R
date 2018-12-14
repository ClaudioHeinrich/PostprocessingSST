## Where I compare the full model to non statistical alternatives
rm(list = ls())

library(PostProcessing)
library(data.table)
library(quantreg)
setwd("~/PostClimDataNoBackup/SFE/")
path_out = "~/"
print_figs = FALSE

load("./Fc_201812/DT_merged.RData")

DT_train_y = DT_t2m
setkey(DT_train_y, Lat,Lon,month,year)

Lat_nordic = c(58,63)
Lon_nordic = c(5,11)

DT_nordic = DT_train_y[between(Lat, Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2])]
model_names = c("obs","ecmwf","dwd","cmcc","mf","ukmo")
for(n in 1:length(model_names)){
    
    nm = paste0(model_names[n], "_t2m_climatology")
    nm_var = paste0(model_names[n], "_t2m")
    DT_nordic[, eval(nm) := median(get(nm_var)), .(Lat, Lon, month)]
}

for(n in 1:length(model_names)){
    nm = paste0(model_names[n], "_t2m_anomaly")
    nm_var1 = paste0(model_names[n], "_t2m")
    nm_var2 = paste0(model_names[n], "_t2m_climatology")
    DT_nordic[, eval(nm) := get(nm_var1) - get(nm_var2)]
}

    f = paste0("obs_t2m_anomaly ~ ",paste0(model_names[-1], "_t2m_anomaly", collapse=" + "))
    mod = lm(as.formula(f), data = DT_train_y)
    DT_y[, model := obs_t2m_climatology + predict(mod, newdata = DT_y)]

    for(m in months){
        mod_month = lm(as.formula(f),
                       data = DT_train_y[month == m])
        DT_y[month == m, model_month := obs_t2m_climatology + predict(mod_month, newdata = DT_y[month == m])]
    }
    nms_anamoly = paste0(model_names[-1],"_t2m_anomaly")
    DT_y[, equal_weighted_anamoly := obs_t2m_climatology + rowMeans(.SD), .SDcols = nms_anamoly]
    return(DT_y)
}

DT_y_all = parallel::mclapply(1:length(vintage_years), "helper", mc.cores = 25, mc.silent = FALSE)

gg = function(x,y){return(sqrt(mean( (x - y)^2)))}

mod_names = c("obs_t2m_climatology", "equal_weighted_anamoly", "model", "model_month")
DT_y = rbindlist(DT_y_all)
Scores_ts = DT_y[between(Lat, Lat_nordic[1], Lat_nordic[2]) & between(Lon, Lon_nordic[1], Lon_nordic[2]),
              lapply(.SD,gg,obs_t2m),
              .SDcols = mod_names,
              month]

Scores_all = DT_y[,
                  lapply(.SD,gg,obs_t2m),
                  .SDcols = mod_names,
                  month]

Scores_grid = DT_y[,
              lapply(.SD,gg,obs_t2m),
              .SDcols = mod_names,
              .(month, Lon,Lat)]

Scores_grid[, Rel := model / obs_t2m_climatology]

X11();plot_smooth(Scores_grid[month == 1,.(Lon,Lat,Rel)])
X11();plot_smooth(Scores_grid[month == 1,.(Lon,Lat,model)])

DT_y[,model_residual := obs_t2m - model]
DT_y[,sd_err := sd(model_residual), .(Lon,Lat,month)]

X11();plot_smooth(DT_y[month == 1 & year == 1994, .(Lon, Lat, sd_err)])
