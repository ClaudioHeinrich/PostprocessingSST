## Where I make a netcdf file for export
rm(list = ls())

library(PostProcessing)
library(data.table)
library(ncdf4)

setwd("~/PostClimDataNoBackup/SFE/")

load("./FcNov2018/Forecast_ts.RData")
DT_fit_ts = DT_forecast
DT_fit_ts[,Lon:=lon]
DT_fit_ts[,Lat:=lat]
load("./FcNov2018/Forecast_prect.RData")
DT_fit_prect = DT_forecast
DT_fit_prect[,Lon:=lon]
DT_fit_prect[,Lat:=lat]

ncpath <- "./FcNov2018/"
ncname <- "fcnov2018"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

lon_all = unique(sort(DT_fit_ts[, Lon]))
lat_all = unique(sort(DT_fit_ts[, Lat]))
months_all = c(12, 1, 2, 3, 4)
q_print = c(10,25,33,50,67,75,90)

londim <- ncdim_def("lon", "degrees_east", as.double(lon_all))
latdim <- ncdim_def("lat", "degrees_north", as.double(lat_all))
qdim = ncdim_def("q", "quantile", as.double(q_print / 100))
monthdim = ncdim_def("month", "month", as.double(months_all))
fillvalue <- 1e32

dlname <- "2m air temperature mean forecast"
ts_mean <- ncvar_def("ts_mean",
                     "deg_C",
                     list(londim, latdim, monthdim),
                     fillvalue,
                     dlname,
                     prec="double")

ts_mean_array = array(dim = c(length(lon_all), length(lat_all), length(months_all)))
for(m in 1:length(months_all)){
    A = squarify(DT_fit_ts[month == months_all[m], .(Lon, Lat, model_mean)])
    setkey(A,"Lon","Lat")
    ts_mean_array[,,m] = matrix(A[,model_mean],nrow=length(lon_all),byrow = TRUE)
}

dlname <- "2m air temperature climatology"
ts_climatology <- ncvar_def("ts_climatology",
                            "deg_C",
                            list(londim, latdim, monthdim),
                            fillvalue,
                            dlname,
                            prec="double")

ts_climatology_array = array(dim = c(length(lon_all), length(lat_all), length(months_all)))
for(m in 1:length(months_all)){
    A = squarify(DT_fit_ts[month == months_all[m], .(Lon, Lat, obs_ts_climatology)])
    setkey(A,"Lon","Lat")
    ts_climatology_array[,,m] = matrix(A[,obs_ts_climatology],nrow = length(lon_all), byrow = TRUE)
}

dlname <- "2m air temperature forecasted anamoly"
ts_anamoly <- ncvar_def("ts_anamoly",
                        "deg_C",
                        list(londim, latdim, monthdim),
                        fillvalue,
                        dlname,
                        prec="double")

ts_anamoly_array = array(dim = c(length(lon_all), length(lat_all), length(months_all)))
for(m in 1:length(months_all)){
    A = squarify(DT_fit_ts[month == months_all[m], .(Lon, Lat, "Anamoly" = model_mean - obs_ts_climatology)])
    setkey(A,"Lon","Lat")
    ts_anamoly_array[,,m] = matrix(A[, Anamoly],nrow = length(lon_all), byrow = TRUE)
}


fillvalue <- 1e32
dlname <- "2m air temperature quantiles"
ts_q <- ncvar_def("ts_q",
                     "deg_C",
                     list(londim, latdim, monthdim, qdim),
                     fillvalue,
                     dlname,
                     prec="double")

ts_q_array = array(dim = c(length(lon_all), length(lat_all), length(months_all), length(q_print)))
##for(m in 1:length(months_all)){
##    for(q in 1:length(q_print)){
##        A = squarify(DT_fit_ts[month == months_all[m],
##                              .SD,
##                               .SDcols = c("Lon","Lat",paste0("q_hat_", q_print[q]))])
       
##        setkey(A, "Lon", "Lat")
##        B = matrix(A[, .SD, .SDcols = paste0("q_hat_", q_print[q])][[1]],
##                   nrow = length(lon_all),
##                   byrow = TRUE)
##        ts_q_array[, , m, q] = B
##    }
##}


fillvalue <- 1e32
dlname <- "2m air temperature probability below quantile"
ts_q_below <- ncvar_def("ts_q_below",
                        "deg_C",
                        list(londim, latdim, monthdim, qdim),
                        fillvalue,
                        dlname,
                        prec="double")
ts_q_below_array = array(dim = c(length(lon_all), length(lat_all), length(months_all), length(q_print)))
for(q in 1:length(q_print)){
    nm = paste0("p_below_",q)
    for(m in 1:length(months_all)){
        A = squarify(DT_fit_ts[month == months_all[m], .(Lon, Lat, get(nm))])
        ts_q_below_array[,,m,q] = matrix(A[,V3], nrow = length(lon_all), byrow = TRUE)
    }
}

DT_fit_prect = DT_fit_prect[Lat != 29.5]  ## !!! HACK !!!
dlname <- "Monthly precipitation mean forecast"
prect_mean <- ncvar_def("prect_mean",
                     "deg_C",
                     list(londim, latdim, monthdim),
                     fillvalue,
                     dlname,
                     prec="double")

prect_mean_array = array(dim = c(length(lon_all), length(lat_all), length(months_all)))
for(m in 1:length(months_all)){
    A = squarify(DT_fit_prect[month == months_all[m], .(Lon, Lat, model_mean)])
    setkey(A,"Lon","Lat")
    prect_mean_array[, , m] = matrix(A[, model_mean], nrow = length(lon_all), byrow = TRUE)
}

dlname <- "Monthly precipitation climatology"
prect_climatology <- ncvar_def("prect_climatology",
                               "deg_C",
                               list(londim, latdim, monthdim),
                               fillvalue,
                               dlname,
                               prec="double")

prect_climatology_array = array(dim = c(length(lon_all), length(lat_all), length(months_all)))
for(m in 1:length(months_all)){
    A = squarify(DT_fit_prect[month == months_all[m], .(Lon, Lat, obs_prect_climatology)])
    setkey(A,"Lon","Lat")
    prect_climatology_array[,,m] = matrix(A[, obs_prect_climatology],nrow = length(lon_all), byrow = TRUE)
}

dlname <- "Monthly precipitation forecasted anamoly"
prect_anamoly <- ncvar_def("prect_anamoly",
                           "deg_C",
                           list(londim, latdim, monthdim),
                           fillvalue,
                           dlname,
                           prec="double")

prect_anamoly_array = array(dim = c(length(lon_all), length(lat_all), length(months_all)))
for(m in 1:length(months_all)){
    A = squarify(DT_fit_prect[month == months_all[m], .(Lon, Lat, "anamoly" = model_mean - obs_prect_climatology)])
    setkey(A,"Lon","Lat")
    prect_anamoly_array[,,m] = matrix(A[,anamoly],nrow = length(lon_all), byrow = TRUE)
}

dlname <- "Monthly precipitation quantiles"
prect_q <- ncvar_def("prect_q",
                  "deg_C",
                  list(londim, latdim, monthdim, qdim),
                  fillvalue,
                  dlname,
                  prec="double")

prect_q_array = array(dim = c(length(lon_all),
                              length(lat_all),
                              length(months_all),
                              length(q_print)))

##for(m in 1:length(months_all)){
##    for(q in 1:length(q_print)){
##        A = squarify(DT_fit_prect[month == months_all[m],
##                                  .SD,
##                                  .SDcols = c("Lon", "Lat", paste0("q_hat_", q_print[q]))])
##       
##        setkey(A, "Lon", "Lat")
##        B = matrix(A[, .SD, .SDcols = paste0("q_hat_", q_print[q])][[1]],
##                   nrow = length(lon_all),
##                   byrow = TRUE)
##        prect_q_array[, , m, q] = B
##    }
##}

fillvalue <- 1e32
dlname <- "Month precipitation probability below quantile"
prect_q_below <- ncvar_def("prect_q_below",
                        "deg_C",
                        list(londim, latdim, monthdim, qdim),
                        fillvalue,
                        dlname,
                        prec="double")
prect_q_below_array = array(dim = c(length(lon_all),
                                    length(lat_all),
                                    length(months_all),
                                    length(q_print)))
for(q in 1:length(q_print)){
    nm = paste0("p_below_",q)
    for(m in 1:length(months_all)){
        A = squarify(DT_fit_prect[month == months_all[m], .(Lon, Lat, get(nm))])
        prect_q_below_array[,,m,q] = matrix(A[, V3], nrow = length(lon_all), byrow = TRUE)
    }
}

ncout <- nc_create(ncfname,
                   list(ts_mean, ts_climatology, ts_anamoly, ts_q,ts_q_below,
                        prect_mean, prect_climatology, prect_anamoly, prect_q, prect_q_below),
                   force_v4 = TRUE)

# put variables
ncvar_put(ncout, ts_mean, ts_mean_array)
ncvar_put(ncout, ts_climatology, ts_climatology_array)
ncvar_put(ncout, ts_anamoly, ts_anamoly_array)
ncvar_put(ncout, ts_q, ts_q_array)
ncvar_put(ncout, ts_q_below, ts_q_below_array)
ncvar_put(ncout, prect_mean, prect_mean_array)
ncvar_put(ncout, prect_climatology, prect_climatology_array)
ncvar_put(ncout, prect_anamoly, prect_anamoly_array)
ncvar_put(ncout, prect_q, prect_q_array)
ncvar_put(ncout, prect_q_below, prect_q_below_array)


nc_close(ncout)




