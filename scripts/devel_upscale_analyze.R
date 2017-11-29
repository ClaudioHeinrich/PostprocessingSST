rm(list = ls())

##----- Setup --------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
options(max.print = 1e3)
print_figs = FALSE
##--------------------

##------- Load --------
in.dir = "~/PostClimDataNoBackup/SFE/Derived/"
out.dir = "~/PostClimDataNoBackup/SFE/Derived/"
##---------------------

##------ Load Some Data ----
load(paste0(in.dir, "/senorge2_upscale.RData"))
load(paste0(in.dir, "/senorge2_gcfs1_map.RData"))
##-------------------------

bias_correct = function(Y,X)
{
  N = length(Y)
  E = Y - X
  X_tilde = rep(NA,N)
  X_tilde[-1] = X[-1] + tail(cumsum(E) - E,-1)/2:N
  return(X_tilde)
}

##----- Bias Correct ---
setkey(dt_senorge_upscale, GCFS1_id, month, year)
dt_senorge_upscale[,Ens_tilde:= bias_correct(temp, Ens_bar),
                   .(GCFS1_id,month)]
##----------------------

##--- Leave one out ----
dt_senorge_upscale[,Ens_leave_out:= Ens_bar + ( sum(temp - Ens_bar) - (temp - Ens_bar) )/(.N - 1),.(GCFS1_id, month)]
dt_senorge_upscale[,Climatology_leave_out:= (sum(temp) - temp)/(.N - 1),.(GCFS1_id, month)]
##--------------------

##--- Model a bit ----
BETA = list()
for(m in 1:12)
{
  B = dt_senorge_upscale[month == m]
  if(dim(B)[1] > 0)
  {
    yy = unique(B[,year])
    B_y = list()
    for(y in 3:length(yy))
    {
      M = B[year < yy[y],.(temp, Ens_tilde)]
      mod = lm("temp ~ Ens_tilde", data = M)
      B_y[[y]] = data.table(month = m, year = yy[y], alpha = mod$coef[1], beta = mod$coef[2])
    }
    BETA[[m]] = rbindlist(B_y)
  }
}
BETA_all= rbindlist(BETA)
##-----------------------

setkey(BETA_all, "year","month")
setkey(dt_senorge_upscale, "year", "month", GCFS1_id)
dt_senorge_upscale = merge(dt_senorge_upscale, BETA_all)
dt_senorge_upscale[,T_hat:=alpha + beta * Ens_tilde]
setkey(dt_senorge_upscale, GCFS1_id, month, year)
dt_senorge_upscale[,temp_climatology:= (cumsum(temp) - temp ) / (year - min(year)), .(GCFS1_id, month)]


id = 3174
m = 11
P = dt_senorge_upscale[GCFS1_id == id][month == m][,.(year, temp, T_hat, Ens_bar,Ens_tilde, temp_climatology)]
rr = range(P[,.(temp,T_hat,Ens_bar, Ens_tilde, temp_climatology)], na.rm=TRUE)
plot(P[,.(year,temp)], type="l", ylim = rr, id, main = paste0("ID: ",id," month: ",m))
points(P[,.(year,temp)], pch = 20)
lines(P[,.(year,T_hat)], col="red")
points(P[,.(year,T_hat)], col="red", pch= 20)
lines(P[,.(year,Ens_bar)],col="blue")
points(P[,.(year,Ens_bar)],col="blue", pch = 20)
lines(P[,.(year,Ens_tilde)],col="purple")
points(P[,.(year,Ens_tilde)],col="purple", pch = 20)
lines(P[,.(year, temp_climatology)], col="green")
points(P[,.(year, temp_climatology)], col="green", pch = 20)
legend("topright", lty = 1, col=c("black","red","blue","purple", "green"), legend = c("Obs","Model","Bias Correct","Raw","Climate"), pch = 20)


dt_senorge_upscale[,SE_raw:= (temp - Ens_bar)^2]
dt_senorge_upscale[,SE_tilde:= (temp - Ens_tilde)^2]
dt_senorge_upscale[,SE_climatology:= (temp - temp_climatology)^2]
dt_senorge_upscale[,SE_model := (temp - T_hat)^2]
dt_senorge_upscale[,SE_ens_out:= (temp - Ens_leave_out)^2]
dt_senorge_upscale[,SE_cl_out := (temp - Climatology_leave_out)^2]

dt_senorge_upscale[,AE_raw:= abs(temp - Ens_bar)]
dt_senorge_upscale[,AE_tilde:= abs(temp - Ens_tilde)]
dt_senorge_upscale[,AE_climatology:= abs(temp - temp_climatology)]
dt_senorge_upscale[,AE_model := abs(temp - T_hat)]

Results_MSE = dt_senorge_upscale[,lapply(.SD,mean, na.rm = TRUE),
                   .(year, month),
                   .SDcols = paste0("SE_",c("raw","tilde","climatology","model","ens_out","cl_out"))]

Results_MAE = dt_senorge_upscale[,lapply(.SD,mean, na.rm = TRUE),
                   .(year, month),
                   .SDcols = paste0("AE_",c("raw","tilde","climatology","model"))]

Results_MSE[month == 1,.(mean(SE_ens_out), mean(SE_cl_out), mean(SE_climatology,na.rm=TRUE), mean(SE_tilde)),year]
y_start = 2005
for(m in unique(Results_MSE$month))
{
  A = Results_MSE[month == m][,lapply(.

}
rr = range(Results_MSE[,!c("year","month"),with=FALSE])

plot(Results_MSE[,.(year*12 + month,SE_raw)], ylim = rr, xlab="year", ylab="MSE", type="l", col="purple")
points(Results_MSE[,.(year,SE_raw)], pch = 20, col="purple")

lines(Results_MSE[,.(year,SE_tilde)], col="blue")
points(Results_MSE[,.(year,SE_tilde)], pch = 20, col="blue")

lines(Results_MSE[,.(year,SE_climatology)], col="green")
points(Results_MSE[,.(year,SE_climatology)], pch = 20, col="green")

lines(Results_MSE[,.(year,SE_model)], col="red")
points(Results_MSE[,.(year,SE_model)], pch = 20, col="red")

legend("topleft", lty = 1, col=c("red","blue","purple", "green"), legend = c("Model","Bias Correct","Raw","Climate"), pch = 20)


##---- Diagnostic Plot -------
w_id = unique(dt_senorge_grid[,ens_id_1])
cols = 1:length(w_id)
if(print_figs){pdf("./figures/grid_ids.pdf")}else{X11()}
plot(dt_ens_grid[w_id,.(Lon,Lat)], pch = 20, col=cols, xlim = range(dt_senorge_grid[,Lon]), ylim = range(dt_senorge_grid[,Lat]))
for(j in 1:length(w_id))
{
  w_j = dt_senorge_grid[ens_id_1 == w_id[j], .(Lon,Lat)]
  points(w_j, pch =".", col=cols[j])
}
map("world", add = TRUE)
if(print_figs)dev.off()
##-----------------------------
