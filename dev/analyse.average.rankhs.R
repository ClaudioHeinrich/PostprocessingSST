 

########## analyse multivariate rank histograms ############

# get data


rm(list = ls())



setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO/lv" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")



load(file = paste0(save_dir,"setup.RData"))



load(paste0(PCA_dir,"fc_mc.RData"))
load(paste0(PCA_dir,"fc_ac.RData"))
load(paste0(GS_dir,"fc.RData"))
load(paste0(ECC_dir,"fc.RData"))

# just keep what you need


rm(DT)

PCA_fc_ac = PCA_fc_ac[,-c(paste0('fc',51:500))]
PCA_fc_mc = PCA_fc_mc[,-c(paste0('fc',51:500))]
GS_fc = GS_fc[,-c(paste0('fc',51:500))]

rks_PCA_ac

# show average rank histograms

brks = seq(1,501,length.out = 11)
brks_ECC = 0.5:10.5

par(mfrow = c(2,2),oma = c(0,0,2,0))

hist(rks_PCA_mc[,av.rk.obs],main = 'PCA_mc',breaks = brks)
hist(rks_PCA_ac[,av.rk.obs],main = 'PCA_ac',breaks = brks)
hist(rks_GS[,av.rk.obs],main = 'GS',breaks = brks)
hist(rks_ECC[,av.rk.obs],main = 'ECC', breaks = brks_ECC)
title('average rank histograms by time', outer = TRUE)




par(mfrow = c(2,2),oma = c(0,0,2,0))

plot(rks_PCA_ac[,.(YM/12,av.rk.obs)],main = 'PCA_ac',type = 'l')
plot(rks_PCA_mc[,.(YM/12,av.rk.obs)],main = 'PCA_mc',type = 'l')
plot(rks_GS[,.(YM/12,av.rk.obs)],main = 'GS',type = 'l')
plot(rks_ECC[,.(YM/12,av.rk.obs)],main = 'ECC',type = 'l')
title('average ranks by time', outer = TRUE)



plot(rks_PCA_ac[,.(YM/12,bd.rk.obs)],main = 'PCA_ac')
plot(rks_PCA_mc[,.(YM/12,bd.rk.obs)],main = 'PCA_mc')
plot(rks_GS[,.(YM/12,bd.rk.obs)],main = 'GS')
plot(rks_ECC[,.(YM/12,bd.rk.obs)],main = 'ECC')
title('band depth ranks by time', outer = TRUE)

## have a closer look at a particularly degenerate month to see what is going on ##

y = 2015
m = 3
ym = 12 * y + m

# focus on PCA_mc first



fc_sm = PCA_fc_mc[YM == ym,]
rk_sm = rks_PCA_mc[YM == ym]

# plot residuals

for(i in 1:5)
{
  plot_diagnostic(fc_sm[,.(Lon,Lat,SST_bar - SST_hat)],mn = "obs",rr = c(-3,3))
  plot_diagnostic(fc_sm[,.(Lon,Lat,.SD - SST_hat),.SDcols = paste0('fc',3*(i-1)+1)],mn = paste0("fc ",3*(i-1)+1),rr = c(-3,3))
  plot_diagnostic(fc_sm[,.(Lon,Lat,.SD - SST_hat),.SDcols = paste0('fc',3*(i-1)+2)],mn = paste0("fc ",3*(i-1)+2),rr = c(-3,3))
  plot_diagnostic(fc_sm[,.(Lon,Lat,.SD - SST_hat),.SDcols = paste0('fc',3*(i-1)+3)],mn = paste0("fc ",3*(i-1)+3),rr = c(-3,3))
 }

# plot indicator of simulated residual >= observed residual

for(i in 1:3)
{
  plot_diagnostic(fc_sm[,.(Lon,Lat,SST_bar >= SST_hat )],mn = "obs",rr = c(-3,3))
  plot_diagnostic(fc_sm[,.(Lon,Lat,.SD <= SST_bar),.SDcols = paste0('fc',3*(i-1)+1)],mn = paste0("fc ",3*(i-1)+1),rr = c(-3,3))
  plot_diagnostic(fc_sm[,.(Lon,Lat,.SD <= SST_bar),.SDcols = paste0('fc',3*(i-1)+2)],mn = paste0("fc ",3*(i-1)+2),rr = c(-3,3))
  plot_diagnostic(fc_sm[,.(Lon,Lat,.SD <= SST_bar),.SDcols = paste0('fc',3*(i-1)+3)],mn = paste0("fc ",3*(i-1)+3),rr = c(-3,3))
}


# generate for comparison an uncorrelated forecast with the same marginal standard deviation

fc_uc = fc_sm[,1:12]
test = PCA_fc_mc[YM == ym+12,SD_hat] - fc_sm[,SD_hat]
for(i in 1:50)
{
    fc = data.table(  suppressWarnings(rnorm(n = fc_sm[,.N],mean = fc_sm[,SST_hat],sd = fc_sm[,SD_hat])))
  setnames(fc,paste0('fc',i))
  fc_uc = data.table(fc_uc,fc)
}


# plot indicator of simulated residual >= observed residual for uncorrelated variables

for(i in 1:3)
{
  plot_diagnostic(fc_uc[,.(Lon,Lat,SST_bar >= SST_hat )],mn = "obs",rr = c(-3,3))
  plot_diagnostic(fc_uc[,.(Lon,Lat,.SD <= SST_bar),.SDcols = paste0('fc',3*(i-1)+1)],mn = paste0("fc ",3*(i-1)+1),rr = c(-3,3))
  plot_diagnostic(fc_uc[,.(Lon,Lat,.SD <= SST_bar),.SDcols = paste0('fc',3*(i-1)+2)],mn = paste0("fc ",3*(i-1)+2),rr = c(-3,3))
  plot_diagnostic(fc_uc[,.(Lon,Lat,.SD <= SST_bar),.SDcols = paste0('fc',3*(i-1)+3)],mn = paste0("fc ",3*(i-1)+3),rr = c(-3,3))
}

# get univariate ranks

univ_rks = matrixStats::rowRanks(as.matrix(fc_sm[,.SD,.SDcols = c('SST_bar',paste0('fc',1:50))]))
univ_rks_dt = data.table(univ_rks)
setnames(univ_rks_dt,c('rk.obs',paste0('rk.fc',1:50)))
univ_rks_dt = data.table(fc_sm[,1:12],univ_rks_dt)

par(mfrow = c(1,1))

plot_diagnostic(univ_rks_dt,'rk.obs',rr = c(1,51),mn = 'univariate ranks for PCA forecast')

# do the same for decorrelated forecast

univ_rks_uc = matrixStats::rowRanks(as.matrix(fc_uc[,.SD,.SDcols = c('SST_bar',paste0('fc',1:50))]))
univ_rks_uc = data.table(univ_rks_uc)
setnames(univ_rks_uc,c('rk.obs',paste0('rk.fc',1:50)))
univ_rks_uc = data.table(fc_sm[,1:12],univ_rks_uc)

par(mfrow = c(1,1))

plot_diagnostic(univ_rks_uc,'rk.obs',rr = c(1,51),mn = 'univariate ranks for uncorrelated forecast')


### repeat this for many months ###


Y = 2001:2016

par(mfcol = c(2,2), oma = c(0,0,2,0))

for ( y in Y)
{
m = 4
ym = 12 * y + m

# focus on PCA_mc first



fc_sm = PCA_fc_mc[YM == ym,]
rk_sm = rks_PCA_mc[YM == ym]

# generate for comparison an uncorrelated forecast with the same marginal standard deviation

fc_uc = fc_sm[,1:12]
test = PCA_fc_mc[YM == ym+12,SD_hat] - fc_sm[,SD_hat]
for(i in 1:50)
{
  fc = data.table(  suppressWarnings(rnorm(n = fc_sm[,.N],mean = fc_sm[,SST_hat],sd = fc_sm[,SD_hat])))
  setnames(fc,paste0('fc',i))
  fc_uc = data.table(fc_uc,fc)
}


# get univariate ranks

univ_rks = matrixStats::rowRanks(as.matrix(fc_sm[,.SD,.SDcols = c('SST_bar',paste0('fc',1:50))]))
univ_rks_dt = data.table(univ_rks)
setnames(univ_rks_dt,c('rk.obs',paste0('rk.fc',1:50)))
univ_rks_dt = data.table(fc_sm[,1:12],univ_rks_dt)



plot_diagnostic(univ_rks_dt,'rk.obs',rr = c(1,51),mn =  paste0('PCA_mc, ',m,'/',y))

title('univariate ranks', outer = TRUE)

# get histogram

brks = seq(1,51,length.out = 11)

mean_rank = univ_rks_dt[,mean(rk.obs,na.rm = TRUE)]

hist(univ_rks_dt[,rk.obs],breaks = brks,main = paste0('average rank = ',round(mean_rank,1)))


# do the same for decorrelated forecast

univ_rks_uc = matrixStats::rowRanks(as.matrix(fc_uc[,.SD,.SDcols = c('SST_bar',paste0('fc',1:50))]))
univ_rks_uc = data.table(univ_rks_uc)
setnames(univ_rks_uc,c('rk.obs',paste0('rk.fc',1:50)))
univ_rks_uc = data.table(fc_sm[,1:12],univ_rks_uc)


plot_diagnostic(univ_rks_uc,'rk.obs',rr = c(1,51),mn = paste0('uncorrelated, ',m,'/',y))

# get histogram

brks = seq(1,51,length.out = 11)

mean_rank = univ_rks_uc[,mean(rk.obs,na.rm = TRUE)]

avg_rank = avg.rank(as.matrix(fc_uc[,.SD,.SDcols = c('SST_bar',paste0('fc',1:50))]))

hist(univ_rks_uc[,rk.obs],breaks = brks,main = paste0('mean rank = ',round(mean_rank,1)))

}
