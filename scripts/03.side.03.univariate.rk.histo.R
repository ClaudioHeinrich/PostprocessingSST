###############################################################################

#############  side script 3.3 - comparing univariate methods  ################

###############################################################################

# This script generates univariate rank histogram plots.
# 
# Files generated:
#   
# Plots: 
#   
# Requires previous run of 03.master.var.est.R with the same value of name_abbr as below.

##### setting up ######

rm(list = ls())

setwd("~/NR/SFE")
options(max.print = 1e3)

library(PostProcessing)
library(data.table)

name_abbr = "NAO_2" 

save_dir = paste0("~/PostClimDataNoBackup/SFE/Derived/", name_abbr,"/")

load(file = paste0(save_dir,"setup.RData"))



DT[,SST_bar := trc(SST_bar)]

na_loc = which(DT[,is.na(SST_bar) | is.na(Ens_bar) ])  

n=99
nbreaks = 21


means = DT[-na_loc,][year %in% validation_years & month %in% months,SST_hat]
sds = DT[-na_loc,][year %in% validation_years & month %in% months,SD_hat]
obs = DT[-na_loc,][year %in% validation_years & month %in% months,SST_bar]

fcs = trc(rnorm(n*length(means),mean = means,sd = sds))

rh_mat = matrix(c(obs,fcs),ncol = n+1)
  
rank_mat = t(apply(rh_mat, MARGIN = 1, rank, ties = 'random'))
  
pdf(file=paste0(plot_dir,"rkh.pdf"))
  hist(rank_mat[,1],
       breaks=seq(0, n+1, length.out = nbreaks),
       main = 'univariate rank histogram',
       xlab = "", ylab = "", axes=FALSE, col="gray80", border="gray60")
  abline(a=dim(rank_mat)[1]/(nbreaks-1), b=0, lty=2, col="gray30")
dev.off()


DT[,SD_hat := SD_hat_sv]

save.image(file = paste0(save_dir,"univariate.model.comparison.RData"))
