rm(list = ls())

##-------- Setup ---------
library(SeasonalForecasting)
setwd("~/NR/SFE/")
data.dir = "~/PostClimDataNoBackup/"
options(max.print = 1e3)
print_figs = FALSE
##------------------------

##------ Set up -------
dt = load_combined_wide()
dt[,MeanResid := SST_bar - Ens_bar]
setkey(dt,"grid_id","month","year")
##---------------------

##--- Fit Local Model ----------------
dt[,"MeanResidCum" := (cumsum(MeanResid) - MeanResid) / (year - min(year)),.(grid_id, month)]
dt[,"SST_hat_local":=Ens_bar + MeanResidCum]
for(j in 1:9)
{
  dt[,paste0("SST_hat_",j):=eval(parse(text = paste0(paste0("Ens",j)," + MeanResidCum")))]
}
##-----------------------------------

##------ Compute calibration --------
nms = c(paste0("SST",1:10), paste0("SST_hat_",1:9))
R = apply(dt[!is.na(SST_hat_local) & !is.na(SST_bar), ..nms], 1, "rank")
f = function(x){b = rep(0,19);b[x] = 1;return(b)}
B = apply(R[1:10,],2,f)
A = rowMeans(B)
##-------------------------------------

pdf("./figures/combined_calibration.pdf")

plot(A, xlab = "Bin", ylab="Frequency", pch = 20)
dev.off()
