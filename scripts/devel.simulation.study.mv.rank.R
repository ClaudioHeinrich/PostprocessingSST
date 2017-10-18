rm(list = ls())

##------ Set up ------------
library(SeasonalForecasting)
options(max.print = 1e3)
setwd("~/NR/SFE/")
set.seed(1)
print_figs = FALSE
##--------------------------

##-- Simulation Paramters --
p_obs = 10
p_ens = 9
N = 1e4
p = p_obs + p_ens
sigma_nature = 1
sigma_obs = 1
##--------------------------

mat_rank = function(Ens,Obs)
{
  N = dim(Ens)[1]
  if(length(dim(Obs)) > 1)
  {
    p_obs = dim(Obs)[2]
    R = matrix(NA,N,p_obs)
    for(j in 1:p_obs)
    {
      T = (Obs[,j] > Ens)
      R[,j] = rowSums(T) + 1
    }
  }else{
    T = Obs > Ens
    R = rowSums(T) + 1
  }
  return(R)
}

pit_pretty = function(r)
{
  n = length(r)
  n_bins = max(r)
  freq = rep(0, n_bins)
  for(i in 1:n_bins)
  {
    freq[i] = mean(r == i)
  }
  plot(freq, type="h", ylim = c(0,max(freq)))
}
##--- Simulate Observations --
mu = rnorm(N,0,sigma_nature)
Ens = matrix(rnorm(N * p_ens,0,sigma_nature),N,p_ens)
Obs = matrix(rnorm(N * p_obs,mu,sigma_obs),N,p_obs)
Ens_jitter = Ens + matrix(rnorm(N * p_ens,0,sigma_obs),N,p_ens)
R = mat_rank(Ens,Obs)
R_mean = mat_rank(Ens, rowMeans(Obs))
R_jitter = mat_rank(Ens_jitter, Obs)
pit_pretty(R_jitter[,1])
##----------------------------

Z_ind = matrix(rnorm(N * p),N,p)
Z_over = cbind(matrix(rnorm(N * p_obs,0,1),N,p_obs),
               matrix(rnorm(N * p_ens,0,2),N,p_ens))
Z_under = cbind(matrix(rnorm(N * p_obs,0,2),N,p_obs),
               matrix(rnorm(N * p_ens,0,1),N,p_ens))
S_ind = mv_pit(Z_ind, p_obs)
S_over = mv_pit(Z_over, p_obs)
S_under = mv_pit(Z_under, p_obs)
##----------------------------

##---- Look at mean --------
Z_mean = cbind(rowMeans(Z_over[,1:p_obs]), Z_over[,-(1:p_obs)])
Z_1 = cbind(rowMeans(Z_over[,1]), Z_over[,-(1:p_obs)])
#--------------------------

##------ Plot Basic Result -------
if(print_figs){pdf("./figures/mv_pit_basic.pdf")}else{X11()}
plot(1:p, S_ind, col="blue", pch = 20, ylim =c(0,1))
points(1:p + .01, S_under, col="red", pch = 20)
points(1:p + .02, S_over, col="green", pch = 20)
legend("bottom", col=c("blue","red","green"), pch = 20, legend = c("Calibrated", "Under","Over"))
if(print_figs)dev.off()
##----------------------------------

## ************ Alternative ***********
mu = rnorm(N)
kappa = .01
sigma = .8
Z_obs = matrix(mu,N,p_obs) + matrix(rnorm(N * p_obs,0,sqrt(kappa)), N, p_obs)
Z_ens = matrix(rnorm(N * p_ens,0,sqrt(sigma)), N, p_ens)

Z_true = cbind(mu, Z_ens)
Z_product = cbind(Z_obs, Z_ens)

S = mv_pit(Z_product, p_obs)
