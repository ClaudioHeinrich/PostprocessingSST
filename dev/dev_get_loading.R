rm(list = ls())

p = 10
N = 1e3

beta = rnorm(p)
Z = rnorm(N)
X = matrix(NA,N,p)
for(j in 1:p){
    X[,j] = beta[j] * Z + rnorm(N)
}

Sigma = cov(X)
a = svd(Sigma)
L = a$u %*% diag(sqrt(a$d))

load("/home/alex/PostClimDataNoBackup/SFE/PCACov/CovFU_mon6_obs2.RData")
