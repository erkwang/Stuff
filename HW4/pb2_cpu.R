#STA 250 HOMEWORK 4
#YICHUAN WANG
#PROBLEM 2

#Probit MCMC with cpu

#load necessary packages
library(coda)
library(mvtnorm)
library(MASS)
library(msm)

#CPU function for Probit MCMC
probit_mcmc_cpu = function(
  y,           # vector of length n
  X,           # (n x p) design matrix
  beta_0,      # (p x 1) prior mean
  Sigma_0_inv, # (p x p) prior precision
  niter,       # number of post burnin iterations
  burnin      # number of burnin iterations
){
  #create matrix to store MCMC results for beta
  beta = matrix(0, nrow = length(beta_0), ncol = burnin + niter+1)
  beta[,1] = beta_0
  #posterior variance matrix for beta's
  Sigma.post = ginv(Sigma_0_inv + t(X) %*% X)
  #pre-allocate vector z
  z = numeric(length(y))
  for (i in 2:ncol(beta)){
    #sample z_k's
    for (k in 1:length(y)) {
      if (y[k] == 0) {
        z[k] = rtnorm(1, mean = X[k,] %*% beta[,i-1], upper = 0)
      }
      else {
        z[k] = rtnorm(1, mean = X[k,] %*% beta[,i-1], lower = 0)
      }
    }
    beta[,i] = rmvnorm(n=1, mean = Sigma.post %*% (Sigma_0_inv %*% beta_0+t(X) %*% z),
                       sigma=Sigma.post)
  }
  return(as.mcmc(t(beta[,(burnin+2):ncol(beta)])))
}

cpu_time = list()
estimates = list()

for (j in 1:5) {
dat = read.table(sprintf("~/STA250/Stuff/HW4/data_0%d.txt",j),header = TRUE)
y = dat$y
X = as.matrix(dat[,-1])

cpu_time[[j]] = system.time(expr={
  estimates[[j]] = probit_mcmc_cpu(y, X, rep(0, ncol(X)), diag(rep(1, ncol(X))), 500, 100)
})
}
save(cpu_time, file = "./cpu_time.rda")
save(estimates, file = "./estimates.rda")

