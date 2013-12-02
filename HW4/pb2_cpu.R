#STA 250 HOMEWORK 4
#YICHUAN WANG
#PROBLEM 2

#load necessary packages
library(coda)
library(mvtnorm)
library(MASS)

#set up truncated normal function
rtrunc = function(mu, sigma, interval){
  l = interval[1]
  u = interval[2]
  val = rnorm(1, mu, sigma)
  while (val < l | val > u) {
    val = rnorm(1, mu, sigma)
  }
  val
}

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
  for (i in 2:ncol(beta)){
    #sample z_k's
    z = numeric(length(y))
    for (k in 1:length(y)) {
      if (y[k] == 0) {
        z[k] = rtrunc(mu = X[k,] %*% beta[,i-1], sigma = 1, c(-Inf,0))
      }
      else {
        z[k] = rtrunc(mu = X[k,] %*% beta[,i-1], sigma = 1, c(0,Inf))
      }
    }
    beta[,i] = rmvnorm(n=1, mean = Sigma.post %*% (Sigma_0_inv %*% beta[,i]+t(X) %*% z),
                       sigma=Sigma.post)
  }
  return(as.mcmc(t(beta[,(burnin+2):ncol(beta)])))
}

dat = read.table("~/Desktop/SkyDrive/STA 250/Stuff/HW4/mini_data.txt",header = TRUE)
y = dat$y
X = as.matrix(dat[,-1])

mcmc_cpu_time = system.time(expr={
  betas = probit_mcmc_cpu(y, X, rep(0, ncol(X)), diag(rep(1, ncol(X))), 2000, 100)
})


