#STA 250 HOMEWORK 4
#YICHUAN WANG
#PROBLEM 2

#load necessary packages
library(coda)
library(mvtnorm)
library(MASS)
library(msm)

#load RCUDA and kernel
library(RCUDA)
m = loadModule("rtruncnorm.ptx")
myKernel = m$rtruncnorm_kernel

#function for computing grid size
compute_grid <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
    # if...
    # N = 1,000,000
    # => 1954 blocks of 512 threads will suffice
    # => (62 x 32) grid, (512 x 1 x 1) blocks
    # Fix block dims:
    block_dims <- c(as.integer(sqrt_threads_per_block), as.integer(sqrt_threads_per_block), 1L)
    threads_per_block <- prod(block_dims)
    if (grid_nd==1){
      grid_d1 <- as.integer(max(1L,ceiling(N/threads_per_block)))
      grid_d2 <- 1L
    } else {
      grid_d1 <- as.integer(max(1L, floor(sqrt(N/threads_per_block))))
      grid_d2 <- as.integer(ceiling(N/(grid_d1*threads_per_block)))
    }
    grid_dims <- c(grid_d1, grid_d2, 1L)
    return(list("grid_dims"=grid_dims,"block_dims"=block_dims))
}

#GPU function for Probit MCMC
probit_mcmc_gpu = function(
  y,           # vector of length n
  X,           # (n x p) design matrix
  beta_0,      # (p x 1) prior mean
  Sigma_0_inv, # (p x p) prior precision
  niter,       # number of post burnin iterations
  burnin,      # number of burnin iterations
  block_dim, grid_dim
){
  #create matrix to store MCMC results for beta
  beta = matrix(0, nrow = length(beta_0), ncol = burnin + niter+1)
  beta[,1] = beta_0
  #posterior variance matrix for beta's
  Sigma.post = ginv(Sigma_0_inv + t(X) %*% X)
  for (i in 2:ncol(beta)){
    #sample z_k's
    #pre-allocate vector z
    N = length(y)
  	z = numeric(N)
  	mu = as.numeric(X %*% beta[,i-1])
  	lovec = rep(-Inf, N);lovec[(y == 1)] = 0
  	hivec = rep(Inf, N); hivec[(y == 0)] = 0
    z = .cuda(myKernel, "vals" = z, "n" = N, "mu" = mu, "sigma" = rep(1, N),
		"lo" = lovec, "hi" = hivec, "rng_a" = 2L, "rng_b" = 3L, "rng_c" = as.integer(2+i),
		blockDim = block_dim, gridDim = grid_dim, outputs = "vals")
    beta[,i] = rmvnorm(n=1, mean = Sigma.post %*% (Sigma_0_inv %*% beta_0+t(X) %*% z),
                       sigma=Sigma.post)
  }
  return(as.mcmc(t(beta[,(burnin+2):ncol(beta)])))
}

gpu_time = list()
gpu_estimates = list()

for (j in 1:4) {
dat = read.table(sprintf("./data_0%d.txt",j),header = TRUE)
y = dat$y
X = as.matrix(dat[,-1])
gridsize = compute_grid(length(y))

gpu_time[[j]] = system.time(expr={
  gpu_estimates[[j]] = probit_mcmc_gpu(y, X, rep(0, ncol(X)), diag(rep(1, ncol(X))), 500, 100,gridsize[[2]],gridsize[[1]])
})
}
save(gpu_time, file = "./gpu_time.rda")
save(gpu_estimates, file = "./gpu_estimates.rda")