#function for obtaining blocksize and gridsize

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

#load RCUDA
library(RCUDA)
m = loadModule("rtruncnorm.ptx")
myKernel = m$rtruncnorm_kernel
values = list()
gpu_time = list()
for (k in 1:8) {
N = as.integer(10^k)
gridsize = compute_grid(N)
x = numeric(N)
gpu_time[[k]] = system.time({
values[[k]] = .cuda(myKernel, "vals" = x, "n" = N, "mu" = rep(2, N), "sigma" = rep(1, N),
	"lo" = rep(0, N), "hi" = rep(1.5, N), "rng_a" = 2L, "rng_b" = 3L, "rng_c" = 4L,
	blockDim = gridsize[[2]], gridDim = gridsize[[1]], outputs = "vals")
})
}
gpu_results = list(values, gpu_time)
save(gpu_results, file = "gpu_results.RDA")
