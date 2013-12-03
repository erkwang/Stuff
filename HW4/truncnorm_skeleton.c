#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{

__global__ void 
rtruncnorm_kernel(float *vals, int n, 
                  float *mu, float *sigma, 
                  float *lo, float *hi,
                  int mu_len, int sigma_len,
                  int lo_len, int hi_len,
                  int maxtries)
{
    // Usual block/thread indexing...
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
    int idx = myblock * blocksize + subthread;

    // Setup the RNG:
	curandState rng;
    curand_init(rng_a+idx*rng_b, rng_c, 0, &rng);
    
    // Sample:
	if (__device__ int isfinite(lo[idx]) && __device__ int isfinite(hi[idx])) {
		// sample from both finite
	}
	if (!(__device__ int isfinite(lo[idx]))) {
		// sample from truncated norm -b to infinity, then reverse sign
		float mu_neg = -hi[idx] - mu[idx];
		float alpha = (mu_neg + sqrt(mu_neg^2 + 4))/2;
		float expo_rand = log(1 - curand_uniform(&rng))/(-alpha);
		float z = mu_neg + expo_rand;
		if (mu_neg < alpha) float psi_z = exp(-(alpha - z)^2/2) else float psi_z = exp(-(mu_neg - alpha)^2/2 - (alpha - z)^2/2)
		float u = curand_uniform(&rng)
		while (u >= z) {
			mu_neg = -hi[idx] - mu[idx];
			alpha = (mu_neg + sqrt(mu_neg^2 + 4))/2;
			expo_rand = log(1 - curand_uniform(&rng))/(-alpha);
			z = mu_neg + expo_rand;
			if (mu_neg < alpha) psi_z = exp(-(alpha - z)^2/2) else psi_z = exp(-(mu_neg - alpha)^2/2 - (alpha - z)^2/2)
			u = curand_uniform(&rng)
		}
		val[idx] = -z
	}
	if (!(__device__ int isfinite(hi[idx]))) {
		// sample from truncated norm a to infinity
		float mu_neg = lo[idx] - mu[idx];
		float alpha = (mu_neg + sqrt(mu_neg^2 + 4))/2;
		float expo_rand = log(1 - curand_uniform(&rng))/(-alpha);
		float z = mu_neg + expo_rand;
		if (mu_neg < alpha) float psi_z = exp(-(alpha - z)^2/2) else float psi_z = exp(-(mu_neg - alpha)^2/2 - (alpha - z)^2/2)
		float u = curand_uniform(&rng)
		while (u >= z) {
			mu_neg = lo[idx] - mu[idx];
			alpha = (mu_neg + sqrt(mu_neg^2 + 4))/2;
			expo_rand = log(1 - curand_uniform(&rng))/(-alpha);
			z = mu_neg + expo_rand;
			if (mu_neg < alpha) psi_z = exp(-(alpha - z)^2/2) else psi_z = exp(-(mu_neg - alpha)^2/2 - (alpha - z)^2/2)
			u = curand_uniform(&rng)
		}
		val[idx] = z
	}
	
	
    return;
}

} // END extern "C"

