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
                  int rng_a, int rng_b,	int rng_c)
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
    if (idx < n) {
	if (isfinite(lo[idx]) && isfinite(hi[idx])) {
		
		// sample from both finite
		float mu_neg = (lo[idx] - mu[idx])/sigma[idx];
		float mu_pos = (hi[idx] - mu[idx])/sigma[idx];
		float z = mu_neg + curand_uniform(&rng)*(mu_pos-mu_neg);
		float psi_z = expf(-z*z/2);
		if (mu_neg > 0) float psi_z = expf(-(mu_neg*mu_neg - z*z)/2);
		if (mu_pos < 0) float psi_z = expf(-(mu_pos*mu_pos - z*z)/2);
		float u = curand_uniform(&rng);
		while (u >= psi_z) {
			z = mu_neg + curand_uniform(&rng)*(mu_pos-mu_neg);
			psi_z = expf(-z*z/2);
			if (mu_neg > 0) psi_z = expf(-(mu_neg*mu_neg - z*z)/2);
			if (mu_pos < 0) psi_z = expf(-(mu_pos*mu_pos - z*z)/2);
			u = curand_uniform(&rng);
		}
		vals[idx] = sigma[idx]*z+mu[idx];
	}
	if (!isfinite(lo[idx])) {
		// sample from truncated norm -b to infinity, then reverse sign
		float mu_neg = (-hi[idx] - mu[idx])/sigma[idx];
		float alpha = (mu_neg + sqrtf(mu_neg*mu_neg + 4))/2;
		float expo_rand = logf(1 - curand_uniform(&rng))/(-alpha);
		float z = mu_neg + expo_rand;
		float psi_z = expf(-(mu_neg - alpha)*(mu_neg - alpha)/2 - (alpha - z)*(alpha - z)/2);
		if (mu_neg < alpha) float psi_z = expf(-(alpha - z)*(alpha - z)/2);
		float u = curand_uniform(&rng);
		while (u >= psi_z) {
			expo_rand = logf(1 - curand_uniform(&rng))/(-alpha);
			z = mu_neg + expo_rand;
			psi_z = expf(-(mu_neg - alpha)*(mu_neg - alpha)/2 - (alpha - z)*(alpha - z)/2);
			if (mu_neg < alpha) psi_z = expf(-(alpha - z)*(alpha - z)/2);
			u = curand_uniform(&rng);
		}
		vals[idx] = -(sigma[idx]*z+mu[idx]);
	}
	if (!isfinite(hi[idx])) {
		// sample from truncated norm a to infinity
		float mu_neg = (lo[idx] - mu[idx])/sigma[idx];
		float alpha = (mu_neg + sqrtf(mu_neg*mu_neg + 4))/2;
		float expo_rand = logf(1 - curand_uniform(&rng))/(-alpha);
		float z = mu_neg + expo_rand;
		float psi_z = expf(-(mu_neg - alpha)*(mu_neg - alpha)/2 - (alpha - z)*(alpha - z)/2);
		if (mu_neg < alpha) float psi_z = expf(-(alpha - z)*(alpha - z)/2);
		float u = curand_uniform(&rng);
		while (u >= psi_z) {
			expo_rand = logf(1 - curand_uniform(&rng))/(-alpha);
			z = mu_neg + expo_rand;
			psi_z = expf(-(mu_neg - alpha)*(mu_neg - alpha)/2 - (alpha - z)*(alpha - z)/2);
			if (mu_neg < alpha) psi_z = expf(-(alpha - z)*(alpha - z)/2);
			u = curand_uniform(&rng);
		}
		vals[idx] = sigma[idx]*z+mu[idx];
	}
	}
	return;
}

} // END extern "C"

