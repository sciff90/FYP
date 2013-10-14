/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. B=mexFunction(A).
 *
 * Copyright 2012 The MathWorks, Inc.
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <thrust/random.h>
/*
 * Device code
 */
__device__
unsigned int hash(unsigned int a)
{
	a = (a+0x7ed55d16) + (a<<12);
	a = (a^0xc761c23c) ^ (a>>19);
	a = (a+0x165667b1) + (a<<5);
	a = (a+0xd3a2646c) ^ (a<<9);
	a = (a+0xfd7046c5) + (a<<3);
	a = (a^0xb55a4f09) ^ (a>>16);
	return a;
}
int __device__ gridmap(int const row, int const col,int const height)
{
	return row + col*height;
}
void __global__ mcmc(double const * const u,double const * const y,
                         double * const theta,
                         int const order,int const chain_length,int const threadsPerBlock)
{
    /* Calculate the global linear index, assuming a 1-d grid. */
    int const globalIdx = blockDim.x * blockIdx.x + threadIdx.x;
    
    	int M = chain_length*threadsPerBlock;
    	if(globalIdx<M)
    	{	
    		unsigned int seed_normal = hash(globalIdx);
		thrust::default_random_engine rng_normal(seed_normal);
		thrust::random::experimental::normal_distribution<double> dist_norm(0, 1);
		
    		for(int ii=0;ii<chain_length;ii++)
    			theta[gridmap(chain_length*globalIdx+ii,1,M)] = dist_norm(rng_normal);
    	}	
    	
    
    
}

/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    /* Declare all variables.*/
    mxGPUArray const *u_m;
    mxGPUArray const *y_m;
    double const *order_m;
    double const *chain_length_m;

    mxGPUArray *theta_m;
    
    double const *u;
    double const *y;
    double *theta;    
 

    /* Choose a reasonably sized number of threads for the block. */
    int const threadsPerBlock = 256;
    int blocksPerGrid;
    
    order_m =  (double *)mxGetData(prhs[2]); 
    int order = (int)*order_m;  
    mexPrintf("order = %d\n",order);
    
    chain_length_m =  (double *)mxGetData(prhs[3]); 
    int chain_length = (int)*chain_length_m;  
    mexPrintf("chain_length = %d\n",chain_length);

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();
    
    
    u_m = mxGPUCreateFromMxArray(prhs[0]);
    y_m = mxGPUCreateFromMxArray(prhs[1]);

    u = (double const *)(mxGPUGetDataReadOnly(u_m));
    y = (double const *)(mxGPUGetDataReadOnly(y_m));

    /* Create a GPUArray to hold the result and get its underlying pointer. */
    mwSize const ndims = 2;
    mwSize const dim[] = {chain_length*threadsPerBlock,0,order,0};
    mxClassID const cid = mxDOUBLE_CLASS;
    mxComplexity const cxx = mxREAL;
    mxGPUInitialize const init0  =  MX_GPU_INITIALIZE_VALUES;    
    

    theta_m = mxGPUCreateGPUArray(ndims,dim,cid,cxx,init0);
    theta = (double *)(mxGPUGetData(theta_m));
    
    blocksPerGrid = 2;
    mcmc<<<blocksPerGrid, threadsPerBlock>>>(y,u,theta,order,chain_length,threadsPerBlock);
    
    
    

    /* Wrap the result up as a MATLAB gpuArray for return. */
    plhs[0] = mxGPUCreateMxArrayOnGPU(theta_m);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    mxGPUDestroyGPUArray(y_m);
    mxGPUDestroyGPUArray(theta_m);
}
