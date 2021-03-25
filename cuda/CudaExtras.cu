// Basic cuda functions callable from C/C++ code
#include <cuda.h>

extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size){
	cudaMalloc(address,size);
    cudaMemset(*address,0,size);
}

extern "C" void dvc_CopyToDevice(void* dest, void* source, size_t size){
	cudaMemcpy(dest,source,size,cudaMemcpyHostToDevice);
}


extern "C" void dvc_CopyToHost(void* dest, void* source, size_t size){
	cudaMemcpy(dest,source,size,cudaMemcpyDeviceToHost);
}

extern "C" void dvc_Barrier(){
	cudaDeviceSynchronize();
}
/*
#if __CUDA_ARCH__ < 600 
__device__ double atomicAdd(double* address, double val) { 
unsigned long long int* address_as_ull = (unsigned long long int*)address; unsigned long long int old = *address_as_ull, assumed;
 do { 
      assumed = old; 
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed))); 
      // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) 
    } 
    while (assumed != old); return __longlong_as_double(old); 
} 

#endif
*/