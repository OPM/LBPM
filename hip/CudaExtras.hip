// Basic hip functions callable from C/C++ code
#include "hip/hip_runtime.h"

extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size){
	hipMalloc(address,size);
    hipMemset(*address,0,size);
}

extern "C" void dvc_CopyToDevice(void* dest, void* source, size_t size){
	hipMemcpy(dest,source,size,hipMemcpyHostToDevice);
}


extern "C" void dvc_CopyToHost(void* dest, void* source, size_t size){
	hipMemcpy(dest,source,size,hipMemcpyDeviceToHost);
}

extern "C" void dvc_Barrier(){
	hipDeviceSynchronize();
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
