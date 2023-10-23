
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

