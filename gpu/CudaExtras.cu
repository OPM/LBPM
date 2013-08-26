// Basic cuda functions callable from C/C++ code
#include <cuda.h>

extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size){
	cudaMalloc(address,size);
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
