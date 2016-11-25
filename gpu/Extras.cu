// Basic cuda functions callable from C/C++ code
#include <cuda.h>
#include <stdio.h>

extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size){
       	cudaMalloc(address,size);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
	   printf("Error in cudaMalloc: %s \n",cudaGetErrorString(err));
	}	
}

extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size){
	cudaMemcpy(dest,source,size,cudaMemcpyHostToDevice);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
	   printf("Error in cudaMemcpy (host->device): %s \n",cudaGetErrorString(err));
	}
}


extern "C" void ScaLBL_CopyToHost(void* dest, const void* source, size_t size){
	cudaMemcpy(dest,source,size,cudaMemcpyDeviceToHost);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
	   printf("Error in cudaMemcpy (device->host): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_DeviceBarrier(){
	cudaDeviceSynchronize();
}
