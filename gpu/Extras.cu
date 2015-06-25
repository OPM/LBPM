// Basic cuda functions callable from C/C++ code
#include <cuda.h>

extern "C" void AllocateDeviceMemory(void** address, size_t size){
       	cudaMalloc(address,size);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
	   printf("Error in cudaMalloc\n");
	}	
}

extern "C" void CopyToDevice(void* dest, void* source, size_t size){
	cudaMemcpy(dest,source,size,cudaMemcpyHostToDevice);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
	   printf("Error in cudaMemcpy host to device \n");
	}
}


extern "C" void CopyToHost(void* dest, void* source, size_t size){
	cudaMemcpy(dest,source,size,cudaMemcpyDeviceToHost);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
	   printf("Error in cudaMemcpy device to host \n");
	}
}

extern "C" void DeviceBarrier(){
	cudaDeviceSynchronize();
}
