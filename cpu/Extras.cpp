// Basic cuda functions callable from C/C++ code
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size){
	//cudaMalloc(address,size);
	(*address) = malloc(size);
	
	if (*address==NULL){
		printf("Memory allocation failed! \n");
	}
}

extern "C" void dvc_CopyToDevice(void* dest, void* source, size_t size){
//	cudaMemcpy(dest,source,size,cudaMemcpyHostToDevice);
	memcpy(dest, source, size);
}


extern "C" void dvc_CopyToHost(void* dest, void* source, size_t size){
//	cudaMemcpy(dest,source,size,cudaMemcpyDeviceToHost);
	memcpy(dest, source, size);
}

extern "C" void dvc_Barrier(){
//	cudaDeviceSynchronize();
}
