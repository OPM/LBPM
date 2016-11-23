// Basic cuda functions callable from C/C++ code
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern "C" void ScaLBL_ScaLBL_ScaLBL_AllocateDeviceMemory(void** address, size_t size){
	//cudaMalloc(address,size);
	(*address) = malloc(size);
    memset(*address,0,size);
	
	if (*address==NULL){
		printf("Memory allocation failed! \n");
	}
}

extern "C" void ScaLBL_ScaLBL_CopyToDevice(void* dest, const void* source, size_t size){
//	cudaMemcpy(dest,source,size,cudaMemcpyHostToDevice);
	memcpy(dest, source, size);
}


extern "C" void ScaLBL_CopyToHost(void* dest, const void* source, size_t size){
//	cudaMemcpy(dest,source,size,cudaMemcpyDeviceToHost);
	memcpy(dest, source, size);
}

extern "C" void ScaLBL_DeviceBarrier(){
//	cudaDeviceSynchronize();
}
