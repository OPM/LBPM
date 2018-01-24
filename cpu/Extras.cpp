// Basic cuda functions callable from C/C++ code
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mm_malloc.h>

extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size){
	//cudaMalloc(address,size);
	(*address) = _mm_malloc(size,64);
	memset(*address,0,size);
	
	if (*address==NULL){
		printf("Memory allocation failed! \n");
	}
}

extern "C" void ScaLBL_FreeDeviceMemory(void* pointer){
	_mm_free(pointer);
}

extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size){
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
