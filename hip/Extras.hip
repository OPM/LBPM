// Basic hip functions callable from C/C++ code
#include "hip/hip_runtime.h"
#include <stdio.h>

extern "C" int ScaLBL_SetDevice(int rank){
	int n_devices; 
	//int local_rank = atoi(getenv("MV2_COMM_WORLD_LOCAL_RANK"));
	hipGetDeviceCount(&n_devices); 
	//int device = local_rank % n_devices; 
	int device = rank % n_devices; 
	hipSetDevice(device); 
	if (rank < n_devices) printf("MPI rank=%i will use GPU ID %i / %i \n",rank,device,n_devices);
	return device;
}

extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size){
	hipMalloc(address,size);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("Error in hipMalloc: %s \n",hipGetErrorString(err));
	}	
}

extern "C" void ScaLBL_FreeDeviceMemory(void* pointer){
       hipFree(pointer);
}

extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size){
	hipMemcpy(dest,source,size,hipMemcpyHostToDevice);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
	   printf("Error in hipMemcpy (host->device): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_AllocateZeroCopy(void** address, size_t size){
	//hipMallocHost(address,size);
	hipMalloc(address,size);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
		printf("Error in hipMallocHost: %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_CopyToZeroCopy(void* dest, const void* source, size_t size){
        hipMemcpy(dest,source,size,hipMemcpyHostToDevice);
        hipError_t err = hipGetLastError();
        //memcpy(dest, source, size);

}

extern "C" void ScaLBL_CopyToHost(void* dest, const void* source, size_t size){
	hipMemcpy(dest,source,size,hipMemcpyDeviceToHost);
	hipError_t err = hipGetLastError();
	if (hipSuccess != err){
	   printf("Error in hipMemcpy (device->host): %s \n",hipGetErrorString(err));
	}
}

extern "C" void ScaLBL_DeviceBarrier(){
	hipDeviceSynchronize();
}
