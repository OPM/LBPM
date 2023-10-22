// Basic cuda functions callable from C/C++ code
#include <cuda.h>
#include <stdio.h>

extern "C" int ScaLBL_SetDevice(int rank){
	int n_devices; 
	//int local_rank = atoi(getenv("MV2_COMM_WORLD_LOCAL_RANK"));
	cudaGetDeviceCount(&n_devices); 
	//int device = local_rank % n_devices; 
	int device = rank % n_devices; 
	cudaSetDevice(device); 
	if (rank < n_devices) printf("MPI rank=%i will use GPU ID %i / %i \n",rank,device,n_devices);
	return device;
}

extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size){
	cudaMalloc(address,size);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("Error in cudaMalloc: %s \n",cudaGetErrorString(err));
	}	
}

extern "C" void ScaLBL_FreeDeviceMemory(void* pointer){
       cudaFree(pointer);
}

extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size){
	cudaMemcpy(dest,source,size,cudaMemcpyHostToDevice);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
	   printf("Error in cudaMemcpy (host->device): %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_AllocateZeroCopy(void** address, size_t size){
	//cudaMallocHost(address,size);
	cudaMalloc(address,size);
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("Error in cudaMallocHost: %s \n",cudaGetErrorString(err));
	}
}

extern "C" void ScaLBL_CopyToZeroCopy(void* dest, const void* source, size_t size){
        cudaMemcpy(dest,source,size,cudaMemcpyHostToDevice);
        cudaError_t err = cudaGetLastError();
        //memcpy(dest, source, size);

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
