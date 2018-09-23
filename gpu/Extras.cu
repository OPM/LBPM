/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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
