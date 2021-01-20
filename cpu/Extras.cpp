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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mm_malloc.h>

extern "C" int ScaLBL_SetDevice(int rank){
	return 0;
}

extern "C" void ScaLBL_AllocateZeroCopy(void** address, size_t size){
	//cudaMalloc(address,size);
	(*address) = _mm_malloc(size,64);
	memset(*address,0,size);
	
	if (*address==NULL){
		printf("Memory allocation failed! \n");
	}
}

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

extern "C" void ScaLBL_CopyToZeroCopy(void* dest, const void* source, size_t size){
//	cudaMemcpy(dest,source,size,cudaMemcpyDeviceToHost);
	memcpy(dest, source, size);
}

extern "C" void ScaLBL_DeviceBarrier(){
//	cudaDeviceSynchronize();
}
