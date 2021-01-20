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

extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size){
	cudaMalloc(address,size);
    cudaMemset(*address,0,size);
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
/*
#if __CUDA_ARCH__ < 600 
__device__ double atomicAdd(double* address, double val) { 
unsigned long long int* address_as_ull = (unsigned long long int*)address; unsigned long long int old = *address_as_ull, assumed;
 do { 
      assumed = old; 
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed))); 
      // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) 
    } 
    while (assumed != old); return __longlong_as_double(old); 
} 

#endif
*/