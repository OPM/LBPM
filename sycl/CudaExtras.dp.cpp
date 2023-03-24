// Basic cuda functions callable from C/C++ code
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size){
        *address = (void *)sycl::malloc_device(size, dpct::get_default_queue());
    dpct::get_default_queue().memset(*address, 0, size).wait();
}

extern "C" void dvc_CopyToDevice(void* dest, void* source, size_t size){
        dpct::get_default_queue().memcpy(dest, source, size).wait();
}


extern "C" void dvc_CopyToHost(void* dest, void* source, size_t size){
        dpct::get_default_queue().memcpy(dest, source, size).wait();
}

extern "C" void dvc_Barrier(){
        dpct::get_current_device().queues_wait_and_throw();
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