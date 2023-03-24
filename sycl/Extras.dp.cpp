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
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <stdio.h>

extern "C" int ScaLBL_SetDevice(int rank){
	int n_devices; 
	//int local_rank = atoi(getenv("MV2_COMM_WORLD_LOCAL_RANK"));
        n_devices = dpct::dev_mgr::instance().device_count();
        //int device = local_rank % n_devices; 
	int device = rank % n_devices;
        /*
        DPCT1093:208: The "device" may not be the best XPU device. Adjust the
        selected device if needed.
        */
        dpct::select_device(device);
        if (rank < n_devices) printf("MPI rank=%i will use GPU ID %i / %i \n",rank,device,n_devices);
	return device;
}

extern "C" void ScaLBL_AllocateDeviceMemory(void** address, size_t size){
        *address = (void *)sycl::malloc_device(size, dpct::get_default_queue());
        /*
        DPCT1010:209: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_FreeDeviceMemory(void* pointer){
       sycl::free(pointer, dpct::get_default_queue());
}

extern "C" void ScaLBL_CopyToDevice(void* dest, const void* source, size_t size){
        dpct::get_default_queue().memcpy(dest, source, size).wait();
        /*
        DPCT1010:211: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_AllocateZeroCopy(void** address, size_t size){
	//cudaMallocHost(address,size);
        *address = (void *)sycl::malloc_device(size, dpct::get_default_queue());
        /*
        DPCT1010:213: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_CopyToZeroCopy(void* dest, const void* source, size_t size){
        dpct::get_default_queue().memcpy(dest, source, size).wait();
        /*
        DPCT1010:215: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
        //memcpy(dest, source, size);

}

extern "C" void ScaLBL_CopyToHost(void* dest, const void* source, size_t size){
        int errorStatus=0;
        // BUILD_FIX_AFTER_MIGRATION : Commenting below line , instead add try catch block. Fix TestBubbleDFH test case.
        //dpct::get_default_queue().memcpy(dest, source, size).wait();        
        /*
        DPCT1010:216: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        try {
        errorStatus=(dpct::get_default_queue().memcpy(dest, source, size).wait(),0);
        }
        catch (sycl::exception const &exc) {
        printf("Error in cudaMemcpy (device->host): %s \n",exc.what());
        }
}

extern "C" void ScaLBL_DeviceBarrier(){
        dpct::get_current_device().queues_wait_and_throw();
}
