SET( HIP_SEPERABLE_COMPILATION ON )
SET_SOURCE_FILES_PROPERTIES( BGK.hip Color.hip CudaExtras.hip D3Q19.hip D3Q7.hip dfh.hip Extras.hip MRT.hip PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1 )
HIP_ADD_LIBRARY( lbpm-hip BGK.hip Color.hip CudaExtras.hip D3Q19.hip D3Q7.hip dfh.hip Extras.hip MRT.hip SHARED HIPCC_OPTIONS ${HIP_HIPCC_OPTIONS} HCC_OPTIONS ${HIP_HCC_OPTIONS} NVCC_OPTIONS ${HIP_NVCC_OPTIONS} ${HIP_NVCC_FLAGS} )
TARGET_LINK_LIBRARIES( lbpm-hip /opt/rocm-3.3.0/lib/libhip_hcc.so )
TARGET_LINK_LIBRARIES( lbpm-wia lbpm-hip )
ADD_DEPENDENCIES( lbpm-hip copy-include )



