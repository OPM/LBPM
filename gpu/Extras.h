extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size);

extern "C" void dvc_CopyToDevice(void* dest, void* source, size_t size);

extern "C" void dvc_CopyToHost(void* dest, void* source, size_t size);

extern "C" void dvc_Barrier();
