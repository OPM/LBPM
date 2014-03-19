extern "C" void AllocateDeviceMemory(void** address, size_t size);

extern "C" void CopyToDevice(void* dest, void* source, size_t size);

extern "C" void CopyToHost(void* dest, void* source, size_t size);

extern "C" void DeviceBarrier();
