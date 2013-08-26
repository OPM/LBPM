CUDA_FLAGS=-arch sm_20

bin/ColorLBM:gpu/lb2_Color_mpi.cpp lib/libcuColor.a lib/libcuD3Q19.a lib/libcuD3Q7.a lib/libcuExtras.a
	mkdir -p bin
	mpicxx -O3 -o bin/ColorLBM gpu/lb2_Color_mpi.cpp -lcuColor -lcuD3Q19 -lcuD3Q7 -lcuExtras -Llib

#bin/gpuMRT:gpu/lb1_MRT.cu lib/libcuMRT.a lib/libcuD3Q19.a
#	mkdir -p bin
#	nvcc -O3 -o bin/gpuMRT $(CUDA_FLAGS) gpu/lb1_MRT.cu -lcuMRT -lcuD3Q19 -Llib

#bin/gpuColor:gpu/lb2_Color.cu lib/libcuColor.a lib/libcuD3Q19.a
#	mkdir -p bin
#	nvcc -o bin/gpuColor $(CUDA_FLAGS) gpu/lb2_Color.cu -lcuColor -lcuD3Q19 -Llib

lib/libcuExtras.a: gpu/CudaExtras.cu
	mkdir -p lib
	nvcc -lib $(CUDA_FLAGS) gpu/CudaExtras.cu -o lib/libcuExtras.a

#lib/libcuMRT.a: gpu/MRT.cu
#	mkdir -p lib
#	nvcc -lib $(CUDA_FLAGS) gpu/MRT.cu -o lib/libcuMRT.a

lib/libcuD3Q7.a: gpu/D3Q7.cu
	mkdir -p lib
	nvcc -lib $(CUDA_FLAGS) gpu/D3Q7.cu -o lib/libcuD3Q7.a

lib/libcuD3Q19.a: gpu/D3Q19.cu
	mkdir -p lib
	nvcc -lib $(CUDA_FLAGS) gpu/D3Q19.cu -o lib/libcuD3Q19.a

lib/libcuColor.a: gpu/Color.cu
	mkdir -p lib
	nvcc -lib $(CUDA_FLAGS) gpu/Color.cu -o lib/libcuColor.a

clean:
	rm bin/*
	rm lib/*