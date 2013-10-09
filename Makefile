CUDA_FLAGS=-arch sm_35

bin/Color-WIA:gpu/lb2_Color_wia_mpi.o lib/libcuColor.a lib/libcuD3Q19.a lib/libcuD3Q7.a lib/libcuExtras.a
	mkdir -p bin
	CC -O3 -o bin/Color-WIA gpu/lb2_Color_wia_mpi.o -lcuColor -lcuD3Q19 -lcuD3Q7 -lcuExtras -Llib -Iinclude


bin/Color-WIA-CBUB:gpu/lb2_Color_wia_mpi_cbub.o lib/libcuColor.a lib/libcuD3Q19.a lib/libcuD3Q7.a lib/libcuExtras.a
	mkdir -p bin
	CC -O3 -D CBUB -o bin/Color-WIA-CBUB gpu/lb2_Color_wia_mpi_cbub.o -lcuColor -lcuD3Q19 -lcuD3Q7 -lcuExtras -Llib -Iinclude


bin/ColorLBM:gpu/lb2_Color_mpi.cpp lib/libcuColor.a lib/libcuD3Q19.a lib/libcuD3Q7.a lib/libcuExtras.a
	mkdir -p bin
	CC -O3 -o bin/ColorLBM gpu/lb2_Color_mpi.cpp -lcuColor -lcuD3Q19 -lcuD3Q7 -lcuExtras -Llib

#bin/gpuMRT:gpu/lb1_MRT.cu lib/libcuMRT.a lib/libcuD3Q19.a
#	mkdir -p bin
#	nvcc -O3 -o bin/gpuMRT $(CUDA_FLAGS) gpu/lb1_MRT.cu -lcuMRT -lcuD3Q19 -Llib

#bin/gpuColor:gpu/lb2_Color.cu lib/libcuColor.a lib/libcuD3Q19.a
#	mkdir -p bin
#	nvcc -o bin/gpuColor $(CUDA_FLAGS) gpu/lb2_Color.cu -lcuColor -lcuD3Q19 -Llib

gpu/lb2_Color_wia_mpi.o:gpu/lb2_Color_wia_mpi.cpp
	CC -c -o gpu/lb2_Color_wia_mpi.o gpu/lb2_Color_wia_mpi.cpp -Iinclude
 

gpu/lb2_Color_wia_mpi_cbub.o:gpu/lb2_Color_wia_mpi.cpp
	CC -c -DCBUB -o gpu/lb2_Color_wia_mpi_cbub.o gpu/lb2_Color_wia_mpi.cpp -Iinclude

#pmmc/pmmc.o:pmmc/pmmc.cpp
#	CC -c -o pmmc/pmmc.o pmmc/pmmc.cpp -Iinclude

#include/Array.o:pmmc/Array.cpp
#	CC -c -o include/Array.o include/Array.cpp -Iinclude

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
	rm gpu/*.o
	rm bin/*
	rm lib/*
