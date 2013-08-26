CXX=mpicxx
FLAGS=-O3

ColorLBM-cpu:D3Q19.o D3Q7.o Color.o lb2_Color_mpi.o
	$(CXX) $(FLAGS) -o ColorLBM-cpu lb2_Color_mpi.o D3Q19.o D3Q7.o Color.o

D3Q19.o:D3Q19.cpp
	$(CXX) $(FLAGS) -c -o D3Q19.o D3Q19.cpp

D3Q7.o:D3Q7.cpp
	$(CXX) $(FLAGS) -c -o D3Q7.o D3Q7.cpp

Color.o:Color.cpp
	$(CXX) $(FLAGS) -c -o Color.o Color.cpp

lb2_Color_mpi.o:lb2_Color_mpi.cpp
	$(CXX) $(FLAGS) -c -o lb2_Color_mpi.o lb2_Color_mpi.cpp

#MRT-MPI.o:lb1_MRT_mpi.cpp
#	$(CXX) -c -o MRT-MPI.o lb1_MRT_mpi.cpp

#MRT-MPI:D3Q19.o MRT-MPI.o
#	$(CXX) -o MRT-MPI D3Q19.o MRT-MPI.o

clean:
	rm *.o
#	rm bin/*
#	rm lib/*