CXX=mpicxx
FLAGS=-O3
INC=../include

all:ColorLBM-cpu ColorLBM-CBUB

ColorLBM-cpu:Extras.o D3Q19.o D3Q7.o Color.o lb2_Color_wia_mpi.o
	$(CXX) $(FLAGS) -I$(INC) -o ColorLBM-cpu lb2_Color_wia_mpi.o D3Q19.o D3Q7.o Color.o Extras.o

ColorLBM-CBUB:Extras.o D3Q19.o D3Q7.o Color.o lb2_Color_wia_mpi-CBUB.o
	$(CXX) $(FLAGS) -I$(INC) -DCBUB -o ColorLBM-CBUB lb2_Color_wia_mpi-CBUB.o D3Q19.o D3Q7.o Color.o Extras.o

ColorLBM-BUB:Extras.o D3Q19.o D3Q7.o Color.o lb2_Color_wia_mpi_bubble.o
	$(CXX) $(FLAGS) -I$(INC) -DCBUB -o ColorLBM-BUB lb2_Color_wia_mpi_bubble.o D3Q19.o D3Q7.o Color.o Extras.o

D3Q19.o:D3Q19.cpp
	$(CXX) $(FLAGS) -c -o D3Q19.o D3Q19.cpp

D3Q7.o:D3Q7.cpp
	$(CXX) $(FLAGS) -c -o D3Q7.o D3Q7.cpp

Extras.o:Extras.cpp
	$(CXX) $(FLAGS) -c -o Extras.o Extras.cpp

Color.o:Color.cpp
	$(CXX) $(FLAGS) -c -o Color.o Color.cpp

lb2_Color_wia_mpi.o:lb2_Color_wia_mpi.cpp
	$(CXX) $(FLAGS) -I$(INC) -c -o lb2_Color_wia_mpi.o lb2_Color_wia_mpi.cpp

lb2_Color_wia_mpi-CBUB.o:lb2_Color_wia_mpi.cpp
	$(CXX) $(FLAGS) -I$(INC) -DCBUB -c -o lb2_Color_wia_mpi-CBUB.o lb2_Color_wia_mpi.cpp

lb2_Color_wia_mpi_bubble.o:lb2_Color_wia_mpi_bubble.cpp
	$(CXX) $(FLAGS) -I$(INC) -DCBUB -c -o lb2_Color_wia_mpi_bubble.o lb2_Color_wia_mpi_bubble.cpp

#MRT-MPI.o:lb1_MRT_mpi.cpp
#	$(CXX) -c -o MRT-MPI.o lb1_MRT_mpi.cpp

#MRT-MPI:D3Q19.o MRT-MPI.o
#	$(CXX) -o MRT-MPI D3Q19.o MRT-MPI.o

clean:
	rm *.o
#	rm bin/*
#	rm lib/*