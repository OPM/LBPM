#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "models/GreyscaleColorModel.h"
//#define WRITE_SURFACES

using namespace std;


int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	{
		// parallel domain size (# of sub-domains)
		int nprocx,nprocy,nprocz;
		int iproc,jproc,kproc;

		if (rank == 0){
			printf("****************************************\n");
			printf("Running Greyscale Two-Phase Calculation \n");
			printf("****************************************\n");
		}
		// Initialize compute device
		int device=ScaLBL_SetDevice(rank);
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		
		ScaLBL_GreyscaleColorModel GreyscaleColor(rank,nprocs,comm);
		auto filename = argv[1];
		GreyscaleColor.ReadParams(filename);
		GreyscaleColor.SetDomain();    // this reads in the domain 
		GreyscaleColor.ReadInput();
		GreyscaleColor.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		GreyscaleColor.Initialize();   // initializing the model will set initial conditions for variables
		GreyscaleColor.Run();	 
		GreyscaleColor.WriteDebug();
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
