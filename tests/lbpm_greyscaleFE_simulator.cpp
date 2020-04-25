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
#include "models/GreyscaleFEModel.h"
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
		
		ScaLBL_GreyscaleFEModel GreyscaleFE(rank,nprocs,comm);
		auto filename = argv[1];
		GreyscaleFE.ReadParams(filename);
		GreyscaleFE.SetDomain();    // this reads in the domain 
		GreyscaleFE.ReadInput();
		GreyscaleFE.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		GreyscaleFE.Initialize();   // initializing the model will set initial conditions for variables
		GreyscaleFE.Run();	 
		//GreyscaleFE.VelocityField();
		GreyscaleFE.WriteDebug();
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
