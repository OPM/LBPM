#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "models/GreyscaleModel.h"
//#define WRITE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */

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
			printf("********************************************************\n");
			printf("Running Greyscale Single Phase Permeability Calculation \n");
			printf("********************************************************\n");
		}
		// Initialize compute device
		int device=ScaLBL_SetDevice(rank);
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		
		ScaLBL_GreyscaleModel Greyscale(rank,nprocs,comm);
		auto filename = argv[1];
		Greyscale.ReadParams(filename);
		Greyscale.SetDomain();    // this reads in the domain 
		Greyscale.ReadInput();
		Greyscale.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		Greyscale.Initialize();   // initializing the model will set initial conditions for variables
		Greyscale.Run();	 
		//Greyscale.VelocityField();
		//Greyscale.WriteDebug();
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
