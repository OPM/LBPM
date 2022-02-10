#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "common/MPI.h"
#include "models/BGKModel.h"
//#define WRITE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */

using namespace std;


int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();
	{
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Single Phase Permeability Calculation \n");
			printf("********************************************************\n");
		}
		// Initialize compute device
		int device=ScaLBL_SetDevice(rank);
        NULL_USE( device );
		ScaLBL_DeviceBarrier();
		comm.barrier();
		
		ScaLBL_BGKModel BGK(rank,nprocs,comm);
		auto filename = argv[1];
		BGK.ReadParams(filename);
		BGK.SetDomain();    // this reads in the domain 
		BGK.ReadInput();
		BGK.Create();       // creating the model will create data structure to match the pore structure and allocate variables
		BGK.Initialize();   // initializing the model will set initial conditions for variables
		BGK.Run();	 
		BGK.VelocityField();
		cout << flush;
	}
    Utilities::shutdown();
}
