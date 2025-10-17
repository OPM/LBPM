#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "common/Utilities.h"
#include "models/ColorModel.h"

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */


//*************************************************************************
// Implementation of Two-Phase Immiscible LBM using CUDA
//*************************************************************************

int main( int argc, char **argv )
{

	// Initialize
        Utilities::startup( argc, argv, true );

	{ // Limit scope so variables that contain communicators will free before MPI_Finialize

		Utilities::MPI comm( MPI_COMM_WORLD );
		int rank   = comm.getRank();
		int nprocs = comm.getSize();

		if ( rank == 0 ) {
			printf( "********************************************************\n" );
			printf( "Running Color LBM Unsteady state \n" );
			printf( "********************************************************\n" );
		}
		// Initialize compute device
		int device = ScaLBL_SetDevice( rank );
		NULL_USE( device );
		ScaLBL_DeviceBarrier();
		comm.barrier();

		PROFILE_ENABLE( 1 );
		PROFILE_SYNCHRONIZE();
		PROFILE_START( "Main" );
		Utilities::setErrorHandlers();

		auto filename = argv[1];
		ScaLBL_ColorModel ColorModel( rank, nprocs, comm );
		ColorModel.ReadParams( filename );
		ColorModel.SetDomain();
		ColorModel.ReadInput();
		ColorModel.Create();     // creating the model will create data structure to match the pore
		// structure and allocate variables
		ColorModel.Initialize(); // initializing the model will set initial conditions for variables

		ColorModel.Run();

	} // Limit scope so variables that contain communicators will free before MPI_Finialize
	cout << flush;
	Utilities::shutdown();
	return 0;
}
