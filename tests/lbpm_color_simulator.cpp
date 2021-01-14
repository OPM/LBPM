#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "models/ColorModel.h"
#include "common/Utilities.h"

//#define WRE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */


//*************************************************************************
// Implementation of Two-Phase Immiscible LBM using CUDA
//*************************************************************************

int main(int argc, char **argv)
{
	  
	// Initialize MPI
  Utilities::startup( argc, argv );
  Utilities::MPI comm( MPI_COMM_WORLD );
  int rank = comm.getRank();
  int nprocs = comm.getSize();
  
  // Load the input database
  auto db = std::make_shared<Database>( argv[1] );

  // Initialize MPI and error handlers
  auto multiple = db->getWithDefault<bool>( "MPI_THREAD_MULTIPLE", true );
  Utilities::startup( argc, argv, multiple );
  Utilities::MPI::changeProfileLevel( 1 );

  { // Limit scope so variables that contain communicators will free before MPI_Finialize

    if (rank == 0){
	    printf("********************************************************\n");
	    printf("Running Color LBM	\n");
	    printf("********************************************************\n");
    }
	// Initialize compute device
	int device=ScaLBL_SetDevice(rank);
    NULL_USE( device );
	ScaLBL_DeviceBarrier();
	comm.barrier();

    PROFILE_ENABLE(1);
    //PROFILE_ENABLE_TRACE();
    //PROFILE_ENABLE_MEMORY();
    PROFILE_SYNCHRONIZE();
    PROFILE_START("Main");
    Utilities::setErrorHandlers();

    auto filename = argv[1];
    ScaLBL_ColorModel ColorModel(rank,nprocs,comm);
    ColorModel.ReadParams(filename);
    ColorModel.SetDomain();    
    ColorModel.ReadInput();    
    ColorModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables
    ColorModel.Initialize();   // initializing the model will set initial conditions for variables
    ColorModel.Run();	       
    //ColorModel.WriteDebug();

    PROFILE_STOP("Main");
    auto file = db->getWithDefault<std::string>( "TimerFile", "lbpm_color_simulator" );
    auto level = db->getWithDefault<int>( "TimerLevel", 1 );
    PROFILE_SAVE(file,level);
    // ****************************************************


  } // Limit scope so variables that contain communicators will free before MPI_Finialize

  Utilities::shutdown();
}


