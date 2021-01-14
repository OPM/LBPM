#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "models/DFHModel.h"

//#define WRE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */

using namespace std;

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
  auto thread_support = Utilities::MPI::queryThreadSupport();
  if ( rank==0 && thread_support != Utilities::MPI::ThreadSupport::MULTIPLE )
  std::cerr << "Warning: Failed to start MPI with necessary thread support, thread support will be disabled" << std::endl;
  { // Limit scope so variables that contain communicators will free before MPI_Finialize

	if (rank == 0){
		printf("********************************************************\n");
		printf("Running Color LBM	\n");
		printf("********************************************************\n");
	}
    PROFILE_ENABLE(1);
    //PROFILE_ENABLE_TRACE();
    //PROFILE_ENABLE_MEMORY();
    PROFILE_SYNCHRONIZE();
    PROFILE_START("Main");
    Utilities::setErrorHandlers();

	auto filename = argv[1];
	ScaLBL_DFHModel DFHModel(rank,nprocs,comm);
	DFHModel.ReadParams(filename);
	DFHModel.SetDomain();    
	DFHModel.ReadInput();    
	DFHModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables
	DFHModel.Initialize();   // initializing the model will set initial conditions for variables
	DFHModel.Run();	       
	DFHModel.WriteDebug();
	
    PROFILE_STOP("Main");
    PROFILE_SAVE("lbpm_color_simulator",1);
	// ****************************************************
	comm.barrier();
  } // Limit scope so variables that contain communicators will free before MPI_Finialize
    Utilities::shutdown();
}


