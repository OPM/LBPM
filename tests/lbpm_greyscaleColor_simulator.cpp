#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "models/GreyscaleColorModel.h"
#include "common/Utilities.h"
//#define WRITE_SURFACES

//*************************************************************************
// Implementation of Greyscale Two-Fluid Color LBM using CUDA
//*************************************************************************

using namespace std;


int main(int argc, char **argv)
{

    // Initialize MPI and error handlers
    Utilities::startup( argc, argv );

    { // Limit scope so variables that contain communicators will free before MPI_Finialize
  
      MPI_Comm comm;
      MPI_Comm_dup(MPI_COMM_WORLD,&comm);
      int rank = comm_rank(comm);
      int nprocs = comm_size(comm);
  
      if (rank == 0){
	      printf("****************************************\n");
	      printf("Running Greyscale Two-Phase Calculation \n");
	      printf("****************************************\n");
      }
      // Initialize compute device
      ScaLBL_SetDevice(rank);
      ScaLBL_DeviceBarrier();
      MPI_Barrier(comm);
  
      PROFILE_ENABLE(1);
      //PROFILE_ENABLE_TRACE();
      //PROFILE_ENABLE_MEMORY();
      PROFILE_SYNCHRONIZE();
      PROFILE_START("Main");
      Utilities::setErrorHandlers();
  
      auto filename = argv[1];
      ScaLBL_GreyscaleColorModel GreyscaleColor(rank,nprocs,comm);
      GreyscaleColor.ReadParams(filename);
      GreyscaleColor.SetDomain();    
      GreyscaleColor.ReadInput();    
      GreyscaleColor.Create();       // creating the model will create data structure to match the pore structure and allocate variables
      GreyscaleColor.Initialize();   // initializing the model will set initial conditions for variables
      GreyscaleColor.Run();	       
      GreyscaleColor.WriteDebug();
  
      PROFILE_STOP("Main");
      PROFILE_SAVE("lbpm_greyscaleColor_simulator",1);
      // ****************************************************
  
      MPI_Barrier(comm);
      MPI_Comm_free(&comm);
  
    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();

}
