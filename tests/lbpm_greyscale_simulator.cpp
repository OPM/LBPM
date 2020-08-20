#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "models/GreyscaleModel.h"
#include "common/Utilities.h"
//#define WRITE_SURFACES

//****************************************************************
// Implementation of Greyscale Single-Fluid LBM using CUDA
//****************************************************************

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
	      printf("********************************************************\n");
	 	  printf("Running Greyscale Single Phase Permeability Calculation \n");
	 	  printf("********************************************************\n");
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
      ScaLBL_GreyscaleModel Greyscale(rank,nprocs,comm);
      Greyscale.ReadParams(filename);
      Greyscale.SetDomain();    
      Greyscale.ReadInput();    
      Greyscale.Create();       // creating the model will create data structure to match the pore structure and allocate variables
      Greyscale.Initialize();   // initializing the model will set initial conditions for variables
      Greyscale.Run();	       
	  Greyscale.VelocityField();
      //Greyscale.WriteDebug();
  
      PROFILE_STOP("Main");
      PROFILE_SAVE("lbpm_greyscale_simulator",1);
      // ****************************************************
  
      MPI_Barrier(comm);
      MPI_Comm_free(&comm);
  
    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();

}
