#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <math.h>

#include "models/PoissonSolver.h"

using namespace std;

//********************************************************
// Test lattice-Boltzmann solver of Poisson equation
//********************************************************

int main(int argc, char **argv)
{
    // Initialize MPI
    int provided_thread_support = -1;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided_thread_support);
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD,&comm);
    int rank = comm_rank(comm);
    int nprocs = comm_size(comm);
    if ( rank==0 && provided_thread_support<MPI_THREAD_MULTIPLE ){
      std::cerr << "Warning: Failed to start MPI with necessary thread support, thread support will be disabled" << std::endl;
    }

    // Limit scope so variables that contain communicators will free before MPI_Finialize
    { 
        if (rank == 0){
            printf("********************************************************\n");
            printf("Running Test for LB-Poisson Solver \n");
            printf("********************************************************\n");
        }
        //PROFILE_ENABLE_TRACE();
        //PROFILE_ENABLE_MEMORY();
        PROFILE_SYNCHRONIZE();
        PROFILE_START("Main");
        Utilities::setErrorHandlers();

        auto filename = argv[1];
        ScaLBL_Poisson PoissonSolver(rank,nprocs,comm); 

        // Initialize LB-Poisson model
        PoissonSolver.ReadParams(filename);
        PoissonSolver.SetDomain();    
        PoissonSolver.ReadInput();    
        PoissonSolver.Create();       
        PoissonSolver.Initialize();   

        //Initialize dummy charge density for test
        PoissonSolver.DummyChargeDensity();   

        PoissonSolver.Run(PoissonSolver.ChargeDensityDummy);
        PoissonSolver.getElectricPotential(1);
        PoissonSolver.getElectricField(1);

        if (rank==0) printf("Maximum timestep is reached and the simulation is completed\n");
        if (rank==0) printf("*************************************************************\n");

        PROFILE_STOP("Main");
        PROFILE_SAVE("TestPoissonSolver",1);
        // ****************************************************
        MPI_Barrier(comm);
    } // Limit scope so variables that contain communicators will free before MPI_Finialize
    MPI_Comm_free(&comm);
    MPI_Finalize();
}


