#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <math.h>

#include "models/PoissonSolver.h"
#include "common/Utilities.h"

using namespace std;

//********************************************************
// Test lattice-Boltzmann solver of Poisson equation
//********************************************************

int main(int argc, char **argv)
{
    // Initialize MPI
    Utilities::startup( argc, argv );
    
    {// Limit scope so variables that contain communicators will free before MPI_Finialize

        MPI_Comm comm;
        MPI_Comm_dup(MPI_COMM_WORLD,&comm);
        int rank = comm_rank(comm);
        int nprocs = comm_size(comm);

        if (rank == 0){
            printf("********************************************************\n");
            printf("Running Test for LB-Poisson Solver \n");
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
        ScaLBL_Poisson PoissonSolver(rank,nprocs,comm); 

        // Initialize LB-Poisson model
        PoissonSolver.ReadParams(filename);
        PoissonSolver.SetDomain();    
        PoissonSolver.ReadInput();    
        PoissonSolver.Create();       
        if (PoissonSolver.TestPeriodic==true){
            PoissonSolver.Initialize(PoissonSolver.TestPeriodicTimeConv);   
        }
        else {
            PoissonSolver.Initialize(0);   
        }

        //Initialize dummy charge density for test
        PoissonSolver.DummyChargeDensity();   

        if (PoissonSolver.TestPeriodic==true){
            if (rank==0) printf("Testing periodic voltage input is enabled. Total test time is %.3g[s], saving data every %.3g[s];
                                 user-specified time resolution is %.3g[s/lt]\n",
                                PoissonSolver.TestPeriodicTime,PoissonSolver.TestPeriodicSaveInterval,PoissonSolver.TestPeriodicTimeConv);
            int timestep = 0;
            while (timestep<(PoissonSolver.TestPeriodicTime/PoissonSolver.TestPeriodicTimeConv)){
                timestep++;
                PoissonSolver.Run(PoissonSolver.ChargeDensityDummy,timestep);
                if (timestep%(PoissonSolver.TestPeriodicSaveInterval/PoissonSolver.TestPeriodicTimeConv)==0){
                    if (rank==0) printf("   Time = %.3g[s]; saving electric potential and field\n",timestep*PoissonSolver.TestPeriodicTimeConv);
                    PoissonSolver.getElectricPotential_debug(timestep*PoissonSolver.TestPeriodicTimeConv);
                    PoissonSolver.getElectricField_debug(timestep*PoissonSolver.TestPeriodicTimeConv);
                }
            }
        }
        else {
            PoissonSolver.Run(PoissonSolver.ChargeDensityDummy,1);
            PoissonSolver.getElectricPotential_debug(1);
            PoissonSolver.getElectricField_debug(1);
        }

        if (rank==0) printf("Maximum timestep is reached and the simulation is completed\n");
        if (rank==0) printf("*************************************************************\n");

        PROFILE_STOP("Main");
        PROFILE_SAVE("TestPoissonSolver",1);
        // ****************************************************
        
        MPI_Barrier(comm);
        MPI_Comm_free(&comm);

    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
}


