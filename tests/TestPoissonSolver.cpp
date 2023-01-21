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
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();    
    {// Limit scope so variables that contain communicators will free before MPI_Finialize

        if (rank == 0){
            printf("********************************************************\n");
            printf("Running Test for LB-Poisson Solver \n");
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
            if (rank==0) printf("Testing periodic voltage input is enabled. Total test time is %.3g[s], saving data every %.3g[s]; user-specified time resolution is %.3g[s/lt]\n",
                                PoissonSolver.TestPeriodicTime,PoissonSolver.TestPeriodicSaveInterval,PoissonSolver.TestPeriodicTimeConv);
            int timestep = 0;
            int timeMax = int(PoissonSolver.TestPeriodicTime/PoissonSolver.TestPeriodicTimeConv);
            int timeSave = int(PoissonSolver.TestPeriodicSaveInterval/PoissonSolver.TestPeriodicTimeConv);
            while (timestep<timeMax){
                timestep++;
                PoissonSolver.Run(PoissonSolver.ChargeDensityDummy,false,timestep);
                if (timestep%timeSave==0){
                    if (rank==0) printf("   Time = %.3g[s]; saving electric potential and field\n",timestep*PoissonSolver.TestPeriodicTimeConv);
                    PoissonSolver.getElectricPotential_debug(timestep);
                    PoissonSolver.getElectricField_debug(timestep);
                }
            }
            PoissonSolver.WriteVis(timestep);
        }
        else {
        	int timestep = 1;
            PoissonSolver.Run(PoissonSolver.ChargeDensityDummy,false,1);
            PoissonSolver.getElectricPotential_debug(1);
            PoissonSolver.getElectricField_debug(1);
            PoissonSolver.WriteVis(timestep);
        }

        if (rank==0) printf("Maximum timestep is reached and the simulation is completed\n");
        if (rank==0) printf("*************************************************************\n");

        PROFILE_STOP("Main");
        PROFILE_SAVE("TestPoissonSolver",1);
        // ****************************************************
        
    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
}


