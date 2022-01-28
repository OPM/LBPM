#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <math.h>

#include "models/IonModel.h"
#include "models/StokesModel.h"
#include "models/PoissonSolver.h"
#include "models/MultiPhysController.h"
#include "common/Utilities.h"

using namespace std;

//***************************************************************************
// Test lattice-Boltzmann Ion Model coupled with Poisson equation
//***************************************************************************

int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int nprocs = comm.getSize();

    { // Limit scope so variables that contain communicators will free before MPI_Finialize


        if (rank == 0){
            printf("********************************************************\n");
            printf("Running Test for LB-Poisson-Ion Coupling \n");
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
        ScaLBL_StokesModel StokesModel(rank,nprocs,comm);
        ScaLBL_IonModel IonModel(rank,nprocs,comm);
        ScaLBL_Poisson PoissonSolver(rank,nprocs,comm); 
        ScaLBL_Multiphys_Controller Study(rank,nprocs,comm);//multiphysics controller coordinating multi-model coupling

        // Load controller information
        Study.ReadParams(filename);

        // Load user input database files for Navier-Stokes and Ion solvers
        StokesModel.ReadParams(filename);
        IonModel.ReadParams(filename);

        // Setup other model specific structures
        StokesModel.SetDomain();    
        StokesModel.ReadInput();    
        StokesModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables

        IonModel.SetDomain();    
        IonModel.ReadInput();    
        IonModel.Create();       

        // Get internal iteration number
        StokesModel.timestepMax = Study.getStokesNumIter_PNP_coupling(StokesModel.time_conv,IonModel.time_conv);
        StokesModel.Initialize();   // initializing the model will set initial conditions for variables

        IonModel.timestepMax = Study.getIonNumIter_PNP_coupling(StokesModel.time_conv,IonModel.time_conv);
        IonModel.Initialize();   

        // Initialize LB-Poisson model
        PoissonSolver.ReadParams(filename);
        PoissonSolver.SetDomain();    
        PoissonSolver.ReadInput();    
        PoissonSolver.Create();       
        PoissonSolver.Initialize(0);   


        int timestep=0;
        double error = 1.0;
        double error_ion    = 1.0;
        double error_stokes = 1.0;
        vector<double>ci_avg_previous{0.0,0.0};//assuming 1:1 solution
        double vel_avg_previous = 0.0;
        while (timestep < Study.timestepMax && error > Study.tolerance){
            
            timestep++;
            PoissonSolver.Run(IonModel.ChargeDensity,StokesModel.UseSlippingVelBC,0);//solve Poisson equtaion to get steady-state electrical potental
            StokesModel.Run_Lite(IonModel.ChargeDensity, PoissonSolver.ElectricField);// Solve the N-S equations to get velocity
            IonModel.Run(StokesModel.Velocity,PoissonSolver.ElectricField); //solve for ion transport and electric potential
            
            timestep++;//AA operations

            if (timestep%Study.analysis_interval==0){
                error_ion = IonModel.CalIonDenConvergence(ci_avg_previous);
                error_stokes = StokesModel.CalVelocityConvergence(vel_avg_previous,IonModel.ChargeDensity,PoissonSolver.ElectricField);
                error = max(error_ion,error_stokes);
            }
        }

        PoissonSolver.getElectricPotential_debug(timestep);
        PoissonSolver.getElectricField_debug(timestep);
        IonModel.getIonConcentration_debug(timestep);
        StokesModel.getVelocity_debug(timestep);

        if (rank==0) printf("Maximum timestep is reached and the simulation is completed\n");
        if (rank==0) printf("*************************************************************\n");

        PROFILE_STOP("Main");
        PROFILE_SAVE("TestPNP_Stokes",1);
        // ****************************************************


    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
}


