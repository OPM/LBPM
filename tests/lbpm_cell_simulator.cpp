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
#include "analysis/ElectroChemistry.h"

using namespace std;

//***************************************************************************
// Test lattice-Boltzmann Ion Model coupled with Poisson equation
//***************************************************************************

int main(int argc, char **argv)
{
    // Initialize MPI and error handlers
    	Utilities::startup( argc, argv );
    	Utilities::MPI comm( MPI_COMM_WORLD );
    	int rank = comm.getRank();
    	int nprocs = comm.getSize();

    { // Limit scope so variables that contain communicators will free before MPI_Finialize

        if (rank == 0){
            printf("********************************************************\n");
            printf("Running LBPM electrokinetic single-fluid solver \n");
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

	    bool SlipBC = false;
	
        // Load controller information
        Study.ReadParams(filename);

        // Load user input database files for Navier-Stokes and Ion solvers
        StokesModel.ReadParams(filename);

        // Setup other model specific structures
        StokesModel.SetDomain();    
        StokesModel.ReadInput();    
        StokesModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables
        comm.barrier();
        if (rank == 0) printf("Stokes model setup complete\n");
	
        IonModel.ReadParams(filename);
        IonModel.SetDomain();    
        IonModel.ReadInput();    
        IonModel.Create();      
        IonModel.SetMembrane();
        comm.barrier();
        if (rank == 0) printf("Ion model setup complete\n");
        fflush(stdout);
        
        // Create analysis object
        ElectroChemistryAnalyzer Analysis(IonModel.Dm);

        // Get internal iteration number
        StokesModel.timestepMax = Study.getStokesNumIter_PNP_coupling(StokesModel.time_conv,IonModel.time_conv);
        StokesModel.Initialize();   // initializing the model will set initial conditions for variables
        comm.barrier();
        if (rank == 0) printf("Stokes model initialized \n");
        fflush(stdout);
	
        IonModel.timestepMax = Study.getIonNumIter_PNP_coupling(StokesModel.time_conv,IonModel.time_conv);
        IonModel.Initialize();   
        IonModel.DummyFluidVelocity();
        comm.barrier();
        if (rank == 0) printf("Ion model initialized \n");
        // Get maximal time converting factor based on Sotkes and Ion solvers
        Study.getTimeConvMax_PNP_coupling(StokesModel.time_conv,IonModel.time_conv);

        // Initialize LB-Poisson model
        PoissonSolver.ReadParams(filename);
        PoissonSolver.SetDomain();    
        PoissonSolver.ReadInput();    
        PoissonSolver.Create();       
        comm.barrier();
        if (rank == 0) printf("Poisson solver created \n");
        fflush(stdout);
        PoissonSolver.Initialize(Study.time_conv_max);   
        comm.barrier();
        if (rank == 0) printf("Poisson solver initialized \n");
        fflush(stdout);

        int timestep=0;
        while (timestep < Study.timestepMax){
            
            timestep++;
            PoissonSolver.Run(IonModel.ChargeDensity,SlipBC,timestep);//solve Poisson equtaion to get steady-state electrical potental
            comm.barrier();
            //if (rank == 0) printf("    Poisson step %i \n",timestep);
            //StokesModel.Run_Lite(IonModel.ChargeDensity, PoissonSolver.ElectricField);// Solve the N-S equations to get velocity
    	    //fflush(stdout);

            IonModel.RunMembrane(IonModel.FluidVelocityDummy,PoissonSolver.ElectricField,PoissonSolver.Psi); //solve for ion transport with membrane
	        comm.barrier();
            //if (rank == 0) printf("    Membrane step %i \n",timestep);
	        //fflush(stdout);

            //timestep++;//AA operations

            if (timestep%Study.analysis_interval==0){
	            //Analysis.Basic(IonModel,PoissonSolver,StokesModel,timestep);
            }
            if (timestep%Study.visualization_interval==0){
            	//Analysis.WriteVis(IonModel,PoissonSolver,StokesModel,Study.db,timestep);
            	// PoissonSolver.getElectricPotential(timestep);
                //PoissonSolver.getElectricField(timestep);
                //IonModel.getIonConcentration(timestep);
                //StokesModel.getVelocity(timestep);
            	PoissonSolver.getElectricPotential_debug(timestep);
                PoissonSolver.getElectricField_debug(timestep);
                IonModel.getIonConcentration_debug(timestep);
            }
        }

        if (rank==0) printf("Save simulation raw data at maximum timestep\n");
    	Analysis.WriteVis(IonModel,PoissonSolver,StokesModel,Study.db,timestep);

        if (rank==0) printf("Maximum timestep is reached and the simulation is completed\n");
        if (rank==0) printf("*************************************************************\n");

        PROFILE_STOP("Main");
        PROFILE_SAVE("lbpm_electrokinetic_SingleFluid_simulator",1);
        // ****************************************************

    } // Limit scope so variables that contain communicators will free before MPI_Finialize

  Utilities::shutdown();
}


