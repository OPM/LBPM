#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <math.h>

#include "models/IonModel.h"
#include "models/PoissonSolver.h"
#include "models/MultiPhysController.h"
#include "common/Utilities.h"
#include "analysis/ElectroChemistry.h"

using namespace std;

int main(int argc, char **argv)
{
    // Initialize MPI and error handlers
    	Utilities::startup( argc, argv );
    	Utilities::MPI comm( MPI_COMM_WORLD );
    	int rank = comm.getRank();
    	int nprocs = comm.getSize();

    { // Limit scope so variables that contain communicators will free before MPI_Finialize

        if (rank == 0){
            printf("*************************************************\n");
            printf("Running LBPM Nernst-Planck solver \n");
            printf("*************************************************\n");
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
        ScaLBL_IonModel IonModel(rank,nprocs,comm);
        ScaLBL_Poisson PoissonSolver(rank,nprocs,comm); 
        ScaLBL_Multiphys_Controller Study(rank,nprocs,comm);//multiphysics controller coordinating multi-model coupling
        
        // Load controller information
        Study.ReadParams(filename);

        // Load user input database files for Ion solvers
        IonModel.ReadParams(filename);
        IonModel.SetDomain();    
        IonModel.ReadInput();    
        IonModel.Create();      
        
        // Create analysis object
        //ElectroChemistryAnalyzer Analysis(IonModel.Dm);


        //IonModel.timestepMax = Study.getIonNumIter_PNP_coupling(StokesModel.time_conv,IonModel.time_conv);
        IonModel.timestepMax = Study.getIonNumIter_NernstPlanck_coupling(IonModel.time_conv);
        IonModel.Initialize();   
        IonModel.DummyFluidVelocity();
        
        // Get maximal time converting factor based on Sotkes and Ion solvers
        //Study.getTimeConvMax_PNP_coupling(StokesModel.time_conv,IonModel.time_conv);
        // Get time conversion factor for the main iteration loop in electrokinetic single fluid simulator
        Study.time_conv_MainLoop = IonModel.timestepMax[0]*IonModel.time_conv[0];

        //----------------------------------- print out for debugging ------------------------------------------//
        if (rank==0){
            for (size_t i=0;i<IonModel.timestepMax.size();i++){
                printf("Main loop time_conv computed from ion %lu: %.5g[s/lt]\n",i+1,IonModel.timestepMax[i]*IonModel.time_conv[i]);
            }
        }
        //----------------------------------- print out for debugging ------------------------------------------//

        // Initialize LB-Poisson model
        PoissonSolver.ReadParams(filename);
        PoissonSolver.SetDomain();    
        PoissonSolver.ReadInput();    
        PoissonSolver.Create();       
        PoissonSolver.Initialize(Study.time_conv_MainLoop);   


        if (rank == 0){
            printf("********************************************************\n");
            printf("Key Summary of LBPM electrokinetic single-fluid solver \n");
            printf("   1. Max LB Timestep: %i [lt]\n", Study.timestepMax);
            printf("   2. Time conversion factor per LB Timestep: %.6g [sec/lt]\n",Study.time_conv_MainLoop);
            printf("   3. Max Physical Time: %.6g [sec]\n",Study.timestepMax*Study.time_conv_MainLoop);
            printf("********************************************************\n");
        }

        int timestep=0;
        while (timestep < Study.timestepMax){
            
            timestep++;
            PoissonSolver.Run(IonModel.ChargeDensity,false,timestep);//solve Poisson equtaion to get steady-state electrical potental
            //StokesModel.Run_Lite(IonModel.ChargeDensity, PoissonSolver.ElectricField);// Solve the N-S equations to get velocity
            IonModel.Run(IonModel.FluidVelocityDummy,PoissonSolver.ElectricField); //solve for ion transport and electric potential
            
            //if (timestep%Study.analysis_interval==0){
            //	Analysis.Basic(IonModel,PoissonSolver,StokesModel,timestep);
            //}
            if (timestep%Study.visualization_interval==0){
                //Analysis.WriteVis(IonModel,PoissonSolver,StokesModel,Study.db,timestep);
            	PoissonSolver.getElectricPotential_debug(timestep);
                PoissonSolver.getElectricField_debug(timestep);
                IonModel.getIonConcentration_debug(timestep);
                //StokesModel.getVelocity(timestep);
            	
            }
        }

        //if (rank==0) printf("Save simulation raw data at maximum timestep\n");
    	//Analysis.WriteVis(IonModel,PoissonSolver,StokesModel,Study.db,timestep);

        if (rank==0) printf("Maximum LB timestep = %i is reached and the simulation is completed\n",Study.timestepMax);
        if (rank==0) printf("*************************************************************\n");

        PROFILE_STOP("Main");
        PROFILE_SAVE("lbpm_nernst_planck_simulator",1);
        // ****************************************************

    } // Limit scope so variables that contain communicators will free before MPI_Finialize

  Utilities::shutdown();
}


