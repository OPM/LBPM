#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <math.h>

#include "models/StokesModel.h"
#include "models/IonModel.h"
#include "models/PoissonSolver.h"
#include "models/MultiPhysController.h"

using namespace std;

//***************************************************************************
// Implementation of Multiphysics simulator using lattice-Boltzmann method
//***************************************************************************

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
            printf("Running Electrokinetic LBM Simulator \n");
            printf("********************************************************\n");
        }
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

        // Initialize LB Navier-Stokes model
        StokesModel.ReadParams(filename,Study.num_iter_Stokes);
        StokesModel.SetDomain();    
        StokesModel.ReadInput();    
        StokesModel.Create();       // creating the model will create data structure to match the pore structure and allocate variables
        StokesModel.Initialize();   // initializing the model will set initial conditions for variables

        // Initialize LB-Ion model
        IonModel.ReadParams(filename,Study.num_iter_Ion,Study.num_iter_Stokes,StokesModel.time_conv);
        IonModel.SetDomain();    
        IonModel.ReadInput();    
        IonModel.Create();       
        IonModel.Initialize();   

        // Initialize LB-Poisson model
        PoissonSolver.ReadParams(filename);
        PoissonSolver.SetDomain();    
        PoissonSolver.ReadInput();    
        PoissonSolver.Create();       
        PoissonSolver.Initialize();   

        int timestep=0;
        while (timestep < Study.timestepMax){
            
            timestep++;
            PoissonSolver.Run(IonModel.ChargeDensity);//solve Poisson equtaion to get steady-state electrical potental
            StokesModel.Run_Lite(IonModel.ChargeDensity, PoissonSolver.ElectricField);// Solve the N-S equations to get velocity
            IonModel.Run(StokesModel.Velocity,PoissonSolver.ElectricField); //solve for ion transport and electric potential
            
            
            timestep++;//AA operations
            //--------------------------------------------
            //potentially leave analysis module for future
            //--------------------------------------------
        }

        //StokesModel.WriteDebug();

        PROFILE_STOP("Main");
        PROFILE_SAVE("lbpm_electrokinetic_simulator",1);
        // ****************************************************
        MPI_Barrier(comm);
    } // Limit scope so variables that contain communicators will free before MPI_Finialize
    MPI_Comm_free(&comm);
    MPI_Finalize();
}


