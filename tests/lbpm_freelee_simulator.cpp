#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "common/Utilities.h"
#include "models/FreeLeeModel.h"
#include "analysis/FreeEnergy.h"

//*******************************************************************
// Implementation of Free-Energy Two-Phase LBM (Lee model)
//*******************************************************************

int main( int argc, char **argv )
{
	  
    // Initialize
    Utilities::startup( argc, argv );

    // Load the input database
    auto db = std::make_shared<Database>( argv[1] );
  
    { // Limit scope so variables that contain communicators will free before MPI_Finialize

        Utilities::MPI comm( MPI_COMM_WORLD );
        int rank = comm.getRank();
        int nprocs = comm.getSize();

        if (rank == 0){
            printf("********************************************************\n");
            printf("Running Free Energy Lee LBM	\n");
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
        ScaLBL_FreeLeeModel LeeModel( rank,nprocs,comm );
        LeeModel.ReadParams( filename );
        LeeModel.SetDomain();    
        LeeModel.ReadInput();    
        LeeModel.Create_TwoFluid();       
        
        FreeEnergyAnalyzer Analysis(LeeModel.Dm);
        
        LeeModel.Initialize_TwoFluid();   
        
        /*** RUN MAIN TIMESTEPS HERE ************/
        double MLUPS=0.0;
        int timestep = 0;
        int visualization_time = LeeModel.timestepMax;
    	if (LeeModel.vis_db->keyExists( "visualization_interval" )){
    		visualization_time = LeeModel.vis_db->getScalar<int>( "visualization_interval" );
    		timestep += visualization_time;
    	}
        while (LeeModel.timestep < LeeModel.timestepMax){
        	MLUPS = LeeModel.Run_TwoFluid(timestep);
        	if (rank==0) printf("Lattice update rate (per MPI process)= %f MLUPS \n", MLUPS);
        	Analysis.WriteVis(LeeModel,LeeModel.db, timestep);
        	timestep += visualization_time;
        }
        LeeModel.WriteDebug_TwoFluid();
    	if (rank==0) printf("********************************************************\n");
    	if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
    	MLUPS *= nprocs;
    	if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
    	if (rank==0) printf("********************************************************\n");
    	// ************************************************************************
    	
        PROFILE_STOP("Main");
        auto file = db->getWithDefault<std::string>( "TimerFile", "lbpm_freelee_simulator" );
        auto level = db->getWithDefault<int>( "TimerLevel", 1 );
	NULL_USE(level);
        PROFILE_SAVE( file,level );
        // ****************************************************


    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
    return 0;
}
