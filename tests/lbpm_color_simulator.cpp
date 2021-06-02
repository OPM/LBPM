#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "common/Utilities.h"
#include "models/ColorModel.h"

//#define WRE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */


//*************************************************************************
// Implementation of Two-Phase Immiscible LBM using CUDA
//*************************************************************************

int main( int argc, char **argv )
{

    // Initialize
    Utilities::startup( argc, argv );

    { // Limit scope so variables that contain communicators will free before MPI_Finialize

        Utilities::MPI comm( MPI_COMM_WORLD );
        int rank   = comm.getRank();
        int nprocs = comm.getSize();
        std::string SimulationMode = "production";
        // Load the input database
        auto db = std::make_shared<Database>( argv[1] );
        if (argc > 2) {
        	SimulationMode = "development";
        }

        if ( rank == 0 ) {
            printf( "********************************************************\n" );
            printf( "Running Color LBM	\n" );
            printf( "********************************************************\n" );
            if (SimulationMode == "development")
            	printf("**** DEVELOPMENT MODE ENABLED *************\n");
        }
        // Initialize compute device
        int device = ScaLBL_SetDevice( rank );
        NULL_USE( device );
        ScaLBL_DeviceBarrier();
        comm.barrier();

        PROFILE_ENABLE( 1 );
        // PROFILE_ENABLE_TRACE();
        // PROFILE_ENABLE_MEMORY();
        PROFILE_SYNCHRONIZE();
        PROFILE_START( "Main" );
        Utilities::setErrorHandlers();

        auto filename = argv[1];
        ScaLBL_ColorModel ColorModel( rank, nprocs, comm );
        ColorModel.ReadParams( filename );
        ColorModel.SetDomain();
        ColorModel.ReadInput();
        ColorModel.Create();     // creating the model will create data structure to match the pore
                                 // structure and allocate variables
        ColorModel.Initialize(); // initializing the model will set initial conditions for variables
       
        if (SimulationMode == "development"){
            double MLUPS=0.0;
            int timestep = 0;
            int analysis_interval = ColorModel.timestepMax;
        	if (ColorModel.analysis_db->keyExists( "" )){
        		analysis_interval = ColorModel.analysis_db->getScalar<int>( "analysis_interval" );
        	}
        	FlowAdaptor Adapt(ColorModel);
        	runAnalysis analysis(ColorModel);
            while (ColorModel.timestep < ColorModel.timestepMax){
            	timestep += analysis_interval;
            	MLUPS = ColorModel.Run(timestep);
            	if (rank==0) printf("Lattice update rate (per MPI process)= %f MLUPS \n", MLUPS);
            	
            	Adapt.MoveInterface(ColorModel);
            }
            ColorModel.WriteDebug();
        }            	//Analysis.WriteVis(LeeModel,LeeModel.db, timestep);

        else
        	ColorModel.Run();        

        PROFILE_STOP( "Main" );
        auto file  = db->getWithDefault<std::string>( "TimerFile", "lbpm_color_simulator" );
        auto level = db->getWithDefault<int>( "TimerLevel", 1 );
        PROFILE_SAVE( file, level );
        // ****************************************************


    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
    return 0;
}
