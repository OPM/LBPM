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
        	/* flow adaptor keys to control */
        	int SKIP_TIMESTEPS = 0;
        	int MAX_STEADY_TIME = 1000000;
        	double FRACTIONAL_FLOW_INCREMENT = 0.05
        	if (ColorModel.db->keyExists( "FlowAdaptor" )){
        		auto flow_db = ColorModel.db->getDatabase( "FlowAdaptor" );
        		MAX_STEADY_TIME = flow_db->getWithDefault<int>( "max_steady_timesteps", 1000000 );
        		SKIP_TIMESTEPS = flow_db->getWithDefault<int>( "skip_timesteps", 100000 );
        		FRACTIONAL_FLOW_INCREMENT = flow_db->getWithDefault<double>( "fractional_flow_increment", 0.05);
        	}
        	int ANALYSIS_INTERVAL = ColorModel.timestepMax;
        	if (ColorModel.analysis_db->keyExists( "analysis_interval" )){
        		ANALYSIS_INTERVAL = ColorModel.analysis_db->getScalar<int>( "analysis_interval" );
        	}
        	/* Launch the simulation */
        	FlowAdaptor Adapt(ColorModel);
        	runAnalysis analysis(ColorModel);
        	while (ColorModel.timestep < ColorModel.timestepMax){
            	/* this will run steady points */
            	timestep += MAX_STEADY_TIME;
            	MLUPS = ColorModel.Run(timestep);
            	if (rank==0) printf("Lattice update rate (per MPI process)= %f MLUPS \n", MLUPS);
            	
            	/* update the fractional flow by adding mass */
            	int skip_time = 0;
            	timestep = ColorModel.timestep;
    			double SaturationChange = 0.0;
    			double volB = M.Averages->gwb.V; 
    			double volA = M.Averages->gnb.V; 
    			double initialSaturation = volB/(volA + volB);
            	while (skip_time < SKIP_TIMESTEPS && fabs(SaturationChange) < fabs(FRACTIONAL_FLOW_INCREMENT) ){
            		timestep += ANALYSIS_INTERVAL;
                	Adapt.UpdateFractionalFlow(ColorModel);
                	MLUPS = ColorModel.Run(timestep);
        			double volB = M.Averages->gwb.V; 
        			double volA = M.Averages->gnb.V;
        			SaturationChange = volB/(volA + volB) - initialSaturation;
            		skip_time += ANALYSIS_INTERVAL;
            	}
            	if (rank==0) printf("  *********************************************************************  \n");
            	if (rank==0) printf("   Updated fractional flow with saturation change = %f  \n", SaturationChange);
            	if (rank==0) printf("  *********************************************************************  \n");

            	/* apply timestep skipping algorithm to accelerate steady-state */
            	/* skip_time = 0;
            	timestep = ColorModel.timestep;
            	while (skip_time < SKIP_TIMESTEPS){
            		timestep += ANALYSIS_INTERVAL;
                	MLUPS = ColorModel.Run(timestep);
            		Adapt.MoveInterface(ColorModel);
            		skip_time += ANALYSIS_INTERVAL;
            	}
            	*/
            }
            ColorModel.WriteDebug();
        }

        else
        	ColorModel.Run();        

        PROFILE_STOP( "Main" );
        auto file  = db->getWithDefault<std::string>( "TimerFile", "lbpm_color_simulator" );
        auto level = db->getWithDefault<int>( "TimerLevel", 1 );
	NULL_USE(level);
        PROFILE_SAVE( file, level );
        // ****************************************************


    } // Limit scope so variables that contain communicators will free before MPI_Finialize

    Utilities::shutdown();
    return 0;
}
