#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "common/Utilities.h"
#include "models/ColorModel.h"

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
        Utilities::startup( argc, argv, true );

	{ // Limit scope so variables that contain communicators will free before MPI_Finialize

		Utilities::MPI comm( MPI_COMM_WORLD );
		int rank   = comm.getRank();
		int nprocs = comm.getSize();
		std::string SimulationMode = "production";
		// Load the input database
		//auto db = std::make_shared<Database>( argv[1] );
		//if (argc > 2) {
		//	SimulationMode = "legacy";
		//}

		if ( rank == 0 ) {
			printf( "********************************************************\n" );
			printf( "Running Color LBM	\n" );
			printf( "********************************************************\n" );
			if (SimulationMode == "legacy")
				printf("**** LEGACY MODE ENABLED *************\n");
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

		if (SimulationMode == "legacy"){
			ColorModel.Run();        
		}
		else {
			double MLUPS=0.0;
			int timestep = 0;
			bool ContinueSimulation = true;
			
			/* Variables for simulation protocols */
			auto PROTOCOL = ColorModel.color_db->getWithDefault<std::string>( "protocol", "default" );
			/* image sequence protocol */
			int IMAGE_INDEX = 0;
			int IMAGE_COUNT = 0;
			std::vector<std::string> ImageList;
			/* flow adaptor keys to control behavior */			
			int SKIP_TIMESTEPS = 0;
			int MAX_STEADY_TIME = 1000000;
			double ENDPOINT_THRESHOLD = 0.1;
			double FRACTIONAL_FLOW_INCREMENT = 0.0; // this will skip the flow adaptor if not enabled
			double SEED_WATER = 0.0;
			if (ColorModel.db->keyExists( "FlowAdaptor" )){
				auto flow_db = ColorModel.db->getDatabase( "FlowAdaptor" );
				MAX_STEADY_TIME = flow_db->getWithDefault<int>( "max_steady_timesteps", 1000000 );
				SKIP_TIMESTEPS = flow_db->getWithDefault<int>( "skip_timesteps", 50000 );
				ENDPOINT_THRESHOLD = flow_db->getWithDefault<double>( "endpoint_threshold", 0.1);
				/* protocol specific key values */
				if (PROTOCOL == "image sequence" || PROTOCOL == "core flooding")
					SKIP_TIMESTEPS = 0;
				if (PROTOCOL == "fractional flow")
					FRACTIONAL_FLOW_INCREMENT = flow_db->getWithDefault<double>( "fractional_flow_increment", 0.05);
				if (PROTOCOL == "seed water"){
					SEED_WATER = flow_db->getWithDefault<double>( "seed_water", 0.01);
					FRACTIONAL_FLOW_INCREMENT = flow_db->getWithDefault<double>( "fractional_flow_increment", 0.05);
				}
			}
			/* analysis keys*/
			int ANALYSIS_INTERVAL = ColorModel.timestepMax;
			if (ColorModel.analysis_db->keyExists( "analysis_interval" )){
				ANALYSIS_INTERVAL = ColorModel.analysis_db->getScalar<int>( "analysis_interval" );
			}
			/* Launch the simulation */
			FlowAdaptor Adapt(ColorModel);
			runAnalysis analysis(ColorModel);
			while (ContinueSimulation){
				/* this will run steady points */
				if (PROTOCOL == "fractional flow" || PROTOCOL == "seed water" || PROTOCOL == "shell aggregation" || PROTOCOL == "image sequence" )
					timestep += MAX_STEADY_TIME;
				else 
					timestep += ColorModel.timestepMax;
				/* Run the simulation timesteps*/
				MLUPS = ColorModel.Run(timestep);
				if (rank==0) printf("Lattice update rate (per MPI process)= %f MLUPS \n", MLUPS);
				if (ColorModel.timestep >= ColorModel.timestepMax){
					ContinueSimulation = false;
				}
				
				/* Load a new image if image sequence is specified */
				if (PROTOCOL == "image sequence"){
					IMAGE_INDEX++;
					if (IMAGE_INDEX < IMAGE_COUNT){
						std::string next_image = ImageList[IMAGE_INDEX];
						if (rank==0) printf("***Loading next image in sequence (%i) ***\n",IMAGE_INDEX);
						ColorModel.color_db->putScalar<int>("image_index",IMAGE_INDEX);
						Adapt.ImageInit(ColorModel, next_image);
					}
					else{
						if (rank==0) printf("Finished simulating image sequence \n");
						ColorModel.timestep =  ColorModel.timestepMax;
						ContinueSimulation = false;
					}
				}
				/*********************************************************/
				/* update the fluid configuration with the flow adapter */
				int skip_time = 0;
				timestep = ColorModel.timestep;
				/* get the averaged flow measures computed internally for the last simulation point*/
				double SaturationChange = 0.0;
				double volB = ColorModel.Averages->gwb.V; 
				double volA = ColorModel.Averages->gnb.V; 
				double initialSaturation = volB/(volA + volB);
				double vA_x = ColorModel.Averages->gnb.Px/ColorModel.Averages->gnb.M; 
				double vA_y = ColorModel.Averages->gnb.Py/ColorModel.Averages->gnb.M; 
				double vA_z = ColorModel.Averages->gnb.Pz/ColorModel.Averages->gnb.M; 
				double vB_x = ColorModel.Averages->gwb.Px/ColorModel.Averages->gwb.M; 
				double vB_y = ColorModel.Averages->gwb.Py/ColorModel.Averages->gwb.M; 
				double vB_z = ColorModel.Averages->gwb.Pz/ColorModel.Averages->gwb.M; 			
				double speedA = sqrt(vA_x*vA_x + vA_y*vA_y + vA_z*vA_z);
				double speedB = sqrt(vB_x*vB_x + vB_y*vB_y + vB_z*vB_z);
				/* stop simulation if previous point was sufficiently close to the endpoint*/
				if (volA*speedA < ENDPOINT_THRESHOLD*volB*speedB) ContinueSimulation = false;
				if (ContinueSimulation && SKIP_TIMESTEPS > 0 ){
					while (skip_time < SKIP_TIMESTEPS && fabs(SaturationChange) < fabs(FRACTIONAL_FLOW_INCREMENT) ){
						timestep += ANALYSIS_INTERVAL;
						if (PROTOCOL == "fractional flow")	{							
							Adapt.UpdateFractionalFlow(ColorModel);
						}
						else if (PROTOCOL == "shell aggregation"){
							double target_volume_change = FRACTIONAL_FLOW_INCREMENT*initialSaturation - SaturationChange;
							Adapt.ShellAggregation(ColorModel,target_volume_change);
						}
						else if (PROTOCOL == "seed water"){
							Adapt.SeedPhaseField(ColorModel,SEED_WATER);
						}
						/* Run some LBM timesteps to let the system relax a bit */
						MLUPS = ColorModel.Run(timestep);
						/* Recompute the volume fraction now that the system has adjusted */
						double volB = ColorModel.Averages->gwb.V; 
						double volA = ColorModel.Averages->gnb.V;
						SaturationChange = volB/(volA + volB) - initialSaturation;
						skip_time += ANALYSIS_INTERVAL;
					}
					if (rank==0) printf("  *********************************************************************  \n");
					if (rank==0) printf("   Updated fractional flow with saturation change = %f  \n", SaturationChange);
					if (rank==0) printf("   Used protocol = %s  \n", PROTOCOL.c_str());
					if (rank==0) printf("  *********************************************************************  \n");
				}
				/*********************************************************/
				if (rank==0) printf("   (flatten density field)  \n");
				if (PROTOCOL == "fractional flow")	{							
					Adapt.Flatten(ColorModel);
				}
			}
		}
		/*
		PROFILE_STOP( "Main" );
		auto file  = db->getWithDefault<std::string>( "TimerFile", "lbpm_color_simulator" );
		auto level = db->getWithDefault<int>( "TimerLevel", 1 );
		NULL_USE(level);
		PROFILE_SAVE( file, level );
		*/
		// ****************************************************


	} // Limit scope so variables that contain communicators will free before MPI_Finialize
	cout << flush;
	Utilities::shutdown();
	return 0;
}
