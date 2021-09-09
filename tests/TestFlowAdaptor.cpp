/* Test the flow adaptor class with color model */
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "common/Utilities.h"
#include "models/ColorModel.h"

inline void InitializeBubble(ScaLBL_ColorModel &ColorModel, double BubbleRadius){
	// initialize a bubble
	int i,j,k,n;
	int rank = ColorModel.Dm->rank();
	int nprocx = ColorModel.Dm->nprocx();
	int nprocy = ColorModel.Dm->nprocy();
	int nprocz = ColorModel.Dm->nprocz();
	int Nx = ColorModel.Dm->Nx;
	int Ny = ColorModel.Dm->Ny;
	int Nz = ColorModel.Dm->Nz;
	if (rank == 0) cout << "Setting up bubble..." << endl;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				ColorModel.Averages->SDs(i,j,k) = 100.f;
			}
		}
	}
	// Initialize the bubble
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nx + i;
				double iglobal= double(i+(Nx-2)*ColorModel.Dm->iproc())-double((Nx-2)*nprocx)*0.5;
				double jglobal= double(j+(Ny-2)*ColorModel.Dm->jproc())-double((Ny-2)*nprocy)*0.5;
				double kglobal= double(k+(Nz-2)*ColorModel.Dm->kproc())-double((Nz-2)*nprocz)*0.5;
				// Initialize phase position field for parallel bubble test
				if ((iglobal*iglobal)+(jglobal*jglobal)+(kglobal*kglobal) < BubbleRadius*BubbleRadius){
					ColorModel.Mask->id[n] = 2;
					ColorModel.Mask->id[n] = 2;
				}
				else{
					ColorModel.Mask->id[n]=1;
					ColorModel.Mask->id[n]=1;
				}
				ColorModel.id[n] = ColorModel.Mask->id[n];
				ColorModel.Dm->id[n] = ColorModel.Mask->id[n];
			}
		}
	}
	
	FILE *OUTFILE;
	char LocalRankFilename[40];
	sprintf(LocalRankFilename,"Bubble.%05i.raw",rank);
	OUTFILE = fopen(LocalRankFilename,"wb");
	fwrite(ColorModel.id,1,Nx*Ny*Nz,OUTFILE);
	fclose(OUTFILE);
	// initialize the phase indicator field
}


int main( int argc, char **argv )
{

	// Initialize
	Utilities::startup( argc, argv );

	{ // Limit scope so variables that contain communicators will free before MPI_Finialize

		Utilities::MPI comm( MPI_COMM_WORLD );
		int rank   = comm.getRank();
		int nprocs = comm.getSize();
		
		// Initialize compute device
		int device = ScaLBL_SetDevice( rank );
		NULL_USE( device );
		ScaLBL_DeviceBarrier();
		comm.barrier();
		
		Utilities::setErrorHandlers();
		
		if ( rank == 0 ) {
			printf( "********************************************************\n" );
			printf( "Test Flow Adaptor with Color LBM	\n" );
			printf( "********************************************************\n" );
		}
		
		/*  Populate the input database */
		auto domain_db = std::make_shared<Database>();
		auto color_db = std::make_shared<Database>();
		auto vis_db = std::make_shared<Database>();
		auto flow_db = std::make_shared<Database>();
		auto analysis_db = std::make_shared<Database>();
		auto db = std::make_shared<Database>();

		domain_db->putVector<int>( "nproc", { 1, 1, 1 } );
		domain_db->putVector<int>( "N", { 40, 40, 40 } );
		domain_db->putVector<int>( "n", { 40, 40, 40 } );
		domain_db->putScalar<int>( "BC", 0 );
		
		flow_db->putScalar<double>("mass_fraction_factor",0.01);
		flow_db->putScalar<int>("min_steady_timesteps",400);
		flow_db->putScalar<int>("max_steady_timesteps",400);
		flow_db->putScalar<int>("skiptimesteps",100);
		
		color_db->putScalar<double>("alpha",0.01);
		color_db->putScalar<int>("timestepMax",2000);
		color_db->putVector<int>( "ComponentLabels", { 0, -1 } );
		color_db->putVector<double>( "ComponentAffinity", { 0.0, 0.0 } );

		db->putDatabase("Color",color_db);
		db->putDatabase("Domain",domain_db);
		db->putDatabase("FlowAdaptor",flow_db);
		db->putDatabase("Visualization",vis_db);
		db->putDatabase("Analysis",analysis_db);

		ScaLBL_ColorModel ColorModel( rank, nprocs, comm );
		ColorModel.color_db = color_db;
		ColorModel.domain_db = domain_db;
		ColorModel.vis_db = vis_db;
		ColorModel.analysis_db = analysis_db;
		ColorModel.db = db;
		
		//ColorModel.ReadParams( filename );
		ColorModel.SetDomain();
		//ColorModel.ReadInput();
		double radius=15.5;
		InitializeBubble(ColorModel,radius);
		
		ColorModel.Create();     // creating the model will create data structure to match the pore
		// structure and allocate variables
		ColorModel.Initialize(); // initializing the model will set initial conditions for variables

		double MLUPS=0.0;
		int timestep = 0;
		bool ContinueSimulation = true;

		/* Variables for simulation protocols */
		//auto PROTOCOL = ColorModel.color_db->getWithDefault<std::string>( "protocol", "default" );
		std::string PROTOCOL = "fractional flow";
		/* image sequence protocol */
		int IMAGE_INDEX = 0;
		int IMAGE_COUNT = 0;
		std::vector<std::string> ImageList;
		/* flow adaptor keys to control behavior */			
		int SKIP_TIMESTEPS = 0;
		int MAX_STEADY_TIME = 1000000;
		double SEED_WATER = 0.01;
		double ENDPOINT_THRESHOLD = 0.1;
		double FRACTIONAL_FLOW_INCREMENT = 0.0; // this will skip the flow adaptor if not enabled

		MAX_STEADY_TIME = flow_db->getWithDefault<int>( "max_steady_timesteps", 1000000 );
		SKIP_TIMESTEPS = flow_db->getWithDefault<int>( "skip_timesteps", 50000 );
		ENDPOINT_THRESHOLD = flow_db->getWithDefault<double>( "endpoint_threshold", 0.1);
		/* protocol specific key values */
		if (PROTOCOL == "fractional flow")
			FRACTIONAL_FLOW_INCREMENT = flow_db->getWithDefault<double>( "fractional_flow_increment", 0.05);
		else if (PROTOCOL == "seed water")
			SEED_WATER = flow_db->getWithDefault<double>( "seed_water", 0.01);

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
		}		
		// ****************************************************


	} // Limit scope so variables that contain communicators will free before MPI_Finialize

	Utilities::shutdown();
	return 0;
}