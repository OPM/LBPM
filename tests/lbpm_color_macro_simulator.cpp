#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/TwoPhase.h"
#include "common/MPI_Helpers.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"
#include "lbpm_color_simulator.h"

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2018
 */

using namespace std;

//*************************************************************************
// Implementation of Two-Phase Immiscible LBM 
//*************************************************************************

int main(int argc, char **argv)
{
  // Initialize MPI
	int provided_thread_support = -1;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided_thread_support);
	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD,&comm);
	int rank = comm_rank(comm);
	int nprocs = comm_size(comm);
	if ( rank==0 && provided_thread_support<MPI_THREAD_MULTIPLE )
		std::cerr << "Warning: Failed to start MPI with necessary thread support, thread support will be disabled" << std::endl;
	{ // Limit scope so variables that contain communicators will free before MPI_Finialize

		// parallel domain size (# of sub-domains)
		int nprocx,nprocy,nprocz;
		int iproc,jproc,kproc;

		MPI_Request req1[18],req2[18];
		MPI_Status stat1[18],stat2[18];

		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Color LBM	\n");
			printf("********************************************************\n");
		}

		PROFILE_ENABLE(1);
		//PROFILE_ENABLE_TRACE();
		//PROFILE_ENABLE_MEMORY();
		PROFILE_SYNCHRONIZE();
		PROFILE_START("Main");
		Utilities::setErrorHandlers();

		// Variables that specify the computational domain  
		string FILENAME;
		int Nx,Ny,Nz;		// local sub-domain size
		double Lx,Ly,Lz;	// Domain length
		double D = 1.0;		// reference length for non-dimensionalization
		// Color Model parameters
		int timestepMax;
		double tau1, tau2, rho1,rho2;
		double Fx,Fy,Fz,tol,err;
		double alpha, beta;
		int BoundaryCondition;
		int InitialCondition;
		//	bool pBC,Restart;
		int i,j,k;
		double din, dout, flux;
		double inletA,inletB,outletA,outletB;
		inletA=1.f;
		inletB=0.f;
		outletA=0.f;
		outletB=1.f;
		typeBC=4;
		flux = 10.f;
		dout=1.f;

		int RESTART_INTERVAL=20000;
		//int ANALYSIS_)INTERVAL=1000;	
		int BLOB_ANALYSIS_INTERVAL=1000;
		int timestep = -1;

		if (rank==0){
			//.............................................................
			//		READ SIMULATION PARMAETERS FROM INPUT FILE
			//.............................................................
			ifstream input("Color.in");
			if (input.is_open()){
				// Line 1: model parameters (tau, alpha, beta, das, dbs)
				input >> tau1;			// Viscosity non-wetting
				input >> tau2;			// Viscosity wetting
				input >> rho1;			// density non-wetting
				input >> rho2;			// density wetting
				input >> alpha;			// Surface Tension parameter
				input >> beta;			// Width of the interface
				// Line 2:  External force components (Fx,Fy, Fz)
				input >> Fx;
				input >> Fy;
				input >> Fz;
				// Line 4: Pressure Boundary conditions
				input >> InitialCondition;
				input >> BoundaryCondition;
				input >> din;
				input >> dout;
				// Line 5: time-stepping criteria
				input >> timestepMax;		// max no. of timesteps
				input >> RESTART_INTERVAL;	// restart interval
				input >> tol;		      	// error tolerance
				// Line 6: Analysis options
				input >> BLOB_ANALYSIS_INTERVAL; // interval to analyze blob states
				//.............................................................
			}
			else{
				// Set default values
				// Print warning
				printf("WARNING: No input file provided (Color.in is missing)! Default parameters will be used. \n");
				tau1 = tau2 = 1.0;
				rho1 = rho2 = 1.0;
				alpha=0.005;
				beta= 0.9;
				Fx = Fy = Fz = 0.0;
				InitialCondition=0;
				BoundaryCondition=0;
				din=dout=1.0;
				timestepMax=0;
			}

			//.......................................................................
			// Reading the domain information file
			//.......................................................................
			ifstream domain("Domain.in");
			if (input.is_open()){
				domain >> nprocx;
				domain >> nprocy;
				domain >> nprocz;
				domain >> Nx;
				domain >> Ny;
				domain >> Nz;
				domain >> nspheres;
				domain >> Lx;
				domain >> Ly;
				domain >> Lz;
				//.......................................................................
			}
			else{
				// Set default values
				// Print warning
				printf("WARNING: No input file provided (Domain.in is missing)! Default parameters will be used. \n");
				nprocx=nprocy=nprocz=1;
				Nx=Ny=Nz=10;
				nspheres=0;
				Lx=Ly=Lz=1.0;
			}
		}
		// **************************************************************
		// Broadcast simulation parameters from rank 0 to all other procs
		MPI_Barrier(comm);
		//.................................................
		MPI_Bcast(&tau1,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&tau2,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&rho1,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&rho2,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&alpha,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&beta,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&BoundaryCondition,1,MPI_INT,0,comm);
		MPI_Bcast(&InitialCondition,1,MPI_INT,0,comm);
		MPI_Bcast(&din,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&dout,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&timestepMax,1,MPI_INT,0,comm);
		MPI_Bcast(&RESTART_INTERVAL,1,MPI_INT,0,comm);
		MPI_Bcast(&tol,1,MPI_DOUBLE,0,comm);
		// Computational domain
		MPI_Bcast(&Nx,1,MPI_INT,0,comm);
		MPI_Bcast(&Ny,1,MPI_INT,0,comm);
		MPI_Bcast(&Nz,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
		MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
		MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
		//.................................................

		double flux = 0.f;
		if (BoundaryCondition==4) flux = din*rho1; // mass flux must adjust for density (see formulation for details)

		// Get the rank info
		const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

		MPI_Barrier(comm);

		if (nprocs != nprocx*nprocy*nprocz){
			printf("nprocx =  %i \n",nprocx);
			printf("nprocy =  %i \n",nprocy);
			printf("nprocz =  %i \n",nprocz);
			INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
		}

		if (rank==0){
			printf("********************************************************\n");
			printf("tau (non-wetting) = %f \n", tau1);
			printf("tau (wetting) = %f \n", tau2);
			printf("density (non-wetting) = %f \n", rho1);
			printf("density (wetting) = %f \n", rho2);
			printf("alpha = %f \n", alpha);		
			printf("beta = %f \n", beta);
			printf("gamma_{wn} = %f \n", 5.796*alpha);
			printf("Force(x) = %f \n", Fx);
			printf("Force(y) = %f \n", Fy);
			printf("Force(z) = %f \n", Fz);
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
			if (BoundaryCondition==0) printf("Periodic boundary conditions will applied \n");
			if (BoundaryCondition==1) printf("Pressure boundary conditions will be applied \n");
			if (BoundaryCondition==2) printf("Velocity boundary conditions will be applied \n");
			if (BoundaryCondition==3) printf("Dynamic pressure boundary conditions will be applied \n");
			if (BoundaryCondition==4) printf("Average flux boundary conditions will be applied \n");
			if (InitialCondition==0) printf("Initial conditions assigned from phase ID file \n");
			if (InitialCondition==1) printf("Initial conditions assigned from restart file \n");
			printf("********************************************************\n");
		}

		// Initialized domain and averaging framework for Two-Phase Flow
		bool pBC,velBC;
		if (BoundaryCondition==1 || BoundaryCondition==3 || BoundaryCondition == 4)
			pBC=true;
		else						pBC=false;
		if (BoundaryCondition==2)	velBC=true;
		else						velBC=false;

		bool Restart;
		if (InitialCondition==1)    Restart=true;
		else 						Restart=false;
		NULL_USE(pBC);  NULL_USE(velBC);

		// Full domain used for averaging (do not use mask for analysis)
		Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);
		for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = 1;
		std::shared_ptr<TwoPhase> Averages( new TwoPhase(Dm) );
		//   TwoPhase Averages(Dm);
		Dm.CommInit(comm);

		// Mask that excludes the solid phase
		Domain Mask(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);
		MPI_Barrier(comm);

		Nx+=2; Ny+=2; Nz += 2;
		int N = Nx*Ny*Nz;
		int dist_mem_size = N*sizeof(double);

		//.......................................................................
		if (rank == 0)	printf("Read input media... \n");
		//.......................................................................

		//.......................................................................
		// Filenames used
		char LocalRankString[8];
		char LocalRankFilename[40];
		char LocalRestartFile[40];
		char tmpstr[10];
		sprintf(LocalRankString,"%05d",rank);
		sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
		sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);

		//	printf("Local File Name =  %s \n",LocalRankFilename);
		// .......... READ THE INPUT FILE .......................................
		//	char value;
		char *id;
		id = new char[N];
		double sum_local;
		double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
		if (BoundaryCondition > 0) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
		double porosity, pore_vol;
		//...........................................................................
		if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;

		//.......................................................................
		// Read the signed distance
		sprintf(LocalRankString,"%05d",rank);
		sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
		ReadBinaryFile(LocalRankFilename, Averages->SDs.data(), N);
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;

		//.......................................................................
		// Assign the phase ID field based on the signed distance
		//.......................................................................

		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					int n = k*Nx*Ny+j*Nx+i;
					id[n] = 0;
				}
			}
		}
		double sum=0;
		pore_vol = 0.0;
		for ( k=0;k<Nz;k++){
			for ( j=0;j<Ny;j++){
				for ( i=0;i<Nx;i++){
					int n = k*Nx*Ny+j*Nx+i;
					if (Averages->SDs(n) > 0.0){
						id[n] = 2;	
					}
					// compute the porosity (actual interface location used)
					if (Averages->SDs(n) > 0.0){
						sum++;	
					}
				}
			}
		}

		if (rank==0) printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
		sprintf(LocalRankFilename,"ID.%05i",rank);
		size_t readID;
		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
		readID=fread(id,1,N,IDFILE);
		if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading ID (rank=%i) \n",rank);
		fclose(IDFILE);
		
		//.......................................................................
		// Compute the media porosity, assign phase labels and solid composition
		//.......................................................................
		double sum;
		double sum_local=0.0, porosity;
		int Np=0;  // number of local pore nodes
		double *PhaseLabel;
		PhaseLabel = new double[N];
		Dm.AssignComponentLabels(PhaseLabel);
		//.......................................................................
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (id[n] > 0){
						sum_local+=1.0;
						Np++;
					}
				}
			}
		}
		MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
		porosity = sum*iVol_global;
		if (rank==0) printf("Media porosity = %f \n",porosity);

		//.........................................................
		// don't perform computations at the eight corners
		id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
		id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
		//.........................................................

		// Initialize communication structures in averaging domain
		for (i=0; i<Mask.Nx*Mask.Ny*Mask.Nz; i++) Mask.id[i] = id[i];
		Mask.CommInit(comm);

		//...........................................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device (will use optimized layout)
		ScaLBL_Communicator ScaLBL_Comm(Mask);
		//Create a second communicator based on the regular data layout
		ScaLBL_Communicator ScaLBL_Comm_Regular(Mask);
		
		if (rank==0)	printf ("Set up memory efficient layout \n");
		int neighborSize=18*Np*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Np];
		ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Dm.id,Np);
		MPI_Barrier(comm);

		//...........................................................................
		//				MAIN  VARIABLES ALLOCATED HERE
		//...........................................................................
		// LBM variables
		if (rank==0)	printf ("Allocating distributions \n");
		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);
		if (rank==0)	printf ("Allocating distributions \n");

		int *NeighborList;
		int *dvcMap;
		//		double *f_even,*f_odd;
		double *fq, *Aq, *Bq;
		double *Den, *Phi;
		double *ColorGrad;
		double *Vel;
		double *Pressure;
		
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);

		ScaLBL_AllocateDeviceMemory((void **) &fq, 19*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Aq, 7*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Bq, 7*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Nx*Ny*Nz);		
		ScaLBL_AllocateDeviceMemory((void **) &Pressure, sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &Vel, 3*sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*sizeof(double)*Np);
		
		//...........................................................................
		// Update GPU data structures
		if (rank==0)	printf ("Setting up device map and neighbor list \n");
		int *TmpMap;
		TmpMap=new int[Np];
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx < 0))
						TmpMap[idx] = k*Nx*Ny+j*Nx+i;
				}
			}
		}
		//for (int idx=0; idx<Np; idx++) printf("Map=%i\n",TmpMap[idx]);

		ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		delete [] TmpMap;
		
		// copy the neighbor list 
		ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
		// initialize phi based on PhaseLabel (include solid component labels)
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
		//...........................................................................

		if (rank==0)	printf ("Initializing distributions \n");
		// Initialize the phase field and variables
		ScaLBL_D3Q19_Init(fq, Np);
		if (rank==0)	printf ("Initializing phase field \n");
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, Np);
		if (Dm.kproc==0){
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
		}
		if (Dm.kproc == nprocz-1){
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
		}
		

		if (Restart == true){
			if (rank==0){
				printf("Reading restart file! \n");
				ifstream restart("Restart.txt");
				if (restart.is_open()){
					restart  >> timestep;
					printf("Restarting from timestep =%i \n",timestep);
				}
				else{
					printf("WARNING:No Restart.txt file, setting timestep=0 \n");
					timestep=0;
				}
			}
			MPI_Bcast(&timestep,1,MPI_INT,0,comm);

			// Read in the restart file to CPU buffers
			double *cDen = new double[2*N];
			double *cDistEven = new double[10*N];
			double *cDistOdd = new double[9*N];
			ReadCheckpoint(LocalRestartFile, cDen, cDistEven, cDistOdd, N);
			// Copy the restart data to the GPU
			ScaLBL_CopyToDevice(f_even,cDistEven,10*N*sizeof(double));
			ScaLBL_CopyToDevice(f_odd,cDistOdd,9*N*sizeof(double));
			ScaLBL_CopyToDevice(Den,cDen,2*N*sizeof(double));
			ScaLBL_DeviceBarrier();
			delete [] cDen;
			delete [] cDistEven;
			delete [] cDistOdd;
			MPI_Barrier(comm);
		}

		//.......................................................................
		// Once phase has been initialized, map solid to account for 'smeared' interface
		//for (i=0; i<N; i++)	Averages.SDs(i) -= (1.0);
		// Make sure the id match for the two domains
		for (i=0; i<N; i++)	Dm.id[i] = Mask.id[i];
		//.......................................................................
		// Finalize setup for averaging domain
		Averages->UpdateSolid();
		//.......................................................................		
		ScaLBL_D3Q19_Pressure(fq,Pressure,Np);
    	ScaLBL_D3Q19_Momentum(fq,Vel,Np);
		//...........................................................................
		// Copy the phase indicator field for the earlier timestep
		ScaLBL_DeviceBarrier();
		ScaLBL_CopyToHost(Averages->Phase_tplus.data(),Phi,N*sizeof(double));
		//...........................................................................
		// Copy the data for for the analysis timestep
		//...........................................................................
		// Copy the phase from the GPU -> CPU
		//...........................................................................
		ScaLBL_DeviceBarrier();
		ScaLBL_CopyToHost(Averages->Phase.data(),Phi,N*sizeof(double));
		ScaLBL_CopyToHost(Averages->Press.data(),Pressure,Np*sizeof(double));
		ScaLBL_CopyToHost(Averages->Vel_x.data(),&Velocity[0],Np*sizeof(double));
		ScaLBL_CopyToHost(Averages->Vel_y.data(),&Velocity[N],Np*sizeof(double));
		ScaLBL_CopyToHost(Averages->Vel_z.data(),&Velocity[2*N],Np*sizeof(double));
	
		double *TmpDat;
		TmpDat = new double [Np];
		ScaLBL_CopyToHost(&TmpDat[0],&Pressure[0], Np*sizeof(double));
		ScaLBL_Comm.RegularLayout(Map,TmpDat,Averages->Press.data());
		//...........................................................................

		if (rank==0) printf("********************************************************\n");
		if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);

		//.......create and start timer............
		double starttime,stoptime,cputime;
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		starttime = MPI_Wtime();
		//.........................................

		err = 1.0; 	
		double sat_w_previous = 1.01; // slightly impossible value!
		if (rank==0) printf("Begin timesteps: error tolerance is %f \n", tol);
		// Create the thread pool
		int N_threads = 4;
		if ( provided_thread_support < MPI_THREAD_MULTIPLE )
			N_threads = 0;
		if ( N_threads > 0 ) {
			// Set the affinity
			int N_procs = ThreadPool::getNumberOfProcessors();
			std::vector<int> procs(N_procs);
			for (int i=0; i<N_procs; i++)
				procs[i] = i;
			ThreadPool::setProcessAffinity(procs);
		}
		ThreadPool tpool(N_threads);

		// Create the MeshDataStruct
		fillHalo<double> fillData(Dm.Comm,Dm.rank_info,Nx-2,Ny-2,Nz-2,1,1,1,0,1);
		std::vector<IO::MeshDataStruct> meshData(1);
		meshData[0].meshName = "domain";
		meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(Dm.rank_info,Nx-2,Ny-2,Nz-2,Lx,Ly,Lz) );
		std::shared_ptr<IO::Variable> PhaseVar( new IO::Variable() );
		std::shared_ptr<IO::Variable> PressVar( new IO::Variable() );
		std::shared_ptr<IO::Variable> SignDistVar( new IO::Variable() );
		std::shared_ptr<IO::Variable> BlobIDVar( new IO::Variable() );
		PhaseVar->name = "phase";
		PhaseVar->type = IO::VariableType::VolumeVariable;
		PhaseVar->dim = 1;
		PhaseVar->data.resize(Nx-2,Ny-2,Nz-2);
		meshData[0].vars.push_back(PhaseVar);
		PressVar->name = "Pressure";
		PressVar->type = IO::VariableType::VolumeVariable;
		PressVar->dim = 1;
		PressVar->data.resize(Nx-2,Ny-2,Nz-2);
		meshData[0].vars.push_back(PressVar);
		SignDistVar->name = "SignDist";
		SignDistVar->type = IO::VariableType::VolumeVariable;
		SignDistVar->dim = 1;
		SignDistVar->data.resize(Nx-2,Ny-2,Nz-2);
		meshData[0].vars.push_back(SignDistVar);
		BlobIDVar->name = "BlobID";
		BlobIDVar->type = IO::VariableType::VolumeVariable;
		BlobIDVar->dim = 1;
		BlobIDVar->data.resize(Nx-2,Ny-2,Nz-2);
		meshData[0].vars.push_back(BlobIDVar);

		//************ MAIN ITERATION LOOP ***************************************/
		PROFILE_START("Loop");
		BlobIDstruct last_ids, last_index;
		BlobIDList last_id_map;
		writeIDMap(ID_map_struct(),0,id_map_filename);
		AnalysisWaitIdStruct work_ids;
		while (timestep < timestepMax && err > tol ) {
			//if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
			PROFILE_START("Update");
			// *************ODD TIMESTEP*************
			timestep++;
			// Compute the Phase indicator field
			// Read for Aq, Bq happens in this routine (requires communication)
			ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
			ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
			ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);

			// Halo exchange for phase field
			ScaLBL_Comm_Regular.SendHalo(Phi);
			ScaLBL_Comm_Regular.RecvHalo(Phi);
			
			// Perform the collision operation
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
			ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			// Set BCs
			if (typeBC > 0){
				ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
				ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
			}
			if (typeBC == 3){
				ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			if (typeBC == 4){
				din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

			// *************EVEN TIMESTEP*************
			timestep++;
			// Compute the Phase indicator field
			ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
			ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
			ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);

			// Halo exchange for phase field
			ScaLBL_Comm_Regular.SendHalo(Phi);
			ScaLBL_Comm_Regular.RecvHalo(Phi);

			// Perform the collision operation
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FORM NORMAL
			ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm.next, Np, Np);
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			// Set boundary conditions
			if (typeBC > 0){
				ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
				ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
			}
			if (typeBC == 3){
				ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			else if (typeBC == 4){
				din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Vel, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			//************************************************************************
			
			MPI_Barrier(comm);
			PROFILE_STOP("Update");

			// Run the analysis
			run_analysis(timestep,RESTART_INTERVAL,rank_info,*Averages,last_ids,last_index,last_id_map,
					Nx,Ny,Nz,pBC,beta,err,Phi,Pressure,Velocity,ID,f_even,f_odd,Den,
					LocalRestartFile,meshData,fillData,tpool,work_ids);

			// Save the timers
			if ( timestep%50==0 )
				PROFILE_SAVE("lbpm_color_simulator",1);
		}
		tpool.wait_pool_finished();
		PROFILE_STOP("Loop");
		//************************************************************************
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		stoptime = MPI_Wtime();
		if (rank==0) printf("-------------------------------------------------------------------\n");
		// Compute the walltime per timestep
		cputime = (stoptime - starttime)/timestep;
		// Performance obtained from each node
		double MLUPS = double(Nx*Ny*Nz)/cputime/1000000;

		if (rank==0) printf("********************************************************\n");
		if (rank==0) printf("CPU time = %f \n", cputime);
		if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
		MLUPS *= nprocs;
		if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
		if (rank==0) printf("********************************************************\n");

		// ************************************************************************

		PROFILE_STOP("Main");
		PROFILE_SAVE("lbpm_color_simulator",1);
		// ****************************************************
		MPI_Barrier(comm);
	} // Limit scope so variables that contain communicators will free before MPI_Finialize
	MPI_Comm_free(&comm);
	MPI_Finalize();
}


