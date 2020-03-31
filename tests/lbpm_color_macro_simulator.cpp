#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "analysis/runAnalysis.h"
#include "common/MPI_Helpers.h"
#include "ProfilerApp.h"
#include "threadpool/thread_pool.h"

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
		// Initialize compute device
		//		int device=ScaLBL_SetDevice(rank);
		//printf("Using GPU ID %i for rank %i \n",device,rank);
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);

		PROFILE_ENABLE(1);
		//PROFILE_ENABLE_TRACE();
		//PROFILE_ENABLE_MEMORY();
		PROFILE_SYNCHRONIZE();
		PROFILE_START("Main");
		Utilities::setErrorHandlers();

		int ANALYSIS_INTERVAL = 1000;
		int BLOBID_INTERVAL = 1000;
        std::string analysis_method = "independent";
		if (argc >= 3) {
			ANALYSIS_INTERVAL = atoi(argv[1]);
			BLOBID_INTERVAL   = atoi(argv[2]);
		}
		if (argc >= 4)
			analysis_method = std::string(argv[3]);

		// Variables that specify the computational domain  
		string FILENAME;
		int Nx,Ny,Nz,Np;		// local sub-domain size
		double Lx,Ly,Lz;	// Domain length
		double D = 1.0;		// reference length for non-dimensionalization
		// Color Model parameters
		int timestepMax;
		double tauA, tauB, rhoA,rhoB;
		double Fx,Fy,Fz,tol,err;
		double alpha, beta;
		int BoundaryCondition;
		int InitialCondition;
		//	bool pBC,Restart;
		int i,j,k,n;
		double din, dout, flux;
		double inletA,inletB,outletA,outletB;
		inletA=1.f;
		inletB=0.f;
		outletA=0.f;
		outletB=1.f;
		flux = 10.f;
		dout=1.f;

		int RESTART_INTERVAL=20000;
		//int ANALYSIS_)INTERVAL=1000;	
		int BLOB_ANALYSIS_INTERVAL=1000;
		int timestep = 6;

		if (rank==0){
			//.............................................................
			//		READ SIMULATION PARMAETERS FROM INPUT FILE
			//.............................................................
			ifstream input("Color.in");
			if (input.is_open()){
				// Line 1: model parameters (tau, alpha, beta, das, dbs)
				input >> tauA;			// Viscosity non-wetting
				input >> tauB;			// Viscosity wetting
				input >> rhoA;			// density non-wetting
				input >> rhoB;			// density wetting
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
				tauA = tauB = 1.0;
				rhoA = rhoB = 1.0;
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
				Lx=Ly=Lz=1.0;
			}
		}
		// **************************************************************
		// Broadcast simulation parameters from rank 0 to all other procs
		MPI_Barrier(comm);
		//.................................................
		MPI_Bcast(&tauA,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&tauB,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&rhoA,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&rhoB,1,MPI_DOUBLE,0,comm);
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
		MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
		//.................................................

		flux = 0.f;
		if (BoundaryCondition==4) flux = din*rhoA; // mass flux must adjust for density (see formulation for details

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
			printf("tau (non-wetting) = %f \n", tauA);
			printf("tau (wetting) = %f \n", tauB);
			printf("density (non-wetting) = %f \n", rhoA);
			printf("density (wetting) = %f \n", rhoB);
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
		Dm.CommInit();

		// Mask that excludes the solid phase
		Domain Mask(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);
		MPI_Barrier(comm);

		Nx+=2; Ny+=2; Nz += 2;
		int N = Nx*Ny*Nz;
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
		double sum, sum_local;
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

		if (rank==0) printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
		sprintf(LocalRankFilename,"ID.%05i",rank);
		size_t readID;
		FILE *IDFILE = fopen(LocalRankFilename,"rb");
		if (IDFILE==NULL) ERROR("lbpm_color_simulator: Error opening file: ID.xxxxx");
		readID=fread(id,1,N,IDFILE);
		if (readID != size_t(N)) printf("lbpm_color_simulator: Error reading ID (rank=%i) \n",rank);
		fclose(IDFILE);
		
		// Read id from restart
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
			FILE *RESTART = fopen(LocalRestartFile,"rb");
			if (IDFILE==NULL) ERROR("lbpm_color_simulator: Error opening file: Restart.xxxxx");
			readID=fread(id,1,N,RESTART);
			if (readID != size_t(N)) printf("lbpm_color_simulator: Error reading Restart (rank=%i) \n",rank);
			fclose(RESTART);
			/*
			// Read in the restart file to CPU buffers
			double *cDen = new double[2*Np];
			double *cfq = new double[19*Np];
			ReadCheckpoint(LocalRestartFile, cDen, cfq, Np);
			// Copy the restart data to the GPU
			ScaLBL_CopyToDevice(fq,cfq,19*Np*sizeof(double));
			ScaLBL_CopyToDevice(Den,cDen,2*Np*sizeof(double));
			ScaLBL_DeviceBarrier();
			delete [] cDen;
			delete [] cfq;
			*/
			MPI_Barrier(comm);
		}
		
		fflush(stdout);
		//.......................................................................
		// Compute the media porosity, assign phase labels and solid composition
		//.......................................................................
		sum_local=0.0;
		Np=0;  // number of local pore nodes
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
		// If external boundary conditions are applied remove solid
		if (BoundaryCondition >  0  && Dm.kproc() == 0){
			for (k=0; k<3; k++){
				for (j=0;j<Ny;j++){
					for (i=0;i<Nx;i++){
						int n = k*Nx*Ny+j*Nx+i;
						//id[n] = 1;
						Averages->SDs(n) = max(Averages->SDs(n),1.0*(2.5-k));
					}					
				}
			}
		}
		if (BoundaryCondition >  0  && Dm.kproc() == nprocz-1){
			for (k=Nz-3; k<Nz; k++){
				for (j=0;j<Ny;j++){
					for (i=0;i<Nx;i++){
						int n = k*Nx*Ny+j*Nx+i;
						//id[n] = 2;
						Averages->SDs(n) = max(Averages->SDs(n),1.0*(k-Nz+2.5));
					}					
				}
			}
		}
		//.........................................................
		// don't perform computations at the eight corners
		id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
		id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
		//.........................................................

		// Initialize communication structures in averaging domain
		for (i=0; i<Mask.Nx*Mask.Ny*Mask.Nz; i++) Mask.id[i] = id[i];
		Mask.CommInit(comm);
		double *PhaseLabel;
		PhaseLabel = new double[N];
		Mask.AssignComponentLabels(PhaseLabel);
		
		fflush(stdout);
		//...........................................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device (will use optimized layout)
		ScaLBL_Communicator ScaLBL_Comm(Mask);
		//Create a second communicator based on the regular data layout
		ScaLBL_Communicator ScaLBL_Comm_Regular(Mask);
		
		int Npad=Np+32;
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Npad];
		Np = ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Mask.id,Np);
		if (rank==0)	printf ("Set up memory efficient layout Npad=%i, Np=%i \n",Npad,Np);
		MPI_Barrier(comm);
		//...........................................................................
		//				MAIN  VARIABLES ALLOCATED HERE
		//...........................................................................
		// LBM variables
		if (rank==0)	printf ("Allocating distributions \n");
		//......................device distributions.................................
		fflush(stdout);
		int dist_mem_size = Np*sizeof(double);
		int neighborSize=18*(Np*sizeof(int));

		int *NeighborList;
		int *dvcMap;
		double *fq, *Aq, *Bq;
		double *Den, *Phi;
		double *ColorGrad;
		double *Velocity;
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
		ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*sizeof(double)*Np);
		
		//...........................................................................
		// Update GPU data structures
		if (rank==0)	printf ("Setting up device map and neighbor list \n");
		fflush(stdout);
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
		// check that TmpMap is valid
		for (int idx=0; idx<ScaLBL_Comm.last_interior; idx++){
			if (idx == ScaLBL_Comm.next) idx = ScaLBL_Comm.first_interior;
			int n = TmpMap[idx];
			if (n > Nx*Ny*Nz){
				printf("Bad value! idx=%i \n");
				TmpMap[idx] = Nx*Ny*Nz-1;
			}
		}
		ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		delete [] TmpMap;
		
		// copy the neighbor list 
		ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
		// initialize phi based on PhaseLabel (include solid component labels)
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
		//...........................................................................

		if (rank==0)	printf ("Initializing distributions \n");
		ScaLBL_D3Q19_Init(fq, Np);
		if (rank==0)	printf ("Initializing phase field \n");
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, 0, ScaLBL_Comm.next, Np);
		ScaLBL_PhaseField_Init(dvcMap, Phi, Den, Aq, Bq, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);

		if (BoundaryCondition >0 ){
			if (Dm.kproc()==0){
				ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,0);
				ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,1);
				ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,2);
			}
			if (Dm.kproc() == nprocz-1){
				ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-1);
				ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-2);
				ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,Nz-3);
			}
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
		ScaLBL_D3Q19_Momentum(fq,Velocity,Np);
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
		ScaLBL_Comm.RegularLayout(Map,Pressure,Averages->Press);
		ScaLBL_Comm.RegularLayout(Map,&Velocity[0],Averages->Vel_x);
		ScaLBL_Comm.RegularLayout(Map,&Velocity[Np],Averages->Vel_y);
		ScaLBL_Comm.RegularLayout(Map,&Velocity[2*Np],Averages->Vel_z);
		//...........................................................................

		if (rank==0){
			printf("********************************************************\n");
			printf("No. of timesteps: %i \n", timestepMax);
			fflush(stdout);
		}

		//.......create and start timer............
		double starttime,stoptime,cputime;
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		starttime = MPI_Wtime();
		//.........................................

		err = 1.0; 	
		double sat_w_previous = 1.01; // slightly impossible value!
		if (rank==0) printf("Begin timesteps: error tolerance is %f \n", tol);
		if (rank==0){
		  printf("Analysis intervals: (restart) %i, (TCAT) %i, (blobtracking) %i \n",RESTART_INTERVAL,ANALYSIS_INTERVAL,BLOBID_INTERVAL);
		}

		//************ MAIN ITERATION LOOP ***************************************/
		PROFILE_START("Loop");
        std::shared_ptr<Database> analysis_db;
        runAnalysis analysis( analysis_db, rank_info, ScaLBL_Comm, Dm, Np, pBC, beta, Map );
        analysis.createThreads( analysis_method, 4 );
		while (timestep < timestepMax && err > tol ) {
			//if ( rank==0 ) { printf("Running timestep %i (%i MB)\n",timestep+1,(int)(Utilities::getMemoryUsage()/1048576)); }
			PROFILE_START("Update");
			// *************ODD TIMESTEP*************
			timestep++;
			// Compute the Phase indicator field
			// Read for Aq, Bq happens in this routine (requires communication)
			ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
			ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
			ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
			ScaLBL_D3Q7_AAodd_PhaseField(NeighborList, dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);
			
			// Perform the collision operation
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FROM NORMAL
			// Halo exchange for phase field
			ScaLBL_Comm_Regular.SendHalo(Phi);

			ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
			ScaLBL_Comm_Regular.RecvHalo(Phi);
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			// Set BCs
			if (BoundaryCondition > 0){
				ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
				ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
			}
			if (BoundaryCondition == 3){
				ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			if (BoundaryCondition == 4){
				din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			ScaLBL_D3Q19_AAodd_Color(NeighborList, dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);

			// *************EVEN TIMESTEP*************
			timestep++;
			// Compute the Phase indicator field
			ScaLBL_Comm.BiSendD3Q7AA(Aq,Bq); //READ FROM NORMAL
			ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
			ScaLBL_Comm.BiRecvD3Q7AA(Aq,Bq); //WRITE INTO OPPOSITE
			ScaLBL_D3Q7_AAeven_PhaseField(dvcMap, Aq, Bq, Den, Phi, 0, ScaLBL_Comm.next, Np);

			// Perform the collision operation
			ScaLBL_Comm.SendD3Q19AA(fq); //READ FORM NORMAL
			// Halo exchange for phase field
			ScaLBL_Comm_Regular.SendHalo(Phi);
			ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz,  Nx, Nx*Ny, ScaLBL_Comm.first_interior, ScaLBL_Comm.last_interior, Np);
			ScaLBL_Comm_Regular.RecvHalo(Phi);
			ScaLBL_Comm.RecvD3Q19AA(fq); //WRITE INTO OPPOSITE
			// Set boundary conditions
			if (BoundaryCondition > 0){
				ScaLBL_Comm.Color_BC_z(dvcMap, Phi, Den, inletA, inletB);
				ScaLBL_Comm.Color_BC_Z(dvcMap, Phi, Den, outletA, outletB);
			}
			if (BoundaryCondition == 3){
				ScaLBL_Comm.D3Q19_Pressure_BC_z(NeighborList, fq, din, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			else if (BoundaryCondition == 4){
				din = ScaLBL_Comm.D3Q19_Flux_BC_z(NeighborList, fq, flux, timestep);
				ScaLBL_Comm.D3Q19_Pressure_BC_Z(NeighborList, fq, dout, timestep);
			}
			ScaLBL_D3Q19_AAeven_Color(dvcMap, fq, Aq, Bq, Den, Phi, Velocity, rhoA, rhoB, tauA, tauB,
					alpha, beta, Fx, Fy, Fz, Nx, Nx*Ny, 0, ScaLBL_Comm.next, Np);
			ScaLBL_DeviceBarrier(); MPI_Barrier(comm);
			//************************************************************************
			
			MPI_Barrier(comm);
			PROFILE_STOP("Update");

			// Run the analysis
            analysis.run( timestep, *Averages, Phi, Pressure, Velocity, fq, Den );

		}
        analysis.finish();
		PROFILE_STOP("Loop");
		PROFILE_SAVE("lbpm_color_simulator",1);
		//************************************************************************
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		stoptime = MPI_Wtime();
		if (rank==0) printf("-------------------------------------------------------------------\n");
		// Compute the walltime per timestep
		cputime = (stoptime - starttime)/timestep;
		// Performance obtained from each node
		double MLUPS = double(Np)/cputime/1000000;

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


