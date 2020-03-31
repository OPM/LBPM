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

#include "lbpm_nonnewtonian_simulator.h"

//#define WRITE_SURFACES

/*
 * Simulator for single-phase non-newtonian flow
 * James E. McClure 2013-2014 & Christopher P. Fowler 2017
 */

using namespace std;

//*************************************************************************
// Implementation of Steady State Single-Phase LBM for permeability measurement
//*************************************************************************
inline void PackID(int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

inline void UnpackID(int *list, int count, char *recvbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}

//***************************************************************************************

inline void ZeroHalo(double *Data, int Nx, int Ny, int Nz)
{
	int i,j,k,n;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			i=0;
			n = k*Nx*Ny+j*Nx+i;
			Data[2*n] = 0.0;
			Data[2*n+1] = 0.0;
			i=Nx-1;
			n = k*Nx*Ny+j*Nx+i;
			Data[2*n] = 0.0;
			Data[2*n+1] = 0.0;
		}
	}

	for (k=0;k<Nz;k++){
		for (i=0;i<Nx;i++){
			j=0;
			n = k*Nx*Ny+j*Nx+i;
			Data[2*n] = 0.0;
			Data[2*n+1] = 0.0;
			j=Ny-1;
			n = k*Nx*Ny+j*Nx+i;
			Data[2*n] = 0.0;
			Data[2*n+1] = 0.0;
		}
	}

	for (j=0;j<Ny;j++){
		for (i=0;i<Nx;i++){
			k=0;
			n = k*Nx*Ny+j*Nx+i;
			Data[2*n] = 0.0;
			Data[2*n+1] = 0.0;
			k=Nz-1;
			n = k*Nx*Ny+j*Nx+i;
			Data[2*n] = 0.0;
			Data[2*n+1] = 0.0;
		}
	}
}
//***************************************************************************************


int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	//MPI_Init(&argc,&argv);

	/*
	 * Definitely seems to be an issue - let's hope James gets back to me...
	 */
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
		//*****************************************
		// MPI ranks for all 18 neighbors
		//**********************************
		//		int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
		//		int rank_xy,rank_XY,rank_xY,rank_Xy;
		//		int rank_xz,rank_XZ,rank_xZ,rank_Xz;
		//		int rank_yz,rank_YZ,rank_yZ,rank_Yz;
		//**********************************
		MPI_Request req1[18],req2[18];
		MPI_Status stat1[18],stat2[18];

		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Single Phase Non-Newtonian Calculation \n");
			printf("********************************************************\n");
		}

		// Variables that specify the computational domain
		string FILENAME;
		int Nx,Ny,Nz;		// local sub-domain size
		int nspheres;		// number of spheres in the packing
		double Lx,Ly,Lz;	// Domain length
		double D = 1.0;		// reference length for non-dimensionalization
		// Color Model parameters
		int timestepMax, interval;
		double tau,Fx,Fy,Fz,tol,err;
		double din,dout;
		bool pBC,Restart;
		int i,j,k,n;


		/*
		 * Analysis flags
		 */
		int RESTART_INTERVAL=20000;
		int BLOB_ANALYSIS_INTERVAL=1000;
		int timestep = -1;

		/*
		 *  Read file data in
		 *
		 */

		if (rank==0){
			//.............................................................
			//		READ SIMULATION PARMAETERS FROM INPUT FILE
			//.............................................................
			ifstream input("Permeability.in");
			// Line 1: model parameters (tau, alpha, beta, das, dbs)
			input >> tau;			// Viscosity parameter
			// Line 2: External force components (Fx,Fy, Fz)
			input >> Fx;
			input >> Fy;
			input >> Fz;
			// Line 3: Pressure Boundary conditions
			input >> Restart;
			input >> pBC;
			input >> din;
			input >> dout;
			// Line 4: time-stepping criteria
			input >> timestepMax;		// max no. of timesteps
			input >> interval;			// restart interval
			input >> tol;				// error tolerance
			//.............................................................

			//.......................................................................
			// Reading the domain information file
			//.......................................................................
			ifstream domain("Domain.in");
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

		/*
		 * Broadcast
		 */

		// **************************************************************
		// Broadcast simulation parameters from rank 0 to all other procs
		MPI_Barrier(comm);
		//.................................................
		MPI_Bcast(&tau,1,MPI_DOUBLE,0,comm);
		//MPI_Bcast(&pBC,1,MPI_LOGICAL,0,comm);
		//	MPI_Bcast(&Restart,1,MPI_LOGICAL,0,comm);
		MPI_Bcast(&din,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&dout,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&timestepMax,1,MPI_INT,0,comm);
		MPI_Bcast(&interval,1,MPI_INT,0,comm);
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
		MPI_Barrier(comm);

		//?
		RESTART_INTERVAL=interval;
		// **************************************************************
		// **************************************************************

		/*
		 * Set up rank info struct
		 */

		const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

		MPI_Barrier(comm);

		/*
		 * Set up the relaxation rates and STATIC VISCOSITY
		 */

		double rlxA = 1.f/tau;
		double rlxB = 8.f*(2.f-rlxA)/(8.f-rlxA);
		double viscosity=(tau-0.5)/3.0;

		/*
		 *  Debug block 1
		 */
		printf("\npBC=%d (an int) \n",pBC);
		printf("viscosity=%f\n",viscosity);

		/*
		 * Check processor counts
		 */

		if (nprocs != nprocx*nprocy*nprocz){
			printf("nprocx =  %i \n",nprocx);
			printf("nprocy =  %i \n",nprocy);
			printf("nprocz =  %i \n",nprocz);
			INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
		}

		/*
		 * Display what we've got thus far
		 */

		if (rank==0){
			printf("********************************************************\n");
			printf("tau = %f \n", tau);
			printf("Force(x) = %f \n", Fx);
			printf("Force(y) = %f \n", Fy);
			printf("Force(z) = %f \n", Fz);
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("Process grid = %i x %i x %i\n",nprocx,nprocy,nprocz);
			printf("********************************************************\n");
		}

		/*
		 *  Initialized domain and averaging framework for Two-Phase flow
		 */

		// not needed right now

		// Initialized domain and averaging framework for Two-Phase Flow
		int BC=pBC;

		printf("BC=pBC=%d (an int)\n",BC);

		Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);                    /*    1  */
		for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = 1;
		std::shared_ptr<TwoPhase> Averages( new TwoPhase(Dm) );
		Dm.CommInit();																/*   2 */

		Domain Mask(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);

		//		TwoPhase Averages(Dm);
		//
		//		InitializeRanks( rank, nprocx, nprocy, nprocz, iproc, jproc, kproc,
		//				rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z,
		//				rank_xy, rank_XY, rank_xY, rank_Xy, rank_xz, rank_XZ, rank_xZ, rank_Xz,
		//				rank_yz, rank_YZ, rank_yZ, rank_Yz );

		MPI_Barrier(comm);
		Nx+=2; Ny+=2; Nz+=2;
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

		/*
		 * Debug block 2
		 */
		printf("LocalRankString=%s\n",LocalRankString);
		printf("LocalRestartFile=%s\n",LocalRestartFile);
		printf("LocalRankFilename=%s\n",LocalRankFilename);

		// .......... READ THE INPUT FILE .......................................
		//	char value;
		char *id;
		id = new char[N];
		int sum = 0;
		double sum_local;
		double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
		if (pBC) {
			printf("tripped if (pBC) statement\n");
			iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
		}
		double porosity, pore_vol;


		/*
		 * Debug block 3
		 */
		printf("iVol_global=%f\n",iVol_global);

		//...........................................................................
		if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;

		//.......................................................................
		sprintf(LocalRankString,"%05d",rank);
		//	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
		//	WriteLocalSolidID(LocalRankFilename, id, N);
		sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
		ReadBinaryFile(LocalRankFilename, Averages->SDs.data(), N);
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;                                               /*    3      */

		//.......................................................................
		// Assign the phase ID field based on the signed distance
		//.......................................................................
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 0;
				}
			}
		}
		sum=0;
		pore_vol = 0.0;
		for ( k=1;k<Nz-1;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (Averages->SDs(n) > 0.0){
						id[n] = 2;
					}
					// compute the porosity (actual interface location used)
					if (Averages->SDs(n) > 0.0){
						sum++;
					}
				}
			}
		}																							/*   4     */

		/*
		 *  Initialize from segmented data 															     5
		 *
		 */

		// not needed at the moment


		/*
		 * Debug block 4
		 */
		printf("sum=%d\n",sum);

		// Set up kstart, kfinish so that the reservoirs are excluded from averaging
		int kstart,kfinish;
		kstart = 1;
		kfinish = Nz-1;
		if (pBC && kproc==0)		kstart = 4;
		if (pBC && kproc==nprocz-1)	kfinish = Nz-4;

		// Compute the pore volume
		sum_local = 0.0;
		for ( k=kstart;k<kfinish;k++){
			for ( j=1;j<Ny-1;j++){
				for ( i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (id[n] > 0){
						sum_local += 1.0;
					}
				}
			}
		}

		MPI_Allreduce(&sum_local,&pore_vol,1,MPI_DOUBLE,MPI_SUM,comm);					/*    6     */
		//MPI_Allreduce(&sum_local,&porosity,1,MPI_DOUBLE,MPI_SUM,comm);
		porosity = pore_vol*iVol_global;

		if (rank==0) printf("Media porosity = %f \n",porosity);
		//.........................................................
		// If external boundary conditions are applied remove solid
		if (pBC && kproc == 0){
			printf("Tripped if (pcB && kproc == 0)\n");
			for (k=0; k<3; k++){
				for (j=0;j<Ny;j++){
					for (i=0;i<Nx;i++){
						n = k*Nx*Ny+j*Nx+i;
						id[n] = 1;
						Averages->SDs(n) = max(Averages->SDs(n),1.0*(2.5-k));
					}
				}
			}
		}
		if (pBC && kproc == nprocz-1){
			printf("Tripped if (pcB && kproc == nprocz-1)\n");
			for (k=Nz-3; k<Nz; k++){
				for (j=0;j<Ny;j++){
					for (i=0;i<Nx;i++){
						n = k*Nx*Ny+j*Nx+i;
						id[n] = 2;
						Averages->SDs(n) = max(Averages->SDs(n),1.0*(k-Nz+2.5));
					}
				}
			}
		}

		//.........................................................
		// don't perform computations at the eight corners
		id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
		id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;     /*    7   */
		//.........................................................

		/*
		 *  To use a mask or not - that is the question!
		 */
		// maybe mask?																									/*   8  */

		// Initialize communication structures in averaging domain
		//		for (i=0; i<Dm.Nx*Dm.Ny*Dm.Nz; i++) Dm.id[i] = id[i];
		//		Dm.CommInit();
		for (i=0; i<Mask.Nx*Mask.Ny*Mask.Nz; i++) Mask.id[i] = id[i];
		Mask.CommInit(comm);

		//...........................................................................
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");
		// Create a communicator for the device
		//		ScaLBL_Communicator ScaLBL_Comm(Dm);
		ScaLBL_Communicator ScaLBL_Comm(Mask);																		/*   9   */

		// set reservoirs   (not needed, right?)
		if (pBC > 0){
			for ( k=0;k<Nz;k++){
				for ( j=0;j<Ny;j++){
					for ( i=0;i<Nx;i++){
						int n = k*Nx*Ny+j*Nx+i;
						if (Dm.kproc==0 && k==0)			id[n]=1;
						if (Dm.kproc==0 && k==1)			id[n]=1;
						if (Dm.kproc==nprocz-1 && k==Nz-2)	id[n]=2;
						if (Dm.kproc==nprocz-1 && k==Nz-1)	id[n]=2;
						Mask.id[n] = id[n];
					}
				}
			}
		}																									/* 10 */

		//...........device phase ID.................................................
		if (rank==0)	printf ("Copying phase ID to device \n");
		char *ID;
		ScaLBL_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					int n = k*Nx*Ny+j*Nx+i;
					if (i==0 || i==Nx-1 || j==0 || j==Ny-1 || k==0 || k==Nz-1)	id[n] = 0;
				}
			}
		}
		// Copy to the device
		ScaLBL_CopyToDevice(ID, id, N);
		ScaLBL_DeviceBarrier();
		//...........................................................................                       11
		//...........................................................................
		//				MAIN  VARIABLES ALLOCATED HERE
		//...........................................................................
		// LBM variables
		if (rank==0)	printf ("Allocating distributions \n");
		//......................device distributions.................................
		double *f_even,*f_odd;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
		ScaLBL_AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory
		//...........................................................................
		double *Velocity, *Pressure, *dvcSignDist;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &Pressure, dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &dvcSignDist, dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*dist_mem_size);
		//...........................................................................

		// Copy signed distance for device initialization
		ScaLBL_CopyToDevice(dvcSignDist, Averages->SDs.data(), dist_mem_size);
		//...........................................................................

		int logcount = 0; // number of surface write-outs												/*   12   */

		/*
		 *  Display simulation metrics:
		 */
		if (rank == 0) {
			printf("Displaying simulation metrics... \n");
		}


		//...........................................................................
		//				MAIN  VARIABLES INITIALIZED HERE
		//...........................................................................
		//...........................................................................
		if (rank==0)	printf("Setting the distributions, size = %i\n", N);
		//...........................................................................
		ScaLBL_DeviceBarrier();
		ScaLBL_D3Q19_Init(ID, f_even, f_odd, Nx, Ny, Nz);
		ScaLBL_DeviceBarrier();
		//......................................................................        /*   13   */


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
					timestep=5;
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
			// ScaLBL_CopyToDevice(Den,cDen,2*N*sizeof(double));  /*    Two-phase stuff   */
			ScaLBL_DeviceBarrier();
			delete [] cDen;
			delete [] cDistEven;
			delete [] cDistOdd;
			MPI_Barrier(comm);
		}                                                                                   /*  14 */

//		//......................................................................
//		ScaLBL_D3Q7_Init(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
//		ScaLBL_D3Q7_Init(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
//		ScaLBL_DeviceBarrier();
//		MPI_Barrier(comm);																/*  15  */

		//.......................................................................
		// Once phase has been initialized, map solid to account for 'smeared' interface
		//for (i=0; i<N; i++)	Averages->SDs(i) -= (1.0);
		// Make sure the id match for the two domains
		for (i=0; i<N; i++)	Dm.id[i] = Mask.id[i];
		//.......................................................................	/* 16 */

		//.......................................................................
		// Finalize setup for averaging domain
		Averages->UpdateSolid();													/*    17   */



//		//.......................................................................
//
//		//*************************************************************************
//		// 		Compute the phase indicator field and reset Copy, Den
//		//*************************************************************************
//		ScaLBL_ComputePhaseField(ID, Phi, Den, N);
//		//*************************************************************************
//		ScaLBL_DeviceBarrier();
//		ScaLBL_Comm.SendHalo(Phi);
//		ScaLBL_Comm.RecvHalo(Phi);
//		ScaLBL_DeviceBarrier();
//		MPI_Barrier(comm);
//		//*************************************************************************   /*  18  */


			if (rank==0 && pBC){
				printf("Setting inlet pressure = %f \n", din);
				printf("Setting outlet pressure = %f \n", dout);
			}
			if (pBC && kproc == 0)	{
				ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
			}

			if (pBC && kproc == nprocz-1){
				ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
			}																				/* 19   */

			// will have to fill in stuff for 2 phase later

			//timestepMax = 500;

			if (rank==0) printf("********************************************************\n");
			if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);

			//...........................................................................
			// Copy the data for for the analysis timestep
			//...........................................................................
			// Copy the phase from the GPU -> CPU
			//...........................................................................
			ScaLBL_DeviceBarrier();
			ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
//			ScaLBL_CopyToHost(Averages->Phase.data(),Phi,N*sizeof(double));
			ScaLBL_CopyToHost(Averages->Press.data(),Pressure,N*sizeof(double));
			ScaLBL_CopyToHost(Averages->Vel_x.data(),&Velocity[0],N*sizeof(double));
			ScaLBL_CopyToHost(Averages->Vel_y.data(),&Velocity[N],N*sizeof(double));
			ScaLBL_CopyToHost(Averages->Vel_z.data(),&Velocity[2*N],N*sizeof(double));
			//...........................................................................	/* 20 */

			//.......create and start timer............
			double starttime,stoptime,cputime;
			MPI_Barrier(comm);
			starttime = MPI_Wtime();

			/*
			 *  Create the thread pool
			 *
			 */

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
			ThreadPool tpool(N_threads);											/* 21 */

			printf("N_threads=%d\n",N_threads);

			/*
			 *   Create the MeshDataStruct - seems particularly important
			 */

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
			meshData[0].vars.push_back(BlobIDVar);									/* 22 */











			//.........................................

			double D32,Fo,Re,velocity,err1D,mag_force,vel_prev;
			err = vel_prev = 1.0;
			if (rank==0) printf("Begin timesteps: error tolerance is %f \n", tol);
			//************ MAIN ITERATION LOOP ***************************************/














			BlobIDstruct last_ids, last_index;
			BlobIDList last_id_map;
			writeIDMap(ID_map_struct(),0,id_map_filename);
			AnalysisWaitIdStruct work_ids;








			double beta = 0;
//			double * Phi;
//			double * Den;




			while (timestep < timestepMax && err > tol ){

				//*************************************************************************
				// Fused Color Gradient and Collision
				//*************************************************************************
				ScaLBL_D3Q19_MRT(ID,f_even,f_odd,rlxA,rlxB,Fx,Fy,Fz,Nx,Ny,Nz);
				//*************************************************************************
				// Pack and send the D3Q19 distributions
				ScaLBL_Comm.SendD3Q19(f_even, f_odd);
				//*************************************************************************
				// Swap the distributions for momentum transport
				//*************************************************************************
				ScaLBL_D3Q19_Swap(ID, f_even, f_odd, Nx, Ny, Nz);
				//*************************************************************************
				// Wait for communications to complete and unpack the distributions
				ScaLBL_Comm.RecvD3Q19(f_even, f_odd);
				//*************************************************************************

				if (pBC && kproc == 0)	{
					ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
				}

				if (pBC && kproc == nprocz-1){
					ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
				}
				//...................................................................................
				ScaLBL_DeviceBarrier();
				MPI_Barrier(comm);

				// Timestep completed!
				timestep++;

				/// Perform the analysis
				run_analysis(timestep,RESTART_INTERVAL,rank_info,*Averages,last_ids,last_index,last_id_map,
						Nx,Ny,Nz,pBC,err,Pressure,Velocity,ID,f_even,f_odd,
						LocalRestartFile,meshData,fillData,tpool,work_ids);


			}
			//************************************************************************/
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

			NULL_USE(RESTART_INTERVAL);
		}
		MPI_Barrier(comm);
		MPI_Finalize();
	 //****************************************************
}



















// Scrap

// if (rank==0){
// ************* DIMENSIONLESS FORCHEIMER EQUATION *************************
//  Dye, A.L., McClure, J.E., Gray, W.G. and C.T. Miller
//  Description of Non-Darcy Flows in Porous Medium Systems
//  Physical Review E 87 (3), 033012
//  Fo := density*D32^3*(density*force) / (viscosity^2)
//	Re := density*D32*velocity / viscosity
//  Fo = a*Re + b*Re^2
// *************************************************************************
//viscosity = (tau-0.5)*0.333333333333333333;
/*
 * Original formula for D32: D32 = 6.0*(Dm.Volume-Averages.vol_w_global)/Averages.As_global;
 */
//	D32 = 6.0*(Dm.Volume-Averages.vol_w_global)/Averages.As_global;
//		    	printf("Dm.Volume=%f Averages.vol_w_global=%f Averages.As_global=%f \n",Dm.Volume,Averages.vol_w_global,Averages.As_global);
//		    	D32 = 6.0*(Dm.Volume-Averages.vol_w_global);
//				printf("Sauter Mean Diameter = %f \n",D32);
//	mag_force = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
//				Fo = D32*D32*D32*mag_force/viscosity/viscosity;
// .... 1-D flow should be aligned with force ...
//	velocity = vawx*Fx/mag_force + vawy*Fy/mag_force + vawz*Fz/mag_force;
//	err1D = fabs(velocity-sqrt(vawx*vawx+vawy*vawy+vawz*vawz))/velocity;
//.......... Computation of the Reynolds number Re ..............
//				Re = D32*velocity/viscosity;
//	printf("Force: %.5g,%.5g,%.5g \n",Fx,Fy,Fz);
//	printf("Velocity: %.5g,%.5g,%.5g \n",vawx,vawy,vawz);
//	printf("Relative error for 1D representation: %.5g \n",err1D);
//				printf("Dimensionless force: %5g \n", Fo);
//				printf("Reynolds number: %.5g \n", Re);
//				printf("Dimensionless Permeability (k/D^2): %.5g \n", Re/Fo);
// }




/*
 * 		if (timestep%5 == 0){
			//...........................................................................
			// Copy the data for for the analysis timestep
			//...........................................................................
			// Copy the phase from the GPU -> CPU
			//...........................................................................
			ScaLBL_DeviceBarrier();
			ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
			ScaLBL_D3Q19_Velocity(ID,f_even,f_odd,Velocity,Nx,Ny,Nz);
			ScaLBL_CopyToHost(Averages.Press.data(),Pressure,N*sizeof(double));
			ScaLBL_CopyToHost(Averages.Vel_x.data(),&Velocity[0],N*sizeof(double));
			ScaLBL_CopyToHost(Averages.Vel_y.data(),&Velocity[N],N*sizeof(double));
			ScaLBL_CopyToHost(Averages.Vel_z.data(),&Velocity[2*N],N*sizeof(double));

			// Way more work than necessary -- this is just to get the solid interfacial area!!
			Averages.Initialize();
			Averages.UpdateMeshValues();
			Averages.ComputeLocal();
			Averages.Reduce();

			double vawx = -Averages.vaw_global(0);
			double vawy = -Averages.vaw_global(1);
			double vawz = -Averages.vaw_global(2);


		    if (rank==0){
				mag_force = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
				// .... 1-D flow should be aligned with force ...
				velocity = vawx*Fx/mag_force + vawy*Fy/mag_force + vawz*Fz/mag_force;
				err1D = fabs(velocity-sqrt(vawx*vawx+vawy*vawy+vawz*vawz))/velocity;
				//printf("Force: %.5g,%.5g,%.5g \n",Fx,Fy,Fz);
				printf("vel_z=%.5g\n",vawz);
				//printf("Velocity: %.5g,%.5g,%.5g \n",vawx,vawy,vawz);
				printf("Relative error for 1D representation: %.5g \n",err1D);
			}

		}
 */



//		// Initialize two phase flow variables (all wetting phase)
	//		for (k=0;k<Nz;k++){
	//			for (j=0;j<Ny;j++){
	//				for (i=0;i<Nx;i++){
	//					n=k*Nx*Ny+j*Nx+i;
	//					Averages.Phase(i,j,k) = -1.0;
	//					Averages.SDn(i,j,k) = Averages.Phase(i,j,k);
	//					Averages.Phase_tplus(i,j,k) = Averages.SDn(i,j,k);
	//					Averages.Phase_tminus(i,j,k) = Averages.SDn(i,j,k);
	//					Averages.DelPhi(i,j,k) = 0.0;
	//					Averages.Press(i,j,k) = 0.0;
	//					Averages.Vel_x(i,j,k) = 0.0;
	//					Averages.Vel_y(i,j,k) = 0.0;
	//					Averages.Vel_z(i,j,k) = 0.0;
	//				}
	//			}
	//		}
	//
	//		//.......................................................................
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
