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

//#define WRE_SURFACES

/*
 * Simulator for two-phase flow in porous media
 * James E. McClure 2013-2014
 */

using namespace std;

//*************************************************************************
// Implementation of Two-Phase Immiscible LBM using CUDA
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
	unsigned int nBlocks, nthreads;
	int Nx,Ny,Nz;		// local sub-domain size
	int nspheres;		// number of spheres in the packing
	double Lx,Ly,Lz;	// Domain length
	double D = 1.0;		// reference length for non-dimensionalization
	// Color Model parameters
	int timestepMax;
	double tau,Fx,Fy,Fz,tol,err;
	double alpha, beta;
	double das, dbs, phi_s;
	double din,dout;
	double wp_saturation;
	int BoundaryCondition;
	int InitialCondition;
//	bool pBC,Restart;
	int i,j,k;

	// pmmc threshold values
	//double fluid_isovalue,solid_isovalue;
	//fluid_isovalue = 0.0;
	//solid_isovalue = 0.0;
	
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
			input >> tau;			// Viscosity parameter
			input >> alpha;			// Surface Tension parameter
			input >> beta;			// Width of the interface
			input >> phi_s;			// value of phi at the solid surface
			// Line 2: wetting phase saturation to initialize
			input >> wp_saturation;
			// Line 3: External force components (Fx,Fy, Fz)
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
			tau = 1.0;
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
	MPI_Bcast(&tau,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&alpha,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&beta,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&das,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&dbs,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&phi_s,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&wp_saturation,1,MPI_DOUBLE,0,comm);
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

    // Get the rank info
    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);
	MPI_Barrier(comm);
	
	// **************************************************************
	// **************************************************************
	double Ps = -(das-dbs)/(das+dbs);
	double rlxA = 1.f/tau;
	double rlxB = 8.f*(2.f-rlxA)/(8.f-rlxA);
	//double xIntPos = log((1.0+phi_s)/(1.0-phi_s))/(2.0*beta); 	
	
	// Set the density values inside the solid based on the input value phi_s
 	das = (phi_s+1.0)*0.5;
	dbs = 1.0 - das;
	
	if (nprocs != nprocx*nprocy*nprocz){
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
		INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
	}

	if (rank==0){
		printf("********************************************************\n");
		printf("tau = %f \n", tau);
		printf("alpha = %f \n", alpha);		
		printf("beta = %f \n", beta);
		printf("das = %f \n", das);
		printf("dbs = %f \n", dbs);
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
		if (InitialCondition==0) printf("Initial conditions assigned from phase ID file \n");
		if (InitialCondition==1) printf("Initial conditions assigned from restart file \n");
		printf("********************************************************\n");
	}

	// Initialized domain and averaging framework for Two-Phase Flow
	bool pBC,velBC;
	if (BoundaryCondition==1 || BoundaryCondition==3)
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
	Dm.CommInit(comm);

	// Mask that excludes the solid phase
	Domain Mask(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);

	MPI_Barrier(comm);

	Nx+=2; Ny+=2; Nz += 2;
	//Nx = Ny = Nz;	// Cubic domain

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
	int sum = 0;
	double sum_local;
	double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
	if (BoundaryCondition > 0) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
	double porosity, pore_vol;
	//...........................................................................
	if (rank == 0) cout << "Reading in domain from signed distance function..." << endl;

	//.......................................................................
	sprintf(LocalRankString,"%05d",rank);
//	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
//	WriteLocalSolidID(LocalRankFilename, id, N);
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
	sum=0;
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

	/*	for ( k=0;k<Nz;k++){
		for ( j=0;j<Ny;j++){
			for ( i=0;i<Nx;i++){
				int n = k*Nx*Ny+j*Nx+i;
				// The following turns off communication if external BC are being set
				if (BoundaryCondition > 0){
					if (kproc==0 && k==0)			id[n]=0;
					if (kproc==0 && k==1)			id[n]=0;
					if (kproc==nprocz-1 && k==Nz-2)	id[n]=0;
					if (kproc==nprocz-1 && k==Nz-1)	id[n]=0;
				}
			}
		}
	}
	*/
	// Set up kstart, kfinish so that the reservoirs are excluded from averaging
	int kstart,kfinish;
	kstart = 1;
	kfinish = Nz-1;
	if (BoundaryCondition >  0 && Dm.kproc==0)		kstart = 4;
	if (BoundaryCondition >  0 && Dm.kproc==nprocz-1)	kfinish = Nz-4;

	// Compute the pore volume
	sum_local = 0.0;
	for ( k=kstart;k<kfinish;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				int n = k*Nx*Ny+j*Nx+i;
				if (id[n] > 0){
					sum_local += 1.0;
				}
			}
		}
	}
	MPI_Allreduce(&sum_local,&pore_vol,1,MPI_DOUBLE,MPI_SUM,comm);
//	MPI_Allreduce(&sum_local,&porosity,1,MPI_DOUBLE,MPI_SUM,comm);
	porosity = pore_vol*iVol_global;
	if (rank==0) printf("Media porosity = %f \n",porosity);
	//.........................................................
	// If external boundary conditions are applied remove solid
	if (BoundaryCondition >  0  && Dm.kproc == 0){
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
	if (BoundaryCondition >  0  && Dm.kproc == nprocz-1){
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

	//...........................................................................
	if (rank==0)	printf ("Create ScaLBL_Communicator \n");
	// Create a communicator for the device
	ScaLBL_Communicator ScaLBL_Comm(Mask);

	// set reservoirs
	if (BoundaryCondition > 0){
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
	}

	//...........device phase ID.................................................
	if (rank==0)	printf ("Copy phase ID to device \n");
	char *ID;
	ScaLBL_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
	// Don't compute in the halo
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
	//...........................................................................

	//...........................................................................
	//				MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	// LBM variables
	if (rank==0)	printf ("Allocating distributions \n");
	//......................device distributions.................................
	double *f_even,*f_odd;
	double *A_even,*A_odd,*B_even,*B_odd;
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &A_even, 4*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &A_odd, 3*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &B_even, 4*dist_mem_size);	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &B_odd, 3*dist_mem_size);	// Allocate device memory
	//...........................................................................
	double *Phi,*Den;
	double *ColorGrad, *Velocity, *Pressure, *dvcSignDist;
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &Phi, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Pressure, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &dvcSignDist, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*dist_mem_size);
	//...........................................................................

	// Copy signed distance for device initialization
	ScaLBL_CopyToDevice(dvcSignDist, Averages->SDs.data(), dist_mem_size);
	//...........................................................................

	int logcount = 0; // number of surface write-outs
	
	//...........................................................................
	//				MAIN  VARIABLES INITIALIZED HERE
	//...........................................................................
	//...........................................................................
	//...........................................................................
	if (rank==0)	printf("Setting the distributions, size = %i\n", N);
	//...........................................................................
	ScaLBL_DeviceBarrier();
	ScaLBL_D3Q19_Init(ID, f_even, f_odd, Nx, Ny, Nz);
	ScaLBL_Color_Init(ID, Den, Phi, das, dbs, Nx, Ny, Nz);
	ScaLBL_DeviceBarrier();
	//......................................................................

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
		ScaLBL_CopyToDevice(Den,cDen,2*N*sizeof(double));
		ScaLBL_DeviceBarrier();
	    delete [] cDen;
	    delete [] cDistEven;
	    delete [] cDistOdd;
		MPI_Barrier(comm);
	}

	//......................................................................
	ScaLBL_D3Q7_Init(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
	ScaLBL_D3Q7_Init(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
	ScaLBL_DeviceBarrier();
	MPI_Barrier(comm);
	//.......................................................................
	// Once phase has been initialized, map solid to account for 'smeared' interface
	//for (i=0; i<N; i++)	Averages->SDs(i) -= (1.0);
	// Make sure the id match for the two domains
	for (i=0; i<N; i++)	Dm.id[i] = Mask.id[i];
	//.......................................................................
	// Finalize setup for averaging domain
	//Averages->SetupCubes(Mask);
	Averages->UpdateSolid();
	//.......................................................................
	
	//*************************************************************************
	// 		Compute the phase indicator field and reset Copy, Den
	//*************************************************************************
	ScaLBL_ComputePhaseField(ID, Phi, Den, N);
	//*************************************************************************
	ScaLBL_DeviceBarrier();
	ScaLBL_Comm.SendHalo(Phi);
	ScaLBL_Comm.RecvHalo(Phi);
	ScaLBL_DeviceBarrier();
	MPI_Barrier(comm);
	//*************************************************************************

	if (rank==0 && BoundaryCondition==1){
		printf("Setting inlet pressure = %f \n", din);
		printf("Setting outlet pressure = %f \n", dout);
	}
	if (BoundaryCondition==1 && Mask.kproc == 0)	{
		ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
		ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
	}
		
	if (BoundaryCondition==1 && Mask.kproc == nprocz-1){
		ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
		ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
	}

	if (rank==0 && BoundaryCondition==2){
		printf("Setting inlet velocity = %f \n", din);
		printf("Setting outlet velocity = %f \n", dout);
	}
	if (BoundaryCondition==2 && Mask.kproc == 0)	{
		ScaLBL_D3Q19_Velocity_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
		//ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,0);
	}

	if (BoundaryCondition==2 && Mask.kproc == nprocz-1){
		ScaLBL_D3Q19_Velocity_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
		//ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);	
		ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,Nz-1);
	}
	// Set dynamic pressure boundary conditions
	double dp, slope;
	dp = slope = 0.0;
	if (BoundaryCondition==3){
		slope = abs(dout-din)/timestepMax;
		dp = abs(timestep)*slope;
		if (rank==0) printf("Change in pressure / time =%.3e \n",slope);
		// set the initial value
		din  = 1.0+0.5*dp;
		dout = 1.0-0.5*dp;
		// set the initial boundary conditions
		if (Mask.kproc == 0)	{
			ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
			ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		}
		if (Mask.kproc == nprocz-1){
			ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
			ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		}
	}

	ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
	ScaLBL_D3Q19_Velocity(ID,f_even,f_odd,Velocity,Nx,Ny,Nz);

	//...........................................................................
	// Copy the phase indicator field for the earlier timestep
	ScaLBL_DeviceBarrier();
	ScaLBL_CopyToHost(Averages->Phase_tplus.data(),Phi,N*sizeof(double));
	//...........................................................................
	//...........................................................................
	// Copy the data for for the analysis timestep
	//...........................................................................
	// Copy the phase from the GPU -> CPU
	//...........................................................................
	ScaLBL_DeviceBarrier();
	ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
	ScaLBL_CopyToHost(Averages->Phase.data(),Phi,N*sizeof(double));
	ScaLBL_CopyToHost(Averages->Press.data(),Pressure,N*sizeof(double));
	ScaLBL_CopyToHost(Averages->Vel_x.data(),&Velocity[0],N*sizeof(double));
	ScaLBL_CopyToHost(Averages->Vel_y.data(),&Velocity[N],N*sizeof(double));
	ScaLBL_CopyToHost(Averages->Vel_z.data(),&Velocity[2*N],N*sizeof(double));
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

		//*************************************************************************
		// Fused Color Gradient and Collision 
		//*************************************************************************
		ScaLBL_D3Q19_ColorCollide( ID,f_even,f_odd,Phi,ColorGrad,
							 Velocity,Nx,Ny,Nz,rlxA,rlxB,alpha,beta,Fx,Fy,Fz);
		//*************************************************************************

		ScaLBL_DeviceBarrier();
		//*************************************************************************
		// Pack and send the D3Q19 distributions
		ScaLBL_Comm.SendD3Q19(f_even, f_odd);
		//*************************************************************************

		//*************************************************************************
		// 		Carry out the density streaming step for mass transport
		//*************************************************************************
		ScaLBL_D3Q7_ColorCollideMass(ID, A_even, A_odd, B_even, B_odd, Den, Phi,
								ColorGrad, Velocity, beta, N, pBC);
		//*************************************************************************

		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		//*************************************************************************
		// 		Swap the distributions for momentum transport
		//*************************************************************************
		ScaLBL_D3Q19_Swap(ID, f_even, f_odd, Nx, Ny, Nz);
		//*************************************************************************

		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
		//*************************************************************************
		// Wait for communications to complete and unpack the distributions
		ScaLBL_Comm.RecvD3Q19(f_even, f_odd);
		//*************************************************************************

		ScaLBL_DeviceBarrier();
		//*************************************************************************
		// Pack and send the D3Q7 distributions
		ScaLBL_Comm.BiSendD3Q7(A_even, A_odd, B_even, B_odd);
		//*************************************************************************

		ScaLBL_DeviceBarrier();
		ScaLBL_D3Q7_Swap(ID, A_even, A_odd, Nx, Ny, Nz);
		ScaLBL_D3Q7_Swap(ID, B_even, B_odd, Nx, Ny, Nz);

		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);

		//*************************************************************************
		// Wait for communication and unpack the D3Q7 distributions
		ScaLBL_Comm.BiRecvD3Q7(A_even, A_odd, B_even, B_odd);
		//*************************************************************************

		ScaLBL_DeviceBarrier();
		//..................................................................................
		ScaLBL_D3Q7_Density(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
		ScaLBL_D3Q7_Density(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
		//*************************************************************************
		// 		Compute the phase indicator field 
		//*************************************************************************
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);

		ScaLBL_ComputePhaseField(ID, Phi, Den, N);
		//*************************************************************************
		ScaLBL_Comm.SendHalo(Phi);
		ScaLBL_DeviceBarrier();
		ScaLBL_Comm.RecvHalo(Phi);
		//*************************************************************************

		ScaLBL_DeviceBarrier();
		
		// Pressure boundary conditions
		if (BoundaryCondition==1 && Mask.kproc == 0)	{
			ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
			ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		}
		if (BoundaryCondition==1 && Mask.kproc == nprocz-1){
			ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
			ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		}

		// Velocity boundary conditions
		if (BoundaryCondition==2 && Mask.kproc == 0)	{
			ScaLBL_D3Q19_Velocity_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
			ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
			ScaLBL_SetSlice_z(Phi,-1.0,Nx,Ny,Nz,0);
		}
		if (BoundaryCondition==2 && Mask.kproc == nprocz-1){
			ScaLBL_D3Q19_Velocity_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
			ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
			ScaLBL_SetSlice_z(Phi,1.0,Nx,Ny,Nz,Nz-1);
		}

		if (BoundaryCondition==3){
			// Increase the pressure difference
			dp += slope;
			din  = 1.0+0.5*dp;
			dout = 1.0-0.5*dp;
			// set the initial boundary conditions
			if (Mask.kproc == 0)	{
				ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
				ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
			}
			if (Mask.kproc == nprocz-1){
				ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
				ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
			}
		}

		//...................................................................................

		MPI_Barrier(comm);
        PROFILE_STOP("Update");

		// Timestep completed!
		timestep++;
  
        // Run the analysis, blob identification, and write restart files
	//	if (BLOB_ANALYSIS_INTERVAL > 0){
	//if (timestep > 5)
	  run_analysis(timestep,RESTART_INTERVAL,rank_info,*Averages,last_ids,last_index,last_id_map,
            Nx,Ny,Nz,pBC,beta,err,Phi,Pressure,Velocity,ID,f_even,f_odd,Den,
            LocalRestartFile,meshData,fillData,tpool,work_ids);
	/*		}
		else{
		  ComputeMacroscaleAverages(timestep,ANALYSIS_INTERVAL,RESTART_INTERVAL,rank_info,*Averages,
            Nx,Ny,Nz,pBC,beta,err,Phi,Pressure,Velocity,ID,f_even,f_odd,Den,
				  LocalRestartFile,meshData,fillData,tpool,work_ids);
		}
	*/
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
/*	// Perform component averaging and write tcat averages
	Averages->Initialize();
	Averages->ComponentAverages();
	Averages->SortBlobs();
	Averages->PrintComponents(timestep);
	// ************************************************************************

	int NumberComponents_NWP = ComputeGlobalPhaseComponent(Mask.Nx-2,Mask.Ny-2,Mask.Nz-2,Mask.rank_info,Averages->PhaseID,1,Averages->Label_NWP);
	printf("Number of non-wetting phase components: %i \n ",NumberComponents_NWP);
	ScaLBL_DeviceBarrier();
	ScaLBL_CopyToHost(Averages->Phase.data(),Phi,N*sizeof(double));
*/
    
/*	Averages->WriteSurfaces(0);

	sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
	FILE *PHASE;
	PHASE = fopen(LocalRankFilename,"wb");
	fwrite(Averages->SDn.data(),8,N,PHASE);
	fclose(PHASE);
	*/

	/*	sprintf(LocalRankFilename,"%s%s","Pressure.",LocalRankString);
	FILE *PRESS;
	PRESS = fopen(LocalRankFilename,"wb");
	fwrite(Averages->Press.data(),8,N,PRESS);
	fclose(PRESS);

	ScaLBL_CopyToHost(Averages->Phase.data(),Phi,N*sizeof(double));
	double * Grad;
	Grad = new double [3*N];
	ScaLBL_CopyToHost(Grad,ColorGrad,3*N*sizeof(double));
	sprintf(LocalRankFilename,"%s%s","ColorGrad.",LocalRankString);
	FILE *GRAD;
	GRAD = fopen(LocalRankFilename,"wb");
	fwrite(Grad,8,3*N,GRAD);
	fclose(GRAD);
	*/
    PROFILE_STOP("Main");
    PROFILE_SAVE("lbpm_color_simulator",1);
	// ****************************************************
	MPI_Barrier(comm);
  } // Limit scope so variables that contain communicators will free before MPI_Finialize
  MPI_Comm_free(&comm);
  MPI_Finalize();
}


