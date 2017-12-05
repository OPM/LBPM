#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "pmmc.h"
#include "Domain.h"
#include "Extras.h"
#include "D3Q19.h"
#include "D3Q7.h"
#include "Color.h"
#include "common/MPI_Helpers.h"
#include "Communication.h"

#define WRITE_SURFACES

#define MAX_LOCAL_BLOB_COUNT=100

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
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
	int sendtag,recvtag;
	//*****************************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************
	MPI_Request req1[18],req2[18];
	MPI_Status stat1[18],stat2[18];

	if (rank == 0){
		printf("********************************************************\n");
		printf("Running Color LBM	\n");
		printf("********************************************************\n");
	}

	// Variables that specify the computational domain  
	string FILENAME;
	unsigned int nBlocks, nthreads;
	int Nx,Ny,Nz;		// local sub-domain size
	int nspheres;		// number of spheres in the packing
	double Lx,Ly,Lz;	// Domain length
	double D = 1.0;		// reference length for non-dimensionalization
	// Color Model parameters
	int timestepMax, interval;
	double tau,Fx,Fy,Fz,tol;
	double alpha, beta;
	double das, dbs, phi_s;
	double din,dout;
	double wp_saturation;
	bool pBC,Restart;
	int i,j,k,p,n;

	// pmmc threshold values
	double fluid_isovalue,solid_isovalue;
	fluid_isovalue = 0.0;
	solid_isovalue = 0.0;
	nBlocks = 32;
	nthreads = 128;
	
	int RESTART_INTERVAL=20000;
	
	if (rank==0){
		//.............................................................
		//		READ SIMULATION PARMAETERS FROM INPUT FILE
		//.............................................................
		ifstream input("Color.in");
		// Line 1: Name of the phase indicator file (s=0,w=1,n=2)
//		input >> FILENAME;
		// Line 2: domain size (Nx, Ny, Nz)
//		input >> Nz;				// number of nodes (x,y,z)
//		input >> nBlocks;
//		input >> nthreads;
		// Line 3: model parameters (tau, alpha, beta, das, dbs)
		input >> tau;			// Viscosity parameter
		input >> alpha;			// Surface Tension parameter
		input >> beta;			// Width of the interface
		input >> phi_s;			// value of phi at the solid surface
//		input >> das;
//		input >> dbs;
		// Line 4: wetting phase saturation to initialize
		input >> wp_saturation;
		// Line 5: External force components (Fx,Fy, Fz)
		input >> Fx;
		input >> Fy;
		input >> Fz;
		// Line 6: Pressure Boundary conditions
		input >> Restart;
		input >> pBC;
		input >> din;
		input >> dout;
		// Line 7: time-stepping criteria
		input >> timestepMax;		// max no. of timesteps
		input >> interval;			// restart interval
		input >> tol;				// error tolerance
		//.............................................................
		das = 0.1; dbs = 0.9;	// hard coded for density initialization
								// should be OK to remove these parameters
								// they should have no impact with the 
								// current boundary condition
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
	MPI_Bcast(&pBC,1,MPI_LOGICAL,0,comm);
	MPI_Bcast(&Restart,1,MPI_LOGICAL,0,comm);
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
//	MPI_Bcast(&nBlocks,1,MPI_INT,0,comm);
//	MPI_Bcast(&nthreads,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
	//.................................................
	MPI_Barrier(comm);
	
	RESTART_INTERVAL=interval;
	// **************************************************************
	// **************************************************************
	double Ps = -(das-dbs)/(das+dbs);
	double rlxA = 1.f/tau;
	double rlxB = 8.f*(2.f-rlxA)/(8.f-rlxA);
	double xIntPos;
	xIntPos = log((1.0+phi_s)/(1.0-phi_s))/(2.0*beta); 	
	
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
//		printf("das = %f \n", das);
//		printf("dbs = %f \n", dbs);
		printf("Value of phi at solid surface = %f \n", phi_s);
		printf("Distance to phi = 0.0: %f \n", xIntPos);
		printf("gamma_{wn} = %f \n", 5.796*alpha);
//		printf("cos theta_c = %f \n", 1.05332*Ps);
		printf("Force(x) = %f \n", Fx);
		printf("Force(y) = %f \n", Fy);
		printf("Force(z) = %f \n", Fz);
		printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
		printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
		printf("********************************************************\n");
	}

	 InitializeRanks( rank, nprocx, nprocy, nprocz, iproc, jproc, kproc, 
			 	 	 rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z,
			 	 	 rank_xy, rank_XY, rank_xY, rank_Xy, rank_xz, rank_XZ, rank_xZ, rank_Xz,
			 	 	 rank_yz, rank_YZ, rank_yZ, rank_Yz );
	 
	 MPI_Barrier(comm);

	Nz += 2;
	Nx = Ny = Nz;	// Cubic domain

	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

//	unsigned int nBlocks = 32;
//	int nthreads = 128;
	int S = N/nthreads/nBlocks+1;

//	unsigned int nBlocks = N/nthreads + (N%nthreads == 0?0:1);
//	dim3 grid(nBlocks,1,1);

	if (rank==0) printf("Number of blocks = %i \n", nBlocks);
	if (rank==0) printf("Threads per block = %i \n", nthreads);
	if (rank==0) printf("Sweeps per thread = %i \n", S);
	if (rank==0) printf("Number of nodes per side = %i \n", Nx);
	if (rank==0) printf("Total Number of nodes = %i \n", N);
	if (rank==0) printf("********************************************************\n");

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
	double iVol_global = 1.0/(1.0*Nx*Ny*Nz*nprocs);
	double porosity;
/*	//.......................................................................
	ifstream PM(LocalRankFilename,ios::binary);
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				id[n] = 0;
			}
		}
	}
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				PM.read((char *) (&value), sizeof(value));
				n = k*Nx*Ny+j*Nx+i;
				id[n] = value;
				if (value > 0) sum++;
			}
		}
	}
	PM.close();
//	printf("File porosity = %f\n", double(sum)/N);
*/	//...........................................................................
//	double *SignDist;
//	SignDist = new double[N];
	DoubleArray SignDist(Nx,Ny,Nz);
	//.......................................................................
#ifdef CBUB
	// Initializes a constrained bubble test
	double BubbleBot = 20.0;  // How big to make the NWP bubble
	double BubbleTop = 60.0;  // How big to make the NWP bubble
	double TubeRadius = 15.5; // Radius of the capillary tube
	sum=0;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				// Cylindrical capillary tube aligned with the z direction
				SignDist(i,j,k) = TubeRadius-sqrt(1.0*((i-Nx/2)*(i-Nx/2) 
									+ (j-Ny/2)*(j-Ny/2)));
				// Initialize phase positions field
				if (SignDist(i,j,k) < 0.0){
					id[n] = 0;
				}
				else if (k<BubbleBot){
					id[n] = 2;
					sum++;
				}
				else if (k<BubbleTop && rank == 0 && pBC == 0){
					id[n] = 1;
					sum++;
				}
				else{
					id[n] = 2;
					sum++;
				}
			}
		}
	}
	porosity = double(sum)/double(1.0*N);
#else
	// Read in sphere pack
	if (rank==1) printf("nspheres =%i \n",nspheres);
	//.......................................................................
	double *cx,*cy,*cz,*rad;
	cx = new double[nspheres];
	cy = new double[nspheres];
	cz = new double[nspheres];
	rad = new double[nspheres];
	//.......................................................................
	if (rank == 0)	printf("Reading the sphere packing \n");
	if (rank == 0)	ReadSpherePacking(nspheres,cx,cy,cz,rad);
	MPI_Barrier(comm);
	// Broadcast the sphere packing to all processes
	MPI_Bcast(cx,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(cy,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(cz,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(rad,nspheres,MPI_DOUBLE,0,comm);
	//...........................................................................
	MPI_Barrier(comm);
	if (rank == 0) cout << "Domain set." << endl;
	if (rank == 0){
		// Compute the Sauter mean diameter
		double totVol = 0.0;
		double totArea = 0.0;
		// Compute the total volume and area of all spheres
		for (i=0; i<nspheres; i++){
			totVol += 1.3333333333333*3.14159265359*rad[i]*rad[i]*rad[i];
			totArea += 4.0*3.14159265359*rad[i]*rad[i];
		}
		D = 6.0*(Nx-2)*nprocx*totVol / totArea / Lx;
		printf("Sauter Mean Diameter (computed from sphere packing) = %f \n ",D);
	}
	MPI_Bcast(&D,1,MPI_DOUBLE,0,comm);

	//.......................................................................
//	sprintf(LocalRankString,"%05d",rank);
//	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	//.......................................................................
	SignedDistance(SignDist.data,nspheres,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,
					   iproc,jproc,kproc,nprocx,nprocy,nprocz);

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
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (SignDist.data[n] > 0.0){ 
					id[n] = 2;	
					sum++;
				}
			}
		}
	}
	//.........................................................
	// If pressure boundary conditions are applied remove solid
	if (pBC && kproc == 0){
		for (k=0; k<3; k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 1;
				}					
			}
		}
	}
	if (pBC && kproc == nprocz-1){
		for (k=Nz-3; k<Nz; k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					id[n] = 2;
				}					
			}
		}
	}
	//.........................................................
	// don't perform computations at the eight corners
	id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
	id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
	//.........................................................
	sum_local = 1.0*sum;
	MPI_Allreduce(&sum_local,&porosity,1,MPI_DOUBLE,MPI_SUM,comm);
	porosity = porosity*iVol_global;
	if (rank==0) printf("Media porosity = %f \n",porosity);

	// Generate the residual NWP 
	if (rank==0) printf("Initializing with NWP saturation = %f \n",wp_saturation);
	if (!pBC)	GenerateResidual(id,Nx,Ny,Nz,wp_saturation);
#endif

	// Set up MPI communication structurese
	if (rank==0)	printf ("Setting up communication control structures \n");
	//......................................................................................
	// Get the actual D3Q19 communication counts (based on location of solid phase)
	// Discrete velocity set symmetry implies the sendcount = recvcount
	int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
	int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
	int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
	sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
	sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
	sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
	//......................................................................................
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				// Check the phase ID
				if (id[k*Nx*Ny+j*Nx+i] != 0){
					// Counts for the six faces
					if (i==1)	sendCount_x++;
					if (j==1)	sendCount_y++;
					if (k==1)	sendCount_z++;
					if (i==Nx-2)	sendCount_X++;
					if (j==Ny-2)	sendCount_Y++;
					if (k==Nz-2)	sendCount_Z++;
					// Counts for the twelve edges
					if (i==1 && j==1)	sendCount_xy++;
					if (i==1 && j==Ny-2)	sendCount_xY++;
					if (i==Nx-2 && j==1)	sendCount_Xy++;
					if (i==Nx-2 && j==Ny-2)	sendCount_XY++;

					if (i==1 && k==1)	sendCount_xz++;
					if (i==1 && k==Nz-2)	sendCount_xZ++;
					if (i==Nx-2 && k==1)	sendCount_Xz++;
					if (i==Nx-2 && k==Nz-2)	sendCount_XZ++;

					if (j==1 && k==1)	sendCount_yz++;
					if (j==1 && k==Nz-2)	sendCount_yZ++;
					if (j==Ny-2 && k==1)	sendCount_Yz++;
					if (j==Ny-2 && k==Nz-2)	sendCount_YZ++;
				}
			}
		}
	}
	//......................................................................................
	int *sendList_x, *sendList_y, *sendList_z, *sendList_X, *sendList_Y, *sendList_Z;
	int *sendList_xy, *sendList_yz, *sendList_xz, *sendList_Xy, *sendList_Yz, *sendList_xZ;
	int *sendList_xY, *sendList_yZ, *sendList_Xz, *sendList_XY, *sendList_YZ, *sendList_XZ;
	//......................................................................................
	// send buffers
	sendList_x = new int [sendCount_x];
	sendList_y = new int [sendCount_y];
	sendList_z = new int [sendCount_z];
	sendList_X = new int [sendCount_X];
	sendList_Y = new int [sendCount_Y];
	sendList_Z = new int [sendCount_Z];
	sendList_xy = new int [sendCount_xy];
	sendList_yz = new int [sendCount_yz];
	sendList_xz = new int [sendCount_xz];
	sendList_Xy = new int [sendCount_Xy];
	sendList_Yz = new int [sendCount_Yz];
	sendList_xZ = new int [sendCount_xZ];
	sendList_xY = new int [sendCount_xY];
	sendList_yZ = new int [sendCount_yZ];
	sendList_Xz = new int [sendCount_Xz];
	sendList_XY = new int [sendCount_XY];
	sendList_YZ = new int [sendCount_YZ];
	sendList_XZ = new int [sendCount_XZ];
	if (rank==0)	printf ("Preparing the sendlists \n");
	//......................................................................................
	// Populate the send list
	sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
	sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
	sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				// Local value to send
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] != 0){
					// Counts for the six faces
					if (i==1)		sendList_x[sendCount_x++]=n;
					if (j==1)		sendList_y[sendCount_y++]=n;
					if (k==1)		sendList_z[sendCount_z++]=n;
					if (i==Nx-2)	sendList_X[sendCount_X++]=n;
					if (j==Ny-2)	sendList_Y[sendCount_Y++]=n;
					if (k==Nz-2)	sendList_Z[sendCount_Z++]=n;
					// Counts for the twelve edges
					if (i==1 && j==1)		sendList_xy[sendCount_xy++]=n;
					if (i==1 && j==Ny-2)	sendList_xY[sendCount_xY++]=n;
					if (i==Nx-2 && j==1)	sendList_Xy[sendCount_Xy++]=n;
					if (i==Nx-2 && j==Ny-2)	sendList_XY[sendCount_XY++]=n;

					if (i==1 && k==1)		sendList_xz[sendCount_xz++]=n;
					if (i==1 && k==Nz-2)	sendList_xZ[sendCount_xZ++]=n;
					if (i==Nx-2 && k==1)	sendList_Xz[sendCount_Xz++]=n;
					if (i==Nx-2 && k==Nz-2)	sendList_XZ[sendCount_XZ++]=n;

					if (j==1 && k==1)		sendList_yz[sendCount_yz++]=n;
					if (j==1 && k==Nz-2)	sendList_yZ[sendCount_yZ++]=n;
					if (j==Ny-2 && k==1)	sendList_Yz[sendCount_Yz++]=n;
					if (j==Ny-2 && k==Nz-2)	sendList_YZ[sendCount_YZ++]=n;
				}
			}
		}
	}
	MPI_Barrier(comm);
	if (rank==0)	printf ("SendLists are ready on host\n");
	//......................................................................................
	// Use MPI to fill in the recvCounts form the associated processes
	int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
	int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
	int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
	//......................................................................................
	//**********************************************************************************
	// Fill in the recieve counts using MPI
	sendtag = recvtag = 3;
	CommunicateSendRecvCounts( comm, sendtag, recvtag, 
			rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z,
			rank_xy, rank_XY, rank_xY, rank_Xy,
			rank_xz, rank_XZ, rank_xZ, rank_Xz,
			rank_yz, rank_YZ, rank_yZ, rank_Yz,
			sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z,
			sendCount_xy, sendCount_XY, sendCount_xY, sendCount_Xy,
			sendCount_xz, sendCount_XZ, sendCount_xZ, sendCount_Xz,
			sendCount_yz, sendCount_YZ, sendCount_yZ, sendCount_Yz,
			recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z,
			recvCount_xy, recvCount_XY, recvCount_xY, recvCount_Xy,
			recvCount_xz, recvCount_XZ, recvCount_xZ, recvCount_Xz,
			recvCount_yz, recvCount_YZ, recvCount_yZ, recvCount_Yz );
	//**********************************************************************************
	//......................................................................................
	int *recvList_x, *recvList_y, *recvList_z, *recvList_X, *recvList_Y, *recvList_Z;
	int *recvList_xy, *recvList_yz, *recvList_xz, *recvList_Xy, *recvList_Yz, *recvList_xZ;
	int *recvList_xY, *recvList_yZ, *recvList_Xz, *recvList_XY, *recvList_YZ, *recvList_XZ;
	//......................................................................................
	// recv buffers
	recvList_x = new int [recvCount_x];
	recvList_y = new int [recvCount_y];
	recvList_z = new int [recvCount_z];
	recvList_X = new int [recvCount_X];
	recvList_Y = new int [recvCount_Y];
	recvList_Z = new int [recvCount_Z];
	recvList_xy = new int [recvCount_xy];
	recvList_yz = new int [recvCount_yz];
	recvList_xz = new int [recvCount_xz];
	recvList_Xy = new int [recvCount_Xy];
	recvList_Yz = new int [recvCount_Yz];
	recvList_xZ = new int [recvCount_xZ];
	recvList_xY = new int [recvCount_xY];
	recvList_yZ = new int [recvCount_yZ];
	recvList_Xz = new int [recvCount_Xz];
	recvList_XY = new int [recvCount_XY];
	recvList_YZ = new int [recvCount_YZ];
	recvList_XZ = new int [recvCount_XZ];
	//......................................................................................
	//......................................................................................
	// Use MPI to fill in the appropriate values for recvList
	// Fill in the recieve lists using MPI
	sendtag = recvtag = 4;
	CommunicateRecvLists( comm, sendtag, recvtag, 
			sendList_x, sendList_y, sendList_z, sendList_X, sendList_Y, sendList_Z,
			sendList_xy, sendList_XY, sendList_xY, sendList_Xy,
			sendList_xz, sendList_XZ, sendList_xZ, sendList_Xz,
			sendList_yz, sendList_YZ, sendList_yZ, sendList_Yz,
			sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z,
			sendCount_xy, sendCount_XY, sendCount_xY, sendCount_Xy,
			sendCount_xz, sendCount_XZ, sendCount_xZ, sendCount_Xz,
			sendCount_yz, sendCount_YZ, sendCount_yZ, sendCount_Yz,
			recvList_x, recvList_y, recvList_z, recvList_X, recvList_Y, recvList_Z,
			recvList_xy, recvList_XY, recvList_xY, recvList_Xy,
			recvList_xz, recvList_XZ, recvList_xZ, recvList_Xz,
			recvList_yz, recvList_YZ, recvList_yZ, recvList_Yz,
			recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z,
			recvCount_xy, recvCount_XY, recvCount_xY, recvCount_Xy,
			recvCount_xz, recvCount_XZ, recvCount_xZ, recvCount_Xz,
			recvCount_yz, recvCount_YZ, recvCount_yZ, recvCount_Yz,
			rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z, rank_xy, rank_XY, rank_xY,
			rank_Xy, rank_xz, rank_XZ, rank_xZ, rank_Xz, rank_yz, rank_YZ, rank_yZ, rank_Yz );
	//......................................................................................
	for (int idx=0; idx<recvCount_x; idx++)	recvList_x[idx] -= (Nx-2);
	for (int idx=0; idx<recvCount_X; idx++)	recvList_X[idx] += (Nx-2);
	for (int idx=0; idx<recvCount_y; idx++)	recvList_y[idx] -= (Ny-2)*Nx;
	for (int idx=0; idx<recvCount_Y; idx++)	recvList_Y[idx] += (Ny-2)*Nx;
	for (int idx=0; idx<recvCount_z; idx++)	recvList_z[idx] -= (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Z; idx++)	recvList_Z[idx] += (Nz-2)*Nx*Ny;
	
	for (int idx=0; idx<recvCount_xy; idx++)	recvList_xy[idx] -= (Nx-2)+(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_XY; idx++)	recvList_XY[idx] += (Nx-2)+(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_xY; idx++)	recvList_xY[idx] -= (Nx-2)-(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_Xy; idx++)	recvList_Xy[idx] += (Nx-2)-(Ny-2)*Nx;
	
	for (int idx=0; idx<recvCount_xz; idx++)	recvList_xz[idx] -= (Nx-2)+(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_XZ; idx++)	recvList_XZ[idx] += (Nx-2)+(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_xZ; idx++)	recvList_xZ[idx] -= (Nx-2)-(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Xz; idx++)	recvList_Xz[idx] += (Nx-2)-(Nz-2)*Nx*Ny;
	
	for (int idx=0; idx<recvCount_yz; idx++)	recvList_yz[idx] -= (Ny-2)*Nx + (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_YZ; idx++)	recvList_YZ[idx] += (Ny-2)*Nx + (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_yZ; idx++)	recvList_yZ[idx] -= (Ny-2)*Nx - (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Yz; idx++)	recvList_Yz[idx] += (Ny-2)*Nx - (Nz-2)*Nx*Ny;
	//......................................................................................
	double *sendbuf_x, *sendbuf_y, *sendbuf_z, *sendbuf_X, *sendbuf_Y, *sendbuf_Z;
	double *sendbuf_xy, *sendbuf_yz, *sendbuf_xz, *sendbuf_Xy, *sendbuf_Yz, *sendbuf_xZ;
	double *sendbuf_xY, *sendbuf_yZ, *sendbuf_Xz, *sendbuf_XY, *sendbuf_YZ, *sendbuf_XZ;
	double *recvbuf_x, *recvbuf_y, *recvbuf_z, *recvbuf_X, *recvbuf_Y, *recvbuf_Z;
	double *recvbuf_xy, *recvbuf_yz, *recvbuf_xz, *recvbuf_Xy, *recvbuf_Yz, *recvbuf_xZ;
	double *recvbuf_xY, *recvbuf_yZ, *recvbuf_Xz, *recvbuf_XY, *recvbuf_YZ, *recvbuf_XZ;
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_x, 5*sendCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_X, 5*sendCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_y, 5*sendCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Y, 5*sendCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_z, 5*sendCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Z, 5*sendCount_Z*sizeof(double));	// Allocatevoid * memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xy, sendCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xY, sendCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Xy, sendCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_XY, sendCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xz, sendCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_xZ, sendCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Xz, sendCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_XZ, sendCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_yz, sendCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_yZ, sendCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_Yz, sendCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &sendbuf_YZ, sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_x, 5*recvCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_X, 5*recvCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_y, 5*recvCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Y, 5*recvCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_z, 5*recvCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Z, 5*recvCount_Z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xy, recvCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xY, recvCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Xy, recvCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_XY, recvCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xz, recvCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_xZ, recvCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Xz, recvCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_XZ, recvCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_yz, recvCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_yZ, recvCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_Yz, recvCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &recvbuf_YZ, recvCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X, *dvcSendList_Y, *dvcSendList_Z;
	int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy, *dvcSendList_Yz, *dvcSendList_xZ;
	int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY, *dvcSendList_YZ, *dvcSendList_XZ;
	//......................................................................................
	int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X, *dvcRecvList_Y, *dvcRecvList_Z;
	int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy, *dvcRecvList_Yz, *dvcRecvList_xZ;
	int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY, *dvcRecvList_YZ, *dvcRecvList_XZ;
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_x, sendCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_X, sendCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_y, sendCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Y, sendCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_z, sendCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Z, sendCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xy, sendCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xY, sendCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Xy, sendCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_XY, sendCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xz, sendCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_xZ, sendCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Xz, sendCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_XZ, sendCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_yz, sendCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_yZ, sendCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_Yz, sendCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcSendList_YZ, sendCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_x, recvCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_X, recvCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_y, recvCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Y, recvCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_z, recvCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Z, recvCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateDeviceMemory((void **) &dvcRecvList_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	MPI_Barrier(comm);
	if (rank==0)	printf ("Prepare to copy send/recv Lists to device \n");
	ScaLBL_CopyToDevice(dvcSendList_x,sendList_x,sendCount_x*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_X,sendList_X,sendCount_X*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_y,sendList_y,sendCount_y*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Y,sendList_Y,sendCount_Y*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_z,sendList_z,sendCount_z*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Z,sendList_Z,sendCount_Z*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xy,sendList_xy,sendCount_xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_XY,sendList_XY,sendCount_XY*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xY,sendList_xY,sendCount_xY*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Xy,sendList_Xy,sendCount_Xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xz,sendList_xz,sendCount_xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_XZ,sendList_XZ,sendCount_XZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_xZ,sendList_xZ,sendCount_xZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Xz,sendList_Xz,sendCount_Xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_yz,sendList_yz,sendCount_yz*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_YZ,sendList_YZ,sendCount_YZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_yZ,sendList_yZ,sendCount_yZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcSendList_Yz,sendList_Yz,sendCount_Yz*sizeof(int));
	//......................................................................................
	ScaLBL_CopyToDevice(dvcRecvList_x,recvList_x,recvCount_x*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_X,recvList_X,recvCount_X*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_y,recvList_y,recvCount_y*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Y,recvList_Y,recvCount_Y*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_z,recvList_z,recvCount_z*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Z,recvList_Z,recvCount_Z*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xy,recvList_xy,recvCount_xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_XY,recvList_XY,recvCount_XY*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xY,recvList_xY,recvCount_xY*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Xy,recvList_Xy,recvCount_Xy*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xz,recvList_xz,recvCount_xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_XZ,recvList_XZ,recvCount_XZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_xZ,recvList_xZ,recvCount_xZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Xz,recvList_Xz,recvCount_Xz*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_yz,recvList_yz,recvCount_yz*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_YZ,recvList_YZ,recvCount_YZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_yZ,recvList_yZ,recvCount_yZ*sizeof(int));
	ScaLBL_CopyToDevice(dvcRecvList_Yz,recvList_Yz,recvCount_Yz*sizeof(int));
	//......................................................................................
	// Fill in the phase ID from neighboring processors
	char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
	// send buffers
	sendID_x = new char [sendCount_x];
	sendID_y = new char [sendCount_y];
	sendID_z = new char [sendCount_z];
	sendID_X = new char [sendCount_X];
	sendID_Y = new char [sendCount_Y];
	sendID_Z = new char [sendCount_Z];
	sendID_xy = new char [sendCount_xy];
	sendID_yz = new char [sendCount_yz];
	sendID_xz = new char [sendCount_xz];
	sendID_Xy = new char [sendCount_Xy];
	sendID_Yz = new char [sendCount_Yz];
	sendID_xZ = new char [sendCount_xZ];
	sendID_xY = new char [sendCount_xY];
	sendID_yZ = new char [sendCount_yZ];
	sendID_Xz = new char [sendCount_Xz];
	sendID_XY = new char [sendCount_XY];
	sendID_YZ = new char [sendCount_YZ];
	sendID_XZ = new char [sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvID_x = new char [recvCount_x];
	recvID_y = new char [recvCount_y];
	recvID_z = new char [recvCount_z];
	recvID_X = new char [recvCount_X];
	recvID_Y = new char [recvCount_Y];
	recvID_Z = new char [recvCount_Z];
	recvID_xy = new char [recvCount_xy];
	recvID_yz = new char [recvCount_yz];
	recvID_xz = new char [recvCount_xz];
	recvID_Xy = new char [recvCount_Xy];
	recvID_xZ = new char [recvCount_xZ];
	recvID_xY = new char [recvCount_xY];
	recvID_yZ = new char [recvCount_yZ];
	recvID_Yz = new char [recvCount_Yz];
	recvID_Xz = new char [recvCount_Xz];
	recvID_XY = new char [recvCount_XY];
	recvID_YZ = new char [recvCount_YZ];
	recvID_XZ = new char [recvCount_XZ];
	//......................................................................................
	sendtag = recvtag = 7;
	PackID(sendList_x, sendCount_x ,sendID_x, id);
	PackID(sendList_X, sendCount_X ,sendID_X, id);
	PackID(sendList_y, sendCount_y ,sendID_y, id);
	PackID(sendList_Y, sendCount_Y ,sendID_Y, id);
	PackID(sendList_z, sendCount_z ,sendID_z, id);
	PackID(sendList_Z, sendCount_Z ,sendID_Z, id);
	PackID(sendList_xy, sendCount_xy ,sendID_xy, id);
	PackID(sendList_Xy, sendCount_Xy ,sendID_Xy, id);
	PackID(sendList_xY, sendCount_xY ,sendID_xY, id);
	PackID(sendList_XY, sendCount_XY ,sendID_XY, id);
	PackID(sendList_xz, sendCount_xz ,sendID_xz, id);
	PackID(sendList_Xz, sendCount_Xz ,sendID_Xz, id);
	PackID(sendList_xZ, sendCount_xZ ,sendID_xZ, id);
	PackID(sendList_XZ, sendCount_XZ ,sendID_XZ, id);
	PackID(sendList_yz, sendCount_yz ,sendID_yz, id);
	PackID(sendList_Yz, sendCount_Yz ,sendID_Yz, id);
	PackID(sendList_yZ, sendCount_yZ ,sendID_yZ, id);
	PackID(sendList_YZ, sendCount_YZ ,sendID_YZ, id);
	//......................................................................................
	MPI_Sendrecv(sendID_x,sendCount_x,MPI_CHAR,rank_x,sendtag,
			recvID_X,recvCount_X,MPI_CHAR,rank_X,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_X,sendCount_X,MPI_CHAR,rank_X,sendtag,
			recvID_x,recvCount_x,MPI_CHAR,rank_x,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_y,sendCount_y,MPI_CHAR,rank_y,sendtag,
			recvID_Y,recvCount_Y,MPI_CHAR,rank_Y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Y,sendCount_Y,MPI_CHAR,rank_Y,sendtag,
			recvID_y,recvCount_y,MPI_CHAR,rank_y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_z,sendCount_z,MPI_CHAR,rank_z,sendtag,
			recvID_Z,recvCount_Z,MPI_CHAR,rank_Z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Z,sendCount_Z,MPI_CHAR,rank_Z,sendtag,
			recvID_z,recvCount_z,MPI_CHAR,rank_z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xy,sendCount_xy,MPI_CHAR,rank_xy,sendtag,
			recvID_XY,recvCount_XY,MPI_CHAR,rank_XY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XY,sendCount_XY,MPI_CHAR,rank_XY,sendtag,
			recvID_xy,recvCount_xy,MPI_CHAR,rank_xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xy,sendCount_Xy,MPI_CHAR,rank_Xy,sendtag,
			recvID_xY,recvCount_xY,MPI_CHAR,rank_xY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xY,sendCount_xY,MPI_CHAR,rank_xY,sendtag,
			recvID_Xy,recvCount_Xy,MPI_CHAR,rank_Xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xz,sendCount_xz,MPI_CHAR,rank_xz,sendtag,
			recvID_XZ,recvCount_XZ,MPI_CHAR,rank_XZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XZ,sendCount_XZ,MPI_CHAR,rank_XZ,sendtag,
			recvID_xz,recvCount_xz,MPI_CHAR,rank_xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xz,sendCount_Xz,MPI_CHAR,rank_Xz,sendtag,
			recvID_xZ,recvCount_xZ,MPI_CHAR,rank_xZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xZ,sendCount_xZ,MPI_CHAR,rank_xZ,sendtag,
			recvID_Xz,recvCount_Xz,MPI_CHAR,rank_Xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yz,sendCount_yz,MPI_CHAR,rank_yz,sendtag,
			recvID_YZ,recvCount_YZ,MPI_CHAR,rank_YZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_YZ,sendCount_YZ,MPI_CHAR,rank_YZ,sendtag,
			recvID_yz,recvCount_yz,MPI_CHAR,rank_yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Yz,sendCount_Yz,MPI_CHAR,rank_Yz,sendtag,
			recvID_yZ,recvCount_yZ,MPI_CHAR,rank_yZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yZ,sendCount_yZ,MPI_CHAR,rank_yZ,sendtag,
			recvID_Yz,recvCount_Yz,MPI_CHAR,rank_Yz,recvtag,comm,MPI_STATUS_IGNORE);
	//......................................................................................
	UnpackID(recvList_x, recvCount_x ,recvID_x, id);
	UnpackID(recvList_X, recvCount_X ,recvID_X, id);
	UnpackID(recvList_y, recvCount_y ,recvID_y, id);
	UnpackID(recvList_Y, recvCount_Y ,recvID_Y, id);
	UnpackID(recvList_z, recvCount_z ,recvID_z, id);
	UnpackID(recvList_Z, recvCount_Z ,recvID_Z, id);
	UnpackID(recvList_xy, recvCount_xy ,recvID_xy, id);
	UnpackID(recvList_Xy, recvCount_Xy ,recvID_Xy, id);
	UnpackID(recvList_xY, recvCount_xY ,recvID_xY, id);
	UnpackID(recvList_XY, recvCount_XY ,recvID_XY, id);
	UnpackID(recvList_xz, recvCount_xz ,recvID_xz, id);
	UnpackID(recvList_Xz, recvCount_Xz ,recvID_Xz, id);
	UnpackID(recvList_xZ, recvCount_xZ ,recvID_xZ, id);
	UnpackID(recvList_XZ, recvCount_XZ ,recvID_XZ, id);
	UnpackID(recvList_yz, recvCount_yz ,recvID_yz, id);
	UnpackID(recvList_Yz, recvCount_Yz ,recvID_Yz, id);
	UnpackID(recvList_yZ, recvCount_yZ ,recvID_yZ, id);
	UnpackID(recvList_YZ, recvCount_YZ ,recvID_YZ, id);
	//......................................................................................
	// Fill in the phase MeshData from neighboring processors
	double *sendMeshData_x, *sendMeshData_y, *sendMeshData_z, *sendMeshData_X, *sendMeshData_Y, *sendMeshData_Z;
	double *sendMeshData_xy, *sendMeshData_yz, *sendMeshData_xz, *sendMeshData_Xy, *sendMeshData_Yz, *sendMeshData_xZ;
	double *sendMeshData_xY, *sendMeshData_yZ, *sendMeshData_Xz, *sendMeshData_XY, *sendMeshData_YZ, *sendMeshData_XZ;
	double *recvMeshData_x, *recvMeshData_y, *recvMeshData_z, *recvMeshData_X, *recvMeshData_Y, *recvMeshData_Z;
	double *recvMeshData_xy, *recvMeshData_yz, *recvMeshData_xz, *recvMeshData_Xy, *recvMeshData_Yz, *recvMeshData_xZ;
	double *recvMeshData_xY, *recvMeshData_yZ, *recvMeshData_Xz, *recvMeshData_XY, *recvMeshData_YZ, *recvMeshData_XZ;
	// send buffers
	sendMeshData_x = new double [sendCount_x];
	sendMeshData_y = new double [sendCount_y];
	sendMeshData_z = new double [sendCount_z];
	sendMeshData_X = new double [sendCount_X];
	sendMeshData_Y = new double [sendCount_Y];
	sendMeshData_Z = new double [sendCount_Z];
	sendMeshData_xy = new double [sendCount_xy];
	sendMeshData_yz = new double [sendCount_yz];
	sendMeshData_xz = new double [sendCount_xz];
	sendMeshData_Xy = new double [sendCount_Xy];
	sendMeshData_Yz = new double [sendCount_Yz];
	sendMeshData_xZ = new double [sendCount_xZ];
	sendMeshData_xY = new double [sendCount_xY];
	sendMeshData_yZ = new double [sendCount_yZ];
	sendMeshData_Xz = new double [sendCount_Xz];
	sendMeshData_XY = new double [sendCount_XY];
	sendMeshData_YZ = new double [sendCount_YZ];
	sendMeshData_XZ = new double [sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvMeshData_x = new double [recvCount_x];
	recvMeshData_y = new double [recvCount_y];
	recvMeshData_z = new double [recvCount_z];
	recvMeshData_X = new double [recvCount_X];
	recvMeshData_Y = new double [recvCount_Y];
	recvMeshData_Z = new double [recvCount_Z];
	recvMeshData_xy = new double [recvCount_xy];
	recvMeshData_yz = new double [recvCount_yz];
	recvMeshData_xz = new double [recvCount_xz];
	recvMeshData_Xy = new double [recvCount_Xy];
	recvMeshData_xZ = new double [recvCount_xZ];
	recvMeshData_xY = new double [recvCount_xY];
	recvMeshData_yZ = new double [recvCount_yZ];
	recvMeshData_Yz = new double [recvCount_Yz];
	recvMeshData_Xz = new double [recvCount_Xz];
	recvMeshData_XY = new double [recvCount_XY];
	recvMeshData_YZ = new double [recvCount_YZ];
	recvMeshData_XZ = new double [recvCount_XZ];
	if (rank==0)	printf ("Devices are ready to communicate. \n");
	MPI_Barrier(comm);

	//...........device phase ID.................................................
	if (rank==0)	printf ("Copying phase ID to device \n");
	char *ID;
	ScaLBL_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
	// Copy to the device
	ScaLBL_CopyToDevice(ID, id, N);
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
//	double *Copy;
	double *ColorGrad, *Velocity, *Pressure, *dvcSignDist;
	//...........................................................................
	ScaLBL_AllocateDeviceMemory((void **) &Phi, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Pressure, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &dvcSignDist, dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*dist_mem_size);
	ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*dist_mem_size);
	// Copy signed distance for device initialization
	ScaLBL_CopyToDevice(dvcSignDist, SignDist.data, dist_mem_size);
	//...........................................................................
	// Phase indicator (in array form as needed by PMMC algorithm)
	DoubleArray Phase(Nx,Ny,Nz);

	// Extra copies of phi needed to compute time derivatives on CPU
	DoubleArray Phase_tminus(Nx,Ny,Nz);
	DoubleArray Phase_tplus(Nx,Ny,Nz);
	DoubleArray dPdt(Nx,Ny,Nz);

	//copies of data needed to perform checkpointing from cpu
	double *cDen, *cDistEven, *cDistOdd;
	cDen = new double[2*N];
	cDistEven = new double[10*N];
	cDistOdd = new double[9*N];

	// data needed to perform CPU-based averaging
//	double *Vel;
//	Vel = new double[3*N];		// fluid velocity
//	Press = new double[N];		// fluid pressure

	DoubleArray Press(Nx,Ny,Nz);
	DoubleArray MeanCurvature(Nx,Ny,Nz);
	DoubleArray GaussCurvature(Nx,Ny,Nz);
	DoubleArray SignDist_x(Nx,Ny,Nz);		// Gradient of the signed distance
	DoubleArray SignDist_y(Nx,Ny,Nz);
	DoubleArray SignDist_z(Nx,Ny,Nz);
	DoubleArray Phase_x(Nx,Ny,Nz);			// Gradient of the phase indicator field
	DoubleArray Phase_y(Nx,Ny,Nz);
	DoubleArray Phase_z(Nx,Ny,Nz);
	DoubleArray Vel_x(Nx,Ny,Nz);			// Gradient of the phase indicator field
	DoubleArray Vel_y(Nx,Ny,Nz);
	DoubleArray Vel_z(Nx,Ny,Nz);
	/*****************************************************************
	 VARIABLES FOR THE PMMC ALGORITHM
	 ****************************************************************** */
	//...........................................................................
	// Averaging variables
	//...........................................................................
	// local averages (to each MPI process)
	double trimdist=1.0; 						// pixel distance to trim surface for specified averages
	double awn,ans,aws,lwns,nwp_volume;
	double As;
	double vol_w, vol_n;						// volumes the exclude the interfacial region
	double sat_w;
	double pan,paw;								// local phase averaged pressure
//	double vx_w,vy_w,vz_w,vx_n,vy_n,vz_n;  		// local phase averaged velocity
	// Global averages (all processes)
	double vol_w_global, vol_n_global;			// volumes the exclude the interfacial region
	double awn_global,ans_global,aws_global;
	double lwns_global;
	double efawns,efawns_global;				// averaged contact angle
	double Jwn,Jwn_global;						// average mean curavture - wn interface
	double Kwn,Kwn_global;						// average Gaussian curavture - wn interface
	double trawn,trawn_global;					// trimmed interfacial area
	double trJwn,trJwn_global;					// trimmed interfacial area	
	DoubleArray van(3);
	DoubleArray vaw(3);
	DoubleArray vawn(3);
	DoubleArray Gwn(6);
	DoubleArray Gns(6);
	DoubleArray Gws(6);
	
	// Local blob averages
	double blob_awn,blob_ans,blob_aws,blob_lwns,blob_nwp_volume;
	double blob_As;
	double blob_vol_w, blob_vol_n;						// volumes the exclude the interfacial region
	double blob_pan,blob_paw;								// local phase averaged pressure
	double blob_efawns;				// averaged contact angle
	double blob_Jwn;						// average mean curavture - wn interface
	double blob_Kwn;						// average Gaussian curavture - wn interface
	double blob_trawn;					// trimmed interfacial area
	double blob_trJwn;					// trimmed interfacial area	
	DoubleArray blob_van(3);
	DoubleArray blob_vaw(3);
	DoubleArray blob_vawn(3);
	DoubleArray blob_Gwn(6);
	DoubleArray blob_Gns(6);
	DoubleArray blob_Gws(6);
	
	
	double nwp_volume_global;					// volume for the wetting phase (for saturation)
//	double p_n_global,p_w_global;				// global phase averaged pressure
//	double vx_w_global,vy_w_global,vz_w_global;	// global phase averaged velocity
//	double vx_n_global,vy_n_global,vz_n_global;	// global phase averaged velocity
	double As_global;
	double dEs,dAwn,dAns;						// Global surface energy (calculated by rank=0)
	double pan_global,paw_global;								// local phase averaged pressure
	DoubleArray van_global(3);
	DoubleArray vaw_global(3);
	DoubleArray vawn_global(3);
	DoubleArray Gwn_global(6);
	DoubleArray Gns_global(6);
	DoubleArray Gws_global(6);
	//...........................................................................

	//	bool add=1;			// Set to false if any corners contain nw-phase ( F > fluid_isovalue)
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	DoubleArray CubeValues(2,2,2);
//	int count_in=0,count_out=0;
//	int nodx,nody,nodz;
	// initialize lists for vertices for surfaces, common line
	DTMutableList<Point> nw_pts(20);
	DTMutableList<Point> ns_pts(20);
	DTMutableList<Point> ws_pts(20);
	DTMutableList<Point> nws_pts(20);
	// initialize triangle lists for surfaces
	IntArray nw_tris(3,20);
	IntArray ns_tris(3,20);
	IntArray ws_tris(3,20);
	// initialize list for line segments
	IntArray nws_seg(2,20);
	DTMutableList<Point> tmp(20);
	DoubleArray Values(20);
	DoubleArray ContactAngle(20);
	DoubleArray Curvature(20);
	DoubleArray DistValues(20);
	DoubleArray InterfaceSpeed(20);
	DoubleArray NormalVector(60);
	
	//	IntArray store;
	
	int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0;
	int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
	
//	double s,s1,s2,s3;		// Triangle sides (lengths)
	Point A,B,C,P;
//	double area;
	
	// Initialize arrays for local solid surface
	DTMutableList<Point> local_sol_pts(20);
	int n_local_sol_pts = 0;
	IntArray local_sol_tris(3,18);
	int n_local_sol_tris;
	DoubleArray values(20);
	DTMutableList<Point> local_nws_pts(20);
	int n_local_nws_pts;
	
	//int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg;
	int c;
	//int newton_steps = 0;
	//...........................................................................
	int ncubes = (Nx-2)*(Ny-2)*(Nz-2);	// Exclude the "upper" halo
//	IntArray cubeList(3,ncubes);
	IntArray LocalBlobID(Nx,Ny,Nz);
	IntArray LocalBlobCubeList(3,ncubes);				// store indices for blobs (cubes)
	IntArray temp(3,ncubes);							// temporary storage array
	IntArray Blobs(MAX_LOCAL_BLOB_COUNT);				// number of nodes in each blob
	IntArray LocalBlobTable(20,MAX_LOCAL_BLOB_COUNT); 			// hold the neighboring ranks
	DoubleArray LocalBlobAverages(30,MAX_LOCAL_BLOB_COUNT); 	// hold the local blob averages (tcat)
	int nc=0;
	//...........................................................................

	//...........................................................................
	// Grids used to pack faces on the GPU for MPI
	int faceGrid,edgeGrid,packThreads;
	packThreads=512;
	edgeGrid=1;
	faceGrid=Nx*Ny/packThreads;
	//...........................................................................

	int logcount = 0; // number of surface write-outs
	
	//...........................................................................
	//				MAIN  VARIABLES INITIALIZED HERE
	//...........................................................................
	//...........................................................................
	if (rank==0)	printf("Setting the distributions, size = %i\n", N);
	//...........................................................................
	ScaLBL_D3Q19_Init(ID, f_even, f_odd, Nx, Ny, Nz);
	//......................................................................
//	ScaLBL_Color_InitDistance(ID, Copy, Phi, SignDist.data, das, dbs, beta, xIntPos, Nx, Ny, Nz, S);
	ScaLBL_Color_InitDistance(ID, Den, Phi, dvcSignDist, das, dbs, beta, xIntPos, Nx, Ny, Nz);
	//......................................................................
	//.......................................................................
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	WriteLocalSolidID(LocalRankFilename, id, N);
	sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
	WriteLocalSolidDistance(LocalRankFilename, SignDist.data, N);
	//.......................................................................
	if (Restart == true){
		if (rank==0) printf("Reading restart file! \n");
		// Read in the restart file to CPU buffers
		ReadCheckpoint(LocalRestartFile, cDen, cDistEven, cDistOdd, N);
		// Copy the restart data to the GPU
		ScaLBL_CopyToDevice(f_even,cDistEven,10*N*sizeof(double));
		ScaLBL_CopyToDevice(f_odd,cDistOdd,9*N*sizeof(double));
		ScaLBL_CopyToDevice(Den,cDen,2*N*sizeof(double));
		ScaLBL_DeviceBarrier();
		MPI_Barrier(comm);
	}
	// Set up the cube list (very regular in this case due to lack of blob-ID)
	// Set up kstart, kfinish so that the reservoirs are excluded from averaging
	int kstart,kfinish;
	kstart = 1;
	kfinish = Nz-1;
	nc = 0;
	for (k=0;k<Nz-1;k++){
		for (j=0;j<Ny-1;j++){
			for (i=0;i<Nx-1;i++){
				if ( LocalBlobID(i,j,k) == -1 ){
					if ( Phase(i,j,k) > 0.0 ){
						if ( SignDist(i,j,k) > 0.0 ){
							// node i,j,k is in the porespace
							Blobs(nblobs) = ComputeBlob(LocalBlobCubeList,nblobs,nc,LocalBlobID,Phase,SignDist,0.0,0.0,i,j,k,temp);
							nblobs++;
							INSIST(nblobs < MAX_LOCAL_BLOB_COUNT,"Not enough to store the local blobs!");

						}
					}
				}
			}
		}
	}
	int count_in=0,count_out=0;
	int nodx,nody,nodz;
	for (k=0;k<Nz-1;k++){
		for (j=0;j<Ny-1;j++){
			for (i=0;i<Nx-1;i++){
				// Loop over cube corners
				add=1;				// initialize to true - add unless corner occupied by nw-phase
				for (p=0;p<8;p++){
					nodx=i+cube[p][0];
					nody=j+cube[p][1];
					nodz=k+cube[p][2];
					if ( indicator(nodx,nody,nodz) > -1 ){
						// corner occupied by nw-phase  -> do not add
						add = 0;
					}
				}
				if ( add == 1 ){
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					count_in++;
				}
				else { count_out++; }
			}
		}
	}
	b(nblobs) = count_in;
	nblobs++;
	ScaLBL_D3Q7_Init(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
	ScaLBL_D3Q7_Init(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
	// Once phase has been initialized, map solid to account for 'smeared' interface
	//......................................................................
	for (i=0; i<N; i++)	SignDist.data[i] -= (1.0); // 
	
	//...........................................................................
	// Gradient of the Signed Distance function
	//...........................................................................
	pmmc_MeshGradient(SignDist,SignDist_x,SignDist_y,SignDist_z,Nx,Ny,Nz);
	//...........................................................................
	CommunicateMeshHalo(SignDist_x, comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
			sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
			sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
			sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
			sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
			sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
			recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
			recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
			recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
			recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
			recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
			recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
			rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
			rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
	//...........................................................................
	CommunicateMeshHalo(SignDist_y, comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
			sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
			sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
			sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
			sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
			sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
			recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
			recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
			recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
			recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
			recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
			recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
			rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
			rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
	//...........................................................................
	CommunicateMeshHalo(SignDist_z, comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
			sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
			sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
			sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
			sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
			sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
			recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
			recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
			recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
			recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
			recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
			recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
			rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
			rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
	//...........................................................................
	//......................................................................

	//*************************************************************************
	// 		Compute the phase indicator field and reset Copy, Den
	//*************************************************************************
//	ScaLBL_ComputePhaseField(ID, Phi, Copy, Den, N, S);
	ScaLBL_ComputePhaseField(ID, Phi, Den, N);
	//*************************************************************************
	//...................................................................................
	ScaLBL_Scalar_Pack(dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_xy, sendCount_xy,sendbuf_xy, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_xY, sendCount_xY,sendbuf_xY, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_Xy, sendCount_Xy,sendbuf_Xy, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_XY, sendCount_XY,sendbuf_XY, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_xz, sendCount_xz,sendbuf_xz, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_xZ, sendCount_xZ,sendbuf_xZ, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_Xz, sendCount_Xz,sendbuf_Xz, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_XZ, sendCount_XZ,sendbuf_XZ, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_yz, sendCount_yz,sendbuf_yz, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_yZ, sendCount_yZ,sendbuf_yZ, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_Yz, sendCount_Yz,sendbuf_Yz, Phi, N);
	ScaLBL_Scalar_Pack(dvcSendList_YZ, sendCount_YZ,sendbuf_YZ, Phi, N);
	
	ScaLBL_DeviceBarrier();
	
	//...................................................................................
	// Send / Recv all the phase indcator field values
	//...................................................................................
	MPI_Isend(sendbuf_x, sendCount_x,MPI_DOUBLE,rank_x,sendtag,comm,&req1[0]);
	MPI_Irecv(recvbuf_X, recvCount_X,MPI_DOUBLE,rank_X,recvtag,comm,&req2[0]);
	MPI_Isend(sendbuf_X, sendCount_X,MPI_DOUBLE,rank_X,sendtag,comm,&req1[1]);
	MPI_Irecv(recvbuf_x, recvCount_x,MPI_DOUBLE,rank_x,recvtag,comm,&req2[1]);
	MPI_Isend(sendbuf_y, sendCount_y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[2]);
	MPI_Irecv(recvbuf_Y, recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[2]);
	MPI_Isend(sendbuf_Y, sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[3]);
	MPI_Irecv(recvbuf_y, recvCount_y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[3]);
	MPI_Isend(sendbuf_z, sendCount_z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[4]);
	MPI_Irecv(recvbuf_Z, recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[4]);
	MPI_Isend(sendbuf_Z, sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[5]);
	MPI_Irecv(recvbuf_z, recvCount_z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[5]);
	MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,comm,&req1[6]);
	MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,comm,&req2[6]);
	MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,comm,&req1[7]);
	MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,comm,&req2[7]);
	MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,comm,&req1[8]);
	MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,comm,&req2[8]);
	MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,comm,&req1[9]);
	MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,comm,&req2[9]);
	MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,comm,&req1[10]);
	MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,comm,&req2[10]);
	MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,comm,&req1[11]);
	MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,comm,&req2[11]);
	MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,comm,&req1[12]);
	MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,comm,&req2[12]);
	MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,comm,&req1[13]);
	MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,comm,&req2[13]);
	MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,comm,&req1[14]);
	MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,comm,&req2[14]);
	MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,comm,&req1[15]);
	MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,comm,&req2[15]);
	MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,comm,&req1[16]);
	MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,comm,&req2[16]);
	MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,comm,&req1[17]);
	MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,comm,&req2[17]);
	//...................................................................................
	//...................................................................................
	// Wait for completion of Indicator Field communication
	//...................................................................................
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	ScaLBL_DeviceBarrier();
	//...................................................................................
	//...................................................................................
/*		ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
	ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
	ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
	ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
	ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
	ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
*/		
	ScaLBL_Scalar_Unpack(dvcRecvList_x, recvCount_x,recvbuf_x, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_y, recvCount_y,recvbuf_y, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_z, recvCount_z,recvbuf_z, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_X, recvCount_X,recvbuf_X, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Y, recvCount_Y,recvbuf_Y, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Z, recvCount_Z,recvbuf_Z, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xy, recvCount_xy,recvbuf_xy, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xY, recvCount_xY,recvbuf_xY, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XY, recvCount_XY,recvbuf_XY, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xz, recvCount_xz,recvbuf_xz, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yz, recvCount_yz,recvbuf_yz, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, Phi, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, Phi, N);
	//...................................................................................

	if (rank==0 && pBC){
		printf("Setting inlet pressure = %f \n", din);
		printf("Setting outlet pressure = %f \n", dout);
	}
	if (pBC && kproc == 0)	{
		ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
		ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
	}
		
	if (pBC && kproc == nprocz-1){
		ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
		ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
	}

	//...........................................................................
	// Copy the phase indicator field for the earlier timestep
	ScaLBL_DeviceBarrier();
	ScaLBL_CopyToHost(Phase_tplus.data,Phi,N*sizeof(double));
	//...........................................................................
	//...........................................................................
	// Copy the data for for the analysis timestep
	//...........................................................................
	// Copy the phase from the GPU -> CPU
	//...........................................................................
	ScaLBL_DeviceBarrier();
	ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
	ScaLBL_CopyToHost(Phase.data,Phi,N*sizeof(double));
	ScaLBL_CopyToHost(Press.data,Pressure,N*sizeof(double));
	ScaLBL_CopyToHost(Vel_x.data,&Velocity[0],N*sizeof(double));
	ScaLBL_CopyToHost(Vel_y.data,&Velocity[N],N*sizeof(double));
	ScaLBL_CopyToHost(Vel_z.data,&Velocity[2*N],N*sizeof(double));
	MPI_Barrier(comm);
	//...........................................................................
	
	int timestep = 0;
	if (rank==0) printf("********************************************************\n");
	if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);

	//.......create a stream for the LB calculation.......
//	cudaStream_t stream;
//	cudaStreamCreate(&stream);
	
	//.......create and start timer............
	double starttime,stoptime,cputime;
	MPI_Barrier(comm);
	starttime = MPI_Wtime();
	//.........................................
	
	sendtag = recvtag = 5;
	FILE *TIMELOG;
	if (rank==0){
		TIMELOG= fopen("timelog.tcat","a+");
		if (fseek(TIMELOG,0,SEEK_SET) == fseek(TIMELOG,0,SEEK_CUR)){
			// If timelog is empty, write a short header to list the averages
			fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
			fprintf(TIMELOG,"timestep dEs ");								// Timestep, Change in Surface Energy
			fprintf(TIMELOG,"sw pw pn awn ans aws Jwn Kwn lwns efawns ");	// Scalar averages
			fprintf(TIMELOG,"vw[x, y, z] vn[x, y, z] vwn[x, y, z]");		// Velocity averages
			fprintf(TIMELOG,"Gwn [xx, yy, zz, xy, xz, yz] ");				// Orientation tensors
			fprintf(TIMELOG,"Gws [xx, yy, zz, xy, xz, yz] ");
			fprintf(TIMELOG,"Gns [xx, yy, zz, xy, xz, yz] ");
			fprintf(TIMELOG,"trJwn trawn \n");								// trimmed curvature for wn surface
			fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
		}
	}

	//************ MAIN ITERATION LOOP ***************************************/
	while (timestep < timestepMax){

		//*************************************************************************
		// Fused Color Gradient and Collision 
		//*************************************************************************
		ScaLBL_D3Q19_ColorCollide( ID,f_even,f_odd,Phi,ColorGrad,
							 Velocity,Nx,Ny,Nz,rlxA,rlxB,alpha,beta,Fx,Fy,Fz);
		//*************************************************************************

		//...................................................................................
		ScaLBL_D3Q19_Pack(1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
		ScaLBL_D3Q19_Pack(4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		ScaLBL_D3Q19_Pack(5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		ScaLBL_D3Q19_Pack(6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		ScaLBL_D3Q19_Pack(7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		//...Packing for X face(1,7,9,11,13)................................
		ScaLBL_D3Q19_Pack(0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
		ScaLBL_D3Q19_Pack(3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		ScaLBL_D3Q19_Pack(4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		ScaLBL_D3Q19_Pack(5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		ScaLBL_D3Q19_Pack(6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		//...Packing for y face(4,8,9,16,18).................................
		ScaLBL_D3Q19_Pack(2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
		ScaLBL_D3Q19_Pack(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		ScaLBL_D3Q19_Pack(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
		ScaLBL_D3Q19_Pack(8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		ScaLBL_D3Q19_Pack(9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		//...Packing for Y face(3,7,10,15,17).................................
		ScaLBL_D3Q19_Pack(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
		ScaLBL_D3Q19_Pack(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		ScaLBL_D3Q19_Pack(5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
		ScaLBL_D3Q19_Pack(7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		ScaLBL_D3Q19_Pack(8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		//...Packing for z face(6,12,13,16,17)................................
		ScaLBL_D3Q19_Pack(3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
		ScaLBL_D3Q19_Pack(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		ScaLBL_D3Q19_Pack(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		ScaLBL_D3Q19_Pack(8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		ScaLBL_D3Q19_Pack(8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		//...Packing for Z face(5,11,14,15,18)................................
		ScaLBL_D3Q19_Pack(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
		ScaLBL_D3Q19_Pack(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		ScaLBL_D3Q19_Pack(7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		ScaLBL_D3Q19_Pack(7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		ScaLBL_D3Q19_Pack(9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		//...Pack the xy edge (8)................................
		ScaLBL_D3Q19_Pack(4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
		//...Pack the Xy edge (9)................................
		ScaLBL_D3Q19_Pack(4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
		//...Pack the xY edge (10)................................
		ScaLBL_D3Q19_Pack(5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
		//...Pack the XY edge (7)................................
		ScaLBL_D3Q19_Pack(3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
		//...Pack the xz edge (12)................................
		ScaLBL_D3Q19_Pack(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
		//...Pack the xZ edge (14)................................
		ScaLBL_D3Q19_Pack(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
		//...Pack the Xz edge (13)................................
		ScaLBL_D3Q19_Pack(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
		//...Pack the XZ edge (11)................................
		ScaLBL_D3Q19_Pack(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
		//...Pack the xz edge (12)................................
		//...Pack the yz edge (16)................................
		ScaLBL_D3Q19_Pack(8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
		//...Pack the yZ edge (18)................................
		ScaLBL_D3Q19_Pack(9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
		//...Pack the Yz edge (17)................................
		ScaLBL_D3Q19_Pack(8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
		//...Pack the YZ edge (15)................................
		ScaLBL_D3Q19_Pack(7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
		//...................................................................................

		//...................................................................................
		// Send all the distributions
		MPI_Isend(sendbuf_x, 5*sendCount_x,MPI_DOUBLE,rank_x,sendtag,comm,&req1[0]);
		MPI_Irecv(recvbuf_X, 5*recvCount_X,MPI_DOUBLE,rank_X,recvtag,comm,&req2[0]);
		MPI_Isend(sendbuf_X, 5*sendCount_X,MPI_DOUBLE,rank_X,sendtag,comm,&req1[1]);
		MPI_Irecv(recvbuf_x, 5*recvCount_x,MPI_DOUBLE,rank_x,recvtag,comm,&req2[1]);
		MPI_Isend(sendbuf_y, 5*sendCount_y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[2]);
		MPI_Irecv(recvbuf_Y, 5*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[2]);
		MPI_Isend(sendbuf_Y, 5*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[3]);
		MPI_Irecv(recvbuf_y, 5*recvCount_y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[3]);
		MPI_Isend(sendbuf_z, 5*sendCount_z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[4]);
		MPI_Irecv(recvbuf_Z, 5*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[4]);
		MPI_Isend(sendbuf_Z, 5*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[5]);
		MPI_Irecv(recvbuf_z, 5*recvCount_z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[5]);
		MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,comm,&req1[6]);
		MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,comm,&req2[6]);
		MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,comm,&req1[7]);
		MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,comm,&req2[7]);
		MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,comm,&req1[8]);
		MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,comm,&req2[8]);
		MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,comm,&req1[9]);
		MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,comm,&req2[9]);
		MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,comm,&req1[10]);
		MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,comm,&req2[10]);
		MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,comm,&req1[11]);
		MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,comm,&req2[11]);
		MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,comm,&req1[12]);
		MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,comm,&req2[12]);
		MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,comm,&req1[13]);
		MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,comm,&req2[13]);
		MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,comm,&req1[14]);
		MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,comm,&req2[14]);
		MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,comm,&req1[15]);
		MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,comm,&req2[15]);
		MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,comm,&req1[16]);
		MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,comm,&req2[16]);
		MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,comm,&req1[17]);
		MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,comm,&req2[17]);
		//...................................................................................

		//*************************************************************************
		// 		Carry out the density streaming step for mass transport
		//*************************************************************************
//		DensityStreamD3Q7(ID, Den, Copy, Phi, ColorGrad, Velocity, beta, Nx, Ny, Nz, pBC, S);
		//*************************************************************************
		ScaLBL_D3Q7_ColorCollideMass(ID, A_even, A_odd, B_even, B_odd, Den, Phi,
								ColorGrad, Velocity, beta, N, pBC);
			

		//*************************************************************************
		// 		Swap the distributions for momentum transport
		//*************************************************************************
		ScaLBL_D3Q19_Swap(ID, f_even, f_odd, Nx, Ny, Nz);
		//*************************************************************************

		//...................................................................................
		// Wait for completion of D3Q19 communication
		MPI_Waitall(18,req1,stat1);
		MPI_Waitall(18,req2,stat2);

		//...................................................................................
		// Unpack the distributions on the device
		//...................................................................................
		//...Map recieve list for the X face: q=2,8,10,12,13 .................................
		ScaLBL_D3Q19_Unpack(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(3,-1,-1,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(4,-1,1,0,dvcRecvList_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(5,-1,0,-1,dvcRecvList_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(6,-1,0,1,dvcRecvList_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the x face: q=1,7,9,11,13..................................
		ScaLBL_D3Q19_Unpack(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(4,1,1,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(5,1,-1,0,dvcRecvList_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(6,1,0,1,dvcRecvList_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(7,1,0,-1,dvcRecvList_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
		ScaLBL_D3Q19_Unpack(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(3,-1,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(5,1,-1,0,dvcRecvList_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(7,0,-1,-1,dvcRecvList_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(8,0,-1,1,dvcRecvList_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
		ScaLBL_D3Q19_Unpack(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(4,1,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(4,-1,1,0,dvcRecvList_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(8,0,1,1,dvcRecvList_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(9,0,1,-1,dvcRecvList_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
		ScaLBL_D3Q19_Unpack(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(5,-1,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(7,1,0,-1,dvcRecvList_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(7,0,-1,-1,dvcRecvList_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(9,0,1,-1,dvcRecvList_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
		ScaLBL_D3Q19_Unpack(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(6,1,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(6,-1,0,1,dvcRecvList_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(8,0,1,1,dvcRecvList_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(8,0,-1,1,dvcRecvList_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		//..................................................................................
		//...Map recieve list for the xy edge <<<8)................................
		ScaLBL_D3Q19_Unpack(3,-1,-1,0,dvcRecvList_XY,0,recvCount_XY,recvbuf_XY,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xy edge <<<9)................................
		ScaLBL_D3Q19_Unpack(5,1,-1,0,dvcRecvList_xY,0,recvCount_xY,recvbuf_xY,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xY edge <<<10)................................
		ScaLBL_D3Q19_Unpack(4,-1,1,0,dvcRecvList_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the XY edge <<<7)................................
		ScaLBL_D3Q19_Unpack(4,1,1,0,dvcRecvList_xy,0,recvCount_xy,recvbuf_xy,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xz edge <<<12)................................
		ScaLBL_D3Q19_Unpack(5,-1,0,-1,dvcRecvList_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the xZ edge <<<14)................................
		ScaLBL_D3Q19_Unpack(6,-1,0,1,dvcRecvList_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xz edge <<<13)................................
		ScaLBL_D3Q19_Unpack(7,1,0,-1,dvcRecvList_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the XZ edge <<<11)................................
		ScaLBL_D3Q19_Unpack(6,1,0,1,dvcRecvList_xz,0,recvCount_xz,recvbuf_xz,f_even,Nx,Ny,Nz);
		//...Map recieve list for the yz edge <<<16)................................
		ScaLBL_D3Q19_Unpack(7,0,-1,-1,dvcRecvList_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the yZ edge <<<18)................................
		ScaLBL_D3Q19_Unpack(8,0,-1,1,dvcRecvList_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Yz edge <<<17)................................
		ScaLBL_D3Q19_Unpack(9,0,1,-1,dvcRecvList_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the YZ edge <<<15)................................
		ScaLBL_D3Q19_Unpack(8,0,1,1,dvcRecvList_yz,0,recvCount_yz,recvbuf_yz,f_even,Nx,Ny,Nz);
		//...................................................................................

		//...................................................................................
		ScaLBL_D3Q19_Pack(1,dvcSendList_x,0,sendCount_x,sendbuf_x,A_even,N);
		ScaLBL_D3Q19_Pack(1,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,B_even,N);
		//...Packing for X face(1,7,9,11,13)................................
		ScaLBL_D3Q19_Pack(0,dvcSendList_X,0,sendCount_X,sendbuf_X,A_odd,N);
		ScaLBL_D3Q19_Pack(0,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,B_odd,N);
		//...Packing for y face(4,8,9,16,18).................................
		ScaLBL_D3Q19_Pack(2,dvcSendList_y,0,sendCount_y,sendbuf_y,A_even,N);
		ScaLBL_D3Q19_Pack(2,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,B_even,N);
		//...Packing for Y face(3,7,10,15,17).................................
		ScaLBL_D3Q19_Pack(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,A_odd,N);
		ScaLBL_D3Q19_Pack(1,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,B_odd,N);
		//...Packing for z face(6,12,13,16,17)................................
		ScaLBL_D3Q19_Pack(3,dvcSendList_z,0,sendCount_z,sendbuf_z,A_even,N);
		ScaLBL_D3Q19_Pack(3,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,B_even,N);
		//...Packing for Z face(5,11,14,15,18)................................
		ScaLBL_D3Q19_Pack(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,A_odd,N);
		ScaLBL_D3Q19_Pack(2,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,B_odd,N);
		//...................................................................................

		//...................................................................................
		// Send all the distributions
		MPI_Isend(sendbuf_x, 2*sendCount_x,MPI_DOUBLE,rank_x,sendtag,comm,&req1[0]);
		MPI_Irecv(recvbuf_X, 2*recvCount_X,MPI_DOUBLE,rank_X,recvtag,comm,&req2[0]);
		MPI_Isend(sendbuf_X, 2*sendCount_X,MPI_DOUBLE,rank_X,sendtag,comm,&req1[1]);
		MPI_Irecv(recvbuf_x, 2*recvCount_x,MPI_DOUBLE,rank_x,recvtag,comm,&req2[1]);
		MPI_Isend(sendbuf_y, 2*sendCount_y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[2]);
		MPI_Irecv(recvbuf_Y, 2*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[2]);
		MPI_Isend(sendbuf_Y, 2*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[3]);
		MPI_Irecv(recvbuf_y, 2*recvCount_y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[3]);
		MPI_Isend(sendbuf_z, 2*sendCount_z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[4]);
		MPI_Irecv(recvbuf_Z, 2*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[4]);
		MPI_Isend(sendbuf_Z, 2*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[5]);
		MPI_Irecv(recvbuf_z, 2*recvCount_z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[5]);
		//...................................................................................
		
		ScaLBL_D3Q7_Swap(ID, A_even, A_odd, Nx, Ny, Nz);
		ScaLBL_D3Q7_Swap(ID, B_even, B_odd, Nx, Ny, Nz);
		
		//...................................................................................
		// Wait for completion of D3Q19 communication
		MPI_Waitall(6,req1,stat1);
		MPI_Waitall(6,req2,stat2);
		//...................................................................................
		// Unpack the distributions on the device
		//...................................................................................
		//...Map recieve list for the X face: q=2,8,10,12,13 .................................
		ScaLBL_D3Q19_Unpack(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,A_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(0,-1,0,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,B_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the x face: q=1,7,9,11,13..................................
		ScaLBL_D3Q19_Unpack(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,A_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(1,1,0,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,B_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
		ScaLBL_D3Q19_Unpack(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,A_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(1,0,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,B_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
		ScaLBL_D3Q19_Unpack(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,A_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(2,0,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,B_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
		ScaLBL_D3Q19_Unpack(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,A_odd,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(2,0,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,B_odd,Nx,Ny,Nz);
		//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
		ScaLBL_D3Q19_Unpack(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,A_even,Nx,Ny,Nz);
		ScaLBL_D3Q19_Unpack(3,0,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,B_even,Nx,Ny,Nz);
		//..................................................................................

		//..................................................................................
		ScaLBL_D3Q7_Density(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
		ScaLBL_D3Q7_Density(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
		
		//*************************************************************************
		// 		Compute the phase indicator field 
		//*************************************************************************
//		ScaLBL_ComputePhaseField(ID, Phi, Copy, Den, N);
		ScaLBL_ComputePhaseField(ID, Phi, Den, N);
		//*************************************************************************

		//...................................................................................
		ScaLBL_Scalar_Pack(dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_xy, sendCount_xy,sendbuf_xy, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_xY, sendCount_xY,sendbuf_xY, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_Xy, sendCount_Xy,sendbuf_Xy, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_XY, sendCount_XY,sendbuf_XY, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_xz, sendCount_xz,sendbuf_xz, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_xZ, sendCount_xZ,sendbuf_xZ, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_Xz, sendCount_Xz,sendbuf_Xz, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_XZ, sendCount_XZ,sendbuf_XZ, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_yz, sendCount_yz,sendbuf_yz, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_yZ, sendCount_yZ,sendbuf_yZ, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_Yz, sendCount_Yz,sendbuf_Yz, Phi, N);
		ScaLBL_Scalar_Pack(dvcSendList_YZ, sendCount_YZ,sendbuf_YZ, Phi, N);
		//...................................................................................
		// Send / Recv all the phase indcator field values
		//...................................................................................
		MPI_Isend(sendbuf_x, sendCount_x,MPI_DOUBLE,rank_x,sendtag,comm,&req1[0]);
		MPI_Irecv(recvbuf_X, recvCount_X,MPI_DOUBLE,rank_X,recvtag,comm,&req2[0]);
		MPI_Isend(sendbuf_X, sendCount_X,MPI_DOUBLE,rank_X,sendtag,comm,&req1[1]);
		MPI_Irecv(recvbuf_x, recvCount_x,MPI_DOUBLE,rank_x,recvtag,comm,&req2[1]);
		MPI_Isend(sendbuf_y, sendCount_y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[2]);
		MPI_Irecv(recvbuf_Y, recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[2]);
		MPI_Isend(sendbuf_Y, sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[3]);
		MPI_Irecv(recvbuf_y, recvCount_y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[3]);
		MPI_Isend(sendbuf_z, sendCount_z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[4]);
		MPI_Irecv(recvbuf_Z, recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[4]);
		MPI_Isend(sendbuf_Z, sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[5]);
		MPI_Irecv(recvbuf_z, recvCount_z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[5]);
		MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,comm,&req1[6]);
		MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,comm,&req2[6]);
		MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,comm,&req1[7]);
		MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,comm,&req2[7]);
		MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,comm,&req1[8]);
		MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,comm,&req2[8]);
		MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,comm,&req1[9]);
		MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,comm,&req2[9]);
		MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,comm,&req1[10]);
		MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,comm,&req2[10]);
		MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,comm,&req1[11]);
		MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,comm,&req2[11]);
		MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,comm,&req1[12]);
		MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,comm,&req2[12]);
		MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,comm,&req1[13]);
		MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,comm,&req2[13]);
		MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,comm,&req1[14]);
		MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,comm,&req2[14]);
		MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,comm,&req1[15]);
		MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,comm,&req2[15]);
		MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,comm,&req1[16]);
		MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,comm,&req2[16]);
		MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,comm,&req1[17]);
		MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,comm,&req2[17]);
		//...................................................................................
		//...................................................................................
		// Wait for completion of Indicator Field communication
		//...................................................................................
		MPI_Waitall(18,req1,stat1);
		MPI_Waitall(18,req2,stat2);
		ScaLBL_DeviceBarrier();
		//...................................................................................
		//...................................................................................
/*		ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
		ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
		ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
		ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
		ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
		ScaLBL_Scalar_Unpack(faceGrid, packThreads, dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
*/		
		ScaLBL_Scalar_Unpack(dvcRecvList_x, recvCount_x,recvbuf_x, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_y, recvCount_y,recvbuf_y, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_z, recvCount_z,recvbuf_z, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_X, recvCount_X,recvbuf_X, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_Y, recvCount_Y,recvbuf_Y, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_Z, recvCount_Z,recvbuf_Z, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_xy, recvCount_xy,recvbuf_xy, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_xY, recvCount_xY,recvbuf_xY, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_XY, recvCount_XY,recvbuf_XY, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_xz, recvCount_xz,recvbuf_xz, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_yz, recvCount_yz,recvbuf_yz, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, Phi, N);
		ScaLBL_Scalar_Unpack(dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, Phi, N);
		//...................................................................................
	
		
		if (pBC && kproc == 0)	{
			ScaLBL_D3Q19_Pressure_BC_z(f_even,f_odd,din,Nx,Ny,Nz);
	        ScaLBL_D3Q19_Velocity(ID,f_even,f_odd,Velocity,Nx,Ny,Nz);
			ScaLBL_Color_BC_z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		}
			
		if (pBC && kproc == nprocz-1){
			ScaLBL_D3Q19_Pressure_BC_Z(f_even,f_odd,dout,Nx,Ny,Nz,Nx*Ny*(Nz-2));
	        ScaLBL_D3Q19_Velocity(ID,f_even,f_odd,Velocity,Nx,Ny,Nz);
			ScaLBL_Color_BC_Z(Phi,Den,A_even,A_odd,B_even,B_odd,Nx,Ny,Nz);
		}
		
		//...................................................................................

		MPI_Barrier(comm);

		// Timestep completed!
		timestep++;
		//...................................................................
		if (timestep%1000 == 995){
			//...........................................................................
			// Copy the phase indicator field for the earlier timestep
			ScaLBL_DeviceBarrier();
			ScaLBL_CopyToHost(Phase_tplus.data,Phi,N*sizeof(double));
			//...........................................................................
		}
		if (timestep%1000 == 0){
			//...........................................................................
			// Copy the data for for the analysis timestep
			//...........................................................................
			// Copy the phase from the GPU -> CPU
			//...........................................................................
			ScaLBL_DeviceBarrier();
			ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
			ScaLBL_CopyToHost(Phase.data,Phi,N*sizeof(double));
			ScaLBL_CopyToHost(Press.data,Pressure,N*sizeof(double));
			ScaLBL_CopyToHost(Vel_x.data,&Velocity[0],N*sizeof(double));
			ScaLBL_CopyToHost(Vel_y.data,&Velocity[N],N*sizeof(double));
			ScaLBL_CopyToHost(Vel_z.data,&Velocity[2*N],N*sizeof(double));
			MPI_Barrier(comm);
		}
		if (timestep%1000 == 5){
			//...........................................................................
			// Copy the phase indicator field for the later timestep
			ScaLBL_DeviceBarrier();
			ScaLBL_CopyToHost(Phase_tminus.data,Phi,N*sizeof(double));
			//...........................................................................
			// Calculate the time derivative of the phase indicator field
			for (n=0; n<N; n++)	dPdt(n) = 0.1*(Phase_tplus(n) - Phase_tminus(n));
			//...........................................................................

			//...........................................................................
			CommunicateMeshHalo(Phase, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			
			// Compute the gradients of the phase indicator and signed distance fields
			pmmc_MeshGradient(Phase,Phase_x,Phase_y,Phase_z,Nx,Ny,Nz);
//			pmmc_MeshGradient(SignDist,SignDist_x,SignDist_y,SignDist_z,Nx,Ny,Nz);
			//...........................................................................
			// Compute the mesh curvature of the phase indicator field
			pmmc_MeshCurvature(Phase, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
			//...........................................................................
			// Fill in the halo region for the mesh gradients and curvature
			//...........................................................................
			// Pressure
			//...........................................................................
			CommunicateMeshHalo(Vel_x, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			// Velocity
			//...........................................................................
			CommunicateMeshHalo(Press, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			//...........................................................................
			CommunicateMeshHalo(Vel_y, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			//...........................................................................
			CommunicateMeshHalo(Vel_z, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			// Mean Curvature
			//...........................................................................
			CommunicateMeshHalo(MeanCurvature, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			//...........................................................................
			// Gaussian Curvature
			//...........................................................................
			CommunicateMeshHalo(GaussCurvature, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			// Gradient of the phase indicator field
			//...........................................................................
			CommunicateMeshHalo(Phase_x, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			CommunicateMeshHalo(Phase_y, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			CommunicateMeshHalo(Phase_z, comm,
					sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
					sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
					sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
					recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
					recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
					recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
					sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
					sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
					sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
					sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
					sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
					sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
					recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
					recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
					recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
					recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
					recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
					recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
					rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
					rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
			//...........................................................................
			nc = 0;	nblobs = 0;
			for (k=0;k<Nz-1;k++){
				for (j=0;j<Ny-1;j++){
					for (i=0;i<Nx-1;i++){
						if ( LocalBlobID(i,j,k) == -1 ){
							if ( Phase(i,j,k) > 0.0 ){
								if ( SignDist(i,j,k) > 0.0 ){
									// node i,j,k is in the porespace
									Blobs(nblobs) = ComputeBlob(LocalBlobCubeList,nblobs,nc,LocalBlobID,Phase,SignDist,0.0,0.0,i,j,k,temp);
									nblobs++;
									INSIST(nblobs < MAX_LOCAL_BLOB_COUNT,"Not enough to store the local blobs!");
								}
							}
						}
					}
				}
			}
			count_in=0,count_out=0;
			for (k=0;k<Nz-1;k++){
				for (j=0;j<Ny-1;j++){
					for (i=0;i<Nx-1;i++){
						// Loop over cube corners
						add=1;				// initialize to true - add unless corner occupied by nw-phase
						for (p=0;p<8;p++){
							nodx=i+cube[p][0];
							nody=j+cube[p][1];
							nodz=k+cube[p][2];
							if ( indicator(nodx,nody,nodz) > -1 ){
								// corner occupied by nw-phase  -> do not add
								add = 0;
							}
						}
						if ( add == 1 ){
							blobs(0,ncubes) = i;
							blobs(1,ncubes) = j;
							blobs(2,ncubes) = k;
							ncubes++;
							count_in++;
						}
						else { count_out++; }
					}
				}
			}
			b(nblobs) = count_in;
			nblobs++;
			//...........................................................................
			// Fill in the connected components graph in LocalBlobTable
			//...........................................................................

			//...........................................................................
			// Compute areas using porous medium marching cubes algorithm
			// McClure, Adalsteinsson, et al. (2007)
			//...........................................................................
			awn = aws = ans = lwns = 0.0;
			nwp_volume = 0.0;
			As = 0.0;
			
			// Compute phase averages
			pan = paw = 0.0;
			vaw(0) = vaw(1) = vaw(2) = 0.0;
			van(0) = van(1) = van(2) = 0.0;
			vawn(0) = vawn(1) = vawn(2) = 0.0;
			Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
			Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
			Gws(0) = Gws(1) = Gws(2) = 0.0;
			Gws(3) = Gws(4) = Gws(5) = 0.0;
			Gns(0) = Gns(1) = Gns(2) = 0.0;
			Gns(3) = Gns(4) = Gns(5) = 0.0;
			vol_w = vol_n =0.0;
			Jwn = Kwn = efawns = 0.0;
			trJwn = trawn = 0.0;
			
			start = 0;
			for (a=0;a<nblobs;a++){
				finish = start+cubes_in_blob(a);
				for (c=start;c<finish;c++){
					// Get cube from the list
					i = LocalBlobCubeList(0,c);
					j = LocalBlobCubeList(1,c);
					k = LocalBlobCubeList(2,c);

					//...........................................................................
					// Initialize local blob averages
					blob_awn = blob_aws = blob_ans = blob_lwns = 0.0;
					blob_nwp_volume = 0.0;
					blob_As = 0.0;
					blob_pan = blob_paw = 0.0;
					blob_vaw(0) = blob_vaw(1) = blob_vaw(2) = 0.0;
					blob_van(0) = blob_van(1) = blob_van(2) = 0.0;
					blob_vawn(0) = blob_vawn(1) = blob_vawn(2) = 0.0;
					blob_Gwn(0) = blob_Gwn(1) = blob_Gwn(2) = 0.0;
					blob_Gwn(3) = blob_Gwn(4) = blob_Gwn(5) = 0.0;
					blob_Gws(0) = blob_Gws(1) = blob_Gws(2) = 0.0;
					blob_Gws(3) = blob_Gws(4) = blob_Gws(5) = 0.0;
					blob_Gns(0) = blob_Gns(1) = blob_Gns(2) = 0.0;
					blob_Gns(3) = blob_Gns(4) = blob_Gns(5) = 0.0;
					blob_vol_w = blob_vol_n =0.0;
					blob_Jwn = blob_Kwn = blob_efawns = 0.0;
					blob_trJwn = blob_trawn = 0.0;
					//...........................................................................

					//...........................................................................
					// Compute volume averages
					for (p=0;p<8;p++){
						if ( SignDist(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){

							// 1-D index for this cube corner
							n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;

							// Compute the non-wetting phase volume contribution
							if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 )
								blob_nwp_volume += 0.125;

							// volume averages over the non-wetting phase
							if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0.99 ){
								// volume the excludes the interfacial region
								blob_vol_n += 0.125;
								// pressure
								blob_pan += 0.125*Press.data[n];
								// velocity
								blob_van(0) += 0.125*Vel_x.data[n];
								blob_van(1) += 0.125*Vel_y.data[n];
								blob_van(2) += 0.125*Vel_z.data[n];
							}

							// volume averages over the wetting phase
							if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) < -0.99 ){
								// volume the excludes the interfacial region
								blob_vol_w += 0.125;
								// pressure
								blob_paw += 0.125*Press.data[n];
								// velocity
								blob_vaw(0) += 0.125*Vel_x.data[n];
								blob_vaw(1) += 0.125*Vel_y.data[n];
								blob_vaw(2) += 0.125*Vel_z.data[n];
							}
						}
					}

					//...........................................................................
					// Construct the interfaces and common curve
					pmmc_ConstructLocalCube(SignDist, Phase, solid_isovalue, fluid_isovalue,
							nw_pts, nw_tris, values, ns_pts, ns_tris, ws_pts, ws_tris,
							local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
							n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
							n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
							i, j, k, Nx, Ny, Nz);

					// Integrate the contact angle
					blob_efawns += pmmc_CubeContactAngle(CubeValues,Values,Phase_x,Phase_y,Phase_z,SignDist_x,SignDist_y,SignDist_z,
							local_nws_pts,i,j,k,n_local_nws_pts);

					// Integrate the mean curvature
					blob_Jwn    += pmmc_CubeSurfaceInterpValue(CubeValues,MeanCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
					blob_Kwn    += pmmc_CubeSurfaceInterpValue(CubeValues,GaussCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);

					// Integrate the trimmed mean curvature (hard-coded to use a distance of 4 pixels)
					pmmc_CubeTrimSurfaceInterpValues(CubeValues,MeanCurvature,SignDist,nw_pts,nw_tris,Values,DistValues,
							i,j,k,n_nw_pts,n_nw_tris,trimdist,blob_trawn,blob_trJwn);

					// Compute the normal speed of the interface
					pmmc_InterfaceSpeed(dPdt, Phase_x, Phase_y, Phase_z, CubeValues, nw_pts, nw_tris,
							NormalVector, InterfaceSpeed, blob_vawn, i, j, k, n_nw_pts, n_nw_tris);

					blob_As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);

					// Compute the surface orientation and the interfacial area
					blob_awn += pmmc_CubeSurfaceOrientation(blob_Gwn,nw_pts,nw_tris,n_nw_tris);
					blob_ans += pmmc_CubeSurfaceOrientation(blob_Gns,ns_pts,ns_tris,n_ns_tris);
					blob_aws += pmmc_CubeSurfaceOrientation(blob_Gws,ws_pts,ws_tris,n_ws_tris);
					blob_lwns +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
					//...........................................................................
				}
				// Save blob averages to the table
				LocalBlobAverages(0,a) = blob_vol_n;
				LocalBlobAverages(1,a) = blob_pan;
				LocalBlobAverages(2,a) = blob_awn;
				LocalBlobAverages(3,a) = blob_ans;
				LocalBlobAverages(4,a) = blob_lwns;
				LocalBlobAverages(5,a) = blob_Jwn/blob_awn;
				LocalBlobAverages(6,a) = blob_trJwn/blob_trawn;
				LocalBlobAverages(7,a) = blob_Kwn/blob_awn;
				LocalBlobAverages(8,a) = blob_efawns/blob_lwns;
				LocalBlobAverages(9,a) = blob_van(0);
				LocalBlobAverages(10,a) = blob_van(1);
				LocalBlobAverages(11,a) = blob_van(2);
				LocalBlobAverages(12,a) = blob_vawn(0);
				LocalBlobAverages(13,a) = blob_vawn(1);
				LocalBlobAverages(14,a) = blob_vawn(2);	
				LocalBlobAverages(15,a) = blob_Gwn(0);
				LocalBlobAverages(16,a) = blob_Gwn(1);
				LocalBlobAverages(17,a) = blob_Gwn(2);
				LocalBlobAverages(18,a) = blob_Gwn(3);
				LocalBlobAverages(19,a) = blob_Gwn(4);
				LocalBlobAverages(20,a) = blob_Gwn(5);
				LocalBlobAverages(21,a) = blob_Gns(0);
				LocalBlobAverages(22,a) = blob_Gns(1);
				LocalBlobAverages(23,a) = blob_Gns(2);
				LocalBlobAverages(24,a) = blob_Gns(3);
				LocalBlobAverages(25,a) = blob_Gns(4);
				LocalBlobAverages(26,a) = blob_Gns(5);
				
				// Update the local averages from the blob averages
				awn += blob_awn;
				aws += blob_aws;
				ans += blob_ans;
				lwns += blob_lwns;
				nwp_volume += blob_nwp_volume;
				As += blob_As;
				pan += blob_pan;
				paw += blob_paw;
				vaw(0) += blob_vaw(0);
				vaw(1) += blob_vaw(1);
				vaw(2) += blob_vaw(2);
				van(0) += blob_van(0);
				van(1) += blob_van(1);
				van(2) += blob_van(2);
				vawn(0) += blob_vawn(0);
				vawn(1) += blob_vawn(1);
				vawn(2) += blob_vawn(2);
				Gwn(0) += blob_Gwn(0);
				Gwn(1) += blob_Gwn(1);
				Gwn(2) += blob_Gwn(2);
				Gwn(3) += blob_Gwn(3);
				Gwn(4) += blob_Gwn(4); 
				Gwn(5) += blob_Gwn(5);
				Gws(0) += blob_Gws(0);
				Gws(1) += blob_Gws(1);
				Gws(2) += blob_Gws(2);
				Gws(3) += blob_Gws(3);
				Gws(4) += blob_Gws(4);
				Gws(5) += blob_Gws(5);
				Gns(0) += blob_Gns(0);
				Gns(1) += blob_Gns(1);
				Gns(2) += blob_Gns(2);
				Gns(3) += blob_Gns(3);
				Gns(4) += blob_Gns(4); 
				Gns(5) += blob_Gns(5);
				vol_w += blob_vol_w;
				vol_n += blob_vol_n;
				Jwn += blob_Jwn;
				Kwn += blob_Kwn;
				efawns += blob_efawns;
				trJwn += blob_trJwn;
				trawn += blob_trawn;
				
				
			}
			
			//...........................................................................
			MPI_Barrier(comm);
			MPI_Allreduce(&nwp_volume,&nwp_volume_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&awn,&awn_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&ans,&ans_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&aws,&aws_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&lwns,&lwns_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&As,&As_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&Jwn,&Jwn_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&Kwn,&Kwn_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&efawns,&efawns_global,1,MPI_DOUBLE,MPI_SUM,comm);
			// Phase averages
			MPI_Allreduce(&vol_w,&vol_w_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&vol_n,&vol_n_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&paw,&paw_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&pan,&pan_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&vaw(0),&vaw_global(0),3,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&van(0),&van_global(0),3,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&vawn(0),&vawn_global(0),3,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&Gwn(0),&Gwn_global(0),6,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&Gns(0),&Gns_global(0),6,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&Gws(0),&Gws_global(0),6,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&trawn,&trawn_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&trJwn,&trJwn_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Barrier(comm);
			//.........................................................................
			// Compute the change in the total surface energy based on the defined interval
			// See McClure, Prins and Miller (2014) 
			//.........................................................................
			dAwn += awn_global;
			dAns += ans_global;
			dEs = 6.01603*alpha*(dAwn + 1.05332*Ps*dAns);
			dAwn = -awn_global;		// Get ready for the next analysis interval
			dAns = -ans_global;
			
			// Normalize the phase averages 
			// (density of both components = 1.0)
			paw_global = paw_global / vol_w_global;
			vaw_global(0) = vaw_global(0) / vol_w_global;
			vaw_global(1) = vaw_global(1) / vol_w_global;
			vaw_global(2) = vaw_global(2) / vol_w_global;
			pan_global = pan_global / vol_n_global;
			van_global(0) = van_global(0) / vol_n_global;
			van_global(1) = van_global(1) / vol_n_global;
			van_global(2) = van_global(2) / vol_n_global;
			
			// Normalize surface averages by the interfacial area
			Jwn_global /= awn_global;
			Kwn_global /= awn_global;
			efawns_global /= lwns_global;
			
			if (trawn_global > 0.0)	trJwn_global /= trawn_global;

			if (awn_global > 0.0)	for (i=0; i<3; i++)		vawn_global(i) /= awn_global;
			if (awn_global > 0.0)	for (i=0; i<6; i++)		Gwn_global(i) /= awn_global;
			if (ans_global > 0.0)	for (i=0; i<6; i++)		Gns_global(i) /= ans_global;
			if (aws_global > 0.0)	for (i=0; i<6; i++)		Gws_global(i) /= aws_global;

			sat_w = 1.0 - nwp_volume_global*iVol_global/porosity;
			// Compute the specific interfacial areas and common line length (dimensionless per unit volume)
			awn_global = awn_global*iVol_global*D;
			ans_global = ans_global*iVol_global*D;
			aws_global = aws_global*iVol_global*D;
			dEs = dEs*iVol_global*D;
			lwns_global = lwns_global*iVol_global*D*D;
			
			//.........................................................................
			if (rank==0){

				fprintf(TIMELOG,"%i %.5g ",timestep-5,dEs);										// change in surface energy
				fprintf(TIMELOG,"%.5g %.5g %.5g ",sat_w,paw_global,pan_global);					// saturation and pressure
				fprintf(TIMELOG,"%.5g %.5g %.5g ",awn_global,ans_global,aws_global);				// interfacial areas
				fprintf(TIMELOG,"%.5g %5g ",Jwn_global, Kwn_global);								// curvature of wn interface
				fprintf(TIMELOG,"%.5g ",lwns_global);											// common curve length
				fprintf(TIMELOG,"%.5g ",efawns_global);											// average contact angle
				fprintf(TIMELOG,"%.5g %.5g %.5g ",vaw_global(0),vaw_global(1),vaw_global(2));	// average velocity of w phase
				fprintf(TIMELOG,"%.5g %.5g %.5g ",van_global(0),van_global(1),van_global(2));	// average velocity of n phase
				fprintf(TIMELOG,"%.5g %.5g %.5g ",vawn_global(0),vawn_global(1),vawn_global(2));	// velocity of wn interface
				fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
						Gwn_global(0),Gwn_global(1),Gwn_global(2),Gwn_global(3),Gwn_global(4),Gwn_global(5));	// orientation of wn interface
				fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
						Gns_global(0),Gns_global(1),Gns_global(2),Gns_global(3),Gns_global(4),Gns_global(5));	// orientation of ns interface
				fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
						Gws_global(0),Gws_global(1),Gws_global(2),Gws_global(3),Gws_global(4),Gws_global(5));	// orientation of ws interface
				fprintf(TIMELOG,"%.5g %5g \n",trawn_global, trJwn_global);						// Trimmed curvature
				fflush(TIMELOG);
			}
			
		}
		
		if (timestep%RESTART_INTERVAL == 0){
			// Copy the data to the CPU
			ScaLBL_CopyToHost(cDistEven,f_even,10*N*sizeof(double));
			ScaLBL_CopyToHost(cDistOdd,f_odd,9*N*sizeof(double));
			ScaLBL_CopyToHost(cDen,Den,2*N*sizeof(double));
			// Read in the restart file to CPU buffers
			WriteCheckpoint(LocalRestartFile, cDen, cDistEven, cDistOdd, N);

#ifdef WRITE_SURFACES
			sprintf(tmpstr,"vis%03d",logcount);
			if (rank==0){
				mkdir(tmpstr,0777);
			}
			MPI_Barrier(comm);
			
			FILE *WN_TRIS;
			sprintf(LocalRankFilename,"%s/%s%s",tmpstr,"wn-tris.",LocalRankString);
			WN_TRIS = fopen(LocalRankFilename,"wb");

			FILE *NS_TRIS;
			sprintf(LocalRankFilename,"%s/%s%s",tmpstr,"ns-tris.",LocalRankString);
			NS_TRIS = fopen(LocalRankFilename,"wb");

			FILE *WS_TRIS;
			sprintf(LocalRankFilename,"%s/%s%s",tmpstr,"ws-tris.",LocalRankString);
			WS_TRIS = fopen(LocalRankFilename,"wb");

			FILE *WNS_PTS;
			sprintf(LocalRankFilename,"%s/%s%s",tmpstr,"wns-crv.",LocalRankString);
			WNS_PTS = fopen(LocalRankFilename,"wb");

			for (c=0;c<ncubes;c++){
				// Get cube from the list
				i = LocalBlobCubeList(0,c);
				j = LocalBlobCubeList(1,c);
				k = LocalBlobCubeList(2,c);
				//...........................................................................
				// Construct the interfaces and common curve
				pmmc_ConstructLocalCube(SignDist, Phase, solid_isovalue, fluid_isovalue,
						nw_pts, nw_tris, values, ns_pts, ns_tris, ws_pts, ws_tris,
						local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
						n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
						n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
						i, j, k, Nx, Ny, Nz);

				//.......................................................................................
				// Write the triangle lists to text file
				for (int r=0;r<n_nw_tris;r++){
					A = nw_pts(nw_tris(0,r));
					B = nw_pts(nw_tris(1,r));
					C = nw_pts(nw_tris(2,r));
					// compare the trianlge orientation against the color gradient
					// Orientation of the triangle
					double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
					double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
					double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

					double normal_x = Phase_x(i,j,k);
					double normal_y = Phase_y(i,j,k);
					double normal_z = Phase_z(i,j,k);

					// If the normals don't point in the same direction, flip the orientation of the triangle
					// Right hand rule for triangle orientation is used to determine rendering for most software
					if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
						P = A;
						A = C;
						C = P;
					}
					// Remap the points
					A.x += 1.0*iproc*(Nx-2);
					A.y += 1.0*jproc*(Nx-2);
					A.z += 1.0*kproc*(Nx-2);
					B.x += 1.0*iproc*(Nx-2);
					B.y += 1.0*jproc*(Nx-2);
					B.z += 1.0*kproc*(Nx-2);
					C.x += 1.0*iproc*(Nx-2);
					C.y += 1.0*jproc*(Nx-2);
					C.z += 1.0*kproc*(Nx-2);
					// write the triangle
					fwrite(&A.x,sizeof(A.x),1,WN_TRIS);
					fwrite(&A.y,sizeof(A.y),1,WN_TRIS);
					fwrite(&A.z,sizeof(A.z),1,WN_TRIS);
					fwrite(&B.x,sizeof(B.x),1,WN_TRIS);
					fwrite(&B.y,sizeof(B.y),1,WN_TRIS);
					fwrite(&B.z,sizeof(B.z),1,WN_TRIS);
					fwrite(&C.x,sizeof(C.x),1,WN_TRIS);
					fwrite(&C.y,sizeof(C.y),1,WN_TRIS);
					fwrite(&C.z,sizeof(C.z),1,WN_TRIS);
				}		
				for (int r=0;r<n_ws_tris;r++){
					A = ws_pts(ws_tris(0,r));
					B = ws_pts(ws_tris(1,r));
					C = ws_pts(ws_tris(2,r));
					// Remap the points
					A.x += 1.0*iproc*(Nx-2);
					A.y += 1.0*jproc*(Nx-2);
					A.z += 1.0*kproc*(Nx-2);
					B.x += 1.0*iproc*(Nx-2);
					B.y += 1.0*jproc*(Nx-2);
					B.z += 1.0*kproc*(Nx-2);
					C.x += 1.0*iproc*(Nx-2);
					C.y += 1.0*jproc*(Nx-2);
					C.z += 1.0*kproc*(Nx-2);
					// write the triangle
					fwrite(&A.x,sizeof(A.x),1,WS_TRIS);
					fwrite(&A.y,sizeof(A.y),1,WS_TRIS);
					fwrite(&A.z,sizeof(A.z),1,WS_TRIS);
					fwrite(&B.x,sizeof(B.x),1,WS_TRIS);
					fwrite(&B.y,sizeof(B.y),1,WS_TRIS);
					fwrite(&B.z,sizeof(B.z),1,WS_TRIS);
					fwrite(&C.x,sizeof(C.x),1,WS_TRIS);
					fwrite(&C.y,sizeof(C.y),1,WS_TRIS);
					fwrite(&C.z,sizeof(C.z),1,WS_TRIS);	
				}
				for (int r=0;r<n_ns_tris;r++){
					A = ns_pts(ns_tris(0,r));
					B = ns_pts(ns_tris(1,r));
					C = ns_pts(ns_tris(2,r));
					// Remap the points
					A.x += 1.0*iproc*(Nx-2);
					A.y += 1.0*jproc*(Nx-2);
					A.z += 1.0*kproc*(Nx-2);
					B.x += 1.0*iproc*(Nx-2);
					B.y += 1.0*jproc*(Nx-2);
					B.z += 1.0*kproc*(Nx-2);
					C.x += 1.0*iproc*(Nx-2);
					C.y += 1.0*jproc*(Nx-2);
					C.z += 1.0*kproc*(Nx-2);
					// write the triangle
					fwrite(&A.x,sizeof(A.x),1,NS_TRIS);
					fwrite(&A.y,sizeof(A.y),1,NS_TRIS);
					fwrite(&A.z,sizeof(A.z),1,NS_TRIS);
					fwrite(&B.x,sizeof(B.x),1,NS_TRIS);
					fwrite(&B.y,sizeof(B.y),1,NS_TRIS);
					fwrite(&B.z,sizeof(B.z),1,NS_TRIS);
					fwrite(&C.x,sizeof(C.x),1,NS_TRIS);
					fwrite(&C.y,sizeof(C.y),1,NS_TRIS);
					fwrite(&C.z,sizeof(C.z),1,NS_TRIS);
				}
				for (int p=0; p < n_nws_pts; p++){
					P = nws_pts(p);
				//	fprintf(WNS_PTS,"%f %f %f \n",P.x, P.y, P.z);
					fwrite(&P.x,sizeof(P.x),1,WNS_PTS);
					fwrite(&P.y,sizeof(P.y),1,WNS_PTS);
					fwrite(&P.z,sizeof(P.z),1,WNS_PTS);
				}
			}
			fclose(WN_TRIS);
			fclose(NS_TRIS);
			fclose(WS_TRIS);
			fclose(WNS_PTS);
			logcount++;
#endif 
		}

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
	
	//	printf("Local File Name =  %s \n",LocalRankFilename);
	FILE *FINALSTATE;
	if (rank==0){
		fclose(TIMELOG);
		FINALSTATE= fopen("finalstate.tcat","a");
		fprintf(FINALSTATE,"%i %.5g ",timestep-5,dEs);										// change in surface energy
		fprintf(FINALSTATE,"%.5g %.5g %.5g ",sat_w,paw_global,pan_global);					// saturation and pressure
		fprintf(FINALSTATE,"%.5g %.5g %.5g ",awn_global,ans_global,aws_global);				// interfacial areas
		fprintf(FINALSTATE,"%.5g %5g ",Jwn_global, Kwn_global);								// curvature of wn interface
		fprintf(FINALSTATE,"%.5g ",lwns_global);											// common curve length
		fprintf(FINALSTATE,"%.5g ",efawns_global);											// average contact angle
		fprintf(FINALSTATE,"%.5g %.5g %.5g ",vaw_global(0),vaw_global(1),vaw_global(2));	// average velocity of w phase
		fprintf(FINALSTATE,"%.5g %.5g %.5g ",van_global(0),van_global(1),van_global(2));	// average velocity of n phase
		fprintf(FINALSTATE,"%.5g %.5g %.5g ",vawn_global(0),vawn_global(1),vawn_global(2));	// velocity of wn interface
		fprintf(FINALSTATE,"%.5g %.5g %.5g %.5g %.5g %.5g ",
				Gwn_global(0),Gwn_global(1),Gwn_global(2),Gwn_global(3),Gwn_global(4),Gwn_global(5));	// orientation of wn interface
		fprintf(FINALSTATE,"%.5g %.5g %.5g %.5g %.5g %.5g ",
				Gns_global(0),Gns_global(1),Gns_global(2),Gns_global(3),Gns_global(4),Gns_global(5));	// orientation of ns interface
		fprintf(FINALSTATE,"%.5g %.5g %.5g %.5g %.5g %.5g ",
				Gws_global(0),Gws_global(1),Gws_global(2),Gws_global(3),Gws_global(4),Gws_global(5));	// orientation of ws interface
		fprintf(FINALSTATE,"%.5g %5g \n",trawn_global, trJwn_global);						// Trimmed curvature
		fclose(FINALSTATE);
	}
	
	//************************************************************************/
	// Write out the phase indicator field 
	//************************************************************************/
	//	printf("Local File Name =  %s \n",LocalRankFilename);

	
//#ifdef WriteOutput	
	ScaLBL_CopyToHost(Phase.data,Phi,N*sizeof(double));
	sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
	FILE *PHASE;
	PHASE = fopen(LocalRankFilename,"wb");
	fwrite(Phase.data,8,N,PHASE);
//	fwrite(MeanCurvature.data,8,N,PHASE);
	fclose(PHASE);
//#endif
	ScaLBL_D3Q19_Pressure(ID,f_even,f_odd,Pressure,Nx,Ny,Nz);
	ScaLBL_CopyToHost(Press.data,Pressure,N*sizeof(double));
	CommunicateMeshHalo(Press, comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			sendList_x,sendList_y,sendList_z,sendList_X,sendList_Y,sendList_Z,
			sendList_xy,sendList_XY,sendList_xY,sendList_Xy,sendList_xz,sendList_XZ,
			sendList_xZ,sendList_Xz,sendList_yz,sendList_YZ,sendList_yZ,sendList_Yz,
			sendCount_x,sendCount_y,sendCount_z,sendCount_X,sendCount_Y,sendCount_Z,
			sendCount_xy,sendCount_XY,sendCount_xY,sendCount_Xy,sendCount_xz,sendCount_XZ,
			sendCount_xZ,sendCount_Xz,sendCount_yz,sendCount_YZ,sendCount_yZ,sendCount_Yz,
			recvList_x,recvList_y,recvList_z,recvList_X,recvList_Y,recvList_Z,
			recvList_xy,recvList_XY,recvList_xY,recvList_Xy,recvList_xz,recvList_XZ,
			recvList_xZ,recvList_Xz,recvList_yz,recvList_YZ,recvList_yZ,recvList_Yz,
			recvCount_x,recvCount_y,recvCount_z,recvCount_X,recvCount_Y,recvCount_Z,
			recvCount_xy,recvCount_XY,recvCount_xY,recvCount_Xy,recvCount_xz,recvCount_XZ,
			recvCount_xZ,recvCount_Xz,recvCount_yz,recvCount_YZ,recvCount_yZ,recvCount_Yz,
			rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z,rank_xy,rank_XY,rank_xY,
			rank_Xy,rank_xz,rank_XZ,rank_xZ,rank_Xz,rank_yz,rank_YZ,rank_yZ,rank_Yz);
	sprintf(LocalRankFilename,"%s%s","Pressure.",LocalRankString);
	FILE *PRESS;
	PRESS = fopen(LocalRankFilename,"wb");
	fwrite(Press.data,8,N,PRESS);
	fclose(PRESS);
	
/*	sprintf(LocalRankFilename,"%s%s","dPdt.",LocalRankString);
	FILE *SPEED;
	SPEED = fopen(LocalRankFilename,"wb");
	fwrite(dPdt.data,8,N,SPEED);
	fclose(SPEED);
	
	sprintf(LocalRankFilename,"%s%s","DenA.",LocalRankString);
	FILE *DENA;
	DENA = fopen(LocalRankFilename,"wb");
	fwrite(&Den[0],8,N,DENA);
	fclose(DENA);
	
	sprintf(LocalRankFilename,"%s%s","DenB.",LocalRankString);
	FILE *DENB;
	DENB = fopen(LocalRankFilename,"wb");
	fwrite(&Den[N],8,N,DENB);
	fclose(DENB);
	
	sprintf(LocalRankFilename,"%s%s","GradMag.",LocalRankString);
	FILE *GRAD;
	GRAD = fopen(LocalRankFilename,"wb");
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				Phase(i,j,k) = sqrt(Phase_x(i,j,k)*Phase_x(i,j,k) + Phase_y(i,j,k)*Phase_y(i,j,k) + Phase_z(i,j,k)*Phase_z(i,j,k));
			}
		}
	}
	fwrite(Phase.data,8,N,GRAD);
	fclose(GRAD);
*/	
/*	double *DensityValues;
	DensityValues = new double [2*N];
	ScaLBL_CopyToHost(DensityValues,Copy,2*N*sizeof(double));
	FILE *PHASE;
	PHASE = fopen(LocalRankFilename,"wb");
	fwrite(DensityValues,8,2*N,PHASE);
	fclose(PHASE);
*/	//************************************************************************/

	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
