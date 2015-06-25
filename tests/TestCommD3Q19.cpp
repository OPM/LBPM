//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Domain.h"
#include "D3Q19.h"
#include "D3Q7.h"
#include "Array.h"
#include "Extras.h"
#include "common/MPI_Helpers.h"

using namespace std;


extern void GlobalFlipInitD3Q19(double *dist_even, double *dist_odd, int Nx, int Ny, int Nz, 
								int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Set of Discrete velocities for the D3Q19 Model
	static int D3Q19[18][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
	{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
	{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};
	
	int q,i,j,k,n,N;
	int Cqx,Cqy,Cqz; // Discrete velocity
	int x,y,z;		// Global indices
	int xn,yn,zn; 	// Global indices of neighbor 
	int X,Y,Z;		// Global size
	X = Nx*nprocx;
	Y = Ny*nprocy;
	Z = Nz*nprocz;
	N = (Nx+2)*(Ny+2)*(Nz+2);	// size of the array including halo
	for (k=0; k<Nz; k++){ 
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				
				n = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;
				
				// Get the 'global' index
				x = iproc*Nx+i;
				y = jproc*Ny+j;
				z = kproc*Nz+k;
				for (q=0; q<9; q++){
					// Odd distribution
					Cqx = D3Q19[2*q][0];
					Cqy = D3Q19[2*q][1];
					Cqz = D3Q19[2*q][2];
					xn = x - Cqx;
					yn = y - Cqy;
					zn = z - Cqz;
					if (xn < 0) xn += nprocx*Nx;
					if (yn < 0) yn += nprocy*Ny;
					if (zn < 0) zn += nprocz*Nz;
					if (!(xn < nprocx*Nx)) xn -= nprocx*Nx;
					if (!(yn < nprocy*Ny)) yn -= nprocy*Ny;
					if (!(zn < nprocz*Nz)) zn -= nprocz*Nz;	
					
					dist_even[(q+1)*N+n] = (zn*X*Y+yn*X+xn) + (2*q+1)*0.01;
					
					// Odd distribution
					xn = x + Cqx;
					yn = y + Cqy;
					zn = z + Cqz;
					if (xn < 0) xn += nprocx*Nx;
					if (yn < 0) yn += nprocy*Ny;
					if (zn < 0) zn += nprocz*Nz;
					if (!(xn < nprocx*Nx)) xn -= nprocx*Nx;
					if (!(yn < nprocy*Ny)) yn -= nprocy*Ny;
					if (!(zn < nprocz*Nz)) zn -= nprocz*Nz;
					
					dist_odd[q*N+n] =  (zn*X*Y+yn*X+xn) + 2*(q+1)*0.01;
				
				}
			}
		}
	}
	
}

extern int GlobalCheckDebugDist(double *dist_even, double *dist_odd, int Nx, int Ny, int Nz, 
		int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{

	int returnValue = 0;
	int q,i,j,k,n,N;
	int Cqx,Cqy,Cqz; // Discrete velocity
	int x,y,z;		// Global indices
	int xn,yn,zn; 	// Global indices of neighbor 
	int X,Y,Z;		// Global size
	X = Nx*nprocx;
	Y = Ny*nprocy;
	Z = Nz*nprocz;
	N = (Nx+2)*(Ny+2)*(Nz+2);	// size of the array including halo
	for (k=0; k<Nz; k++){ 
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){

				n = (k+1)*(Nx+2)*(Ny+2) + (j+1)*(Nx+2) + i+1;

				// Get the 'global' index
				x = iproc*Nx+i;
				y = jproc*Ny+j;
				z = kproc*Nz+k;
				for (q=0; q<9; q++){

					if (dist_even[(q+1)*N+n] != (z*X*Y+y*X+x) + 2*(q+1)*0.01){
						printf("******************************************\n");
						printf("error in even distribution q = %i \n", 2*(q+1));
						printf("i,j,k= %i, %i, %i \n", x,y,z);
						printf("dist = %5.2f \n", dist_even[(q+1)*N+n]);
						printf("n= %i \n",z*X*Y+y*X+x);
						returnValue++;
					}


					if (dist_odd[q*N+n] !=  (z*X*Y+y*X+x) + (2*q+1)*0.01){
						printf("******************************************\n");
						printf("error in odd distribution q = %i \n", 2*q+1);
						printf("i,j,k= %i, %i, %i \n", x,y,z);
						printf("dist = %5.2f \n", dist_odd[q*N+n]);
						printf("n= %i \n",z*X*Y+y*X+x);
						returnValue++;
					}
				}
			}
		}
	}
	return returnValue;
}

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
int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	int rank,nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
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
		printf("Running Unit Test for D3Q19 MPI Communication	\n");
		printf("********************************************************\n");
	}

	// BGK Model parameters
	string FILENAME;
	unsigned int nBlocks, nthreads;
	int timestepMax, interval;
	double tau,Fx,Fy,Fz,tol;
	// Domain variables
	double Lx,Ly,Lz;
	int nspheres;
	int Nx,Ny,Nz;
	int i,j,k,n;

	if (rank==0){
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

		// **************************************************************
	}
	// **************************************************************
	// Broadcast simulation parameters from rank 0 to all other procs
	MPI_Barrier(MPI_COMM_WORLD);
	//.................................................
	MPI_Bcast(&Nx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ny,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Nz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nBlocks,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nthreads,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&timestepMax,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(&Nx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ny,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Nz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocy,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nprocz,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nspheres,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//.................................................
	MPI_Barrier(MPI_COMM_WORLD);
	// **************************************************************
	// **************************************************************

	if (nprocs != nprocx*nprocy*nprocz){
		printf("Fatal error in processor number! \n");
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
		abort();		
	}

	if (rank==0){
		printf("********************************************************\n");
		printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
		printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
		printf("********************************************************\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	kproc = rank/(nprocx*nprocy);
	jproc = (rank-nprocx*nprocy*kproc)/nprocx;
	iproc = rank-nprocx*nprocy*kproc-nprocz*jproc;

	int BoundaryCondition=0;
	Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);

	Nx += 2;
	Ny += 2;
	Nz += 2;
	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);
	
	//.......................................................................
	// Assign the phase ID field
	//.......................................................................
	if (rank==0) printf("Assigning phase ID from file \n");
	char LocalRankString[8];
	sprintf(LocalRankString,"%05d",rank);
	char LocalRankFilename[40];
	sprintf(LocalRankFilename,"ID.%05i",rank);
	
	char *id;
	id = new char[Nx*Ny*Nz];

	if (rank==0) printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
	FILE *IDFILE = fopen(LocalRankFilename,"rb");
	if (IDFILE==NULL) ERROR("Error opening file: ID.xxxxx");
	fread(id,1,N,IDFILE);
	fclose(IDFILE);
	// Setup the domain
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				Dm.id[n] = id[n];
			}
		}
	}
	Dm.CommInit(MPI_COMM_WORLD);

	//.......................................................................
	// Compute the media porosity
	//.......................................................................
	double sum;
	double sum_local=0.0, porosity, iVol_global;
	iVol_global = 1.0/Nx/Ny/Nz/nprocx/nprocy/nprocz;
	char component = 0; // solid phase
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] == component){
					sum_local+=1.0;
				}
			}
		}
	}
	MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	porosity = 1.0-sum*iVol_global;
	if (rank==0) printf("Media porosity = %f \n",porosity);
	//.......................................................................

	//...........................................................................
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) cout << "Domain set." << endl;
	//...........................................................................

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
	MPI_Barrier(MPI_COMM_WORLD);
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
	CommunicateSendRecvCounts( MPI_COMM_WORLD, sendtag, recvtag, 
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
	CommunicateRecvLists( MPI_COMM_WORLD, sendtag, recvtag, 
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
	AllocateDeviceMemory((void **) &sendbuf_x, 5*sendCount_x*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_X, 5*sendCount_X*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_y, 5*sendCount_y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Y, 5*sendCount_Y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_z, 5*sendCount_z*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Z, 5*sendCount_Z*sizeof(double));	// Allocatevoid * memory
	AllocateDeviceMemory((void **) &sendbuf_xy, sendCount_xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_xY, sendCount_xY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Xy, sendCount_Xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_XY, sendCount_XY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_xz, sendCount_xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_xZ, sendCount_xZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Xz, sendCount_Xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_XZ, sendCount_XZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_yz, sendCount_yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_yZ, sendCount_yZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_Yz, sendCount_Yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &sendbuf_YZ, sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	AllocateDeviceMemory((void **) &recvbuf_x, 5*recvCount_x*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_X, 5*recvCount_X*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_y, 5*recvCount_y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Y, 5*recvCount_Y*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_z, 5*recvCount_z*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Z, 5*recvCount_Z*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xy, recvCount_xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xY, recvCount_xY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Xy, recvCount_Xy*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_XY, recvCount_XY*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xz, recvCount_xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_xZ, recvCount_xZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Xz, recvCount_Xz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_XZ, recvCount_XZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_yz, recvCount_yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_yZ, recvCount_yZ*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_Yz, recvCount_Yz*sizeof(double));	// Allocate device memory
	AllocateDeviceMemory((void **) &recvbuf_YZ, recvCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X, *dvcSendList_Y, *dvcSendList_Z;
	int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy, *dvcSendList_Yz, *dvcSendList_xZ;
	int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY, *dvcSendList_YZ, *dvcSendList_XZ;
	//......................................................................................
	int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X, *dvcRecvList_Y, *dvcRecvList_Z;
	int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy, *dvcRecvList_Yz, *dvcRecvList_xZ;
	int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY, *dvcRecvList_YZ, *dvcRecvList_XZ;
	//......................................................................................
	AllocateDeviceMemory((void **) &dvcSendList_x, sendCount_x*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_X, sendCount_X*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_y, sendCount_y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Y, sendCount_Y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_z, sendCount_z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Z, sendCount_Z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xy, sendCount_xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xY, sendCount_xY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Xy, sendCount_Xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_XY, sendCount_XY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xz, sendCount_xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_xZ, sendCount_xZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Xz, sendCount_Xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_XZ, sendCount_XZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_yz, sendCount_yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_yZ, sendCount_yZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_Yz, sendCount_Yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcSendList_YZ, sendCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	AllocateDeviceMemory((void **) &dvcRecvList_x, recvCount_x*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_X, recvCount_X*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_y, recvCount_y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Y, recvCount_Y*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_z, recvCount_z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Z, recvCount_Z*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	AllocateDeviceMemory((void **) &dvcRecvList_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank==0)	printf ("Prepare to copy send/recv Lists to device \n");
	CopyToDevice(dvcSendList_x,sendList_x,sendCount_x*sizeof(int));
	CopyToDevice(dvcSendList_X,sendList_X,sendCount_X*sizeof(int));
	CopyToDevice(dvcSendList_y,sendList_y,sendCount_y*sizeof(int));
	CopyToDevice(dvcSendList_Y,sendList_Y,sendCount_Y*sizeof(int));
	CopyToDevice(dvcSendList_z,sendList_z,sendCount_z*sizeof(int));
	CopyToDevice(dvcSendList_Z,sendList_Z,sendCount_Z*sizeof(int));
	CopyToDevice(dvcSendList_xy,sendList_xy,sendCount_xy*sizeof(int));
	CopyToDevice(dvcSendList_XY,sendList_XY,sendCount_XY*sizeof(int));
	CopyToDevice(dvcSendList_xY,sendList_xY,sendCount_xY*sizeof(int));
	CopyToDevice(dvcSendList_Xy,sendList_Xy,sendCount_Xy*sizeof(int));
	CopyToDevice(dvcSendList_xz,sendList_xz,sendCount_xz*sizeof(int));
	CopyToDevice(dvcSendList_XZ,sendList_XZ,sendCount_XZ*sizeof(int));
	CopyToDevice(dvcSendList_xZ,sendList_xZ,sendCount_xZ*sizeof(int));
	CopyToDevice(dvcSendList_Xz,sendList_Xz,sendCount_Xz*sizeof(int));
	CopyToDevice(dvcSendList_yz,sendList_yz,sendCount_yz*sizeof(int));
	CopyToDevice(dvcSendList_YZ,sendList_YZ,sendCount_YZ*sizeof(int));
	CopyToDevice(dvcSendList_yZ,sendList_yZ,sendCount_yZ*sizeof(int));
	CopyToDevice(dvcSendList_Yz,sendList_Yz,sendCount_Yz*sizeof(int));
	//......................................................................................
	CopyToDevice(dvcRecvList_x,recvList_x,recvCount_x*sizeof(int));
	CopyToDevice(dvcRecvList_X,recvList_X,recvCount_X*sizeof(int));
	CopyToDevice(dvcRecvList_y,recvList_y,recvCount_y*sizeof(int));
	CopyToDevice(dvcRecvList_Y,recvList_Y,recvCount_Y*sizeof(int));
	CopyToDevice(dvcRecvList_z,recvList_z,recvCount_z*sizeof(int));
	CopyToDevice(dvcRecvList_Z,recvList_Z,recvCount_Z*sizeof(int));
	CopyToDevice(dvcRecvList_xy,recvList_xy,recvCount_xy*sizeof(int));
	CopyToDevice(dvcRecvList_XY,recvList_XY,recvCount_XY*sizeof(int));
	CopyToDevice(dvcRecvList_xY,recvList_xY,recvCount_xY*sizeof(int));
	CopyToDevice(dvcRecvList_Xy,recvList_Xy,recvCount_Xy*sizeof(int));
	CopyToDevice(dvcRecvList_xz,recvList_xz,recvCount_xz*sizeof(int));
	CopyToDevice(dvcRecvList_XZ,recvList_XZ,recvCount_XZ*sizeof(int));
	CopyToDevice(dvcRecvList_xZ,recvList_xZ,recvCount_xZ*sizeof(int));
	CopyToDevice(dvcRecvList_Xz,recvList_Xz,recvCount_Xz*sizeof(int));
	CopyToDevice(dvcRecvList_yz,recvList_yz,recvCount_yz*sizeof(int));
	CopyToDevice(dvcRecvList_YZ,recvList_YZ,recvCount_YZ*sizeof(int));
	CopyToDevice(dvcRecvList_yZ,recvList_yZ,recvCount_yZ*sizeof(int));
	CopyToDevice(dvcRecvList_Yz,recvList_Yz,recvCount_Yz*sizeof(int));
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
			recvID_X,recvCount_X,MPI_CHAR,rank_X,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_X,sendCount_X,MPI_CHAR,rank_X,sendtag,
			recvID_x,recvCount_x,MPI_CHAR,rank_x,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_y,sendCount_y,MPI_CHAR,rank_y,sendtag,
			recvID_Y,recvCount_Y,MPI_CHAR,rank_Y,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Y,sendCount_Y,MPI_CHAR,rank_Y,sendtag,
			recvID_y,recvCount_y,MPI_CHAR,rank_y,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_z,sendCount_z,MPI_CHAR,rank_z,sendtag,
			recvID_Z,recvCount_Z,MPI_CHAR,rank_Z,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Z,sendCount_Z,MPI_CHAR,rank_Z,sendtag,
			recvID_z,recvCount_z,MPI_CHAR,rank_z,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xy,sendCount_xy,MPI_CHAR,rank_xy,sendtag,
			recvID_XY,recvCount_XY,MPI_CHAR,rank_XY,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XY,sendCount_XY,MPI_CHAR,rank_XY,sendtag,
			recvID_xy,recvCount_xy,MPI_CHAR,rank_xy,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xy,sendCount_Xy,MPI_CHAR,rank_Xy,sendtag,
			recvID_xY,recvCount_xY,MPI_CHAR,rank_xY,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xY,sendCount_xY,MPI_CHAR,rank_xY,sendtag,
			recvID_Xy,recvCount_Xy,MPI_CHAR,rank_Xy,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xz,sendCount_xz,MPI_CHAR,rank_xz,sendtag,
			recvID_XZ,recvCount_XZ,MPI_CHAR,rank_XZ,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XZ,sendCount_XZ,MPI_CHAR,rank_XZ,sendtag,
			recvID_xz,recvCount_xz,MPI_CHAR,rank_xz,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xz,sendCount_Xz,MPI_CHAR,rank_Xz,sendtag,
			recvID_xZ,recvCount_xZ,MPI_CHAR,rank_xZ,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xZ,sendCount_xZ,MPI_CHAR,rank_xZ,sendtag,
			recvID_Xz,recvCount_Xz,MPI_CHAR,rank_Xz,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yz,sendCount_yz,MPI_CHAR,rank_yz,sendtag,
			recvID_YZ,recvCount_YZ,MPI_CHAR,rank_YZ,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_YZ,sendCount_YZ,MPI_CHAR,rank_YZ,sendtag,
			recvID_yz,recvCount_yz,MPI_CHAR,rank_yz,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Yz,sendCount_Yz,MPI_CHAR,rank_Yz,sendtag,
			recvID_yZ,recvCount_yZ,MPI_CHAR,rank_yZ,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yZ,sendCount_yZ,MPI_CHAR,rank_yZ,sendtag,
			recvID_Yz,recvCount_Yz,MPI_CHAR,rank_Yz,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
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
	MPI_Barrier(MPI_COMM_WORLD);

	//...........device phase ID.................................................
	if (rank==0)	printf ("Copying phase ID to device \n");
	char *ID;
	AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
	// Copy to the device
	CopyToDevice(ID, id, N);
	//...........................................................................

	//...........................................................................
	//				MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	// LBM variables
	if (rank==0)	printf ("Allocating distributions \n");
	//......................device distributions.................................
	double *f_even,*f_odd;
	double *A_even,*A_odd,*B_even,*B_odd;
	double *f_even_host,*f_odd_host;
	//...........................................................................
	AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
	AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory

	// Write the communcation structure into a file for debugging
	char LocalCommFile[40];
	sprintf(LocalCommFile,"%s%s","Comm.",LocalRankString);
	FILE *CommFile;
	CommFile = fopen(LocalCommFile,"w");
	fprintf(CommFile,"rank=%d, ",rank);
	fprintf(CommFile,"i=%d,j=%d,k=%d :",iproc,jproc,kproc);
	fprintf(CommFile,"x=%d, ",rank_x);
	fprintf(CommFile,"X=%d, ",rank_X);
	fprintf(CommFile,"y=%d, ",rank_y);
	fprintf(CommFile,"Y=%d, ",rank_Y);
	fprintf(CommFile,"z=%d, ",rank_z);
	fprintf(CommFile,"Z=%d, ",rank_Z);
	fprintf(CommFile,"xy=%d, ",rank_xy);
	fprintf(CommFile,"XY=%d, ",rank_XY);
	fprintf(CommFile,"xY=%d, ",rank_xY);
	fprintf(CommFile,"Xy=%d, ",rank_Xy);
	fprintf(CommFile,"xz=%d, ",rank_xz);
	fprintf(CommFile,"XZ=%d, ",rank_XZ);
	fprintf(CommFile,"xZ=%d, ",rank_xZ);
	fprintf(CommFile,"Xz=%d, ",rank_Xz);
	fprintf(CommFile,"yz=%d, ",rank_yz);
	fprintf(CommFile,"YZ=%d, ",rank_YZ);
	fprintf(CommFile,"yZ=%d, ",rank_yZ);
	fprintf(CommFile,"Yz=%d, ",rank_Yz);
	fprintf(CommFile,"\n");
	fclose(CommFile);

	if (rank==0)	printf("Setting the distributions, size = : %i\n", N);
	//...........................................................................
	GlobalFlipInitD3Q19(f_even_host, f_odd_host, Nx-2, Ny-2, Nz-2,iproc,jproc,kproc,nprocx,nprocy,nprocz);
	CopyToDevice(f_even, f_even_host, 10*dist_mem_size);
	CopyToDevice(f_odd, f_odd_host, 9*dist_mem_size);
	//...........................................................................

	//...................................................................................
	PackDist(1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
	PackDist(4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist(5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist(6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	PackDist(7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
	//...Packing for X face(1,7,9,11,13)................................
	PackDist(0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	PackDist(6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
	//...Packing for y face(4,8,9,16,18).................................
	PackDist(2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
	PackDist(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	PackDist(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
	PackDist(8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	PackDist(9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
	//...Packing for Y face(3,7,10,15,17).................................
	PackDist(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist(5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
	PackDist(7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	PackDist(8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
	//...Packing for z face(6,12,13,16,17)................................
	PackDist(3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
	PackDist(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	PackDist(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	PackDist(8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
	PackDist(8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
	//...Packing for Z face(5,11,14,15,18)................................
	PackDist(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist(7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	PackDist(7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
	PackDist(9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
	//...Pack the xy edge (8)................................
	PackDist(4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
	//...Pack the Xy edge (9)................................
	PackDist(4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
	//...Pack the xY edge (10)................................
	PackDist(5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
	//...Pack the XY edge (7)................................
	PackDist(3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
	//...Pack the xz edge (12)................................
	PackDist(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
	//...Pack the xZ edge (14)................................
	PackDist(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
	//...Pack the Xz edge (13)................................
	PackDist(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
	//...Pack the XZ edge (11)................................
	PackDist(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
	//...Pack the xz edge (12)................................
	//...Pack the yz edge (16)................................
	PackDist(8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
	//...Pack the yZ edge (18)................................
	PackDist(9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
	//...Pack the Yz edge (17)................................
	PackDist(8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
	//...Pack the YZ edge (15)................................
	PackDist(7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
	//...................................................................................

	//...................................................................................
	// Send all the distributions
	MPI_Isend(sendbuf_x, 5*sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_WORLD,&req1[0]);
	MPI_Irecv(recvbuf_X, 5*recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_WORLD,&req2[0]);
	MPI_Isend(sendbuf_X, 5*sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_WORLD,&req1[1]);
	MPI_Irecv(recvbuf_x, 5*recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_WORLD,&req2[1]);
	MPI_Isend(sendbuf_y, 5*sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_WORLD,&req1[2]);
	MPI_Irecv(recvbuf_Y, 5*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_WORLD,&req2[2]);
	MPI_Isend(sendbuf_Y, 5*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_WORLD,&req1[3]);
	MPI_Irecv(recvbuf_y, 5*recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_WORLD,&req2[3]);
	MPI_Isend(sendbuf_z, 5*sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_WORLD,&req1[4]);
	MPI_Irecv(recvbuf_Z, 5*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_WORLD,&req2[4]);
	MPI_Isend(sendbuf_Z, 5*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_WORLD,&req1[5]);
	MPI_Irecv(recvbuf_z, 5*recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_WORLD,&req2[5]);
	MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,MPI_COMM_WORLD,&req1[6]);
	MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,MPI_COMM_WORLD,&req2[6]);
	MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,MPI_COMM_WORLD,&req1[7]);
	MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,MPI_COMM_WORLD,&req2[7]);
	MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,MPI_COMM_WORLD,&req1[8]);
	MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,MPI_COMM_WORLD,&req2[8]);
	MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,MPI_COMM_WORLD,&req1[9]);
	MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,MPI_COMM_WORLD,&req2[9]);
	MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,MPI_COMM_WORLD,&req1[10]);
	MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,MPI_COMM_WORLD,&req2[10]);
	MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,MPI_COMM_WORLD,&req1[11]);
	MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,MPI_COMM_WORLD,&req2[11]);
	MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,MPI_COMM_WORLD,&req1[12]);
	MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,MPI_COMM_WORLD,&req2[12]);
	MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,MPI_COMM_WORLD,&req1[13]);
	MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,MPI_COMM_WORLD,&req2[13]);
	MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,MPI_COMM_WORLD,&req1[14]);
	MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,MPI_COMM_WORLD,&req2[14]);
	MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,MPI_COMM_WORLD,&req1[15]);
	MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,MPI_COMM_WORLD,&req2[15]);
	MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,MPI_COMM_WORLD,&req1[16]);
	MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,MPI_COMM_WORLD,&req2[16]);
	MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,MPI_COMM_WORLD,&req1[17]);
	MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,MPI_COMM_WORLD,&req2[17]);
	//...................................................................................

	//*************************************************************************
	// 		Swap the distributions for momentum transport
	//*************************************************************************
	SwapD3Q19(ID, f_even, f_odd, Nx, Ny, Nz);
	//*************************************************************************

	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);

	//...................................................................................
	// Unpack the distributions on the device
	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,13 .................................
	UnpackDist(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(3,-1,-1,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(4,-1,1,0,dvcRecvList_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(5,-1,0,-1,dvcRecvList_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	UnpackDist(6,-1,0,1,dvcRecvList_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	UnpackDist(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(4,1,1,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(5,1,-1,0,dvcRecvList_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(6,1,0,1,dvcRecvList_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	UnpackDist(7,1,0,-1,dvcRecvList_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
	UnpackDist(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	UnpackDist(3,-1,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	UnpackDist(5,1,-1,0,dvcRecvList_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,Nx,Ny,Nz);
	UnpackDist(7,0,-1,-1,dvcRecvList_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	UnpackDist(8,0,-1,1,dvcRecvList_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
	UnpackDist(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	UnpackDist(4,1,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	UnpackDist(4,-1,1,0,dvcRecvList_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,Nx,Ny,Nz);
	UnpackDist(8,0,1,1,dvcRecvList_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	UnpackDist(9,0,1,-1,dvcRecvList_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
	UnpackDist(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
	UnpackDist(5,-1,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
	UnpackDist(7,1,0,-1,dvcRecvList_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
	UnpackDist(7,0,-1,-1,dvcRecvList_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
	UnpackDist(9,0,1,-1,dvcRecvList_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
	//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
	UnpackDist(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
	UnpackDist(6,1,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
	UnpackDist(6,-1,0,1,dvcRecvList_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
	UnpackDist(8,0,1,1,dvcRecvList_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
	UnpackDist(8,0,-1,1,dvcRecvList_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
	//..................................................................................
	//...Map recieve list for the xy edge <<<8)................................
	UnpackDist(3,-1,-1,0,dvcRecvList_XY,0,recvCount_XY,recvbuf_XY,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the Xy edge <<<9)................................
	UnpackDist(5,1,-1,0,dvcRecvList_xY,0,recvCount_xY,recvbuf_xY,f_even,Nx,Ny,Nz);
	//...Map recieve list for the xY edge <<<10)................................
	UnpackDist(4,-1,1,0,dvcRecvList_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the XY edge <<<7)................................
	UnpackDist(4,1,1,0,dvcRecvList_xy,0,recvCount_xy,recvbuf_xy,f_even,Nx,Ny,Nz);
	//...Map recieve list for the xz edge <<<12)................................
	UnpackDist(5,-1,0,-1,dvcRecvList_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the xZ edge <<<14)................................
	UnpackDist(6,-1,0,1,dvcRecvList_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the Xz edge <<<13)................................
	UnpackDist(7,1,0,-1,dvcRecvList_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,Nx,Ny,Nz);
	//...Map recieve list for the XZ edge <<<11)................................
	UnpackDist(6,1,0,1,dvcRecvList_xz,0,recvCount_xz,recvbuf_xz,f_even,Nx,Ny,Nz);
	//...Map recieve list for the yz edge <<<16)................................
	UnpackDist(7,0,-1,-1,dvcRecvList_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the yZ edge <<<18)................................
	UnpackDist(8,0,-1,1,dvcRecvList_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,Nx,Ny,Nz);
	//...Map recieve list for the Yz edge <<<17)................................
	UnpackDist(9,0,1,-1,dvcRecvList_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,Nx,Ny,Nz);
	//...Map recieve list for the YZ edge <<<15)................................
	UnpackDist(8,0,1,1,dvcRecvList_yz,0,recvCount_yz,recvbuf_yz,f_even,Nx,Ny,Nz);
	//...................................................................................

	//...........................................................................	

	int check;
	CopyToHost(f_even_host,f_even,10*N*sizeof(double));
	CopyToHost(f_odd_host,f_even,10*N*sizeof(double));
	check =	GlobalCheckDebugDist(f_even_host, f_odd_host, Nx-2, Ny-2, Nz-2,iproc,jproc,kproc,nprocx,nprocy,nprocz);
	//...........................................................................

	int timestep = 0;
	if (rank==0) printf("********************************************************\n");
	if (rank==0)	printf("No. of timesteps for timing: %i \n", 100);

	//.......create and start timer............
	double starttime,stoptime,cputime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();
	//.........................................

	sendtag = recvtag = 5;

	//************ MAIN ITERATION LOOP (timing communications)***************************************/
	while (timestep < 100){

		//...................................................................................
		PackDist(1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
		PackDist(4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist(5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist(6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist(7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		//...Packing for X face(1,7,9,11,13)................................
		PackDist(0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		//...Packing for y face(4,8,9,16,18).................................
		PackDist(2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
		PackDist(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		PackDist(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
		PackDist(8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		PackDist(9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		//...Packing for Y face(3,7,10,15,17).................................
		PackDist(1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist(5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
		PackDist(7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist(8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		//...Packing for z face(6,12,13,16,17)................................
		PackDist(3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
		PackDist(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		PackDist(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		PackDist(8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		PackDist(8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		//...Packing for Z face(5,11,14,15,18)................................
		PackDist(2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist(7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		PackDist(7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist(9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		//...Pack the xy edge (8)................................
		PackDist(4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
		//...Pack the Xy edge (9)................................
		PackDist(4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
		//...Pack the xY edge (10)................................
		PackDist(5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
		//...Pack the XY edge (7)................................
		PackDist(3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
		//...Pack the xz edge (12)................................
		PackDist(6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
		//...Pack the xZ edge (14)................................
		PackDist(7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
		//...Pack the Xz edge (13)................................
		PackDist(6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
		//...Pack the XZ edge (11)................................
		PackDist(5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
		//...Pack the xz edge (12)................................
		//...Pack the yz edge (16)................................
		PackDist(8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
		//...Pack the yZ edge (18)................................
		PackDist(9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
		//...Pack the Yz edge (17)................................
		PackDist(8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
		//...Pack the YZ edge (15)................................
		PackDist(7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
		//...................................................................................

		//...................................................................................
		// Send all the distributions
		MPI_Isend(sendbuf_x, 5*sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_WORLD,&req1[0]);
		MPI_Irecv(recvbuf_X, 5*recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_WORLD,&req2[0]);
		MPI_Isend(sendbuf_X, 5*sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_WORLD,&req1[1]);
		MPI_Irecv(recvbuf_x, 5*recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_WORLD,&req2[1]);
		MPI_Isend(sendbuf_y, 5*sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_WORLD,&req1[2]);
		MPI_Irecv(recvbuf_Y, 5*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_WORLD,&req2[2]);
		MPI_Isend(sendbuf_Y, 5*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_WORLD,&req1[3]);
		MPI_Irecv(recvbuf_y, 5*recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_WORLD,&req2[3]);
		MPI_Isend(sendbuf_z, 5*sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_WORLD,&req1[4]);
		MPI_Irecv(recvbuf_Z, 5*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_WORLD,&req2[4]);
		MPI_Isend(sendbuf_Z, 5*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_WORLD,&req1[5]);
		MPI_Irecv(recvbuf_z, 5*recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_WORLD,&req2[5]);
		MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,MPI_COMM_WORLD,&req1[6]);
		MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,MPI_COMM_WORLD,&req2[6]);
		MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,MPI_COMM_WORLD,&req1[7]);
		MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,MPI_COMM_WORLD,&req2[7]);
		MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,MPI_COMM_WORLD,&req1[8]);
		MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,MPI_COMM_WORLD,&req2[8]);
		MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,MPI_COMM_WORLD,&req1[9]);
		MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,MPI_COMM_WORLD,&req2[9]);
		MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,MPI_COMM_WORLD,&req1[10]);
		MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,MPI_COMM_WORLD,&req2[10]);
		MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,MPI_COMM_WORLD,&req1[11]);
		MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,MPI_COMM_WORLD,&req2[11]);
		MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,MPI_COMM_WORLD,&req1[12]);
		MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,MPI_COMM_WORLD,&req2[12]);
		MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,MPI_COMM_WORLD,&req1[13]);
		MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,MPI_COMM_WORLD,&req2[13]);
		MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,MPI_COMM_WORLD,&req1[14]);
		MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,MPI_COMM_WORLD,&req2[14]);
		MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,MPI_COMM_WORLD,&req1[15]);
		MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,MPI_COMM_WORLD,&req2[15]);
		MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,MPI_COMM_WORLD,&req1[16]);
		MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,MPI_COMM_WORLD,&req2[16]);
		MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,MPI_COMM_WORLD,&req1[17]);
		MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,MPI_COMM_WORLD,&req2[17]);
		//...................................................................................

		//*************************************************************************
		// 		Swap the distributions for momentum transport
		//*************************************************************************
		SwapD3Q19(ID, f_even, f_odd, Nx, Ny, Nz);
		//*************************************************************************

		//...................................................................................
		// Wait for completion of D3Q19 communication
		MPI_Waitall(18,req1,stat1);
		MPI_Waitall(18,req2,stat2);

		//...................................................................................
		// Unpack the distributions on the device
		//...................................................................................
		//...Map recieve list for the X face: q=2,8,10,12,13 .................................
		UnpackDist(0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		UnpackDist(3,-1,-1,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		UnpackDist(4,-1,1,0,dvcRecvList_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		UnpackDist(5,-1,0,-1,dvcRecvList_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		UnpackDist(6,-1,0,1,dvcRecvList_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the x face: q=1,7,9,11,13..................................
		UnpackDist(1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		UnpackDist(4,1,1,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		UnpackDist(5,1,-1,0,dvcRecvList_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		UnpackDist(6,1,0,1,dvcRecvList_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		UnpackDist(7,1,0,-1,dvcRecvList_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
		UnpackDist(1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		UnpackDist(3,-1,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		UnpackDist(5,1,-1,0,dvcRecvList_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,Nx,Ny,Nz);
		UnpackDist(7,0,-1,-1,dvcRecvList_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		UnpackDist(8,0,-1,1,dvcRecvList_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
		UnpackDist(2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		UnpackDist(4,1,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		UnpackDist(4,-1,1,0,dvcRecvList_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,Nx,Ny,Nz);
		UnpackDist(8,0,1,1,dvcRecvList_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		UnpackDist(9,0,1,-1,dvcRecvList_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
		UnpackDist(2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		UnpackDist(5,-1,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		UnpackDist(7,1,0,-1,dvcRecvList_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		UnpackDist(7,0,-1,-1,dvcRecvList_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		UnpackDist(9,0,1,-1,dvcRecvList_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
		UnpackDist(3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		UnpackDist(6,1,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		UnpackDist(6,-1,0,1,dvcRecvList_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		UnpackDist(8,0,1,1,dvcRecvList_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		UnpackDist(8,0,-1,1,dvcRecvList_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		//..................................................................................
		//...Map recieve list for the xy edge <<<8)................................
		UnpackDist(3,-1,-1,0,dvcRecvList_XY,0,recvCount_XY,recvbuf_XY,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xy edge <<<9)................................
		UnpackDist(5,1,-1,0,dvcRecvList_xY,0,recvCount_xY,recvbuf_xY,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xY edge <<<10)................................
		UnpackDist(4,-1,1,0,dvcRecvList_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the XY edge <<<7)................................
		UnpackDist(4,1,1,0,dvcRecvList_xy,0,recvCount_xy,recvbuf_xy,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xz edge <<<12)................................
		UnpackDist(5,-1,0,-1,dvcRecvList_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the xZ edge <<<14)................................
		UnpackDist(6,-1,0,1,dvcRecvList_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xz edge <<<13)................................
		UnpackDist(7,1,0,-1,dvcRecvList_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the XZ edge <<<11)................................
		UnpackDist(6,1,0,1,dvcRecvList_xz,0,recvCount_xz,recvbuf_xz,f_even,Nx,Ny,Nz);
		//...Map recieve list for the yz edge <<<16)................................
		UnpackDist(7,0,-1,-1,dvcRecvList_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the yZ edge <<<18)................................
		UnpackDist(8,0,-1,1,dvcRecvList_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Yz edge <<<17)................................
		UnpackDist(9,0,1,-1,dvcRecvList_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the YZ edge <<<15)................................
		UnpackDist(8,0,1,1,dvcRecvList_yz,0,recvCount_yz,recvbuf_yz,f_even,Nx,Ny,Nz);
		//...................................................................................

		MPI_Barrier(MPI_COMM_WORLD);
		// Iteration completed!
		timestep++;
		//...................................................................
	}
	//************************************************************************/

	stoptime = MPI_Wtime();
	//	cout << "CPU time: " << (stoptime - starttime) << " seconds" << endl;
	cputime = stoptime - starttime;
	//	cout << "Lattice update rate: "<< double(Nx*Ny*Nz*timestep)/cputime/1000000 <<  " MLUPS" << endl;
	double MLUPS = double(Nx*Ny*Nz*timestep)/cputime/1000000;
	if (rank==0) printf("********************************************************\n");
	if (rank==0) printf("CPU time = %f \n", cputime);
	if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
	MLUPS *= nprocs;
	if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
	if (rank==0) printf("********************************************************\n");

	// ****************************************************
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	// ****************************************************

	return check;
}
