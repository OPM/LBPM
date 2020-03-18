#include <stdio.h>
#include <iostream>
#include <fstream>
#include "D3Q19.h"
#include "D3Q7.h"
#include "Color.h"
#include "common/MPI_Helpers.h"

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
		printf("Running Hybrid Implementation of Color LBM	\n");
		printf("********************************************************\n");
	}
	// Color Model parameters
	string FILENAME;
	unsigned int nBlocks, nthreads;
	int Nx,Ny,Nz;
	int timestepMax, interval;
	double tau,Fx,Fy,Fz,tol;
	double alpha, beta;
	double das, dbs;
	double din,dout;
	bool pBC;
	int i,j,k,n;

	if (rank==0){
		//.............................................................
		//		READ SIMULATION PARMAETERS FROM INPUT FILE
		//.............................................................
		ifstream input("Color.in");
		// Line 1: Name of the phase indicator file (s=0,w=1,n=2)
		input >> FILENAME;
		// Line 2: domain size (Nx, Ny, Nz)
		input >> Nz;				// number of nodes (x,y,z)
		input >> nBlocks;
		input >> nthreads;
		// Line 3: model parameters (tau, alpha, beta, das, dbs)
		input >> tau;
		input >> alpha;
		input >> beta;
		input >> das;
		input >> dbs;
		// Line 4: External force components (Fx,Fy, Fz)
		input >> Fx;
		input >> Fy;
		input >> Fz;
		// Line 5: Pressure Boundary conditions
		input >> pBC;
		input >> din;
		input >> dout;
		// Line 6: time-stepping criteria
		input >> timestepMax;			// max no. of timesteps
		input >> interval;			// error interval
		input >> tol;				// error tolerance
		//.............................................................

		ifstream domain("Domain.in");
		domain >> nprocx;
		domain >> nprocy;
		domain >> nprocz;
	}
	// **************************************************************
	// Broadcast simulation parameters from rank 0 to all other procs
	MPI_Barrier(comm);
	//.................................................
	MPI_Bcast(&Nz,1,MPI_INT,0,comm);
	MPI_Bcast(&nBlocks,1,MPI_INT,0,comm);
	MPI_Bcast(&nthreads,1,MPI_INT,0,comm);
	MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&tau,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&alpha,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&beta,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&das,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&dbs,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&pBC,1,MPI_LOGICAL,0,comm);
	MPI_Bcast(&din,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&dout,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&timestepMax,1,MPI_INT,0,comm);
	MPI_Bcast(&interval,1,MPI_INT,0,comm);
	MPI_Bcast(&tol,1,MPI_DOUBLE,0,comm);

	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	//.................................................
	MPI_Barrier(comm);
	// **************************************************************
	// **************************************************************

	double rlxA = 1.f/tau;
	double rlxB = 8.f*(2.f-rlxA)/(8.f-rlxA);

	if (nprocs != nprocx*nprocy*nprocz){
		printf("Fatal error in processor number! \n");
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
	}

	if (rank==0){
		printf("********************************************************\n");
		printf("tau = %f \n", tau);
		printf("alpha = %f \n", alpha);
		printf("beta = %f \n", beta);
		printf("das = %f \n", das);
		printf("dbs = %f \n", dbs);
		printf("Force(x) = %f \n", Fx);
		printf("Force(y) = %f \n", Fy);
		printf("Force(z) = %f \n", Fz);
		printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
		printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
		printf("********************************************************\n");

	}

	MPI_Barrier(comm);
	kproc = rank/(nprocx*nprocy);
	jproc = (rank-nprocx*nprocy*kproc)/nprocx;
	iproc = rank-nprocx*nprocy*kproc-nprocz*jproc;

	//..........................................
	// set up the neighbor ranks
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=0;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_X = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=0;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_x = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Y = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j-=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_y = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=0;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-= nprocy;
	if (!(k<nprocz)) k-= nprocz;
	rank_Z = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=0;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-= nprocy;
	if (!(k<nprocz)) k-= nprocz;
	rank_z = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_XY = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j-=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xy = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j-=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Xy = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=1;
	k+=0;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xY = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=0;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_XZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=0;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i-=1;
	j+=0;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_xZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=1;
	j+=0;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Xz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j+=1;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_YZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j-=1;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_yz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k =kproc;
	i+=0;
	j-=1;
	k+=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_yZ = k*nprocx*nprocy+j*nprocx+i;
	//..........................................
	i=iproc; j=jproc; k=kproc;
	i+=0;
	j+=1;
	k-=1;
	if (i<0)	i+=nprocx;
	if (j<0)	j+=nprocy;
	if (k<0)	k+=nprocz;
	if (!(i<nprocx)) i-= nprocx;
	if (!(j<nprocy)) j-=nprocy;
	if (!(k<nprocz)) k-=nprocz;
	rank_Yz = k*nprocx*nprocy+j*nprocx+i;
	//..........................................

	Nz += 2;
	Nx = Ny = Nz;	// Cubic domain

	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

//	unsigned int nBlocks = 32;
//	int nthreads = 128;
	int S = N/nthreads/nBlocks;

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
	char LocalRankString[8];
	char LocalRankFilename[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
//	printf("Local File Name =  %s \n",LocalRankFilename);
	// .......... READ THE INPUT FILE .......................................
	char value;
	char *id;
	id = new char[N];
	int sum = 0;
//	double porosity;
	//.......................................................................
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
	//...........................................................................
	MPI_Barrier(comm);
	if (rank == 0) cout << "Domain set." << endl;
	//...........................................................................
	// Write the communcation structure into a file for debugging
/*	char LocalCommFile[40];
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
*/	//...........................................................................

	// Set up MPI communication structures
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
	MPI_Send(&sendCount_x,1,MPI_INT,rank_X,sendtag,comm);
	MPI_Recv(&recvCount_X,1,MPI_INT,rank_x,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_X,1,MPI_INT,rank_x,sendtag,comm);
	MPI_Recv(&recvCount_x,1,MPI_INT,rank_X,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_y,1,MPI_INT,rank_Y,sendtag,comm);
	MPI_Recv(&recvCount_Y,1,MPI_INT,rank_y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Y,1,MPI_INT,rank_y,sendtag,comm);
	MPI_Recv(&recvCount_y,1,MPI_INT,rank_Y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_z,1,MPI_INT,rank_Z,sendtag,comm);
	MPI_Recv(&recvCount_Z,1,MPI_INT,rank_z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Z,1,MPI_INT,rank_z,sendtag,comm);
	MPI_Recv(&recvCount_z,1,MPI_INT,rank_Z,recvtag,comm,MPI_STATUS_IGNORE);

	MPI_Send(&sendCount_xy,1,MPI_INT,rank_XY,sendtag,comm);
	MPI_Recv(&recvCount_XY,1,MPI_INT,rank_xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_XY,1,MPI_INT,rank_xy,sendtag,comm);
	MPI_Recv(&recvCount_xy,1,MPI_INT,rank_XY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Xy,1,MPI_INT,rank_xY,sendtag,comm);
	MPI_Recv(&recvCount_xY,1,MPI_INT,rank_Xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_xY,1,MPI_INT,rank_Xy,sendtag,comm);
	MPI_Recv(&recvCount_Xy,1,MPI_INT,rank_xY,recvtag,comm,MPI_STATUS_IGNORE);

	MPI_Send(&sendCount_xz,1,MPI_INT,rank_XZ,sendtag,comm);
	MPI_Recv(&recvCount_XZ,1,MPI_INT,rank_xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_XZ,1,MPI_INT,rank_xz,sendtag,comm);
	MPI_Recv(&recvCount_xz,1,MPI_INT,rank_XZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Xz,1,MPI_INT,rank_xZ,sendtag,comm);
	MPI_Recv(&recvCount_xZ,1,MPI_INT,rank_Xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_xZ,1,MPI_INT,rank_Xz,sendtag,comm);
	MPI_Recv(&recvCount_Xz,1,MPI_INT,rank_xZ,recvtag,comm,MPI_STATUS_IGNORE);

	MPI_Send(&sendCount_yz,1,MPI_INT,rank_YZ,sendtag,comm);
	MPI_Recv(&recvCount_YZ,1,MPI_INT,rank_yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_YZ,1,MPI_INT,rank_yz,sendtag,comm);
	MPI_Recv(&recvCount_yz,1,MPI_INT,rank_YZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_Yz,1,MPI_INT,rank_yZ,sendtag,comm);
	MPI_Recv(&recvCount_yZ,1,MPI_INT,rank_Yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Send(&sendCount_yZ,1,MPI_INT,rank_Yz,sendtag,comm);
	MPI_Recv(&recvCount_Yz,1,MPI_INT,rank_yZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Barrier(comm);
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
	MPI_Isend(sendList_x, sendCount_x,MPI_INT,rank_X,sendtag,comm,&req1[0]);
	MPI_Irecv(recvList_X, recvCount_X,MPI_INT,rank_x,recvtag,comm,&req2[0]);
	MPI_Isend(sendList_X, sendCount_X,MPI_INT,rank_x,sendtag,comm,&req1[1]);
	MPI_Irecv(recvList_x, recvCount_x,MPI_INT,rank_X,recvtag,comm,&req2[1]);
	MPI_Isend(sendList_y, sendCount_y,MPI_INT,rank_Y,sendtag,comm,&req1[2]);
	MPI_Irecv(recvList_Y, recvCount_Y,MPI_INT,rank_y,recvtag,comm,&req2[2]);
	MPI_Isend(sendList_Y, sendCount_Y,MPI_INT,rank_y,sendtag,comm,&req1[3]);
	MPI_Irecv(recvList_y, recvCount_y,MPI_INT,rank_Y,recvtag,comm,&req2[3]);
	MPI_Isend(sendList_z, sendCount_z,MPI_INT,rank_Z,sendtag,comm,&req1[4]);
	MPI_Irecv(recvList_Z, recvCount_Z,MPI_INT,rank_z,recvtag,comm,&req2[4]);
	MPI_Isend(sendList_Z, sendCount_Z,MPI_INT,rank_z,sendtag,comm,&req1[5]);
	MPI_Irecv(recvList_z, recvCount_z,MPI_INT,rank_Z,recvtag,comm,&req2[5]);

	MPI_Isend(sendList_xy, sendCount_xy,MPI_INT,rank_XY,sendtag,comm,&req1[6]);
	MPI_Irecv(recvList_XY, recvCount_XY,MPI_INT,rank_xy,recvtag,comm,&req2[6]);
	MPI_Isend(sendList_XY, sendCount_XY,MPI_INT,rank_xy,sendtag,comm,&req1[7]);
	MPI_Irecv(recvList_xy, recvCount_xy,MPI_INT,rank_XY,recvtag,comm,&req2[7]);
	MPI_Isend(sendList_Xy, sendCount_Xy,MPI_INT,rank_xY,sendtag,comm,&req1[8]);
	MPI_Irecv(recvList_xY, recvCount_xY,MPI_INT,rank_Xy,recvtag,comm,&req2[8]);
	MPI_Isend(sendList_xY, sendCount_xY,MPI_INT,rank_Xy,sendtag,comm,&req1[9]);
	MPI_Irecv(recvList_Xy, recvCount_Xy,MPI_INT,rank_xY,recvtag,comm,&req2[9]);

	MPI_Isend(sendList_xz, sendCount_xz,MPI_INT,rank_XZ,sendtag,comm,&req1[10]);
	MPI_Irecv(recvList_XZ, recvCount_XZ,MPI_INT,rank_xz,recvtag,comm,&req2[10]);
	MPI_Isend(sendList_XZ, sendCount_XZ,MPI_INT,rank_xz,sendtag,comm,&req1[11]);
	MPI_Irecv(recvList_xz, recvCount_xz,MPI_INT,rank_XZ,recvtag,comm,&req2[11]);
	MPI_Isend(sendList_Xz, sendCount_Xz,MPI_INT,rank_xZ,sendtag,comm,&req1[12]);
	MPI_Irecv(recvList_xZ, recvCount_xZ,MPI_INT,rank_Xz,recvtag,comm,&req2[12]);
	MPI_Isend(sendList_xZ, sendCount_xZ,MPI_INT,rank_Xz,sendtag,comm,&req1[13]);
	MPI_Irecv(recvList_Xz, recvCount_Xz,MPI_INT,rank_xZ,recvtag,comm,&req2[13]);

	MPI_Isend(sendList_yz, sendCount_yz,MPI_INT,rank_YZ,sendtag,comm,&req1[14]);
	MPI_Irecv(recvList_YZ, recvCount_YZ,MPI_INT,rank_yz,recvtag,comm,&req2[14]);
	MPI_Isend(sendList_YZ, sendCount_YZ,MPI_INT,rank_yz,sendtag,comm,&req1[15]);
	MPI_Irecv(recvList_yz, recvCount_yz,MPI_INT,rank_YZ,recvtag,comm,&req2[15]);
	MPI_Isend(sendList_Yz, sendCount_Yz,MPI_INT,rank_yZ,sendtag,comm,&req1[16]);
	MPI_Irecv(recvList_yZ, recvCount_yZ,MPI_INT,rank_Yz,recvtag,comm,&req2[16]);
	MPI_Isend(sendList_yZ, sendCount_yZ,MPI_INT,rank_Yz,sendtag,comm,&req1[17]);
	MPI_Irecv(recvList_Yz, recvCount_Yz,MPI_INT,rank_yZ,recvtag,comm,&req2[17]);
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	MPI_Barrier(comm);
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
	sendbuf_x= new double[ 5*sendCount_x];	// Allocate device memory
	sendbuf_X= new double[ 5*sendCount_X];	// Allocate device memory
	sendbuf_y= new double[ 5*sendCount_y];	// Allocate device memory
	sendbuf_Y= new double[ 5*sendCount_Y];	// Allocate device memory
	sendbuf_z= new double[ 5*sendCount_z];	// Allocate device memory
	sendbuf_Z= new double[ 5*sendCount_Z];	// Allocate device memory
	sendbuf_xy= new double[ sendCount_xy];	// Allocate device memory
	sendbuf_xY= new double[ sendCount_xY];	// Allocate device memory
	sendbuf_Xy= new double[ sendCount_Xy];	// Allocate device memory
	sendbuf_XY= new double[ sendCount_XY];	// Allocate device memory
	sendbuf_xz= new double[ sendCount_xz];	// Allocate device memory
	sendbuf_xZ= new double[ sendCount_xZ];	// Allocate device memory
	sendbuf_Xz= new double[ sendCount_Xz];	// Allocate device memory
	sendbuf_XZ= new double[ sendCount_XZ];	// Allocate device memory
	sendbuf_yz= new double[ sendCount_yz];	// Allocate device memory
	sendbuf_yZ= new double[ sendCount_yZ];	// Allocate device memory
	sendbuf_Yz= new double[ sendCount_Yz];	// Allocate device memory
	sendbuf_YZ= new double[ sendCount_YZ];	// Allocate device memory
	//......................................................................................
	recvbuf_x= new double[ 5*recvCount_x];	// Allocate device memory
	recvbuf_X= new double[ 5*recvCount_X];	// Allocate device memory
	recvbuf_y= new double[ 5*recvCount_y];	// Allocate device memory
	recvbuf_Y= new double[ 5*recvCount_Y];	// Allocate device memory
	recvbuf_z= new double[ 5*recvCount_z];	// Allocate device memory
	recvbuf_Z= new double[ 5*recvCount_Z];	// Allocate device memory
	recvbuf_xy= new double[ recvCount_xy];	// Allocate device memory
	recvbuf_xY= new double[ recvCount_xY];	// Allocate device memory
	recvbuf_Xy= new double[ recvCount_Xy];	// Allocate device memory
	recvbuf_XY= new double[ recvCount_XY];	// Allocate device memory
	recvbuf_xz= new double[ recvCount_xz];	// Allocate device memory
	recvbuf_xZ= new double[ recvCount_xZ];	// Allocate device memory
	recvbuf_Xz= new double[ recvCount_Xz];	// Allocate device memory
	recvbuf_XZ= new double[ recvCount_XZ];	// Allocate device memory
	recvbuf_yz= new double[ recvCount_yz];	// Allocate device memory
	recvbuf_yZ= new double[ recvCount_yZ];	// Allocate device memory
	recvbuf_Yz= new double[ recvCount_Yz];	// Allocate device memory
	recvbuf_YZ= new double[ recvCount_YZ];	// Allocate device memory
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
	MPI_Sendrecv(sendID_x,sendCount_x,MPI_CHAR,rank_X,sendtag,
			recvID_X,recvCount_X,MPI_CHAR,rank_x,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_X,sendCount_X,MPI_CHAR,rank_x,sendtag,
			recvID_x,recvCount_x,MPI_CHAR,rank_X,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_y,sendCount_y,MPI_CHAR,rank_Y,sendtag,
			recvID_Y,recvCount_Y,MPI_CHAR,rank_y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Y,sendCount_Y,MPI_CHAR,rank_y,sendtag,
			recvID_y,recvCount_y,MPI_CHAR,rank_Y,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_z,sendCount_z,MPI_CHAR,rank_Z,sendtag,
			recvID_Z,recvCount_Z,MPI_CHAR,rank_z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Z,sendCount_Z,MPI_CHAR,rank_z,sendtag,
			recvID_z,recvCount_z,MPI_CHAR,rank_Z,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xy,sendCount_xy,MPI_CHAR,rank_XY,sendtag,
			recvID_XY,recvCount_XY,MPI_CHAR,rank_xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XY,sendCount_XY,MPI_CHAR,rank_xy,sendtag,
			recvID_xy,recvCount_xy,MPI_CHAR,rank_XY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xy,sendCount_Xy,MPI_CHAR,rank_xY,sendtag,
			recvID_xY,recvCount_xY,MPI_CHAR,rank_Xy,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xY,sendCount_xY,MPI_CHAR,rank_Xy,sendtag,
			recvID_Xy,recvCount_Xy,MPI_CHAR,rank_xY,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xz,sendCount_xz,MPI_CHAR,rank_XZ,sendtag,
			recvID_XZ,recvCount_XZ,MPI_CHAR,rank_xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XZ,sendCount_XZ,MPI_CHAR,rank_xz,sendtag,
			recvID_xz,recvCount_xz,MPI_CHAR,rank_XZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xz,sendCount_Xz,MPI_CHAR,rank_xZ,sendtag,
			recvID_xZ,recvCount_xZ,MPI_CHAR,rank_Xz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xZ,sendCount_xZ,MPI_CHAR,rank_Xz,sendtag,
			recvID_Xz,recvCount_Xz,MPI_CHAR,rank_xZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yz,sendCount_yz,MPI_CHAR,rank_YZ,sendtag,
			recvID_YZ,recvCount_YZ,MPI_CHAR,rank_yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_YZ,sendCount_YZ,MPI_CHAR,rank_yz,sendtag,
			recvID_yz,recvCount_yz,MPI_CHAR,rank_YZ,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Yz,sendCount_Yz,MPI_CHAR,rank_yZ,sendtag,
			recvID_yZ,recvCount_yZ,MPI_CHAR,rank_Yz,recvtag,comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yZ,sendCount_yZ,MPI_CHAR,rank_Yz,sendtag,
			recvID_Yz,recvCount_Yz,MPI_CHAR,rank_yZ,recvtag,comm,MPI_STATUS_IGNORE);
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
	//.....................................................................................
	// Once the ID is saved, free memory allocated to the buffers (no longer needed)
/*	//......................................................................................
	free(sendID_x); free(sendID_X); free(sendID_y); free(sendID_Y); free(sendID_z); free(sendID_Z);
	free(sendID_xy); free(sendID_XY); free(sendID_xY); free(sendID_Xy);
	free(sendID_xz); free(sendID_XZ); free(sendID_xZ); free(sendID_Xz);
	free(sendID_yz); free(sendID_YZ); free(sendID_yZ); free(sendID_Yz);
	free(recvID_x); free(recvID_X); free(recvID_y); free(recvID_Y); free(recvID_z); free(recvID_Z);
	free(recvID_xy); free(recvID_XY); free(recvID_xY); free(recvID_Xy);
	free(recvID_xz); free(recvID_XZ); free(recvID_xZ); free(recvID_Xz);
	free(recvID_yz); free(recvID_YZ); free(recvID_yZ); free(recvID_Yz);
*/	//......................................................................................
	if (rank==0)	printf ("Devices are ready to communicate. \n");
	MPI_Barrier(comm);

	//...........device phase ID.................................................
	if (rank==0)	printf ("Copying phase ID to device \n");
	char *ID;
	ID = new char[N];
	for (n=0; n<N; n++) ID[n] = id[n];
	//..............................................

	if (rank==0)	printf ("Allocating distributions \n");
	//......................device distributions.................................
	double *f_even,*f_odd;
	//...........................................................................
	f_even = new double [10*N];
	f_odd = new double [9*N];
	//...........................................................................

	//...........................................................................
	//				MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	double *Phi,*Den,*Copy;
	double *ColorGrad, *Velocity;
	//...........................................................................
	Phi = new double [N];
	Den = new double [2*N];
	Copy = new double [2*N];
	Velocity = new double [3*N];
	ColorGrad = new double [3*N];
	//...........................................................................

	if (rank==0)	printf("Setting the distributions, size = : %i\n", N);
	//...........................................................................
	InitD3Q19(ID, f_even, f_odd, Nx, Ny, Nz);
	InitDenColor(ID, Den, Phi,  das, dbs, N);
/*	for (n=0; n<N; n++){

		if ( ID[n] == 1){
			Den[2*n] = 1.0;
			Den[2*n+1] = 0.0;
			Phi[n] = 1.0;
		}
		else if ( ID[n] == 2){
			Den[2*n] = 0.0;
			Den[2*n+1] = 1.0;
			Phi[n] = -1.0;
		}
		else{
			Den[2*n] = das;
			Den[2*n+1] = dbs;
			Phi[n] = (das-dbs)/(das+dbs);
		}
	}
*/	//...........................................................................
	//...................................................................................
	//*************************************************************************
	// 		Compute the phase indicator field and reset Copy, Den
	//*************************************************************************
	ComputePhi(ID, Phi, Copy, Den, N);
	//*************************************************************************
	//...................................................................................
	PackValues(sendList_x, sendCount_x,sendbuf_x, Phi, N);
	PackValues(sendList_y, sendCount_y,sendbuf_y, Phi, N);
	PackValues(sendList_z, sendCount_z,sendbuf_z, Phi, N);
	PackValues(sendList_X, sendCount_X,sendbuf_X, Phi, N);
	PackValues(sendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
	PackValues(sendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
	PackValues(sendList_xy, sendCount_xy,sendbuf_xy, Phi, N);
	PackValues(sendList_xY, sendCount_xY,sendbuf_xY, Phi, N);
	PackValues(sendList_Xy, sendCount_Xy,sendbuf_Xy, Phi, N);
	PackValues(sendList_XY, sendCount_XY,sendbuf_XY, Phi, N);
	PackValues(sendList_xz, sendCount_xz,sendbuf_xz, Phi, N);
	PackValues(sendList_xZ, sendCount_xZ,sendbuf_xZ, Phi, N);
	PackValues(sendList_Xz, sendCount_Xz,sendbuf_Xz, Phi, N);
	PackValues(sendList_XZ, sendCount_XZ,sendbuf_XZ, Phi, N);
	PackValues(sendList_yz, sendCount_yz,sendbuf_yz, Phi, N);
	PackValues(sendList_yZ, sendCount_yZ,sendbuf_yZ, Phi, N);
	PackValues(sendList_Yz, sendCount_Yz,sendbuf_Yz, Phi, N);
	PackValues(sendList_YZ, sendCount_YZ,sendbuf_YZ, Phi, N);
	//...................................................................................
	// Send / Recv all the phase indcator field values
	MPI_Isend(sendbuf_x, sendCount_x,MPI_DOUBLE,rank_X,sendtag,comm,&req1[0]);
	MPI_Irecv(recvbuf_X, recvCount_X,MPI_DOUBLE,rank_x,recvtag,comm,&req2[0]);
	MPI_Isend(sendbuf_X, sendCount_X,MPI_DOUBLE,rank_x,sendtag,comm,&req1[1]);
	MPI_Irecv(recvbuf_x, recvCount_x,MPI_DOUBLE,rank_X,recvtag,comm,&req2[1]);
	MPI_Isend(sendbuf_y, sendCount_y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[2]);
	MPI_Irecv(recvbuf_Y, recvCount_Y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[2]);
	MPI_Isend(sendbuf_Y, sendCount_Y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[3]);
	MPI_Irecv(recvbuf_y, recvCount_y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[3]);
	MPI_Isend(sendbuf_z, sendCount_z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[4]);
	MPI_Irecv(recvbuf_Z, recvCount_Z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[4]);
	MPI_Isend(sendbuf_Z, sendCount_Z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[5]);
	MPI_Irecv(recvbuf_z, recvCount_z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[5]);
	MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_XY,sendtag,comm,&req1[6]);
	MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_xy,recvtag,comm,&req2[6]);
	MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_xy,sendtag,comm,&req1[7]);
	MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_XY,recvtag,comm,&req2[7]);
	MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_xY,sendtag,comm,&req1[8]);
	MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_Xy,recvtag,comm,&req2[8]);
	MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_Xy,sendtag,comm,&req1[9]);
	MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_xY,recvtag,comm,&req2[9]);
	MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_XZ,sendtag,comm,&req1[10]);
	MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_xz,recvtag,comm,&req2[10]);
	MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_xz,sendtag,comm,&req1[11]);
	MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_XZ,recvtag,comm,&req2[11]);
	MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_xZ,sendtag,comm,&req1[12]);
	MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_Xz,recvtag,comm,&req2[12]);
	MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_Xz,sendtag,comm,&req1[13]);
	MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_xZ,recvtag,comm,&req2[13]);
	MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_YZ,sendtag,comm,&req1[14]);
	MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_yz,recvtag,comm,&req2[14]);
	MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_yz,sendtag,comm,&req1[15]);
	MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_YZ,recvtag,comm,&req2[15]);
	MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_yZ,sendtag,comm,&req1[16]);
	MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_Yz,recvtag,comm,&req2[16]);
	MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_Yz,sendtag,comm,&req1[17]);
	MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_yZ,recvtag,comm,&req2[17]);
	//...................................................................................
	//...................................................................................
	//...................................................................................
	// Wait for completion of Indicator Field communication
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	//...................................................................................
	//...................................................................................
	UnpackValues(recvList_x, recvCount_x,recvbuf_x, Phi, N);
	UnpackValues(recvList_y, recvCount_y,recvbuf_y, Phi, N);
	UnpackValues(recvList_z, recvCount_z,recvbuf_z, Phi, N);
	UnpackValues(recvList_X, recvCount_X,recvbuf_X, Phi, N);
	UnpackValues(recvList_Y, recvCount_Y,recvbuf_Y, Phi, N);
	UnpackValues(recvList_Z, recvCount_Z,recvbuf_Z, Phi, N);
	UnpackValues(recvList_xy, recvCount_xy,recvbuf_xy, Phi, N);
	UnpackValues(recvList_xY, recvCount_xY,recvbuf_xY, Phi, N);
	UnpackValues(recvList_Xy, recvCount_Xy,recvbuf_Xy, Phi, N);
	UnpackValues(recvList_XY, recvCount_XY,recvbuf_XY, Phi, N);
	UnpackValues(recvList_xz, recvCount_xz,recvbuf_xz, Phi, N);
	UnpackValues(recvList_xZ, recvCount_xZ,recvbuf_xZ, Phi, N);
	UnpackValues(recvList_Xz, recvCount_Xz,recvbuf_Xz, Phi, N);
	UnpackValues(recvList_XZ, recvCount_XZ,recvbuf_XZ, Phi, N);
	UnpackValues(recvList_yz, recvCount_yz,recvbuf_yz, Phi, N);
	UnpackValues(recvList_yZ, recvCount_yZ,recvbuf_yZ, Phi, N);
	UnpackValues(recvList_Yz, recvCount_Yz,recvbuf_Yz, Phi, N);
	UnpackValues(recvList_YZ, recvCount_YZ,recvbuf_YZ, Phi, N);
	//...................................................................................

	int timestep = 0;
	if (rank==0) printf("********************************************************\n");
	if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);

	//.......create and start timer............
	double starttime,stoptime,cputime;
	MPI_Barrier(comm);
	starttime = MPI_Wtime();
	//.........................................

	sendtag = recvtag = 5;

	//************ MAIN ITERATION LOOP ***************************************/
	while (timestep < timestepMax){

		//*************************************************************************
		// 		Compute the color gradient
		//*************************************************************************
		ComputeColorGradient(ID, Phi, ColorGrad, Nx, Ny, Nz);
		//*************************************************************************

		//*************************************************************************
		// 		Perform collision step for the momentum transport
		//*************************************************************************
		ColorCollide(ID, f_even, f_odd, ColorGrad, Velocity, Nz, Ny, Nz,
				rlxA, rlxB,alpha, beta, Fx, Fy, Fz, pBC);
		//*************************************************************************

		//...................................................................................
		PackDist(1,sendList_x,0,sendCount_x,sendbuf_x,f_even,N);
		PackDist(4,sendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist(5,sendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist(6,sendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		PackDist(7,sendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		//...Packing for X face(faceGrid,packThreads,1,7,9,11,13)................................
		PackDist(0,sendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(3,sendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(4,sendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(5,sendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		PackDist(6,sendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		//...Packing for y face(faceGrid,packThreads,4,8,9,16,18).................................
		PackDist(2,sendList_y,0,sendCount_y,sendbuf_y,f_even,N);
		PackDist(4,sendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		PackDist(4,sendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
		PackDist(8,sendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		PackDist(9,sendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		//...Packing for Y face(faceGrid,packThreads,3,7,10,15,17).................................
		PackDist(1,sendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist(3,sendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist(5,sendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
		PackDist(7,sendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		PackDist(8,sendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		//...Packing for z face(faceGrid,packThreads,6,12,13,16,17)................................
		PackDist(3,sendList_z,0,sendCount_z,sendbuf_z,f_even,N);
		PackDist(6,sendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		PackDist(6,sendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		PackDist(8,sendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		PackDist(8,sendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		//...Packing for Z face(faceGrid,packThreads,5,11,14,15,18)................................
		PackDist(2,sendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist(5,sendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist(7,sendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		PackDist(7,sendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		PackDist(9,sendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		//...Pack the xy edge (edgeGrid,packThreads,8)................................
		PackDist(4,sendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
		//...Pack the Xy edge (edgeGrid,packThreads,9)................................
		PackDist(4,sendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
		//...Pack the xY edge (edgeGrid,packThreads,10)................................
		PackDist(5,sendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
		//...Pack the XY edge (edgeGrid,packThreads,7)................................
		PackDist(3,sendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
		//...Pack the xz edge (edgeGrid,packThreads,12)................................
		PackDist(6,sendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
		//...Pack the xZ edge (edgeGrid,packThreads,14)................................
		PackDist(7,sendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
		//...Pack the Xz edge (edgeGrid,packThreads,13)................................
		PackDist(6,sendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
		//...Pack the XZ edge (edgeGrid,packThreads,11)................................
		PackDist(5,sendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
		//...Pack the xz edge (edgeGrid,packThreads,12)................................
		PackDist(6,sendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
		//...Pack the xZ edge (edgeGrid,packThreads,14)................................
		PackDist(7,sendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
		//...Pack the Xz edge (edgeGrid,packThreads,13)................................
		PackDist(6,sendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
		//...Pack the XZ edge (edgeGrid,packThreads,11)................................
		PackDist(5,sendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
		//...Pack the yz edge (edgeGrid,packThreads,16)................................
		PackDist(8,sendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
		//...Pack the yZ edge (edgeGrid,packThreads,18)................................
		PackDist(9,sendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
		//...Pack the Yz edge (edgeGrid,packThreads,17)................................
		PackDist(8,sendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
		//...Pack the YZ edge (edgeGrid,packThreads,15)................................
		PackDist(7,sendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
		//...................................................................................

		//...................................................................................
		// Send all the distributions
		MPI_Isend(sendbuf_x, 5*sendCount_x,MPI_DOUBLE,rank_X,sendtag,comm,&req1[0]);
		MPI_Irecv(recvbuf_X, 5*recvCount_X,MPI_DOUBLE,rank_x,recvtag,comm,&req2[0]);
		MPI_Isend(sendbuf_X, 5*sendCount_X,MPI_DOUBLE,rank_x,sendtag,comm,&req1[1]);
		MPI_Irecv(recvbuf_x, 5*recvCount_x,MPI_DOUBLE,rank_X,recvtag,comm,&req2[1]);
		MPI_Isend(sendbuf_y, 5*sendCount_y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[2]);
		MPI_Irecv(recvbuf_Y, 5*recvCount_Y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[2]);
		MPI_Isend(sendbuf_Y, 5*sendCount_Y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[3]);
		MPI_Irecv(recvbuf_y, 5*recvCount_y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[3]);
		MPI_Isend(sendbuf_z, 5*sendCount_z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[4]);
		MPI_Irecv(recvbuf_Z, 5*recvCount_Z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[4]);
		MPI_Isend(sendbuf_Z, 5*sendCount_Z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[5]);
		MPI_Irecv(recvbuf_z, 5*recvCount_z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[5]);
		MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_XY,sendtag,comm,&req1[6]);
		MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_xy,recvtag,comm,&req2[6]);
		MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_xy,sendtag,comm,&req1[7]);
		MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_XY,recvtag,comm,&req2[7]);
		MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_xY,sendtag,comm,&req1[8]);
		MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_Xy,recvtag,comm,&req2[8]);
		MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_Xy,sendtag,comm,&req1[9]);
		MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_xY,recvtag,comm,&req2[9]);
		MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_XZ,sendtag,comm,&req1[10]);
		MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_xz,recvtag,comm,&req2[10]);
		MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_xz,sendtag,comm,&req1[11]);
		MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_XZ,recvtag,comm,&req2[11]);
		MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_xZ,sendtag,comm,&req1[12]);
		MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_Xz,recvtag,comm,&req2[12]);
		MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_Xz,sendtag,comm,&req1[13]);
		MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_xZ,recvtag,comm,&req2[13]);
		MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_YZ,sendtag,comm,&req1[14]);
		MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_yz,recvtag,comm,&req2[14]);
		MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_yz,sendtag,comm,&req1[15]);
		MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_YZ,recvtag,comm,&req2[15]);
		MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_yZ,sendtag,comm,&req1[16]);
		MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_Yz,recvtag,comm,&req2[16]);
		MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_Yz,sendtag,comm,&req1[17]);
		MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_yZ,recvtag,comm,&req2[17]);
		//...................................................................................

		//*************************************************************************
		// 		Carry out the density streaming step for mass transport
		//*************************************************************************
		DensityStreamD3Q7(ID, Den, Copy, Phi, ColorGrad, Velocity,beta, Nx, Ny, Nz, pBC);
		//*************************************************************************


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
		MapRecvDist(0,-1,0,0,recvList_X,0,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist(3,-1,-1,0,recvList_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist(4,-1,1,0,recvList_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist(5,-1,0,-1,recvList_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		MapRecvDist(6,-1,0,1,recvList_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the x face: q=1,7,9,11,13..................................
		MapRecvDist(1,1,0,0,recvList_x,0,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist(4,1,1,0,recvList_x,recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist(5,1,-1,0,recvList_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist(6,1,0,1,recvList_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		MapRecvDist(7,1,0,-1,recvList_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
		MapRecvDist(1,0,-1,0,recvList_Y,0,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		MapRecvDist(3,-1,-1,0,recvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		MapRecvDist(5,1,-1,0,recvList_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,Nx,Ny,Nz);
		MapRecvDist(7,0,-1,-1,recvList_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		MapRecvDist(8,0,-1,1,recvList_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
		MapRecvDist(2,0,1,0,recvList_y,0,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		MapRecvDist(4,1,1,0,recvList_y,recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		MapRecvDist(4,-1,1,0,recvList_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,Nx,Ny,Nz);
		MapRecvDist(8,0,1,1,recvList_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		MapRecvDist(9,0,1,-1,recvList_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the z face<<<faceGrid,packThreads,6,12,13,16,17)..............................................
		MapRecvDist(2,0,0,-1,recvList_Z,0,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		MapRecvDist(5,-1,0,-1,recvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		MapRecvDist(7,1,0,-1,recvList_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		MapRecvDist(7,0,-1,-1,recvList_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		MapRecvDist(9,0,1,-1,recvList_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		//...Map recieve list for the Z face<<<faceGrid,packThreads,5,11,14,15,18)..............................................
		MapRecvDist(3,0,0,1,recvList_z,0,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		MapRecvDist(6,1,0,1,recvList_z,recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		MapRecvDist(6,-1,0,1,recvList_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		MapRecvDist(8,0,1,1,recvList_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		MapRecvDist(8,0,-1,1,recvList_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		//..................................................................................
		//...Map recieve list for the xy edge <<<edgeGrid,packThreads,8)................................
		MapRecvDist(3,-1,-1,0,recvList_XY,0,recvCount_XY,recvbuf_XY,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xy edge <<<edgeGrid,packThreads,9)................................
		MapRecvDist(5,1,-1,0,recvList_xY,0,recvCount_xY,recvbuf_xY,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xY edge <<<edgeGrid,packThreads,10)................................
		MapRecvDist(4,-1,1,0,recvList_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the XY edge <<<edgeGrid,packThreads,7)................................
		MapRecvDist(4,1,1,0,recvList_xy,0,recvCount_xy,recvbuf_xy,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xz edge <<<edgeGrid,packThreads,12)................................
		MapRecvDist(5,-1,0,-1,recvList_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the xZ edge <<<edgeGrid,packThreads,14)................................
		MapRecvDist(6,-1,0,1,recvList_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xz edge <<<edgeGrid,packThreads,13)................................
		MapRecvDist(7,1,0,-1,recvList_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the XZ edge <<<edgeGrid,packThreads,11)................................
		MapRecvDist(6,1,0,1,recvList_xz,0,recvCount_xz,recvbuf_xz,f_even,Nx,Ny,Nz);
		//...Map recieve list for the yz edge <<<edgeGrid,packThreads,16)................................
		MapRecvDist(7,0,-1,-1,recvList_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the yZ edge <<<edgeGrid,packThreads,18)................................
		MapRecvDist(8,0,-1,1,recvList_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Yz edge <<<edgeGrid,packThreads,17)................................
		MapRecvDist(9,0,1,-1,recvList_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the YZ edge <<<edgeGrid,packThreads,15)................................
		MapRecvDist(8,0,1,1,recvList_yz,0,recvCount_yz,recvbuf_yz,f_even,Nx,Ny,Nz);
		//...................................................................................
		
		//...................................................................................
		PackDenD3Q7(recvList_x,recvCount_x,recvbuf_x,2,Den,N);
		PackDenD3Q7(recvList_y,recvCount_y,recvbuf_y,2,Den,N);
		PackDenD3Q7(recvList_z,recvCount_z,recvbuf_z,2,Den,N);
		PackDenD3Q7(recvList_X,recvCount_X,recvbuf_X,2,Den,N);
		PackDenD3Q7(recvList_Y,recvCount_Y,recvbuf_Y,2,Den,N);
		PackDenD3Q7(recvList_Z,recvCount_Z,recvbuf_Z,2,Den,N);
		//...................................................................................
		//...................................................................................
		// Send all the D3Q7 distributions
		MPI_Isend(recvbuf_x, 2*recvCount_x,MPI_DOUBLE,rank_X,sendtag,comm,&req1[0]);
		MPI_Irecv(sendbuf_X, 2*sendCount_X,MPI_DOUBLE,rank_x,recvtag,comm,&req2[0]);
		MPI_Isend(recvbuf_X, 2*recvCount_X,MPI_DOUBLE,rank_x,sendtag,comm,&req1[1]);
		MPI_Irecv(sendbuf_x, 2*sendCount_x,MPI_DOUBLE,rank_X,recvtag,comm,&req2[1]);
		MPI_Isend(recvbuf_y, 2*recvCount_y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[2]);
		MPI_Irecv(sendbuf_Y, 2*sendCount_Y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[2]);
		MPI_Isend(recvbuf_Y, 2*recvCount_Y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[3]);
		MPI_Irecv(sendbuf_y, 2*sendCount_y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[3]);
		MPI_Isend(recvbuf_z, 2*recvCount_z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[4]);
		MPI_Irecv(sendbuf_Z, 2*sendCount_Z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[4]);
		MPI_Isend(recvbuf_Z, 2*recvCount_Z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[5]);
		MPI_Irecv(sendbuf_z, 2*sendCount_z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[5]);
		//...................................................................................
		//...................................................................................
		// Wait for completion of D3Q7 communication
		MPI_Waitall(6,req1,stat1);
		MPI_Waitall(6,req2,stat2);
		//...................................................................................
		//...................................................................................
		UnpackDenD3Q7(sendList_x,sendCount_x,sendbuf_x,2,Den,N);
		UnpackDenD3Q7(sendList_y,sendCount_y,sendbuf_y,2,Den,N);
		UnpackDenD3Q7(sendList_z,sendCount_z,sendbuf_z,2,Den,N);
		UnpackDenD3Q7(sendList_X,sendCount_X,sendbuf_X,2,Den,N);
		UnpackDenD3Q7(sendList_Y,sendCount_Y,sendbuf_Y,2,Den,N);
		UnpackDenD3Q7(sendList_Z,sendCount_Z,sendbuf_Z,2,Den,N);
		//...................................................................................
		//*************************************************************************
		// 		Compute the phase indicator field and reset Copy, Den
		//*************************************************************************
		ComputePhi(ID, Phi, Copy, Den, N);
		//*************************************************************************
		//...................................................................................
		PackValues(sendList_x, sendCount_x,sendbuf_x, Phi, N);
		PackValues(sendList_y, sendCount_y,sendbuf_y, Phi, N);
		PackValues(sendList_z, sendCount_z,sendbuf_z, Phi, N);
		PackValues(sendList_X, sendCount_X,sendbuf_X, Phi, N);
		PackValues(sendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
		PackValues(sendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
		PackValues(sendList_xy, sendCount_xy,sendbuf_xy, Phi, N);
		PackValues(sendList_xY, sendCount_xY,sendbuf_xY, Phi, N);
		PackValues(sendList_Xy, sendCount_Xy,sendbuf_Xy, Phi, N);
		PackValues(sendList_XY, sendCount_XY,sendbuf_XY, Phi, N);
		PackValues(sendList_xz, sendCount_xz,sendbuf_xz, Phi, N);
		PackValues(sendList_xZ, sendCount_xZ,sendbuf_xZ, Phi, N);
		PackValues(sendList_Xz, sendCount_Xz,sendbuf_Xz, Phi, N);
		PackValues(sendList_XZ, sendCount_XZ,sendbuf_XZ, Phi, N);
		PackValues(sendList_yz, sendCount_yz,sendbuf_yz, Phi, N);
		PackValues(sendList_yZ, sendCount_yZ,sendbuf_yZ, Phi, N);
		PackValues(sendList_Yz, sendCount_Yz,sendbuf_Yz, Phi, N);
		PackValues(sendList_YZ, sendCount_YZ,sendbuf_YZ, Phi, N);
		//...................................................................................
		// Send / Recv all the phase indcator field values
		MPI_Isend(sendbuf_x, sendCount_x,MPI_DOUBLE,rank_X,sendtag,comm,&req1[0]);
		MPI_Irecv(recvbuf_X, recvCount_X,MPI_DOUBLE,rank_x,recvtag,comm,&req2[0]);
		MPI_Isend(sendbuf_X, sendCount_X,MPI_DOUBLE,rank_x,sendtag,comm,&req1[1]);
		MPI_Irecv(recvbuf_x, recvCount_x,MPI_DOUBLE,rank_X,recvtag,comm,&req2[1]);
		MPI_Isend(sendbuf_y, sendCount_y,MPI_DOUBLE,rank_Y,sendtag,comm,&req1[2]);
		MPI_Irecv(recvbuf_Y, recvCount_Y,MPI_DOUBLE,rank_y,recvtag,comm,&req2[2]);
		MPI_Isend(sendbuf_Y, sendCount_Y,MPI_DOUBLE,rank_y,sendtag,comm,&req1[3]);
		MPI_Irecv(recvbuf_y, recvCount_y,MPI_DOUBLE,rank_Y,recvtag,comm,&req2[3]);
		MPI_Isend(sendbuf_z, sendCount_z,MPI_DOUBLE,rank_Z,sendtag,comm,&req1[4]);
		MPI_Irecv(recvbuf_Z, recvCount_Z,MPI_DOUBLE,rank_z,recvtag,comm,&req2[4]);
		MPI_Isend(sendbuf_Z, sendCount_Z,MPI_DOUBLE,rank_z,sendtag,comm,&req1[5]);
		MPI_Irecv(recvbuf_z, recvCount_z,MPI_DOUBLE,rank_Z,recvtag,comm,&req2[5]);
		MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_XY,sendtag,comm,&req1[6]);
		MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_xy,recvtag,comm,&req2[6]);
		MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_xy,sendtag,comm,&req1[7]);
		MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_XY,recvtag,comm,&req2[7]);
		MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_xY,sendtag,comm,&req1[8]);
		MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_Xy,recvtag,comm,&req2[8]);
		MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_Xy,sendtag,comm,&req1[9]);
		MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_xY,recvtag,comm,&req2[9]);
		MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_XZ,sendtag,comm,&req1[10]);
		MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_xz,recvtag,comm,&req2[10]);
		MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_xz,sendtag,comm,&req1[11]);
		MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_XZ,recvtag,comm,&req2[11]);
		MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_xZ,sendtag,comm,&req1[12]);
		MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_Xz,recvtag,comm,&req2[12]);
		MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_Xz,sendtag,comm,&req1[13]);
		MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_xZ,recvtag,comm,&req2[13]);
		MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_YZ,sendtag,comm,&req1[14]);
		MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_yz,recvtag,comm,&req2[14]);
		MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_yz,sendtag,comm,&req1[15]);
		MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_YZ,recvtag,comm,&req2[15]);
		MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_yZ,sendtag,comm,&req1[16]);
		MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_Yz,recvtag,comm,&req2[16]);
		MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_Yz,sendtag,comm,&req1[17]);
		MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_yZ,recvtag,comm,&req2[17]);
		//...................................................................................
		//...................................................................................
		//...................................................................................
		// Wait for completion of Indicator Field communication
		MPI_Waitall(18,req1,stat1);
		MPI_Waitall(18,req2,stat2);
		//...................................................................................
		//...................................................................................
		UnpackValues(recvList_x, recvCount_x,recvbuf_x, Phi, N);
		UnpackValues(recvList_y, recvCount_y,recvbuf_y, Phi, N);
		UnpackValues(recvList_z, recvCount_z,recvbuf_z, Phi, N);
		UnpackValues(recvList_X, recvCount_X,recvbuf_X, Phi, N);
		UnpackValues(recvList_Y, recvCount_Y,recvbuf_Y, Phi, N);
		UnpackValues(recvList_Z, recvCount_Z,recvbuf_Z, Phi, N);
		UnpackValues(recvList_xy, recvCount_xy,recvbuf_xy, Phi, N);
		UnpackValues(recvList_xY, recvCount_xY,recvbuf_xY, Phi, N);
		UnpackValues(recvList_Xy, recvCount_Xy,recvbuf_Xy, Phi, N);
		UnpackValues(recvList_XY, recvCount_XY,recvbuf_XY, Phi, N);
		UnpackValues(recvList_xz, recvCount_xz,recvbuf_xz, Phi, N);
		UnpackValues(recvList_xZ, recvCount_xZ,recvbuf_xZ, Phi, N);
		UnpackValues(recvList_Xz, recvCount_Xz,recvbuf_Xz, Phi, N);
		UnpackValues(recvList_XZ, recvCount_XZ,recvbuf_XZ, Phi, N);
		UnpackValues(recvList_yz, recvCount_yz,recvbuf_yz, Phi, N);
		UnpackValues(recvList_yZ, recvCount_yZ,recvbuf_yZ, Phi, N);
		UnpackValues(recvList_Yz, recvCount_Yz,recvbuf_Yz, Phi, N);
		UnpackValues(recvList_YZ, recvCount_YZ,recvbuf_YZ, Phi, N);
		//...................................................................................

		MPI_Barrier(comm);
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
	
	// Write output files
	sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
//	printf("Local File Name =  %s \n",LocalRankFilename);
	
	FILE *PHASE;
	PHASE = fopen(LocalRankFilename,"wb");
	fwrite(Phi,8,N,PHASE);
	fclose(PHASE);	
	
	double *pressure;
	pressure = new double[N];
	for (n=0; n<N; n++){
		pressure[n] = f_even[n];
		for (int q=0; q<9; q++){
			pressure[n] += f_even[(q+1)*N+n];
			pressure[n] += f_odd[q*N+n];
		}
		if (pressure[n] < 0.f) pressure[n] = 0.f;
	}
	sprintf(LocalRankFilename,"%s%s","Pressure.",LocalRankString);
	FILE *PRESSURE;
	PRESSURE = fopen(LocalRankFilename,"wb");
	fwrite(pressure,8,N,PRESSURE);
	fclose(PRESSURE);
	
	sprintf(LocalRankFilename,"%s%s","ColorGrad.",LocalRankString);
	FILE *COLORGRAD;
	COLORGRAD = fopen(LocalRankFilename,"wb");
	fwrite(ColorGrad,8,3*N,COLORGRAD);
	fclose(COLORGRAD);
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
