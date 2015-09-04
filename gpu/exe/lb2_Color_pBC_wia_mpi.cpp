#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <mpi.h>

#include "pmmc.h"
#include "Domain.h"

using namespace std;

//*************************************************************************
// Functions defined in Color.cu
//*************************************************************************
extern "C" void dvc_InitDenColor( int nblocks, int nthreads, int S,
		char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz);
//*************************************************************************
extern "C" void dvc_ComputeColorGradient(int nBlocks, int nthreads, int S,
		char *ID, double *Phi, double *ColorGrad, int Nx, int Ny, int Nz);
//*************************************************************************
extern "C" void dvc_ColorCollide(int nBlocks, int nthreads, int S,
		char *ID, double *f_even, double *f_odd, double *ColorGrad, double *Velocity,
		double rlxA, double rlxB,double alpha, double beta, double Fx, double Fy, double Fz,
		int Nx, int Ny, int Nz, bool pBC);
//*************************************************************************
extern "C" void dvc_DensityStreamD3Q7(int nBlocks, int nthreads, int S,
		char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
		double beta, int Nx, int Ny, int Nz, bool pBC);
//*************************************************************************
extern "C" void dvc_ComputePhi(int nBlocks, int nthreads, int S,
		char *ID, double *Phi, double *Copy, double *Den, int N);
//*************************************************************************
//*************************************************************************
// Functions defined in D3Q19.cu
//*************************************************************************
extern "C" void dvc_InitD3Q19(int nblocks, int nthreads, int S, char *ID, double *f_even, double *f_odd, int Nx,
							  int Ny, int Nz);
//*************************************************************************
extern "C" void dvc_SwapD3Q19(int nblocks, int nthreads, int S,
		char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);
//*************************************************************************
extern "C" void dvc_PackDist(int grid, int threads, int q, int *SendList, int start,
		int sendCount, double *sendbuf, double *Dist, int N);
//*************************************************************************
extern "C" void dvc_UnpackDist(int grid, int threads, int q, int Cqx, int Cqy, int Cqz, int *RecvList, int start,
		int recvCount, double *recvbuf, double *Dist, int Nx, int Ny, int Nz);
//*************************************************************************
//***************************************************************************************
// Functions defined in D3Q7.cu
//***************************************************************************************
extern "C" void dvc_PackDenD3Q7(int grid, int threads, int *list, int count, double *sendbuf,
		int number, double *Data, int N);
//***************************************************************************************
extern "C" void dvc_UnpackDenD3Q7(int grid, int threads, int *list, int count, double *recvbuf,
		int number, double *Data, int N);
//***************************************************************************************
extern "C" void dvc_PackValues(int grid, int threads, int *list, int count, double *sendbuf,
		double *Data, int N);
//***************************************************************************************
extern "C" void dvc_UnpackValues(int grid, int threads, int *list, int count, double *recvbuf,
		double *Data, int N);
//***************************************************************************************
//*************************************************************************
// Functions defined in CudaExtras.cu
//*************************************************************************
extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size);
//*************************************************************************
extern "C" void dvc_CopyToDevice(void* dest, void* source, size_t size);
//*************************************************************************
extern "C" void dvc_CopyToHost(void* dest, void* source, size_t size);
//*************************************************************************
extern "C" void dvc_Barrier();
//*************************************************************************

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

	// Variables that specify the computational domain  
	string FILENAME;
	unsigned int nBlocks, nthreads;
	int Nx,Ny,Nz;
	int nspheres;
	double Lx,Ly,Lz;
	// Color Model parameters
	int timestepMax, interval;
	double tau,Fx,Fy,Fz,tol;
	double alpha, beta;
	double das, dbs;
	double din,dout;
	double wp_saturation;
	bool pBC;
	int i,j,k,p,q,r,n;

	// pmmc threshold values
	double fluid_isovalue,solid_isovalue;
	fluid_isovalue = 0.0;
	solid_isovalue = 0.0;

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
		// Line 4: wetting phase saturation to initialize
		input >> wp_saturation;
		// Line 5: External force components (Fx,Fy, Fz)
		input >> Fx;
		input >> Fy;
		input >> Fz;
		// Line 6: Pressure Boundary conditions
		input >> pBC;
		input >> din;
		input >> dout;
		// Line 7: time-stepping criteria
		input >> timestepMax;			// max no. of timesteps
		input >> interval;			// error interval
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
	// **************************************************************
	// Broadcast simulation parameters from rank 0 to all other procs
	MPI_Barrier(comm);
	//.................................................
	MPI_Bcast(&tau,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&alpha,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&beta,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&das,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&dbs,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&pBC,1,MPI_LOGICAL,0,comm);
	MPI_Bcast(&din,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&dout,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fy,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Fz,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&timestepMax,1,MPI_INT,0,comm);
	MPI_Bcast(&interval,1,MPI_INT,0,comm);
	MPI_Bcast(&tol,1,MPI_DOUBLE,0,comm);
	// Computational domain
	MPI_Bcast(&Nz,1,MPI_INT,0,comm);
	MPI_Bcast(&nBlocks,1,MPI_INT,0,comm);
	MPI_Bcast(&nthreads,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
	//.................................................
	MPI_Barrier(comm);
	// **************************************************************
	// **************************************************************
	double Ps = -(das-dbs)/(das+dbs);
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
		printf("phi_s = %f \n", Ps);
		printf("gamma_{wn} = %f \n", 6.01603*alpha);
		printf("cos theta_c = %f \n", 1.05332*Ps);
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
	char LocalRankString[8];
	char LocalRankFilename[40];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
//	printf("Local File Name =  %s \n",LocalRankFilename);
	// .......... READ THE INPUT FILE .......................................
//	char value;
	char *id;
	id = new char[N];
	int sum = 0;
//	double porosity;
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
*/
	//...........................................................................
//	double *SignDist;
//	SignDist = new double[N];
	DoubleArray SignDist(Nx,Ny,Nz);
	//.......................................................................
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
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (SignDist.data[n] > 0.0){ 
					id[n] = 1;	
					sum++;
				}
			}
		}
	}
	//......................................................................	
	// Once phase ID has been generated, map solid to account for 'smeared' interface
	//......................................................................	
	for (i=0; i<N; i++)	SignDist.data[i] -= 0.5; // Solid appears half a pixel bigger
	//......................................................................	
	// Generate the residual NWP 
	GenerateResidual(id,Nx,Ny,Nz,wp_saturation);
	//.......................................................................
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	WriteLocalSolidID(LocalRankFilename, id, N);
	sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
	WriteLocalSolidDistance(LocalRankFilename, SignDist.data, N);
	//.......................................................................

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
	MPI_Isend(&sendCount_x, 1,MPI_INT,rank_X,sendtag,comm,&req1[0]);
	MPI_Irecv(&recvCount_X, 1,MPI_INT,rank_x,recvtag,comm,&req2[0]);
	MPI_Isend(&sendCount_X, 1,MPI_INT,rank_x,sendtag,comm,&req1[1]);
	MPI_Irecv(&recvCount_x, 1,MPI_INT,rank_X,recvtag,comm,&req2[1]);
	MPI_Isend(&sendCount_y, 1,MPI_INT,rank_Y,sendtag,comm,&req1[2]);
	MPI_Irecv(&recvCount_Y, 1,MPI_INT,rank_y,recvtag,comm,&req2[2]);
	MPI_Isend(&sendCount_Y, 1,MPI_INT,rank_y,sendtag,comm,&req1[3]);
	MPI_Irecv(&recvCount_y, 1,MPI_INT,rank_Y,recvtag,comm,&req2[3]);
	MPI_Isend(&sendCount_z, 1,MPI_INT,rank_Z,sendtag,comm,&req1[4]);
	MPI_Irecv(&recvCount_Z, 1,MPI_INT,rank_z,recvtag,comm,&req2[4]);
	MPI_Isend(&sendCount_Z, 1,MPI_INT,rank_z,sendtag,comm,&req1[5]);
	MPI_Irecv(&recvCount_z, 1,MPI_INT,rank_Z,recvtag,comm,&req2[5]);

	MPI_Isend(&sendCount_xy, 1,MPI_INT,rank_XY,sendtag,comm,&req1[6]);
	MPI_Irecv(&recvCount_XY, 1,MPI_INT,rank_xy,recvtag,comm,&req2[6]);
	MPI_Isend(&sendCount_XY, 1,MPI_INT,rank_xy,sendtag,comm,&req1[7]);
	MPI_Irecv(&recvCount_xy, 1,MPI_INT,rank_XY,recvtag,comm,&req2[7]);
	MPI_Isend(&sendCount_Xy, 1,MPI_INT,rank_xY,sendtag,comm,&req1[8]);
	MPI_Irecv(&recvCount_xY, 1,MPI_INT,rank_Xy,recvtag,comm,&req2[8]);
	MPI_Isend(&sendCount_xY, 1,MPI_INT,rank_Xy,sendtag,comm,&req1[9]);
	MPI_Irecv(&recvCount_Xy, 1,MPI_INT,rank_xY,recvtag,comm,&req2[9]);

	MPI_Isend(&sendCount_xz, 1,MPI_INT,rank_XZ,sendtag,comm,&req1[10]);
	MPI_Irecv(&recvCount_XZ, 1,MPI_INT,rank_xz,recvtag,comm,&req2[10]);
	MPI_Isend(&sendCount_XZ, 1,MPI_INT,rank_xz,sendtag,comm,&req1[11]);
	MPI_Irecv(&recvCount_xz, 1,MPI_INT,rank_XZ,recvtag,comm,&req2[11]);
	MPI_Isend(&sendCount_Xz, 1,MPI_INT,rank_xZ,sendtag,comm,&req1[12]);
	MPI_Irecv(&recvCount_xZ, 1,MPI_INT,rank_Xz,recvtag,comm,&req2[12]);
	MPI_Isend(&sendCount_xZ, 1,MPI_INT,rank_Xz,sendtag,comm,&req1[13]);
	MPI_Irecv(&recvCount_Xz, 1,MPI_INT,rank_xZ,recvtag,comm,&req2[13]);

	MPI_Isend(&sendCount_yz, 1,MPI_INT,rank_YZ,sendtag,comm,&req1[14]);
	MPI_Irecv(&recvCount_YZ, 1,MPI_INT,rank_yz,recvtag,comm,&req2[14]);
	MPI_Isend(&sendCount_YZ, 1,MPI_INT,rank_yz,sendtag,comm,&req1[15]);
	MPI_Irecv(&recvCount_yz, 1,MPI_INT,rank_YZ,recvtag,comm,&req2[15]);
	MPI_Isend(&sendCount_Yz, 1,MPI_INT,rank_yZ,sendtag,comm,&req1[16]);
	MPI_Irecv(&recvCount_yZ, 1,MPI_INT,rank_Yz,recvtag,comm,&req2[16]);
	MPI_Isend(&sendCount_yZ, 1,MPI_INT,rank_Yz,sendtag,comm,&req1[17]);
	MPI_Irecv(&recvCount_Yz, 1,MPI_INT,rank_yZ,recvtag,comm,&req2[17]);
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	MPI_Barrier(comm);
/*	MPI_Send(&sendCount_x,1,MPI_INT,rank_X,sendtag,comm);
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
*/	//**********************************************************************************
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
	dvc_AllocateDeviceMemory((void **) &sendbuf_x, 5*sendCount_x*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_X, 5*sendCount_X*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_y, 5*sendCount_y*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_Y, 5*sendCount_Y*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_z, 5*sendCount_z*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_Z, 5*sendCount_Z*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_xy, sendCount_xy*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_xY, sendCount_xY*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_Xy, sendCount_Xy*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_XY, sendCount_XY*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_xz, sendCount_xz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_xZ, sendCount_xZ*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_Xz, sendCount_Xz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_XZ, sendCount_XZ*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_yz, sendCount_yz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_yZ, sendCount_yZ*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_Yz, sendCount_Yz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &sendbuf_YZ, sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	dvc_AllocateDeviceMemory((void **) &recvbuf_x, 5*recvCount_x*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_X, 5*recvCount_X*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_y, 5*recvCount_y*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_Y, 5*recvCount_Y*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_z, 5*recvCount_z*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_Z, 5*recvCount_Z*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_xy, recvCount_xy*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_xY, recvCount_xY*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_Xy, recvCount_Xy*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_XY, recvCount_XY*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_xz, recvCount_xz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_xZ, recvCount_xZ*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_Xz, recvCount_Xz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_XZ, recvCount_XZ*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_yz, recvCount_yz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_yZ, recvCount_yZ*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_Yz, recvCount_Yz*sizeof(double));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &recvbuf_YZ, recvCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X, *dvcSendList_Y, *dvcSendList_Z;
	int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy, *dvcSendList_Yz, *dvcSendList_xZ;
	int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY, *dvcSendList_YZ, *dvcSendList_XZ;
	//......................................................................................
	int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X, *dvcRecvList_Y, *dvcRecvList_Z;
	int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy, *dvcRecvList_Yz, *dvcRecvList_xZ;
	int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY, *dvcRecvList_YZ, *dvcRecvList_XZ;
	//......................................................................................
	dvc_AllocateDeviceMemory((void **) &dvcSendList_x, sendCount_x*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_X, sendCount_X*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_y, sendCount_y*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_Y, sendCount_Y*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_z, sendCount_z*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_Z, sendCount_Z*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_xy, sendCount_xy*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_xY, sendCount_xY*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_Xy, sendCount_Xy*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_XY, sendCount_XY*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_xz, sendCount_xz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_xZ, sendCount_xZ*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_Xz, sendCount_Xz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_XZ, sendCount_XZ*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_yz, sendCount_yz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_yZ, sendCount_yZ*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_Yz, sendCount_Yz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcSendList_YZ, sendCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_x, recvCount_x*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_X, recvCount_X*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_y, recvCount_y*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_Y, recvCount_Y*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_z, recvCount_z*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_Z, recvCount_Z*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &dvcRecvList_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	if (rank==0)	printf ("Prepare to copy send/recv Lists to device \n");
	dvc_CopyToDevice(dvcSendList_x,sendList_x,sendCount_x*sizeof(int));
	dvc_CopyToDevice(dvcSendList_X,sendList_X,sendCount_X*sizeof(int));
	dvc_CopyToDevice(dvcSendList_y,sendList_y,sendCount_y*sizeof(int));
	dvc_CopyToDevice(dvcSendList_Y,sendList_Y,sendCount_Y*sizeof(int));
	dvc_CopyToDevice(dvcSendList_z,sendList_z,sendCount_z*sizeof(int));
	dvc_CopyToDevice(dvcSendList_Z,sendList_Z,sendCount_Z*sizeof(int));
	dvc_CopyToDevice(dvcSendList_xy,sendList_xy,sendCount_xy*sizeof(int));
	dvc_CopyToDevice(dvcSendList_XY,sendList_XY,sendCount_XY*sizeof(int));
	dvc_CopyToDevice(dvcSendList_xY,sendList_xY,sendCount_xY*sizeof(int));
	dvc_CopyToDevice(dvcSendList_Xy,sendList_Xy,sendCount_Xy*sizeof(int));
	dvc_CopyToDevice(dvcSendList_xz,sendList_xz,sendCount_xz*sizeof(int));
	dvc_CopyToDevice(dvcSendList_XZ,sendList_XZ,sendCount_XZ*sizeof(int));
	dvc_CopyToDevice(dvcSendList_xZ,sendList_xZ,sendCount_xZ*sizeof(int));
	dvc_CopyToDevice(dvcSendList_Xz,sendList_Xz,sendCount_Xz*sizeof(int));
	dvc_CopyToDevice(dvcSendList_yz,sendList_yz,sendCount_yz*sizeof(int));
	dvc_CopyToDevice(dvcSendList_YZ,sendList_YZ,sendCount_YZ*sizeof(int));
	dvc_CopyToDevice(dvcSendList_yZ,sendList_yZ,sendCount_yZ*sizeof(int));
	dvc_CopyToDevice(dvcSendList_Yz,sendList_Yz,sendCount_Yz*sizeof(int));
	//......................................................................................
	dvc_CopyToDevice(dvcRecvList_x,recvList_x,recvCount_x*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_X,recvList_X,recvCount_X*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_y,recvList_y,recvCount_y*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_Y,recvList_Y,recvCount_Y*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_z,recvList_z,recvCount_z*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_Z,recvList_Z,recvCount_Z*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_xy,recvList_xy,recvCount_xy*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_XY,recvList_XY,recvCount_XY*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_xY,recvList_xY,recvCount_xY*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_Xy,recvList_Xy,recvCount_Xy*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_xz,recvList_xz,recvCount_xz*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_XZ,recvList_XZ,recvCount_XZ*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_xZ,recvList_xZ,recvCount_xZ*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_Xz,recvList_Xz,recvCount_Xz*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_yz,recvList_yz,recvCount_yz*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_YZ,recvList_YZ,recvCount_YZ*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_yZ,recvList_yZ,recvCount_yZ*sizeof(int));
	dvc_CopyToDevice(dvcRecvList_Yz,recvList_Yz,recvCount_Yz*sizeof(int));
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
/*	// Once the ID is saved, free memory allocated to the buffers (no longer needed)
	//......................................................................................
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
	dvc_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
	// Copy to the device
	dvc_CopyToDevice(ID, id, N);
	//...........................................................................

	//...........................................................................
	//				MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	// LBM variables
	if (rank==0)	printf ("Allocating distributions \n");
	//......................device distributions.................................
	double *f_even,*f_odd;
	//...........................................................................
	dvc_AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory
	//...........................................................................
	double *Phi,*Den,*Copy;
	double *ColorGrad, *Velocity;
	//...........................................................................
	dvc_AllocateDeviceMemory((void **) &Phi, dist_mem_size);
	dvc_AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
	dvc_AllocateDeviceMemory((void **) &Copy, 2*dist_mem_size);
	dvc_AllocateDeviceMemory((void **) &Velocity, 3*dist_mem_size);
	dvc_AllocateDeviceMemory((void **) &ColorGrad, 3*dist_mem_size);
	//...........................................................................
	// Copy of Phi that is used on the Host
	double *HostPhi;
	HostPhi = new double[N];
	// Phase indicator
	DoubleArray Phase(Nx,Ny,Nz);

	/* ****************************************************************
	 VARIABLES FOR THE PMMC ALGORITHM
	 ****************************************************************** */
	//...........................................................................
	// Averaging variables
	//...........................................................................
	double awn,ans,aws,lwns,nwp_volume;
	double As;
	double dEs,dAwn,dAns;			 // Global surface energy (calculated by rank=0)
	double awn_global,ans_global,aws_global,lwns_global,nwp_volume_global;	
	double As_global;
//	bool add=1;			// Set to false if any corners contain nw-phase ( F > fluid_isovalue)
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
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
	//	IntArray store;
	
	int n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
	int n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
	
	double s,s1,s2,s3;		// Triangle sides (lengths)
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
	
	int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg;
	int c;
	int newton_steps = 0;
	//...........................................................................
	int ncubes = (Nx-2)*(Ny-2)*(Nz-2);	// Exclude the "upper" halo
	IntArray cubeList(3,ncubes);
	int nc=0;
	//...........................................................................
	// Set up the cube list (very regular in this case due to lack of blob-ID)
	for (k=0; k<Nz-2; k++){
		for (j=0; j<Ny-2; j++){
			for (i=0; i<Nx-2; i++){
				cubeList(0,nc) = i;
				cubeList(1,nc) = j;
				cubeList(2,nc) = k;
				nc++;
			}
		}
	}
	if (nc != ncubes)	fprintf("Basic cubeList error \n");
	
	//...........................................................................
	// Grids used to pack faces on the GPU for MPI
	int faceGrid,edgeGrid,packThreads;
	packThreads=512;
	edgeGrid=1;
	faceGrid=Nx*Ny/packThreads;
	//...........................................................................
	
	//...........................................................................
	//				MAIN  VARIABLES INITIALIZED HERE
	//...........................................................................
	//...........................................................................
	if (rank==0)	printf("Setting the distributions, size = : %i\n", N);
	//...........................................................................
	dvc_InitD3Q19(nBlocks, nthreads, S, ID, f_even, f_odd, Nx, Ny, Nz);
	dvc_InitDenColor(nBlocks, nthreads, S, ID, Copy, Phi,  das, dbs, Nx, Ny, Nz);
	dvc_InitDenColor(nBlocks, nthreads, S, ID, Den, Phi,  das, dbs, Nx, Ny, Nz);
	// Pack the buffers (zeros out the halo region)
	dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_x,recvCount_x,recvbuf_x,2,Den,N);
	dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_y,recvCount_y,recvbuf_y,2,Den,N);
	dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_z,recvCount_z,recvbuf_z,2,Den,N);
	dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_X,recvCount_X,recvbuf_X,2,Den,N);
	dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_Y,recvCount_Y,recvbuf_Y,2,Den,N);
	dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_Z,recvCount_Z,recvbuf_Z,2,Den,N);
	//...................................................................................
	dvc_ComputePhi(nBlocks, nthreads, S,ID, Phi, Copy, Den, N);
	//...........................................................................	
	dvc_PackValues(faceGrid, packThreads, dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_xy, sendCount_xy,sendbuf_xy, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_xY, sendCount_xY,sendbuf_xY, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_Xy, sendCount_Xy,sendbuf_Xy, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_XY, sendCount_XY,sendbuf_XY, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_xz, sendCount_xz,sendbuf_xz, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_xZ, sendCount_xZ,sendbuf_xZ, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_Xz, sendCount_Xz,sendbuf_Xz, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_XZ, sendCount_XZ,sendbuf_XZ, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_yz, sendCount_yz,sendbuf_yz, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_yZ, sendCount_yZ,sendbuf_yZ, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_Yz, sendCount_Yz,sendbuf_Yz, Phi, N);
	dvc_PackValues(faceGrid, packThreads, dvcSendList_YZ, sendCount_YZ,sendbuf_YZ, Phi, N);
	//...................................................................................
	// Send / Recv all the phase indcator field values
	//...................................................................................
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
	// Wait for completion of Indicator Field communication
	//...................................................................................
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	//...................................................................................
	//...................................................................................
	/*		dvc_UnpackValues(faceGrid, packThreads, dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
	 dvc_UnpackValues(faceGrid, packThreads, dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
	 dvc_UnpackValues(faceGrid, packThreads, dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
	 dvc_UnpackValues(faceGrid, packThreads, dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
	 dvc_UnpackValues(faceGrid, packThreads, dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
	 dvc_UnpackValues(faceGrid, packThreads, dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
	 */		
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_x, recvCount_x,recvbuf_x, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_y, recvCount_y,recvbuf_y, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_z, recvCount_z,recvbuf_z, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_X, recvCount_X,recvbuf_X, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Y, recvCount_Y,recvbuf_Y, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Z, recvCount_Z,recvbuf_Z, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xy, recvCount_xy,recvbuf_xy, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xY, recvCount_xY,recvbuf_xY, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_XY, recvCount_XY,recvbuf_XY, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xz, recvCount_xz,recvbuf_xz, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_yz, recvCount_yz,recvbuf_yz, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, Phi, N);
	dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, Phi, N);
	//...................................................................................

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
	if (rank==0) printf("-------------------------------------------------------------------\n");
	if (rank==0) printf("timestep dEs Vn Awn Ans Aws Lwns \n");
	if (rank==0) printf("-------------------------------------------------------------------\n");

	//************ MAIN ITERATION LOOP ***************************************/
	while (timestep < timestepMax){

		//*************************************************************************
		// 		Compute the color gradient
		//*************************************************************************
		dvc_ComputeColorGradient(nBlocks, nthreads, S,
				ID, Phi, ColorGrad, Nx, Ny, Nz);
		//*************************************************************************

		//*************************************************************************
		// 		Perform collision step for the momentum transport
		//*************************************************************************
		dvc_ColorCollide(nBlocks, nthreads, S, ID, f_even, f_odd, ColorGrad, Velocity,
				rlxA, rlxB,alpha, beta, Fx, Fy, Fz, Nx, Ny, Nz, pBC);
		//*************************************************************************

		//...................................................................................
		dvc_PackDist(faceGrid,packThreads,1,dvcSendList_x,0,sendCount_x,sendbuf_x,f_even,N);
		dvc_PackDist(faceGrid,packThreads,4,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		dvc_PackDist(faceGrid,packThreads,5,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		dvc_PackDist(faceGrid,packThreads,6,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		dvc_PackDist(faceGrid,packThreads,7,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,f_even,N);
		//...Packing for X face(faceGrid,packThreads,1,7,9,11,13)................................
		dvc_PackDist(faceGrid,packThreads,0,dvcSendList_X,0,sendCount_X,sendbuf_X,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,3,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,4,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,5,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,6,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,f_odd,N);
		//...Packing for y face(faceGrid,packThreads,4,8,9,16,18).................................
		dvc_PackDist(faceGrid,packThreads,2,dvcSendList_y,0,sendCount_y,sendbuf_y,f_even,N);
		dvc_PackDist(faceGrid,packThreads,4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		dvc_PackDist(faceGrid,packThreads,4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,8,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		dvc_PackDist(faceGrid,packThreads,9,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,f_even,N);
		//...Packing for Y face(faceGrid,packThreads,3,7,10,15,17).................................
		dvc_PackDist(faceGrid,packThreads,1,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,5,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,f_even,N);
		dvc_PackDist(faceGrid,packThreads,7,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,8,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,f_odd,N);
		//...Packing for z face(faceGrid,packThreads,6,12,13,16,17)................................
		dvc_PackDist(faceGrid,packThreads,3,dvcSendList_z,0,sendCount_z,sendbuf_z,f_even,N);
		dvc_PackDist(faceGrid,packThreads,6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		dvc_PackDist(faceGrid,packThreads,6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,8,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,f_even,N);
		dvc_PackDist(faceGrid,packThreads,8,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,f_odd,N);
		//...Packing for Z face(faceGrid,packThreads,5,11,14,15,18)................................
		dvc_PackDist(faceGrid,packThreads,2,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,7,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		dvc_PackDist(faceGrid,packThreads,7,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,f_odd,N);
		dvc_PackDist(faceGrid,packThreads,9,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,f_even,N);
		//...Pack the xy edge (edgeGrid,packThreads,8)................................
		dvc_PackDist(edgeGrid,packThreads,4,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,f_even,N);
		//...Pack the Xy edge (edgeGrid,packThreads,9)................................
		dvc_PackDist(edgeGrid,packThreads,4,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,f_odd,N);
		//...Pack the xY edge (edgeGrid,packThreads,10)................................
		dvc_PackDist(edgeGrid,packThreads,5,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,f_even,N);
		//...Pack the XY edge (edgeGrid,packThreads,7)................................
		dvc_PackDist(edgeGrid,packThreads,3,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,f_odd,N);
		//...Pack the xz edge (edgeGrid,packThreads,12)................................
		dvc_PackDist(edgeGrid,packThreads,6,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,f_even,N);
		//...Pack the xZ edge (edgeGrid,packThreads,14)................................
		dvc_PackDist(edgeGrid,packThreads,7,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,f_even,N);
		//...Pack the Xz edge (edgeGrid,packThreads,13)................................
		dvc_PackDist(edgeGrid,packThreads,6,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,f_odd,N);
		//...Pack the XZ edge (edgeGrid,packThreads,11)................................
		dvc_PackDist(edgeGrid,packThreads,5,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,f_odd,N);
		//...Pack the xz edge (edgeGrid,packThreads,12)................................
		//...Pack the yz edge (edgeGrid,packThreads,16)................................
		dvc_PackDist(edgeGrid,packThreads,8,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,f_even,N);
		//...Pack the yZ edge (edgeGrid,packThreads,18)................................
		dvc_PackDist(edgeGrid,packThreads,9,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,f_even,N);
		//...Pack the Yz edge (edgeGrid,packThreads,17)................................
		dvc_PackDist(edgeGrid,packThreads,8,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,f_odd,N);
		//...Pack the YZ edge (edgeGrid,packThreads,15)................................
		dvc_PackDist(edgeGrid,packThreads,7,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,f_odd,N);
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
		dvc_DensityStreamD3Q7(nBlocks, nthreads, S,
				ID, Den, Copy, Phi, ColorGrad, Velocity, beta, Nx, Ny, Nz, pBC);
		//*************************************************************************

		//*************************************************************************
		// 		Swap the distributions for momentum transport
		//*************************************************************************
		dvc_SwapD3Q19(nBlocks, nthreads, S, ID, f_even, f_odd, Nx, Ny, Nz);
		//*************************************************************************

		//...................................................................................
		// Wait for completion of D3Q19 communication
		MPI_Waitall(18,req1,stat1);
		MPI_Waitall(18,req2,stat2);
		//...................................................................................
		// Unpack the distributions on the device
		//...................................................................................
		//...Map recieve list for the X face: q=2,8,10,12,13 .................................
		dvc_UnpackDist(faceGrid,packThreads,0,-1,0,0,dvcRecvList_X,0,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,3,-1,-1,0,dvcRecvList_X,recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,4,-1,1,0,dvcRecvList_X,2*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,5,-1,0,-1,dvcRecvList_X,3*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,6,-1,0,1,dvcRecvList_X,4*recvCount_X,recvCount_X,recvbuf_X,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the x face: q=1,7,9,11,13..................................
		dvc_UnpackDist(faceGrid,packThreads,1,1,0,0,dvcRecvList_x,0,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,4,1,1,0,dvcRecvList_x,recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,5,1,-1,0,dvcRecvList_x,2*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,6,1,0,1,dvcRecvList_x,3*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,7,1,0,-1,dvcRecvList_x,4*recvCount_x,recvCount_x,recvbuf_x,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
		dvc_UnpackDist(faceGrid,packThreads,1,0,-1,0,dvcRecvList_Y,0,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,3,-1,-1,0,dvcRecvList_Y,recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,5,1,-1,0,dvcRecvList_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,7,0,-1,-1,dvcRecvList_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,8,0,-1,1,dvcRecvList_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,f_odd,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
		dvc_UnpackDist(faceGrid,packThreads,2,0,1,0,dvcRecvList_y,0,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,4,1,1,0,dvcRecvList_y,recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,4,-1,1,0,dvcRecvList_y,2*recvCount_y,recvCount_y,recvbuf_y,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,8,0,1,1,dvcRecvList_y,3*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,9,0,1,-1,dvcRecvList_y,4*recvCount_y,recvCount_y,recvbuf_y,f_even,Nx,Ny,Nz);
		//...................................................................................
		//...Map recieve list for the z face<<<faceGrid,packThreads,6,12,13,16,17)..............................................
		dvc_UnpackDist(faceGrid,packThreads,2,0,0,-1,dvcRecvList_Z,0,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,5,-1,0,-1,dvcRecvList_Z,recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,7,1,0,-1,dvcRecvList_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,7,0,-1,-1,dvcRecvList_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,9,0,1,-1,dvcRecvList_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,f_even,Nx,Ny,Nz);
		//...Map recieve list for the Z face<<<faceGrid,packThreads,5,11,14,15,18)..............................................
		dvc_UnpackDist(faceGrid,packThreads,3,0,0,1,dvcRecvList_z,0,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,6,1,0,1,dvcRecvList_z,recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,6,-1,0,1,dvcRecvList_z,2*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,8,0,1,1,dvcRecvList_z,3*recvCount_z,recvCount_z,recvbuf_z,f_even,Nx,Ny,Nz);
		dvc_UnpackDist(faceGrid,packThreads,8,0,-1,1,dvcRecvList_z,4*recvCount_z,recvCount_z,recvbuf_z,f_odd,Nx,Ny,Nz);
		//..................................................................................
		//...Map recieve list for the xy edge <<<edgeGrid,packThreads,8)................................
		dvc_UnpackDist(edgeGrid,packThreads,3,-1,-1,0,dvcRecvList_XY,0,recvCount_XY,recvbuf_XY,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xy edge <<<edgeGrid,packThreads,9)................................
		dvc_UnpackDist(edgeGrid,packThreads,5,1,-1,0,dvcRecvList_xY,0,recvCount_xY,recvbuf_xY,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xY edge <<<edgeGrid,packThreads,10)................................
		dvc_UnpackDist(edgeGrid,packThreads,4,-1,1,0,dvcRecvList_Xy,0,recvCount_Xy,recvbuf_Xy,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the XY edge <<<edgeGrid,packThreads,7)................................
		dvc_UnpackDist(edgeGrid,packThreads,4,1,1,0,dvcRecvList_xy,0,recvCount_xy,recvbuf_xy,f_even,Nx,Ny,Nz);
		//...Map recieve list for the xz edge <<<edgeGrid,packThreads,12)................................
		dvc_UnpackDist(edgeGrid,packThreads,5,-1,0,-1,dvcRecvList_XZ,0,recvCount_XZ,recvbuf_XZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the xZ edge <<<edgeGrid,packThreads,14)................................
		dvc_UnpackDist(edgeGrid,packThreads,6,-1,0,1,dvcRecvList_Xz,0,recvCount_Xz,recvbuf_Xz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Xz edge <<<edgeGrid,packThreads,13)................................
		dvc_UnpackDist(edgeGrid,packThreads,7,1,0,-1,dvcRecvList_xZ,0,recvCount_xZ,recvbuf_xZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the XZ edge <<<edgeGrid,packThreads,11)................................
		dvc_UnpackDist(edgeGrid,packThreads,6,1,0,1,dvcRecvList_xz,0,recvCount_xz,recvbuf_xz,f_even,Nx,Ny,Nz);
		//...Map recieve list for the yz edge <<<edgeGrid,packThreads,16)................................
		dvc_UnpackDist(edgeGrid,packThreads,7,0,-1,-1,dvcRecvList_YZ,0,recvCount_YZ,recvbuf_YZ,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the yZ edge <<<edgeGrid,packThreads,18)................................
		dvc_UnpackDist(edgeGrid,packThreads,8,0,-1,1,dvcRecvList_Yz,0,recvCount_Yz,recvbuf_Yz,f_odd,Nx,Ny,Nz);
		//...Map recieve list for the Yz edge <<<edgeGrid,packThreads,17)................................
		dvc_UnpackDist(edgeGrid,packThreads,9,0,1,-1,dvcRecvList_yZ,0,recvCount_yZ,recvbuf_yZ,f_even,Nx,Ny,Nz);
		//...Map recieve list for the YZ edge <<<edgeGrid,packThreads,15)................................
		dvc_UnpackDist(edgeGrid,packThreads,8,0,1,1,dvcRecvList_yz,0,recvCount_yz,recvbuf_yz,f_even,Nx,Ny,Nz);
		//...................................................................................

		//...................................................................................
		dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_x,recvCount_x,recvbuf_x,2,Den,N);
		dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_y,recvCount_y,recvbuf_y,2,Den,N);
		dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_z,recvCount_z,recvbuf_z,2,Den,N);
		dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_X,recvCount_X,recvbuf_X,2,Den,N);
		dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_Y,recvCount_Y,recvbuf_Y,2,Den,N);
		dvc_PackDenD3Q7(faceGrid,packThreads,dvcRecvList_Z,recvCount_Z,recvbuf_Z,2,Den,N);
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
		dvc_UnpackDenD3Q7(faceGrid,packThreads,dvcSendList_x,sendCount_x,sendbuf_x,2,Den,N);
		dvc_UnpackDenD3Q7(faceGrid,packThreads,dvcSendList_y,sendCount_y,sendbuf_y,2,Den,N);
		dvc_UnpackDenD3Q7(faceGrid,packThreads,dvcSendList_z,sendCount_z,sendbuf_z,2,Den,N);
		dvc_UnpackDenD3Q7(faceGrid,packThreads,dvcSendList_X,sendCount_X,sendbuf_X,2,Den,N);
		dvc_UnpackDenD3Q7(faceGrid,packThreads,dvcSendList_Y,sendCount_Y,sendbuf_Y,2,Den,N);
		dvc_UnpackDenD3Q7(faceGrid,packThreads,dvcSendList_Z,sendCount_Z,sendbuf_Z,2,Den,N);
		//...................................................................................

		//*************************************************************************
		// 		Compute the phase indicator field and reset Copy, Den
		//*************************************************************************
		dvc_ComputePhi(nBlocks, nthreads, S,ID, Phi, Copy, Den, N);
		//*************************************************************************

		//...................................................................................
		dvc_PackValues(faceGrid, packThreads, dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_xy, sendCount_xy,sendbuf_xy, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_xY, sendCount_xY,sendbuf_xY, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_Xy, sendCount_Xy,sendbuf_Xy, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_XY, sendCount_XY,sendbuf_XY, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_xz, sendCount_xz,sendbuf_xz, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_xZ, sendCount_xZ,sendbuf_xZ, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_Xz, sendCount_Xz,sendbuf_Xz, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_XZ, sendCount_XZ,sendbuf_XZ, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_yz, sendCount_yz,sendbuf_yz, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_yZ, sendCount_yZ,sendbuf_yZ, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_Yz, sendCount_Yz,sendbuf_Yz, Phi, N);
		dvc_PackValues(faceGrid, packThreads, dvcSendList_YZ, sendCount_YZ,sendbuf_YZ, Phi, N);
		//...................................................................................
		// Send / Recv all the phase indcator field values
		//...................................................................................
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
		// Wait for completion of Indicator Field communication
		//...................................................................................
		MPI_Waitall(18,req1,stat1);
		MPI_Waitall(18,req2,stat2);
		//...................................................................................
		//...................................................................................
/*		dvc_UnpackValues(faceGrid, packThreads, dvcSendList_x, sendCount_x,sendbuf_x, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads, dvcSendList_y, sendCount_y,sendbuf_y, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads, dvcSendList_z, sendCount_z,sendbuf_z, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads, dvcSendList_X, sendCount_X,sendbuf_X, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads, dvcSendList_Y, sendCount_Y,sendbuf_Y, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads, dvcSendList_Z, sendCount_Z,sendbuf_Z, Phi, N);
*/		
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_x, recvCount_x,recvbuf_x, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_y, recvCount_y,recvbuf_y, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_z, recvCount_z,recvbuf_z, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_X, recvCount_X,recvbuf_X, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Y, recvCount_Y,recvbuf_Y, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Z, recvCount_Z,recvbuf_Z, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xy, recvCount_xy,recvbuf_xy, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xY, recvCount_xY,recvbuf_xY, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_XY, recvCount_XY,recvbuf_XY, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xz, recvCount_xz,recvbuf_xz, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_yz, recvCount_yz,recvbuf_yz, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, Phi, N);
		dvc_UnpackValues(faceGrid, packThreads,dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, Phi, N);
		//...................................................................................
		MPI_Barrier(comm);
		// Iteration completed!
		timestep++;
		//...................................................................
		
		if (timestep%1000 == 0){
			//...........................................................................
			// Copy the phase from the GPU -> CPU
			//...........................................................................
			dvc_Barrier();
			dvc_CopyToHost(Phase.data,Phi,N*sizeof(double));
			MPI_Barrier(comm);
			//...........................................................................
			// Compute areas using porous medium marching cubes algorithm
			// McClure, Adalsteinsson, et al. (2007)
			//...........................................................................
			awn = aws = ans = lwns = 0.0;
			nwp_volume = 0.0;
			As = 0.0;

			for (c=0;c<ncubes;c++){
				// Get cube from the list
				i = cubeList(0,c);
				j = cubeList(1,c);
				k = cubeList(2,c);
	
				for (p=0;p<8;p++){
					if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 
						&&  SignDist(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
						nwp_volume += 0.125;
					}
				}
				
				// Run PMMC
				n_local_sol_tris = 0;
				n_local_sol_pts = 0;
				n_local_nws_pts = 0;
				
				n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
				n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
				
				n_nw_tris_beg = 0;// n_nw_tris;
				n_ns_tris_beg = 0;//n_ns_tris;
				n_ws_tris_beg = 0;//n_ws_tris;
				
				// if there is a solid phase interface in the grid cell
				if (Interface(SignDist,solid_isovalue,i,j,k) == 1){
					
					/////////////////////////////////////////
					/// CONSTRUCT THE LOCAL SOLID SURFACE ///
					/////////////////////////////////////////
					
					// find the local solid surface
					SOL_SURF(SignDist,0.0,Phase,fluid_isovalue,i,j,k, Nx,Ny,Nz,local_sol_pts,n_local_sol_pts,
							 local_sol_tris,n_local_sol_tris,values);
					
					/////////////////////////////////////////
					//////// TRIM THE SOLID SURFACE /////////
					/////////////////////////////////////////
/*					TRIM(local_sol_pts, n_local_sol_pts, fluid_isovalue,local_sol_tris, n_local_sol_tris,
						 ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts, n_ws_pts,
						 ws_tris, n_ws_tris, values, local_nws_pts, n_local_nws_pts,
						 Phase, SignDist, i, j, k, newton_steps);
*/					
					TRIM(local_sol_pts, n_local_sol_pts, fluid_isovalue,local_sol_tris, n_local_sol_tris,
						 ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts, n_ws_pts,
						 ws_tris, n_ws_tris, values, local_nws_pts, n_local_nws_pts);
					
					/////////////////////////////////////////
					//////// WRITE COMMON LINE POINTS ///////
					////////      TO MAIN ARRAYS      ///////
					/////////////////////////////////////////
					map = n_nws_pts;
					for (p=0; p < n_local_nws_pts; p++){
						nws_pts(n_nws_pts++) = local_nws_pts(p);
					}
					for (q=0; q < n_local_nws_pts-1; q++){
						nws_seg(0,n_nws_seg) = map+q;
						nws_seg(1,n_nws_seg) = map+q+1;
						n_nws_seg++;
					}
					
					/////////////////////////////////////////
					////// CONSTRUCT THE nw SURFACE /////////
					/////////////////////////////////////////
					if ( n_local_nws_pts > 0){
						EDGE(Phase, fluid_isovalue, SignDist, i,j,k, Nx, Ny, Nz, nw_pts, n_nw_pts, nw_tris, n_nw_tris,
							 local_nws_pts, n_local_nws_pts);
					}
					else {
						MC(Phase, fluid_isovalue, SignDist, i,j,k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
					}
				}
				
				/////////////////////////////////////////
				////// CONSTRUCT THE nw SURFACE /////////
				/////////////////////////////////////////
				
				else if (Fluid_Interface(Phase,SignDist,fluid_isovalue,i,j,k) == 1){
					MC(Phase, fluid_isovalue, SignDist, i,j,k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
				}
				//******END OF BLOB PMMC*********************************************

				//*******************************************************************
				// Compute the Interfacial Areas, Common Line length for blob p
				// nw surface
				double temp;
				for (r=0;r<n_nw_tris;r++){
					A = nw_pts(nw_tris(0,r));
					B = nw_pts(nw_tris(1,r));
					C = nw_pts(nw_tris(2,r));
					// Compute length of sides (assume dx=dy=dz)
					s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
					s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
					s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
					s = 0.5*(s1+s2+s3);
					temp = s*(s-s1)*(s-s2)*(s-s3);
					if (temp > 0.0) awn += sqrt(temp);
					
				}
				for (r=0;r<n_ns_tris;r++){
					A = ns_pts(ns_tris(0,r));
					B = ns_pts(ns_tris(1,r));
					C = ns_pts(ns_tris(2,r));
					// Compute length of sides (assume dx=dy=dz)
					s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
					s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
					s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
					s = 0.5*(s1+s2+s3);
					//ans=ans+sqrt(s*(s-s1)*(s-s2)*(s-s3));
					temp = s*(s-s1)*(s-s2)*(s-s3);
					if (temp > 0.0) ans += sqrt(temp);
				}
				for (r=0;r<n_ws_tris;r++){
					A = ws_pts(ws_tris(0,r));
					B = ws_pts(ws_tris(1,r));
					C = ws_pts(ws_tris(2,r));
					// Compute length of sides (assume dx=dy=dz)
					s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
					s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
					s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
					s = 0.5*(s1+s2+s3);
					//aws=aws+sqrt(s*(s-s1)*(s-s2)*(s-s3));
					temp = s*(s-s1)*(s-s2)*(s-s3);
					if (temp > 0.0) aws += sqrt(temp);
				}
				for (r=0;r<n_local_sol_tris;r++){
					A = local_sol_pts(local_sol_tris(0,r));
					B = local_sol_pts(local_sol_tris(1,r));
					C = local_sol_pts(local_sol_tris(2,r));
					// Compute length of sides (assume dx=dy=dz)
					s1 = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
					s2 = sqrt((A.x-C.x)*(A.x-C.x)+(A.y-C.y)*(A.y-C.y)+(A.z-C.z)*(A.z-C.z));
					s3 = sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y)+(B.z-C.z)*(B.z-C.z));
					s = 0.5*(s1+s2+s3);
					//aws=aws+sqrt(s*(s-s1)*(s-s2)*(s-s3));
					temp = s*(s-s1)*(s-s2)*(s-s3);
					if (temp > 0.0) As += sqrt(temp);
				}
				for (p=0; p < n_local_nws_pts-1; p++){
					// Extract the line segment
					A = local_nws_pts(p);
					B = local_nws_pts(p+1);
					// Compute the length of the segment
					s = sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
					// Add the length to the common line 
					lwns += s;
				}
				//*******************************************************************
				// Reset the triangle counts to zero
				n_nw_pts=0,n_ns_pts=0,n_ws_pts=0,n_nws_pts=0, map=0;
				n_nw_tris=0, n_ns_tris=0, n_ws_tris=0, n_nws_seg=0;
				
				n_nw_tris_beg = 0;// n_nw_tris;
			//	n_ns_tris_beg = 0;//n_ns_tris;
			//	n_ws_tris_beg = 0;//n_ws_tris;
			//	n_nws_seg_beg = n_nws_seg;
				//*******************************************************************
			}
			//...........................................................................
			MPI_Barrier(comm);
			MPI_Allreduce(&nwp_volume,&nwp_volume_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&awn,&awn_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&ans,&ans_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&aws,&aws_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&lwns,&lwns_global,1,MPI_DOUBLE,MPI_SUM,comm);
			MPI_Allreduce(&As,&As_global,1,MPI_DOUBLE,MPI_SUM,comm);

			MPI_Barrier(comm);
			//.........................................................................
			// Compute the change in the total surface energy based on the defined interval
			// See McClure, Prins and Miller (2013) 
			//.........................................................................
			dAwn += awn_global;
			dAns += ans_global;
			dEs = 6.01603*alpha*(dAwn + 1.05332*Ps*dAns);
			dAwn = -awn_global;		// Get ready for the next analysis interval
			dAns = -ans_global;
			//.........................................................................
			if (rank==0){
/*				printf("-------------------------------- \n");
				printf("Timestep = %i \n", timestep);
				printf("NWP volume = %f \n", nwp_volume_global);
				printf("Area wn = %f \n", awn_global);
				printf("Area ns = %f \n", ans_global);
				printf("Area ws = %f \n", aws_global);
				printf("Change in surface energy = %f \n", dEs);
				printf("-------------------------------- \n");	
*/
				printf("%i %f %f %f %f %f %f %f \n",timestep,dEs,nwp_volume_global,
						awn_global,ans_global,aws_global, As_global, lwns_global);		

			}
		}
	}
	//************************************************************************/
	dvc_Barrier();
	MPI_Barrier(comm);
	stoptime = MPI_Wtime();
	if (rank==0) printf("-------------------------------------------------------------------\n");
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
	
	//************************************************************************/
	// Write out the phase indicator field 
	//************************************************************************/
	sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
	//	printf("Local File Name =  %s \n",LocalRankFilename);
	dvc_CopyToHost(HostPhi,Phi,N*sizeof(double));
	dvc_CopyToHost(Phase.data,Phi,N*sizeof(double));
	
	FILE *PHASE;
	PHASE = fopen(LocalRankFilename,"wb");
	fwrite(Phase.data,8,N,PHASE);
	fclose(PHASE);

/*	double *DensityValues;
	DensityValues = new double [2*N];
	dvc_CopyToHost(DensityValues,Copy,2*N*sizeof(double));
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
