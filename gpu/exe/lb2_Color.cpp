#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

//*************************************************************************
// Functions defined in Color.cu
//*************************************************************************
extern "C" void dvc_InitDenColor( int nblocks, int nthreads, int S,
		char *ID, double *Den, double *Phi, double das, double dbs, int N);
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

using namespace std;


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

	int rank = 0;
	int nprocs =1;
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
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
		printf("das = %f \n", beta);
		printf("dbs = %f \n", beta);
		printf("Force(x) = %f \n", Fx);
		printf("Force(y) = %f \n", Fy);
		printf("Force(z) = %f \n", Fz);
		printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
		printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
		printf("********************************************************\n");

	}

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

	//...........device phase ID.................................................
	if (rank==0)	printf ("Copying phase ID to device \n");
	char *ID;
	dvc_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
	// Copy to the device
	dvc_CopyToDevice(ID, id, N);
	//...........................................................................

	if (rank==0)	printf ("Allocating distributions \n");
	//......................device distributions.................................
	double *f_even,*f_odd;
	//...........................................................................
	dvc_AllocateDeviceMemory((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
	dvc_AllocateDeviceMemory((void **) &f_odd, 9*dist_mem_size);	// Allocate device memory
	//...........................................................................
	//...........................................................................
	//				MAIN  VARIABLES ALLOCATED HERE
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
	if (rank==0)	printf("Setting the distributions, size = : %i\n", N);
	//...........................................................................
	dvc_InitD3Q19(nBlocks, nthreads, S, ID, f_even, f_odd, Nx, Ny, Nz);
	dvc_InitDenColor(nBlocks, nthreads, S, ID, Den, Phi,  das, dbs, N);
	//...........................................................................
	dvc_ComputePhi(nBlocks, nthreads, S,ID, Phi, Copy, Den, N);
	//...........................................................................
	
	//...........................................................................
	// Grids used to pack faces on the GPU for MPI
	int faceGrid,edgeGrid,packThreads;
	packThreads=512;
	edgeGrid=1;
	faceGrid=Nx*Ny/packThreads;


	int timestep = 0;
	if (rank==0) printf("********************************************************\n");
	if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);

	//.......create a stream for the LB calculation.......
//	cudaStream_t stream;
//	cudaStreamCreate(&stream);

	//.......create and start timer............
	double start,stop;
	double walltime;
	start = clock();


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
		dvc_ColorCollide(nBlocks, nthreads, S,
				ID, f_even, f_odd, ColorGrad, Velocity,
				rlxA, rlxB,alpha, beta, Fx, Fy, Fz, Nx, Ny, Nz, pBC);
		//*************************************************************************

		//*************************************************************************
		// 		Carry out the density streaming step for mass transport
		//*************************************************************************
		dvc_DensityStreamD3Q7(nBlocks, nthreads, S,
				ID, Den, Copy, Phi, ColorGrad, Velocity,beta, Nx, Ny, Nz, pBC);
		//*************************************************************************
		
		//*************************************************************************
		// 		Swap the distributions for momentum transport
		//*************************************************************************
		dvc_SwapD3Q19(nBlocks, nthreads, S, ID, f_even, f_odd, Nx, Ny, Nz);
		//*************************************************************************

		//*************************************************************************
		// 		Compute the phase indicator field and reset Copy, Den
		//*************************************************************************
		dvc_ComputePhi(nBlocks, nthreads, S,ID, Phi, Copy, Den, N);
		//*************************************************************************

		// Iteration completed!
		timestep++;
		
		//...................................................................
	}
	//************************************************************************/
	dvc_Barrier();
	stop = clock();

//	cout << "CPU time: " << (stoptime - starttime) << " seconds" << endl;
	walltime = (stop - start)/CLOCKS_PER_SEC;
//	cout << "Lattice update rate: "<< double(Nx*Ny*Nz*timestep)/cputime/1000000 <<  " MLUPS" << endl;
	double MLUPS = double(Nx*Ny*Nz*timestep)/walltime/1000000;
	if (rank==0) printf("********************************************************\n");
	if (rank==0) printf("CPU time = %f \n", walltime);
	if (rank==0) printf("Lattice update rate (per core)= %f MLUPS \n", MLUPS);
	MLUPS *= nprocs;
	if (rank==0) printf("Lattice update rate (total)= %f MLUPS \n", MLUPS);
	if (rank==0) printf("********************************************************\n");
	
	//************************************************************************/
	// Write out the phase indicator field 
	//************************************************************************/
	sprintf(LocalRankFilename,"%s%s","Phase.",LocalRankString);
	//	printf("Local File Name =  %s \n",LocalRankFilename);
	double *phiOut;
	phiOut = new double[N];
	dvc_CopyToHost(phiOut,Phi,N*sizeof(double));
		
	FILE *PHASE;
	PHASE = fopen(LocalRankFilename,"wb");
	fwrite(phiOut,8,N,PHASE);
	fclose(PHASE);
	//************************************************************************/

}
