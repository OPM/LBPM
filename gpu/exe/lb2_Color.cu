#ifdef useMPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cuda.h>

using namespace std;
//*************************************************************************
// HokieSpeed
//nvcc -Xcompiler -fopenmp -lgomp -O3 -arch sm_20 -o hybridATLKR lb2_ATLKR_hybrid.cu
// -I$VT_MPI_INC -L$VT_MPI_LIB -lmpi
//*************************************************************************

//*************************************************************************
// Implementation of Two-Phase Immiscible LBM using CUDA
//*************************************************************************

//*************************************************************************
extern "C" void dvc_InitD3Q19(int nblocks, int nthreads, int S,
		char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);
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
extern "C" void dvc_AllocateDeviceMemory(void** address, size_t size);
//*************************************************************************
extern "C" void dvc_CopyToDevice(void* dest, void* source, size_t size);
//*************************************************************************
extern "C" void dvc_Barrier();
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

int main(int argc, char *argv[])
{
	
	//********** Initialize MPI ****************
	int numprocs,rank;
#ifdef useMPI
	MPI_Status stat;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm,&numprocs);
	MPI_Comm_rank(comm,&rank);
#else
    MPI_Comm comm = MPI_COMM_WORLD;
	numprocs = 1;
	rank = 0;
#endif
	//******************************************
	
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
	}
#ifdef useMPI
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
	//.................................................
	MPI_Barrier(comm);
	// **************************************************************
#endif
	
	double rlxA = 1.f/tau;
	double rlxB = 8.f*(2.f-rlxA)/(8.f-rlxA);

	if (pBC && rank == 0){
		printf("Assigning presusre boundary conditions \n");
		printf("Inlet density = %f \n", din);
		printf("Outlet density = %f \n", dout);
	}

	if (rank==0){
		printf("....Parameters................\n");
		printf("tau = %f \n", tau);
		printf("alpha = %f \n", alpha);
		printf("beta = %f \n", beta);
		printf("das = %f \n", das);
		printf("dbs = %f \n", dbs);
		printf("Force(x) = %f \n", Fx);
		printf("Force(y) = %f \n", Fy);
		printf("Force(z) = %f \n", Fz);
		printf("Nz = %i \n", Nz);
		printf("timestepMax = %i \n", timestepMax);
		printf("...............................\n");
	}
	
	// Identical cubic sub-domains
	Nx = Ny = Nz;// = 16*s;		// Cubic domain
	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

//	unsigned int nBlocks = 32;
//	int nthreads = 128;
	int S = N/nthreads/nBlocks;
	if (nBlocks*nthreads*S < N)	S++;
//	int S = 1;
	
//	unsigned int nBlocks = N/nthreads + (N%nthreads == 0?0:1);
//	dim3 grid(nBlocks,1,1);
	if (rank==1){
		printf("Number of blocks = %i \n", nBlocks);
		printf("Threads per block = %i \n", nthreads);
		printf("Sweeps per thread = %i \n", S);
		printf("Number of nodes per side = %i \n", Nx);
		printf("Total Number of nodes = %i \n", N);
		printf("...............................\n");
	}
	
	//.......................................................................
	// .......... READ THE INPUT FILE .......................................
	int n;
	char value;
	char *id;
	id = new char[N];	
	int sum = 0;
	// RANK 0 READS THE INPUT FILE	
	if (rank==0){
		printf("Read input media... \n");
		ifstream PM(FILENAME.c_str(),ios::binary);
		for (int k=0;k<Nz;k++){
			for (int j=0;j<Ny;j++){
				for (int i=0;i<Nx;i++){
					PM.read((char *) (&value), sizeof(value));
					n = k*Nx*Ny+j*Nx+i;

					if (value>0){
						if (pBC) value=2; 	// Saturate with NWP
						if (k<8){
							value=1;
						}
					}

					id[n] = value;
					if (value > 0) sum++;
				}
			}
		}
		PM.close();
		printf("File porosity = %f\n", double(sum)/N);
	}
	//......... for pressure BC only............................
	// Void the first / last rows if pressure BC are to be used
	if (pBC){
		for (int k=0;k<Nz;k++){
			for (int j=0;j<Ny;j++){
				for (int i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					if (k<4) id[n] = 1;
					if (k>Nz-5) id[n] = 2;
				}
			}
			// Skip the non-boundary values
			if (k==4)	k=Nz-5;
		}
	}
#ifdef useMPI	//............................................................
	MPI_Barrier(comm);
	MPI_Bcast(&id[0],N,MPI_CHAR,0,comm);
	MPI_Barrier(comm);
#endif
	if (rank == 0) printf("Domain set.\n");
	//...........................................................................

	int SBC;
	int outlet = N-Nx*Ny;
	if (pBC){
		SBC = Nx*Ny/nthreads/nBlocks+1;
		printf("Number of sweeps for inlet / outlet: %i \n", SBC);
	}
	//...........................................................................
	
	//...........................................................................
	//...........device phase ID.................................................
	char *ID;
	cudaMalloc((void **) &ID, N);						// Allocate device memory
	// Copy to the device
	cudaMemcpy(ID, id, N, cudaMemcpyHostToDevice);
	//...........................................................................
	
	//......................device distributions.................................
	double *f_even,*f_odd;
	//...........................................................................
	cudaMalloc((void **) &f_even, 10*dist_mem_size);	// Allocate device memory
	cudaMalloc((void **) &f_odd, 9*dist_mem_size);		// Allocate device memory
//	f_even = new double[10*N];
//	f_odd = new double[9*N];
	//...........................................................................

	//...........................................................................
	//				MAIN  VARIABLES ALLOCATED HERE
	//...........................................................................
	double *Phi,*Den,*Copy;
	double *ColorGrad, *Velocity;
	//...........................................................................
	cudaMalloc((void **) &Phi, dist_mem_size);
	cudaMalloc((void **) &Den, 2*dist_mem_size);
	cudaMalloc((void **) &Copy, 2*dist_mem_size);
	cudaMalloc((void **) &Velocity, 3*dist_mem_size);
	cudaMalloc((void **) &ColorGrad, 3*dist_mem_size);
	//...........................................................................
	
	//...........................................................................
	if (rank==0)	printf("Setting the distributions, size = : %i\n", N);
	//...........................................................................
	dvc_InitD3Q19(nBlocks, nthreads, S, ID, f_even, f_odd, Nx, Ny, Nz);
	dvc_InitDenColor(nBlocks, nthreads, S, ID, Den, Phi,  das, dbs, N);
	//...........................................................................
	dvc_ComputePhi(nBlocks, nthreads, S,ID, Phi, Copy, Den, N);
	//...........................................................................

	int timestep;
//	double starttime,stoptime;
	if (rank==0)	printf("No. of timesteps: %i \n", timestepMax);
	timestep = 0;
	//.......create and start timer............
	cudaEvent_t start, stop;
	float time;
	//.......create a stream for the LB calculation.......
	cudaStream_t stream;
	cudaStreamCreate(&stream);

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord( start, 0 );
	//.........................................
	//************ MAIN TIMESTEP LOOP ***************************************/
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

		dvc_Barrier();
		timestep++;
		//.............................................................................
	}
	//************************************************************************/
	dvc_Barrier();
	//.......... stop and destroy timer.............................
	cudaEventRecord( stop, stream);
	cudaEventSynchronize( stop );

	cudaEventElapsedTime( &time, start, stop );
	printf("CPU time = %f \n", time);

	float MLUPS = 0.001*float(Nx*Ny*Nz)*timestep/time;
	printf("MLUPS = %f \n", MLUPS);

	cudaEventDestroy( start );
	cudaEventDestroy( stop );

	double *Data;
	Data = new double[3*N];

	cudaMemcpy(Data, Phi, dist_mem_size, cudaMemcpyDeviceToHost);

	// Write out the Phase Indicator Field
	FILE *phase;
	phase = fopen("Phase.out","wb");
	fwrite(Data,8,N,phase);
	fclose(phase);

	//....................................................
	// Write out the pressure - (reuse Phi arrays since we're done with those)
//	ComputeDensity<<< grid, nthreads>>> (ID, f_even, f_odd, Phi, Nx, Ny, Nz, S);
//	cudaMemcpy(Data, Phi, dist_mem_size, cudaMemcpyDeviceToHost);
//	FILE *PRESSURE;
//	PRESSURE = fopen("Pressure.out","wb");
//	fwrite(Phi,8,N,PRESSURE);
//	fclose(PRESSURE);
	//....................................................

	// Write out the Color Gradient

	cudaMemcpy(Data, ColorGrad, 3*dist_mem_size, cudaMemcpyDeviceToHost);

	FILE *CG;
	CG = fopen("ColorGrad.out","wb");
	fwrite(Data,8,3*N,CG);
	fclose(CG);
	
	// Write out the Velocity
//	FILE *VEL;
//	VEL = fopen("Velocity.out","wb");
//	fwrite(Velocity,8,3*N,VEL);
//	fclose(VEL);

	// cleanup	
	cudaFree(ID);
	cudaFree(f_even);	cudaFree(f_odd);	
	cudaFree(Velocity);
	cudaFree(Phi);
	
	cudaFree (ColorGrad);
	cudaFree (Den);		cudaFree(Copy);
	cudaFree (Phi);
	free(id);
	
	//***********Finish up!*********************************
#ifdef useMPI
	MPI_Finalize();
#endif
	return 0;
	
}
