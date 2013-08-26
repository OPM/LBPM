#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cuda.h>
//#include <cutil.h>

using namespace std;

//*************************************************************************
extern "C" void dvc_InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx,
							  int Ny, int Nz, int nblocks, int nthreads, int S);
//*************************************************************************
extern "C" void dvc_SwapD3Q19(char *ID, double *f_even, double *f_odd, int Nx,
							  int Ny, int Nz, int nblocks, int nthreads, int S);
//*************************************************************************
extern "C" void dvc_MRT(char *ID, double *f_even, double *f_odd, double rlxA, double rlxB, double Fx, double Fy, double Fz,
		int Nx, int Ny, int Nz, int nblocks, int nthreads, int S);
//*************************************************************************

void Write_Out(double *array, int Nx, int Ny, int Nz){
	int value;
	FILE *output;
	output = fopen("dist.list","w");
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				int index = k*Nx*Ny+j*Nx+i;
				value = int(array[index]);
				fprintf(output, "| %i",value);
			}
			fprintf(output, " | \n");
		}
		fprintf(output,"************************************** \n");	
	}
	fclose(output);
}

//**************************************************************************
// MRT implementation of the LBM using CUDA
//**************************************************************************

int main(void)
{

	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	int device = 1;
	printf("Number of devices = %i \n", deviceCount);
	printf("Current device is = %i \n", device);
	cudaSetDevice(device);
	
	// BGK Model parameters
	string FILENAME;	
	unsigned int nBlocks, nthreads;
	int timestepMax, interval;
	double tau,Fx,Fy,Fz,tol;
	// Domain variables
	int Nx,Ny,Nz;

	ifstream input("MRT.in");
	input >> FILENAME;		// name of the input file
	input >> Nz;			// number of nodes (x,y,z)
	input >> nBlocks;				
	input >> nthreads;				
	input >> tau;				// relaxation time 
	input >> Fx;			// External force components (x,y,z)
	input >> Fy;
	input >> Fz;
	input >> timestepMax;			// max no. of timesteps
	input >> interval;			// error interval
	input >> tol;				// error tolerance
	
	double rlx_setA = 1.f/tau;
	double rlx_setB = 8.f*(2.f-rlx_setA)/(8.f-rlx_setA);
	
	printf("tau = %f \n", tau);
	printf("Set A = %f \n", rlx_setA);
	printf("Set B = %f \n", rlx_setB);
	printf("Force(x) = %f \n", Fx);
	printf("Force(y) = %f \n", Fy);
	printf("Force(z) = %f \n", Fz);

	Nx = Ny = Nz;	// Cubic domain
	
	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);
	
//	unsigned int nBlocks = 32;
//	int nthreads = 128;
	int S = N/nthreads/nBlocks;
	
//	unsigned int nBlocks = N/nthreads + (N%nthreads == 0?0:1);
	dim3 grid(nBlocks,1,1);
		
	printf("Number of blocks = %i \n", nBlocks);
	printf("Threads per block = %i \n", nthreads);
	printf("Sweeps per thread = %i \n", S);
	printf("Number of nodes per side = %i \n", Nx);
	printf("Total Number of nodes = %i \n", N);
	
	//.......................................................................
	printf("Read input media... \n");
	// .......... READ THE INPUT FILE .......................................
	int n;
	char value;
	char *id;
	id = new char[N];	
	int sum = 0;
	double porosity;
	ifstream PM(FILENAME.c_str(),ios::binary);
	for (int k=0;k<Nz;k++){
		for (int j=0;j<Ny;j++){
			for (int i=0;i<Nx;i++){
				PM.read((char *) (&value), sizeof(value));
				n = k*Nx*Ny+j*Nx+i;
				id[n] = value;
				if (value > 0) sum++;
			}
		}
	}
	PM.close();
	printf("File porosity = %f\n", double(sum)/N);
	//.......................................................................
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
	//...........................................................................

	//...........................................................................
//	cudaHostAlloc(&fa,dist_mem_size,cudaHostAllocPortable);
//	cudaHostAlloc(&fb,dist_mem_size,cudaHostAllocPortable);
//	cudaHostRegister(fa,dist_mem_size,cudaHostRegisterPortable);
//	cudaHostRegister(fb,dist_mem_size,cudaHostRegisterPortable);
//	cudaHostRegister(id,N*sizeof(char),cudaHostAllocPortable);

	printf("Setting the distributions, size = : %i\n", N);
	//...........................................................................
//	INITIALIZE <<< grid, nthreads >>>  (ID, f_even, f_odd, Nx, Ny, Nz, S);
	//...........................................................................
	dvc_InitD3Q19(ID,f_even,f_odd,Nx,Ny,Nz,nBlocks,nthreads,S);
	//*************************************************************************

	int timestep = 0;
	printf("No. of timesteps: %i \n", timestepMax);
	
	//.......create a stream for the LB calculation.......
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	
	//.......create and start timer............
	cudaEvent_t start, stop; 
	float time; 
 
	cudaEventCreate(&start); 
	cudaEventCreate(&stop); 
	cudaEventRecord( start, 0 ); 
	//.........................................
	
	//************ MAIN ITERATION LOOP ***************************************/
	while (timestep < timestepMax){
	
		//...................................................................
		//........ Execute the swap kernel (device) .........................
//		SWAP <<< grid, nthreads >>> (ID, f_even, f_odd, Nx, Ny, Nz, S);
		//...................................................................
		dvc_SwapD3Q19(ID,f_even,f_odd,Nx,Ny,Nz,nBlocks,nthreads,S);

		//........ Execute the collision kernel (device) ....................
//		MRT <<< grid, nthreads >>> (ID, f_even, f_odd, Nx, Ny, Nz, S,
//									rlx_setA, rlx_setB, Fx, Fy, Fz);
		//............................................................
		dvc_MRT(ID, f_even, f_odd, rlx_setA, rlx_setB, Fx, Fy, Fz,Nx,Ny,Nz,nBlocks,nthreads,S);
		// Iteration completed!

		timestep++;
		//...................................................................
		
	}
	//************************************************************************/
	
	cudaThreadSynchronize();
	//.......... stop and destroy timer.............................
	cudaEventRecord( stop, stream); 
	cudaEventSynchronize( stop ); 
 
	cudaEventElapsedTime( &time, start, stop ); 
	printf("CPU time = %f \n", time);
	
	float MLUPS = 0.001*float(Nx*Ny*Nz)*timestep/time;
	printf("MLUPS = %f \n", MLUPS);

	cudaStreamDestroy(stream);
	cudaEventDestroy( start ); 
	cudaEventDestroy( stop ); 
	//..............................................................
	
	//..............................................................
	//.........Compute the velocity and copy result to host ........
	double *velocity;
	velocity = new double[3*N];
	//......................device distributions....................................
	double *vel;
	//..............................................................................
	cudaMalloc((void **) &vel, 3*dist_mem_size);	// Allocate device memory
	//..............................................................................
//	Compute_VELOCITY <<< grid, nthreads >>>  (ID, f_even, f_odd, vel, Nx, Ny, Nz, S);
	//..............................................................................
	cudaMemcpy(velocity, vel, 3*dist_mem_size, cudaMemcpyDeviceToHost);
	//..............................................................................

	//............................................................	
	//....Write the z-velocity to test poiseuille flow............
	double vz,vz_avg;	
	vz_avg = 0.0;

	FILE *output;
	output = fopen("velocity.out","w");
	for (int k=0; k<1; k++){
		for (int j=0; j<1; j++){
			for (int i=0; i<Nx; i++){
				int n = k*Nx*Ny+j*Nx+i;
				//.....print value........
				vz = velocity[2*N+n];
				vz_avg += vz;
				fprintf(output, " %e",vz);
			}
		}
	}
	fclose(output);
	
	vz = vz_avg/double(sum);
	printf("Average Velocity = %e\n", vz);


	// cleanup	
	cudaFree(f_even);	cudaFree(f_odd);	cudaFree(vel);	cudaFree(ID);
	free (velocity);	free(id);
	
}
