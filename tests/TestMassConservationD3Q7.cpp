// Unit test to test mass conservation for D3Q7 mass transport LBE
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
		printf("Running Unit Test for D3Q7 Mass Conservation	\n");
		printf("********************************************************\n");
	}
	
	double das=0.0;
	double dbs=1.0;
	double beta=0.95;
	bool pBC=false;
	int i,j,k,n;
	int Nx,Ny,Nz,N;
	Nx=Nz=Ny=100;
	N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);
	
	double *DenOriginal, *DenFinal;
	DenOriginal = new double [2*N];
	DenFinal = new double [2*N];
	
	double *Vel;
	Vel = new double [3*N];
	char *id;
	id = new char [N];
	
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (k<10)		id[n] = 0;
				else if (k>80)	id[n] = 0;
				else if (i<50){
					id[n]=1;
					Vel[n]=0.1;
					Vel[N+n]=0.1;
					Vel[2*N+n]=0.1;
				}
				else {
					id[n]=2;
					Vel[n]=0.1;
					Vel[N+n]=0.1;
					Vel[2*N+n]=0.1;
				}
			}
		}
	}
	

	
	//......................device distributions.................................
	double *A_even,*A_odd,*B_even,*B_odd;
	//...........................................................................
	AllocateDeviceMemory((void **) &A_even, 4*dist_mem_size);	// Allocate device memory
	AllocateDeviceMemory((void **) &A_odd, 3*dist_mem_size);	// Allocate device memory
	AllocateDeviceMemory((void **) &B_even, 4*dist_mem_size);	// Allocate device memory
	AllocateDeviceMemory((void **) &B_odd, 3*dist_mem_size);	// Allocate device memory
	//...........................................................................
	double *Phi,*Den;
	double *ColorGrad, *Velocity;
	//...........................................................................
	AllocateDeviceMemory((void **) &Phi, dist_mem_size);
	AllocateDeviceMemory((void **) &Velocity, 3*dist_mem_size);
	AllocateDeviceMemory((void **) &ColorGrad, 3*dist_mem_size);
	AllocateDeviceMemory((void **) &Den, 2*dist_mem_size);
	//...........device phase ID.................................................
	if (rank==0)	printf ("Copying phase ID to device \n");
	char *ID;
	AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
	// Copy to the device
	CopyToDevice(ID, id, N);
	CopyToDevice(Velocity, Vel, 3*N*sizeof(double));
	//...........................................................................
	
	InitDenColor(ID, Den, Phi, das, dbs, Nx, Ny, Nz);
	CopyToHost(DenOriginal,Den,2*N*sizeof(double));

	//......................................................................
	InitD3Q7(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
	InitD3Q7(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
	ComputePhi(ID, Phi, Den, N);
	ComputeColorGradient(ID,Phi,ColorGrad,Nx,Ny,Nz);
	//..................................................................................

	//*************************************************************************
	// 		Carry out the density streaming step for mass transport
	//*************************************************************************
	MassColorCollideD3Q7(ID, A_even, A_odd, B_even, B_odd, Den, Phi,
							ColorGrad, Velocity, beta, N, pBC);
	//*************************************************************************
	
	//..................................................................................
	ComputeDensityD3Q7(ID, A_even, A_odd, &Den[0], Nx, Ny, Nz);
	ComputeDensityD3Q7(ID, B_even, B_odd, &Den[N], Nx, Ny, Nz);
	//..................................................................................
	ComputePhi(ID, Phi, Den, N);
	ComputeColorGradient(ID,Phi,ColorGrad,Nx,Ny,Nz);
	//..................................................................................

	// Compare and make sure mass is conserved at every lattice site
	bool CleanCheck = true;
	double original,final;
	CopyToHost(DenFinal,Den,2*N*sizeof(double));
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (fabs(DenFinal[n] - DenOriginal[n]) > 1e-15){
					final = DenFinal[n];
					original = DenOriginal[n];
					if (id[n] == 0) printf("Solid phase! \n");
					if (id[n] == 1) printf("Wetting phase! \n");
					if (id[n] == 2) printf("Non-wetting phase! \n");							
					printf("Mass not conserved: WP density, site=%i,%i,%i, original = %f, final = %f \n",i,j,k,original,final);
					CleanCheck=false;
				}
				if (fabs(DenFinal[N+n] - DenOriginal[N+n]) > 1e-15){
					if (id[n] == 0) printf("Solid phase! \n");
					if (id[n] == 1) printf("Wetting phase! \n");
					if (id[n] == 2) printf("Non-wetting phase! \n");
					final = DenFinal[N+n];
					original = DenOriginal[N+n];
					printf("Mass not conserved: NWP density, site=%i,%i,%i, original = %f, final = %f \n",i,j,k,original,final);
					CleanCheck=false;
				}
			}
		}
	}
	if (rank==0) printf("Checking that the correct velocity is retained \n");
	// Swap convention is observed -- velocity is negative
	double *Aeven,*Aodd,*Beven,*Bodd;
	Aeven = new double[4*N];
	Aodd = new double[3*N];
	Beven = new double[4*N];
	Bodd = new double[3*N];
	CopyToHost(Aeven,A_even,4*dist_mem_size);
	CopyToHost(Aodd,A_odd,3*dist_mem_size);
	CopyToHost(Beven,B_even,4*dist_mem_size);
	CopyToHost(Bodd,B_odd,3*dist_mem_size);
	double rho,ux,uy,uz;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				rho = ux = uy = uz = 0.0;
				n = k*Nx*Ny + j*Nx + i;
				if (id[n] != 0){
					rho = Aeven[n]+Aeven[N+n]+Aeven[2*N+n]+Aeven[3*N+n]+Aodd[n]+Aodd[N+n]+Aodd[2*N+n]+
							Beven[n]+Beven[N+n]+Beven[2*N+n]+Beven[3*N+n]+Bodd[n]+Bodd[N+n]+Bodd[2*N+n];
					ux = Aeven[N+n] - Aodd[n] + Beven[N+n] - Bodd[n];
					uy = Aeven[2*N+n] - Aodd[N+n] + Beven[2*N+n] - Bodd[N+n];
					uz = Aeven[3*N+n] - Aodd[2*N+n] + Beven[3*N+n] - Bodd[2*N+n];
					if ( fabs(0.1+ux / rho) > 1e-13 ){
							if (id[n] == 1) printf("Wetting phase! \n");
							if (id[n] == 2) printf("Non-wetting phase! \n");
							final = ux/rho;
							printf("Momentum (x) not conserved, site=%i,%i,%i, final = %f \n",i,j,k,final);
							CleanCheck=false;
					}
					if ( fabs(0.1+uy / rho) > 1e-13 ){
							if (id[n] == 1) printf("Wetting phase! \n");
							if (id[n] == 2) printf("Non-wetting phase! \n");
							final = uy/rho;
							printf("Momentum (y) not conserved, site=%i,%i,%i, final = %f \n",i,j,k,final);
							CleanCheck=false;
					}
					if ( fabs(0.1+uz / rho) > 1e-13 ){
							if (id[n] == 1) printf("Wetting phase! \n");
							if (id[n] == 2) printf("Non-wetting phase! \n");
							final = uz/rho;
							printf("Momentum (z) not conserved, site=%i,%i,%i, final = %f \n",i,j,k,final);
							CleanCheck=false;
					}
				}
			}
		}
	}

	if (CleanCheck){
		if (rank==0) printf("Test passed: mass conservation for D3Q7 \n");
	}
	else {
		if (rank==0) printf("Test failed!: mass conservation for D3Q7 \n");

	}
	
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
