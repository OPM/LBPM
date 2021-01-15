
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>	// Initialize MPI
	Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
	int rank = comm.getRank();
	int nprocs = comm.getSize();
	int check;
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI.h"

using namespace std;


//***************************************************************************************
int main(int argc, char **argv)
{
	//*****************************************
	// ***** MPI STUFF ****************
	//*****************************************
	// Initialize MPI
	Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
	int rank = comm.getRank();
	int nprocs = comm.getSize();
	int check;
	{
		// parallel domain size (# of sub-domains)
		int nprocx,nprocy,nprocz;
		int iproc,jproc,kproc;


		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Color Model: TestColor	\n");
			printf("********************************************************\n");
		}

		// BGK Model parameters
		string FILENAME;
		unsigned int nBlocks, nthreads;
		int timestepMax, interval;
		double Fx,Fy,Fz,tol;
		// Domain variables
		double Lx,Ly,Lz;
		int nspheres;
		int Nx,Ny,Nz;
		int i,j,k,n;
		int dim = 3;
		Nx = Ny = Nz = 32;
		Lx = Ly = Lz = 1.0;
		//if (rank == 0) printf("dim=%d\n",dim);
		int timestep = 0;
		int timesteps = 100;
		int centralNode = 2;


		double iVol_global = 1.0/Nx/Ny/Nz/nprocx/nprocy/nprocz;
		int BoundaryCondition=0;

		Domain Dm(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);

		Nx += 2;
		Ny += 2;
		Nz += 2;
		int N = Nx*Ny*Nz;

		int Np=0;  // number of local pore nodes
		double *PhaseLabel;
		PhaseLabel = new double[N];
		//.......................................................................
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm.id[n]=1;
					Np++;
					// Initialize gradient ColorGrad = (1,2,3)
					double value=double(3*k+2*j+i);
					PhaseLabel[n]= value;
				}
			}
		}
		Dm.CommInit();
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");

		//Create a second communicator based on the regular data layout
		ScaLBL_Communicator ScaLBL_Comm_Regular(Dm);
		ScaLBL_Communicator ScaLBL_Comm(Dm);

		// LBM variables
		if (rank==0)	printf ("Set up the neighborlist \n");

		int neighborSize=18*Np*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Np];

		ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Dm.id,Np,1);
		MPI_Barrier(comm);

		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);
		if (rank==0)	printf ("Allocating distributions \n");

		int *NeighborList;
		int *dvcMap;
		double *Phi;
		double *ColorGrad;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Nx*Ny*Nz);		
		ScaLBL_AllocateDeviceMemory((void **) &ColorGrad, 3*sizeof(double)*Np);
		
		//...........................................................................
		// Update GPU data structures
		if (rank==0)	printf ("Setting up device map and neighbor list \n");
		int *TmpMap;
		TmpMap=new int[Np*sizeof(int)];
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx < 0))
						TmpMap[idx] = k*Nx*Ny+j*Nx+i;
				}
			}
		}
		ScaLBL_CopyToDevice(dvcMap, TmpMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		delete [] TmpMap;
		
		// copy the neighbor list 
		ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
		// initialize phi based on PhaseLabel (include solid component labels)
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
		//...........................................................................

		ScaLBL_D3Q19_Gradient(dvcMap, Phi, ColorGrad, 0, Np, Np, Nx, Ny, Nz);
	
    	double *COLORGRAD;
    	COLORGRAD= new double [3*Np];
    	int SIZE=3*Np*sizeof(double);
    	ScaLBL_CopyToHost(&COLORGRAD[0],&ColorGrad[0],SIZE);


    	double CX,CY,CZ;
    	for (k=1;k<Nz-1;k++){
    		for (j=1;j<Ny-1;j++){
    			for (i=1;i<Nx-1;i++){
    				n = k*Nx*Ny+j*Nx+i;
    				if (Dm.id[n] > 0){
    					int idx = Map(i,j,k);
    					CX=COLORGRAD[idx];
    					CY=COLORGRAD[Np+idx];
    					CZ=COLORGRAD[2*Np+idx];
    					double error=sqrt((CX-1.0)*(CX-1.0)+(CY-2.0)*(CY-2.0)+ (CZ-3.0)*(CZ-3.0));
    					if (error > 1e-8)
    						printf("i,j,k=%i,%i,%i: Color gradient=%f,%f,%f \n",i,j,k,CX,CY,CZ);
    				}
    			}
    		}
    	}

	}
	// ****************************************************
	comm.barrier();
	Utilities::shutdown();
	// ****************************************************
	return check;
}

