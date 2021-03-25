
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/WideHalo.h"
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
	int check=0;
	{

		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Color Model: TestColor	\n");
			printf("********************************************************\n");
		}
		// Domain variables
		int nprocx, nprocy, nprocz;
		double Lx,Ly,Lz;
		int Nx,Ny,Nz;
		int i,j,k,n;
		int dim = 16;
		Lx = Ly = Lz = 1.0;
		int BoundaryCondition=0;

		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		nprocx=nprocy=nprocz=1;
		if (nprocs==1){
			nprocx=nprocy=nprocz=1;
			Nx=Ny=Nz=dim;
			Lx=Ly=Lz=1;
		}
		else if (nprocs==2){
			nprocx=2; nprocy=1;
			nprocz=1;
			Nx=Ny=Nz=dim;
			Nx = dim; Ny = dim; Nz = dim;
			Lx=Ly=Lz=1;
		}
		else if (nprocs==4){
			nprocx=nprocy=2;
			nprocz=1;
			Nx=Ny=Nz=dim;
			Lx=Ly=Lz=1;
		}
		else if (nprocs==8){
			nprocx=nprocy=nprocz=2;
			Nx=Ny=Nz=dim;
			Lx=Ly=Lz=1;
		}
		//.......................................................................
		// **************************************************************

		if (nprocs != nprocx*nprocy*nprocz){
			printf("nprocx =  %i \n",nprocx);
			printf("nprocy =  %i \n",nprocy);
			printf("nprocz =  %i \n",nprocz);
			INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
		}

		if (rank==0){
			printf("********************************************************\n");
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("********************************************************\n");
		}

		comm.barrier();

		std::shared_ptr<Domain> Dm  = std::shared_ptr<Domain>(new Domain(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition));     
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
					Dm->id[n]=1;
					Np++;
					// Initialize gradient ColorGrad = (1,2,3)
					double value=double(3*k+2*j+i);
					PhaseLabel[n]= value;
				}
			}
		}
		Dm->CommInit();
		comm.barrier();
		if (rank == 0) cout << "Domain set." << endl;
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");

		//Create a second communicator based on the regular data layout
		ScaLBL_Communicator ScaLBL_Comm_Regular(Dm);
		ScaLBL_Communicator ScaLBL_Comm(Dm);
		ScaLBLWideHalo_Communicator WideHalo(Dm,2);

		// LBM variables
		if (rank==0)	printf ("Set up the neighborlist \n");

		int neighborSize=18*Np*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Np];

		ScaLBL_Comm.MemoryOptimizedLayoutAA(Map,neighborList,Dm->id.data(),Np,2);
		comm.barrier();

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
		int *WideMap;
		WideMap=new int[Np];
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					int nw = WideHalo.Map(i,j,k);
					int idx = Map(i,j,k);
					WideMap[idx] = nw;
				}
			}
		}
		ScaLBL_CopyToDevice(dvcMap, WideMap, sizeof(int)*Np);
		ScaLBL_DeviceBarrier();
		
		// copy the neighbor list 
		ScaLBL_CopyToDevice(NeighborList, neighborList, neighborSize);
		// initialize phi based on PhaseLabel (include solid component labels)
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					PhaseLabel[n] = 0.2;
				}
			}
		}
		ScaLBL_CopyToDevice(Phi, PhaseLabel, N*sizeof(double));
		//...........................................................................

		int Nxh = Nx+2;
		int Nyh = Ny+2;
		int Nzh = Nz+2;
		ScaLBL_D3Q19_MixedGradient(dvcMap, Phi, ColorGrad, 0, Np, Np, Nxh, Nyh, Nzh);

		double *COLORGRAD;
    	COLORGRAD= new double [3*Np];
    	int SIZE=3*Np*sizeof(double);
    	ScaLBL_CopyToHost(&COLORGRAD[0],&ColorGrad[0],SIZE);

    	double CX,CY,CZ;
    	for (k=1;k<Nz-1;k++){
    		for (j=1;j<Ny-1;j++){
    			for (i=1;i<Nx-1;i++){
    				n = k*Nx*Ny+j*Nx+i;
    				if (Dm->id[n] > 0){
    					int idx = Map(i,j,k);
    					CX=COLORGRAD[idx];
    					CY=COLORGRAD[Np+idx];
    					CZ=COLORGRAD[2*Np+idx];
    					double error=sqrt((CX-1.0)*(CX-1.0)+(CY-2.0)*(CY-2.0)+ (CZ-3.0)*(CZ-3.0));
    					if (error > 1e-8){
    						check++;
    						printf("i,j,k=%i,%i,%i: Color gradient=%f,%f,%f \n",i,j,k,CX,CY,CZ);
    					}
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

