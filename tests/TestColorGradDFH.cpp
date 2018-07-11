
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"

using namespace std;

std::shared_ptr<Database> loadInputs( int nprocs )
{
    auto db = std::make_shared<Database>();
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", { 1, 1, 1 } );
    db->putVector<int>( "n", { 100, 100, 100 } );
    db->putScalar<int>( "nspheres", 1 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
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
	int check=0;
	{
		// parallel domain size (# of sub-domains)
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Color Model: TestColorGradDFH	\n");
			printf("********************************************************\n");
		}

		// BGK Model parameters
		string FILENAME;
		// Domain variables
		int i,j,k,n;
		
	    // Load inputs
	    auto db = loadInputs( nprocs );
	    int Nx = db->getVector<int>( "n" )[0];
	    int Ny = db->getVector<int>( "n" )[1];
	    int Nz = db->getVector<int>( "n" )[2];
	    int nprocx = db->getVector<int>( "nproc" )[0];
	    int nprocy = db->getVector<int>( "nproc" )[1];
	    int nprocz = db->getVector<int>( "nproc" )[2];

	    if (rank==0){
	    	printf("********************************************************\n");
	    	printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
	    	printf("********************************************************\n");
	    }

	    // Get the rank info
	    std::shared_ptr<Domain> Dm(new Domain(db,comm));
		Nx += 2;
		Ny += 2;
		Nz += 2;
		int N = Nx*Ny*Nz;

		double *PhaseLabel;
		PhaseLabel = new double[N];
		//.......................................................................
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm->id[n]=1;
					// Initialize gradient ColorGrad = (1,2,3)
					double value=double(3*k+2*j+i);
					PhaseLabel[n]= value;
				}
			}
		}
		Dm->CommInit();
		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");

		//Create a second communicator based on the regular data layout
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm(new ScaLBL_Communicator(Dm));

		// LBM variables
		if (rank==0)	printf ("Set up the neighborlist \n");
		
		int Np=0;  // number of local pore nodes
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					Np++;
				}
			}
		}
		int Npad=(Np/16 + 2)*16;
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Npad];
		Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Dm->id,Np);
		MPI_Barrier(comm);

		//......................device distributions.................................
		int dist_mem_size = Np*sizeof(double);
		int neighborSize=18*Np*sizeof(int);
		if (rank==0)	printf ("Allocating distributions \n");
		int *NeighborList;
		int *dvcMap;
		double *Phi;
		double *Potential;
		double *ColorGrad;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Np);
		ScaLBL_AllocateDeviceMemory((void **) &Phi, sizeof(double)*Np);		
		ScaLBL_AllocateDeviceMemory((void **) &Potential, sizeof(double)*Np);		
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

		double *PHASE;
		PHASE=new double [Np];
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx < 0)){
						PHASE[idx] = PhaseLabel[k*Nx*Ny+j*Nx+i];
					}
				}
			}
		}
		
		ScaLBL_CopyToDevice(Phi, PHASE, Np*sizeof(double));
		//...........................................................................

		// compute the gradient 
		ScaLBL_D3Q19_Gradient_DFH(neighborList, Phi, ColorGrad, ScaLBL_Comm->first_interior, ScaLBL_Comm->last_interior, Np);
		ScaLBL_Comm->SendHalo(Phi);
		ScaLBL_D3Q19_Gradient_DFH(neighborList, Phi, ColorGrad, 0, ScaLBL_Comm->first_interior, Np);
		ScaLBL_Comm->RecvGrad(Phi,ColorGrad);
		
    	double *COLORGRAD;
    	COLORGRAD= new double [3*Np];
    	int SIZE=3*Np*sizeof(double);
    	ScaLBL_CopyToHost(&COLORGRAD[0],&ColorGrad[0],SIZE);

    	for (k=1;k<Nz-1;k++){
    		for (j=1;j<Ny-1;j++){
    			for (i=1;i<Nx-1;i++){
    				n = k*Nx*Ny+j*Nx+i;
    				if (Dm->id[n] > 0){
    					int idx = Map(i,j,k);
    					printf("%i ",idx);
    				}
    			}
    			printf("\n");
    		}
    		printf("-------\n");
    	}
    	
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
    						printf("i,j,k=%i,%i,%i; idx=%i: Color gradient=%f,%f,%f \n",i,j,k,idx,CX,CY,CZ);
    					/*	for (int q=0; q<18; q++){
    							int nn = neighborList[q*Np+idx]%Np;
    							double value= PHASE[nn];
    							printf("     q=%i, nn=%i, value=%f \n",q,nn,value);
    						}
    						*/
    					}
    				}
    			}
    		}
    	}

	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************

	return check;
}

