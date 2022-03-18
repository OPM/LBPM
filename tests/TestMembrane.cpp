
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/MPI.h"
#include "common/Membrane.h"

using namespace std;

std::shared_ptr<Database> loadInputs( int nprocs )
{
  //auto db = std::make_shared<Database>( "Domain.in" );
    auto db = std::make_shared<Database>();
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", { 1, 1, 1 } );
    db->putVector<int>( "n", { 32, 32, 32 } );
    db->putScalar<int>( "nspheres", 1 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
}

//***************************************************************************************
int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
	int check=0;
	{

		int i,j,k,n;
		
		static int D3Q19[18][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
				{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},
				{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
				{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};

        int rank = comm.getRank();
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running unit test: TestMembrane	\n");
			printf("********************************************************\n");
		}
		
	    // Load inputs
	    auto db = loadInputs( comm.getSize() );
	    int Nx = db->getVector<int>( "n" )[0];
	    int Ny = db->getVector<int>( "n" )[1];
	    int Nz = db->getVector<int>( "n" )[2];
		auto Dm = std::make_shared<Domain>(db,comm);

		Nx += 2;
		Ny += 2;
		Nz += 2;
		int N = Nx*Ny*Nz;
		//.......................................................................
		int Np = 0;
		double distance,radius;
		DoubleArray Distance(Nx,Ny,Nz);
		for (k=0;k<Nz;k++){
			for (j=0;j<Ny;j++){
				for (i=0;i<Nx;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm->id[n] = 1;
					radius = double(Nx)/4;
					distance = sqrt(double((i-0.5*Nx)*(i-0.5*Nx)+ (j-0.5*Ny)*(j-0.5*Ny)+ (k-0.5*Nz)*(k-0.5*Nz)))-radius;
					if (distance < 0.0 ){
						Dm->id[n] = 1;
					}
					Distance(i,j,k) = distance;
					Np++;
				}
			}
		}
		Dm->CommInit();

		// Create a communicator for the device (will use optimized layout)
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm(new ScaLBL_Communicator(Dm));
		//Create a second communicator based on the regular data layout
		std::shared_ptr<ScaLBL_Communicator> ScaLBL_Comm_Regular(new ScaLBL_Communicator(Dm));

		if (rank==0){
			printf("Total domain size = %i \n",N);
			printf("Reduced domain size = %i \n",Np);
		}

		// LBM variables
		if (rank==0)	printf ("Set up the neighborlist \n");
		int Npad=Np+32;
		int neighborSize=18*Npad*sizeof(int);
		int *neighborList;
		IntArray Map(Nx,Ny,Nz);
		neighborList= new int[18*Npad];
		
		//......................device distributions.................................
		int *NeighborList;
		int *dvcMap;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &dvcMap, sizeof(int)*Npad);
	    ScaLBL_CopyToDevice(NeighborList, neighborList, 18*Np*sizeof(int));

		Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Dm->id.data(),Np,1);
		comm.barrier();
		
		double *dist;
		dist = new double [19*Np];
		
		// Check the neighborlist
		printf("Check neighborlist: exterior %i, first interior %i last interior %i \n",ScaLBL_Comm->LastExterior(),ScaLBL_Comm->FirstInterior(),ScaLBL_Comm->LastInterior());
		for (int idx=0; idx<ScaLBL_Comm->LastExterior(); idx++){
			for (int q=0; q<18; q++){
				int nn = neighborList[q*Np+idx]%Np;
				if (nn>Np) printf("neighborlist error (exterior) at q=%i, idx=%i \n",q,idx);
				dist[q*Np + idx] = 0.0;
			}
		}
		for (int idx=ScaLBL_Comm->FirstInterior(); idx<ScaLBL_Comm->LastInterior(); idx++){
			for (int q=0; q<18; q++){
				int nn = neighborList[q*Np+idx]%Np;
				if (nn>Np) printf("neighborlist error (exterior) at q=%i, idx=%i \n",q,idx);
				dist[q*Np + idx] = 0.0;
			}
		}
		
		/* create a membrane data structure */
	    Membrane M(Dm, NeighborList, Np);

	    int MembraneCount = M.Create(Dm, Distance, Map);
		if (rank==0)	printf (" Number of membrane links: %i \n", MembraneCount);

		/* create a tagged array to show where the mebrane is*/
		double *MembraneLinks;
		MembraneLinks = new double [Nx*Ny*Nz];
		for (int n=0; n<Nx*Ny*Nz; n++) {
			MembraneLinks[n] = 0.0;
		}
		for (int mlink=0; mlink<MembraneCount; mlink++){
			int iq = M.membraneLinks[2*mlink];
			int jq = M.membraneLinks[2*mlink+1];
			dist[iq] = -1.0; // set these distributions to non-zero
			dist[jq] = 1.0;
		}
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					int idx = Map(i,j,k);
					double sum = 0.0;
					for (int q=0; q<19; q++){
						sum += dist[q*Np + idx];
					}
					int n = k*Nx*Ny + j*Nx + i;
					MembraneLinks[n] = sum;
					if (sum > 0.f){
						Dm->id[n] = 127;
					}
					if (sum < 0.f){
						Dm->id[n] = 64;
					}
				}
			}
		}
		if (argc > 1)
			Dm->AggregateLabels("membrane.raw");
				 
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
		
		// Create a dummy distribution data structure
		double *fq;
		fq = new double[19*Np];
		if (rank==0)	printf ("Setting up distributions \n");
		for (k=1; k<Nz-1; k++){
			for (j=1; j<Ny-1; j++){
				for (i=1; i<Nx-1; i++){
					int idx=Map(i,j,k);
					if (!(idx<0)){
					  for (int q=0; q<19; q++){
					    fq[q*Np+idx]=k*100.f+j*10.f+i*1.f+0.01*q;
					  }
					}
				}
			}
		}

		// Loop over the distributions for interior lattice sites
		if (rank==0)	printf ("Loop over distributions \n");

		for (int idx=ScaLBL_Comm->first_interior; idx<ScaLBL_Comm->last_interior; idx++){
			n = TmpMap[idx];
			k = n/(Nx*Ny);
			j = (n-Nx*Ny*k)/Nx;
			i = n-Nx*Ny*k-Nx*j;
			for (int q=1; q<19; q++){
				int nn = neighborList[(q-1)*Np+idx];
				double value=fq[nn];
				// 3D index of neighbor
				int iq=i-D3Q19[q-1][0];
				int jq=j-D3Q19[q-1][1];
				int kq=k-D3Q19[q-1][2];
				if (iq==0) iq=1;
				if (jq==0) jq=1;
				if (kq==0) kq=1;
				if (iq==Nx-1) iq=Nx-2;
				if (jq==Ny-1) jq=Ny-2;
				if (kq==Nz-1) kq=Nz-2;
				double check = kq*100.f+jq*10.f+iq*1.f+q*0.01;
				if (value != check)
				  printf("Neighbor q=%i, i=%i,j=%i,k=%i: %f \n",q,iq,jq,kq,value);
			}
		}

		delete [] TmpMap;

	}
    Utilities::shutdown();

	return check;
}

