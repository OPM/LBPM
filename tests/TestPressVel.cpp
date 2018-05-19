
//*************************************************************************
// Lattice Boltzmann Simulator for Single Phase Flow in Porous Media
// James E. McCLure
//*************************************************************************
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "common/ScaLBL.h"
#include "common/MPI_Helpers.h"


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
		if (rank == 0){
			printf("********************************************************\n");
			printf("Running Unit Test: TestPressVel	\n");
			printf("********************************************************\n");
		}
		
		// Domain variables
		int i,j,k,n;
        // Load inputs
		string FILENAME = argv[1];
        // Load inputs
		if (rank==0)	printf("Loading input database \n");
	    auto db = std::make_shared<Database>( FILENAME );
	    auto domain_db = db->getDatabase( "Domain" );
	    int Nx = domain_db->getVector<int>( "n" )[0];
        int Ny = domain_db->getVector<int>( "n" )[1];
        int Nz = domain_db->getVector<int>( "n" )[2];
        int nprocx = domain_db->getVector<int>( "nproc" )[0];
        int nprocy = domain_db->getVector<int>( "nproc" )[1];
        int nprocz = domain_db->getVector<int>( "nproc" )[2];
		if (rank==0){
			printf("********************************************************\n");
			printf("Sub-domain size = %i x %i x %i\n",Nx,Ny,Nz);
			printf("********************************************************\n");
		}

		MPI_Barrier(comm);
		int kproc = rank/(nprocx*nprocy);
		int jproc = (rank-nprocx*nprocy*kproc)/nprocx;
		int iproc = rank-nprocx*nprocy*kproc-nprocz*jproc;

		if (rank == 0) {
			printf("i,j,k proc=%d %d %d \n",iproc,jproc,kproc);
		}
		MPI_Barrier(comm);
		if (rank == 1){
			printf("i,j,k proc=%d %d %d \n",iproc,jproc,kproc);
			printf("\n\n");
		}

		double iVol_global = 1.0/Nx/Ny/Nz/nprocx/nprocy/nprocz;

		std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
		Dm->CommInit();

		Nx += 2;
		Ny += 2;
		Nz += 2;
		int N = Nx*Ny*Nz;

		//.......................................................................
		// Assign the phase ID field
		//.......................................................................
		char LocalRankString[8];
		sprintf(LocalRankString,"%05d",rank);
		char LocalRankFilename[40];
		sprintf(LocalRankFilename,"ID.%05i",rank);

		Dm->CommInit();

		//.......................................................................
		// Compute the media porosity
		//.......................................................................
		double sum;
		double sum_local=0.0, porosity;
		int Np=0;  // number of local pore nodes
		for (k=1;k<Nz-1;k++){
			for (j=1;j<Ny-1;j++){
				for (i=1;i<Nx-1;i++){
					n = k*Nx*Ny+j*Nx+i;
					Dm->id[n] = 1;
					if (Dm->id[n] > 0){
						sum_local+=1.0;
						Np++;
					}
				}
			}
		}
		MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
		porosity = sum*iVol_global;
		if (rank==0) printf("Media porosity = %f \n",porosity);

		MPI_Barrier(comm);
		if (rank == 0) cout << "Domain set." << endl;
		if (rank==0)	printf ("Create ScaLBL_Communicator \n");

		// Create a communicator for the device
		auto ScaLBL_Comm  = std::shared_ptr<ScaLBL_Communicator>(new ScaLBL_Communicator(Dm));
		//...........device phase ID.................................................
		if (rank==0)	printf ("Copying phase ID to device \n");
		char *ID;
		ScaLBL_AllocateDeviceMemory((void **) &ID, N);						// Allocate device memory
		// Copy to the device
		ScaLBL_CopyToDevice(ID, Dm->id, N);
		//...........................................................................

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
		Np = ScaLBL_Comm->MemoryOptimizedLayoutAA(Map,neighborList,Dm->id,Np);
		MPI_Barrier(comm);

		//......................device distributions.................................
		if (rank==0)	printf ("Allocating distributions \n");
		int dist_mem_size = Np*sizeof(double);
		int *NeighborList;
		double * dist;
		double * Velocity;
		//...........................................................................
		ScaLBL_AllocateDeviceMemory((void **) &dist, 19*dist_mem_size);
		ScaLBL_AllocateDeviceMemory((void **) &NeighborList, neighborSize);
		ScaLBL_AllocateDeviceMemory((void **) &Velocity, 3*sizeof(double)*Np);
		ScaLBL_CopyToDevice(NeighborList,     neighborList, neighborSize);
		//...........................................................................

		/*
		 *  AA Algorithm begins here
		 *
		 */
		double *DIST;
		DIST = new double [19*Np];
		double VALUE=0.1;

		for (int n=0; n<Np; n++){
		     //DIST[n]=1.0;
		  // even distributions
		   DIST[Np + n] = 1.0 - VALUE;
		   DIST[2*Np + n] = 1.0 - VALUE;
		   DIST[3*Np + n] = 1.0 - VALUE;
		   DIST[4*Np + n] = 1.0;
		   DIST[5*Np + n] = 1.0;
		   DIST[6*Np + n] = 1.0;
		   DIST[7*Np + n] = 1.0;
		   DIST[8*Np + n] = 1.0;
		   DIST[9*Np + n] = 1.0;
		   // odd distributions
		   DIST[10*Np + n] = 1.0;
		   DIST[11*Np + n] = 1.0;
		   DIST[12*Np + n] = 1.0;
		   DIST[13*Np + n] = 1.0;
		   DIST[14*Np + n] = 1.0;
		   DIST[15*Np + n] = 1.0;
		   DIST[16*Np + n] = 1.0;
		   DIST[17*Np + n] = 1.0;
		   DIST[18*Np + n] = 1.0;
	 	}
		ScaLBL_CopyToDevice(dist, DIST, 19*Np*sizeof(double));	

	   double *Vz;
  	   Vz= new double [Np];
	   size_t SIZE=Np*sizeof(double);
	   ScaLBL_D3Q19_Momentum(dist, Velocity, Np);
	   ScaLBL_CopyToHost(&Vz[0],&Velocity[2],SIZE);
	   // 
	   for (int n=0; n<Np; n++){
	     if (Vz[n] - VALUE > 1e-8){
	       printf("ERROR: site %i, value=%f \n",n,Vz[n]); check = 15;
	     }
	   }
	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
	return check;

}
