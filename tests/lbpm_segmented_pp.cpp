/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"
#include "analysis/TwoPhase.h"
#include "analysis/eikonal.h"

inline void MeanFilter(DoubleArray &Mesh){
	for (int k=1; k<(int)Mesh.size(2)-1; k++){
		for (int j=1; j<(int)Mesh.size(1)-1; j++){
			for (int i=1; i<(int)Mesh.size(0)-1; i++){
				double sum;
				sum=Mesh(i,j,k)+Mesh(i+1,j,k)+Mesh(i-1,j,k)+Mesh(i,j+1,k)+Mesh(i,j-1,k)+
						+Mesh(i,j,k+1)+Mesh(i,j,k-1);
				Mesh(i,j,k) = sum/7.0;
			}
		}
	}
}


int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	{	
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
		double Lx, Ly, Lz;
		int Nx,Ny,Nz;
		int i,j,k,n;
		int BC=0;
		//  char fluidValue,solidValue;

		std::vector<char> solidValues;
		std::vector<char> nwpValues;
		std::string line;

		if (rank==0){
			ifstream domain("Domain.in");
			domain >> nprocx;
			domain >> nprocy;
			domain >> nprocz;
			domain >> nx;
			domain >> ny;
			domain >> nz;
			domain >> nspheres;
			domain >> Lx;
			domain >> Ly;
			domain >> Lz;

		}
		MPI_Barrier(comm);
		// Computational domain
		MPI_Bcast(&nx,1,MPI_INT,0,comm);
		MPI_Bcast(&ny,1,MPI_INT,0,comm);
		MPI_Bcast(&nz,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
		MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
		MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
		MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
		MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
		//.................................................
		MPI_Barrier(comm);

		// Check that the number of processors >= the number of ranks
		if ( rank==0 ) {
			printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
			printf("Number of MPI ranks used: %i \n", nprocs);
			printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
		}
		if ( nprocs < nprocx*nprocy*nprocz ){
			ERROR("Insufficient number of processors");
		}

		char LocalRankFilename[40];

		int N = (nx+2)*(ny+2)*(nz+2);
		Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
		for (n=0; n<N; n++) Dm.id[n]=1;
		Dm.CommInit(comm);

		// Read the phase ID
		size_t readID;
		sprintf(LocalRankFilename,"ID.%05i",rank);
		FILE *ID = fopen(LocalRankFilename,"rb");
		readID=fread(Dm.id,1,N,ID);
		if (readID != size_t(N)) printf("lbpm_segmented_pp: Error reading ID \n");
		fclose(ID);

		// Initialize the domain and communication
		nx+=2; ny+=2; nz+=2;

		char *id;
		id = new char [N];
		TwoPhase Averages(Dm);
		//	DoubleArray Distance(nx,ny,nz);
		//	DoubleArray Phase(nx,ny,nz);

		int count = 0;
		// Solve for the position of the solid phase
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n = k*nx*ny+j*nx+i;
					// Initialize the solid phase
					if (Dm.id[n] == 0)	id[n] = 0;
					else		      	id[n] = 1;
				}
			}
		}
		// Initialize the signed distance function
		for (k=0;k<nz;k++){
			for (j=0;j<ny;j++){
				for (i=0;i<nx;i++){
					n=k*nx*ny+j*nx+i;
					// Initialize distance to +/- 1
					Averages.SDs(i,j,k) = 2.0*double(id[n])-1.0;
				}
			}
		}
		MeanFilter(Averages.SDs);

		double LocalVar, TotalVar;
		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		int Maxtime=10*max(max(Dm.Nx*Dm.nprocx,Dm.Ny*Dm.nprocy),Dm.Nz*Dm.nprocz);
		Maxtime=1000;
		LocalVar = Eikonal(Averages.SDs,id,Dm,Maxtime);

		MPI_Allreduce(&LocalVar,&TotalVar,1,MPI_DOUBLE,MPI_SUM,comm);
		TotalVar /= nprocs;
		if (rank==0) printf("Final variation in signed distance function %f \n",TotalVar);

		sprintf(LocalRankFilename,"SignDist.%05i",rank);
		FILE *DIST = fopen(LocalRankFilename,"wb");
		fwrite(Averages.SDs.data(),8,N,DIST);
		fclose(DIST);

	}
	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;

}
