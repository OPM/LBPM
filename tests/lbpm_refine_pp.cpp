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
#include "common/pmmc.h"

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);

	int InitialWetting;
	double Saturation;
	//	if (argc == 3){
	//sscanf(argv[1],"%lf",&Saturation);
	//sscanf(argv[2],"%d",&InitialWetting);
	Saturation=strtod(argv[1],NULL);
	InitialWetting=atoi(argv[2]);
	if (rank==0){
		printf("Initializing wetting phase saturation of %f \n",Saturation);
		if (InitialWetting == 1)
			printf("Initial connected phase labeled (1) \n");
		else
			printf("Initial connected phase labeled (2) \n");
	}

	if (InitialWetting == 1)	Saturation=1.0-Saturation;
	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
	double Lx, Ly, Lz;
	int i,j,k,n;
	int BC=0;

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

	int BoundaryCondition=0;
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BoundaryCondition);

	nx+=2; ny+=2; nz+=2;
	int N = nx*ny*nz;
	char *id;
	id = new char[N];

	// Define communication sub-domain -- everywhere
	for (int k=0; k<nz; k++){
		for (int j=0; j<ny; j++){
			for (int i=0; i<nx; i++){
				n = k*nx*ny+j*nx+i;
				Dm.id[n] = 1;
			}
		}
	}
	Dm.CommInit(comm);

	DoubleArray SignDist(nx,ny,nz);
	// Read the signed distance from file
	sprintf(LocalRankFilename,"SignDist.%05i",rank);
	FILE *DIST = fopen(LocalRankFilename,"rb");
	size_t ReadSignDist;
	ReadSignDist=fread(SignDist.data(),8,N,DIST);
	if (ReadSignDist != size_t(N)) printf("lbpm_random_pp: Error reading signed distance function (rank=%i)\n",rank);
	fclose(DIST);

	int rnx,rny,rnz;
	rnx=2*(nx-1)+1;
	rny=2*(ny-1)+1;
	rnz=2*(nz-1)+1;
	DoubleArray RefinedSignDist(rnx,rny,rnz);	
	TriLinPoly LocalApprox;
	Point pt;

	int ri,rj,rk,rn; //refined mesh indices
	count = 0;
	for (int k=0; k<nz-1; k++){
		for (int j=0; j<ny-1; j++){
			for (int i=0; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				ri=2*i; 
				rj=2*j;
				rk=2*k;
				// Assign local tri-linear polynomial
				LocalApprox.assign(SignDist,i,j,k); //
				pt.x=1.0*i; pt.y=1.0*j; pt.z=1.0*k;
			}}
		}
	}

	sprintf(LocalRankFilename,"ID.%05i",rank);
	FILE *ID = fopen(LocalRankFilename,"wb");
	fwrite(id,1,N,ID);
	fclose(ID);

	MPI_Barrier(comm);
	MPI_Finalize();
	return 0;
}
