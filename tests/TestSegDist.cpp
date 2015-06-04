// Compute the signed distance from a digitized image 
// Two phases are present
// Phase 1 has value -1
// Phase 2 has value 1
// this code uses the segmented image to generate the signed distance 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <Array.h>
#include <Domain.h>

int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	int i,j,k,n,nn;
	int iproc,jproc,kproc;
	int nx,ny,nz;
	int Nx, Ny, Nz, N;
    int nprocx, nprocy, nprocz, nspheres;
    double Lx, Ly, Lz;
	Nx = Ny = Nz = 50;
	nx = ny = nz = 50;
	N = Nx*Ny*Nz;
	nprocx=nprocy=nprocz=2;
	Lx = Ly = Lz = 1.0;
	int BC=0;

	
	if (nprocs != 8){
		ERROR("TestSegDist: Number of MPI processes must be equal to 8");
	}

    if (nprocx !=2 || nprocz !=2 || nprocy !=2 ){
		ERROR("TestSegDist: MPI process grid must be 2x2x2");
	}

    // Get the rank info
	Domain Dm(nx,ny,nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC);
	Dm.CommInit(MPI_COMM_WORLD);

	nx+=2; ny+=2; nz+=2;
	int count = 0;

	char *id;
	id = new char [N];
	double BubbleRadius = 5;
	// Initialize the bubble
	int x,y,z;
	for (k=1;k<nz-1;k++){
		for (j=1;j<ny-1;j++){
			for (i=1;i<nx-1;i++){
				x = (nx-2)*Dm.iproc+i;
				y = (ny-2)*Dm.jproc+j;
				z = (nz-2)*Dm.kproc+k;
				n = k*nx*ny+j*nx+i;

				// Initialize phase positions
				if ((x-nx+1)*(x-nx+1)+(y-ny+1)*(y-ny+1)+(z-nz+1)*(z-nz+1) < BubbleRadius*BubbleRadius){
					id[n] = 0;
				}
				else{
					id[n]=1;
				}
			}
		}
	}
	
	DoubleArray Distance(nx,ny,nz);
	// Initialize the signed distance function
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n=k*Nx*Ny+j*Nx+i;
				// Initialize distance to +/- 1
				Distance(i,j,k) = 2.0*id[n]-1.0;
			}
		}
	}

	if (rank==0) printf("Nx = %i \n",(int)Distance.size(0));
	if (rank==0) printf("Ny = %i \n",(int)Distance.size(1));
	if (rank==0) printf("Nz = %i \n",(int)Distance.size(2));

	printf("Initialized! Converting to Signed Distance function \n");
	SSO(Distance,id,Dm,1);

	char LocalRankFilename[40];
    sprintf(LocalRankFilename,"Dist.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"wb");
    fwrite(Distance.get(),8,Distance.length(),DIST);
    fclose(DIST);

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;

}
