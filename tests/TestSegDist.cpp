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
#include "common/Array.h"
#include "common/Domain.h"
#include "analysis/eikonal.h"


std::shared_ptr<Database> loadInputs( int nprocs )
{
    INSIST(nprocs==8, "TestSegDist: Number of MPI processes must be equal to 8");
    auto db = std::make_shared<Database>( );
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", { 2, 2, 2 } );
    db->putVector<int>( "n", { 50, 50, 50 } );
    db->putScalar<int>( "nspheres", 0 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
}


//***************************************************************************************
int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
	{
	int i,j,k,n,nn;
	int iproc,jproc,kproc;



    // Load inputs
    auto db = loadInputs( nprocs );
    int Nx = db->getVector<int>( "n" )[0];
    int Ny = db->getVector<int>( "n" )[1];
    int Nz = db->getVector<int>( "n" )[2];
    int nprocx = db->getVector<int>( "nproc" )[0];
    int nprocy = db->getVector<int>( "nproc" )[1];
    int nprocz = db->getVector<int>( "nproc" )[2];


    // Get the rank info
    Domain Dm(db);
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				Dm.id[n] = 1;
			}
		}
	}
	Dm.CommInit(comm);

	int nx = Nx+2;
    int ny = Ny+2;
    int nz = Nz+2;
	int N = nx*ny*nz;
	int count = 0;

	char *id;
	id = new char [N];
	double BubbleRadius = 25;
	// Initialize the bubble
	double x,y,z , Cx,Cy,Cz;

	Cx = 1.0*nx;
	Cy = 1.0*ny;
	Cz = 1.0*nz;

	DoubleArray Distance(nx,ny,nz);
	DoubleArray TrueDist(nx,ny,nz);

	for (k=1;k<nz-1;k++){
		for (j=1;j<ny-1;j++){
			for (i=1;i<nx-1;i++){

				// True signed distance
				x = (nx-2)*Dm.iproc+i-1;
				y = (ny-2)*Dm.jproc+j-1;
				z = (nz-2)*Dm.kproc+k-1;
				TrueDist(i,j,k) = sqrt((x-Cx)*(x-Cx)+(y-Cy)*(y-Cy)+(z-Cz)*(z-Cz)) - BubbleRadius;

				n = k*nx*ny+j*nx+i;

				// Initialize phase positions
				if (TrueDist(i,j,k) < 0.0){
					id[n] = 0;
				}
				else{
					id[n]=1;
				}


			}
		}
	}
	
	// Initialize the signed distance function
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				n=k*nx*ny+j*nx+i;
				// Initialize distance to +/- 1
				Distance(i,j,k) = 1.0*id[n]-0.5;
			}
		}
	}

	if (rank==0) printf("Nx = %i \n",(int)Distance.size(0));
	if (rank==0) printf("Ny = %i \n",(int)Distance.size(1));
	if (rank==0) printf("Nz = %i \n",(int)Distance.size(2));

	MPI_Barrier(comm);
	if (rank==0) printf("Initialized! Converting to Signed Distance function \n");

	double starttime,stoptime,cputime;
	starttime = MPI_Wtime();
	Eikonal(Distance,id,Dm,10);
	stoptime = MPI_Wtime();
	cputime = (stoptime - starttime);

	if (rank==0) printf("Total time: %f seconds \n",cputime);


	double Error=0.0;
	int Count = 0;
	for (k=0;k<nz;k++){
		for (j=0;j<ny;j++){
			for (i=0;i<nx;i++){
				if (fabs(TrueDist(i,j,k)) < 3.0){
					Error += (Distance(i,j,k)-TrueDist(i,j,k))*(Distance(i,j,k)-TrueDist(i,j,k));
					Count += 1;
				}
			}
		}
	}
	Error = sqrt(Error)/(double (Count));
	if (rank==0) printf("Mean error %f \n", Error);


    MPI_Barrier(comm);
    }
    MPI_Finalize();
    return 0;

}
