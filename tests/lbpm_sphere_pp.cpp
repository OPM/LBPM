#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "analysis/pmmc.h"
#include "common/Domain.h"
#include "common/SpherePack.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"

/*
 * Pre-Processor to generate signed distance function from sphere packing
 * to use as an input domain for lattice Boltzmann simulator
 * James E. McClure 2014
 */

using namespace std;

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
		printf("Running Sphere Packing pre-processor for LBPM-WIA	\n");
		printf("********************************************************\n");
	}

	double D = 1.0;		// reference length for non-dimensionalization
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
    int nspheres = domain_db->getScalar<int>( "nsphere" );
    int Lx = domain_db->getVector<double>( "L" )[0];
    int Ly = domain_db->getVector<double>( "L" )[1];
    int Lz = domain_db->getVector<double>( "L" )[2];
    
	int i,j,k,n;

	// **************************************************************
	if (nprocs != nprocx*nprocy*nprocz){
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
		INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
	}

	Nz += 2;
	Nx = Ny = Nz;	// Cubic domain

	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

	if (rank==0) printf("Number of nodes per side = %i \n", Nx);
	if (rank==0) printf("Total Number of nodes = %i \n", N);
	if (rank==0) printf("********************************************************\n");

	//.......................................................................
	if (rank == 0)	printf("Read input media... \n");
	//.......................................................................
	
	//.......................................................................
	// Filenames used
	char LocalRankString[8];
	char LocalRankFilename[40];
	char LocalRestartFile[40];
	char tmpstr[10];
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
	sprintf(LocalRestartFile,"%s%s","Restart.",LocalRankString);
	
//	printf("Local File Name =  %s \n",LocalRankFilename);
	// .......... READ THE INPUT FILE .......................................
//	char value;
	char *id;
	id = new char[N];
	int sum = 0;
	double sum_local;
	double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
	double porosity, pore_vol;
	//...........................................................................
	DoubleArray SignDist(Nx,Ny,Nz);
	//.......................................................................

	// Read in sphere pack
	if (rank==1) printf("nspheres =%i \n",nspheres);
	//.......................................................................
	double *cx,*cy,*cz,*rad;
	cx = new double[nspheres];
	cy = new double[nspheres];
	cz = new double[nspheres];
	rad = new double[nspheres];
	//.......................................................................
	if (rank == 0)	printf("Reading the sphere packing \n");
	if (rank == 0)	ReadSpherePacking(nspheres,cx,cy,cz,rad);
	MPI_Barrier(comm);
	// Broadcast the sphere packing to all processes
	MPI_Bcast(cx,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(cy,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(cz,nspheres,MPI_DOUBLE,0,comm);
	MPI_Bcast(rad,nspheres,MPI_DOUBLE,0,comm);
	//...........................................................................
	MPI_Barrier(comm);
	if (rank == 0) cout << "Domain set." << endl;
	if (rank == 0){
		// Compute the Sauter mean diameter
		double totVol = 0.0;
		double totArea = 0.0;
		// Compute the total volume and area of all spheres
		for (i=0; i<nspheres; i++){
			totVol += 1.3333333333333*3.14159265359*rad[i]*rad[i]*rad[i];
			totArea += 4.0*3.14159265359*rad[i]*rad[i];
		}
		D = 6.0*(Nx-2)*nprocx*totVol / totArea / Lx;
		printf("Sauter Mean Diameter (computed from sphere packing) = %f \n",D);
	}
	MPI_Bcast(&D,1,MPI_DOUBLE,0,comm);

	//.......................................................................
	SignedDistance(SignDist.data(),nspheres,cx,cy,cz,rad,Lx,Ly,Lz,Nx,Ny,Nz,
					   iproc,jproc,kproc,nprocx,nprocy,nprocz);
	//.......................................................................
	// Assign the phase ID field based on the signed distance
	//.......................................................................
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny+j*Nx+i;
				id[n] = 0;
			}
		}
	}
	sum=0;
	pore_vol = 0.0;
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (SignDist(n) > 0.0){ 
					id[n] = 2;	
				}
				// compute the porosity (actual interface location used)
				if (SignDist(n) > 0.0){ 
					sum++;	
				}
			}
		}
	}
	sum_local = 1.0*sum;
	MPI_Allreduce(&sum_local,&porosity,1,MPI_DOUBLE,MPI_SUM,comm);
	porosity = porosity*iVol_global;
	if (rank==0) printf("Media porosity = %f \n",porosity);

	// Compute the pore volume
	sum_local = 0.0;
	for ( k=1;k<Nz-1;k++){
		for ( j=1;j<Ny-1;j++){
			for ( i=1;i<Nx-1;i++){
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] > 0){
					sum_local += 1.0;
				}
			}
		}
	}
	MPI_Allreduce(&sum_local,&pore_vol,1,MPI_DOUBLE,MPI_SUM,comm);
	
	//.........................................................
	// don't perform computations at the eight corners
	id[0] = id[Nx-1] = id[(Ny-1)*Nx] = id[(Ny-1)*Nx + Nx-1] = 0;
	id[(Nz-1)*Nx*Ny] = id[(Nz-1)*Nx*Ny+Nx-1] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx] = id[(Nz-1)*Nx*Ny+(Ny-1)*Nx + Nx-1] = 0;
	//.........................................................

	//.......................................................................
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
	WriteLocalSolidDistance(LocalRankFilename, SignDist.data(), N);
	//......................................................................

	sprintf(LocalRankFilename,"ID.%05i",rank);
	FILE *ID = fopen(LocalRankFilename,"wb");
	fwrite(id,1,N,ID);
	fclose(ID);

	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
