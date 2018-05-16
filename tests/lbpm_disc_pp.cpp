#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "analysis/pmmc.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/MPI_Helpers.h"    // This includes mpi.h
#include "common/SpherePack.h"

/*
 * Pre-Processor to generate signed distance function from disc packing
 * to use as an input domain for lattice Boltzmann simulator
 * James E. McClure 2014
 */

using namespace std;

inline void ReadDiscPacking(int ndiscs, double *List_cx, double *List_cy, double *List_rad)
{
	// Read in the full disc pack
	//...... READ IN THE DISCS...................................
	cout << "Reading the packing file..." << endl;
	FILE *fid = fopen("DiscPack.in","rb");
	INSIST(fid!=NULL,"Error opening DistPack.in");
	//............................................
	char * line = new char[100];
    // We will read until a blank like or end-of-file is reached
	int count = 0;
	while ( !feof(fid) && fgets(line,100,fid)!=NULL ) {
		char* line2 = line;
		List_cx[count] = strtod(line2,&line2);
		List_cy[count] = strtod(line2,&line2);
		List_rad[count] = strtod(line2,&line2);
		count++;
	}
	cout << "Number of discs extracted is: " << count << endl;
    INSIST( count==ndiscs, "Specified number of discs is probably incorrect!" );
	// .............................................................
}

inline void SignedDistanceDiscPack(double *Distance, int ndiscs, double *List_cx, double *List_cy, double *List_rad,
						  double Lx, double Ly, double Lz, int Nx, int Ny, int Nz,
						  int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Use sphere lists to determine which nodes are in porespace
	// Write out binary file for nodes
	int N = Nx*Ny*Nz; 	// Domain size, including the halo
	double hx,hy,hz;
	double x,y,z;
	double cx,cy,cz,r;
	int imin,imax,jmin,jmax,kmin,kmax;
	int p,i,j,k,n;
	//............................................
	double min_x,min_y,min_z;
	double distance;
	//............................................
	// Lattice spacing for the entire domain
	// It should generally be true that hx=hy=hz
	// Otherwise, you will end up with ellipsoids
	hx = Lx/((Nx-2)*nprocx-1);
	hy = Ly/((Ny-2)*nprocy-1);
	hz = Lz/((Nz-2)*nprocz-1);
	//............................................
	// Get maximum and minimum for this domain
	// Halo is included !
	min_x = double(iproc*(Nx-2)-1)*hx;
	min_y = double(jproc*(Ny-2)-1)*hy;
	min_z = double(kproc*(Nz-2)-1)*hz;
    NULL_USE(min_x);
	//............................................

	//............................................
		// Pre-initialize Distance
	for (n=0;n<N;n++){
		Distance[n]=100.0;
	}
	//............................................

	//............................................
	// .........Loop over the spheres.............
	for (p=0;p<ndiscs;p++){
		// Get the sphere from the list, map to local min
		cx = List_cx[p] - min_y;
		cy = List_cy[p] - min_z;
		r = List_rad[p];
		// Check if
		// Range for this sphere in global indexing
		jmin = int ((cx-2*r)/hy);
		jmax = int ((cx+2*r)/hy)+2;
		kmin = int ((cy-2*r)/hz);
		kmax = int ((cy+2*r)/hz)+2;

		// Obviously we have to do something at the edges
		if (jmin<0)		jmin = 0;
		if (jmin>Ny)	jmin = Ny;
		if (jmax<0)		jmax = 0;
		if (jmax>Ny)	jmax = Ny;
		if (kmin<0)		kmin = 0;
		if (kmin>Nz)	kmin = Nz;
		if (kmax<0)		kmax = 0;
		if (kmax>Nz)	kmax = Nz;
		// Loop over the domain for this sphere (may be null)
		for (k=kmin;k<kmax;k++){
			for (j=jmin;j<jmax;j++){
				for (i=0;i<Nx;i++){
					// x,y,z is distance in physical units
					//x = i*hx;
					x = j*hy;
					y = k*hz;
					// if inside disc, set to zero
					// get the position in the list -- x direction assigned the same
					n = k*Nx*Ny+j*Nx+i;
					// Compute the distance
					distance = sqrt((cx-x)*(cx-x)+(cy-y)*(cy-y)) - r;
					// Assign the minimum distance
					if (distance < Distance[n])		Distance[n] = distance;
				}
			}
		}
	}

	// Map the distance to lattice units
	for (n=0; n<N; n++)	Distance[n] = Distance[n]/hx;
}

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
	int nprocx,nprocy,nprocz;
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

	int depth;

	if (rank == 0){
		printf("********************************************************\n");
		printf("Running Disc Packing Pre-Processor for LBPM-WIA	\n");
		printf("********************************************************\n");
	}

	// Variables that specify the computational domain  
	string FILENAME;
	int Nx,Ny,Nz;		// local sub-domain size
	int ndiscs;		// number of spheres in the packing
	double Lx,Ly,Lz;	// Domain length
	double D = 1.0;		// reference length for non-dimensionalization

	int i,j,k,n;

	if (rank==0){
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		ifstream domain("Domain.in");
		domain >> nprocx;
		domain >> nprocy;
		domain >> nprocz;
		domain >> Nx;
		domain >> Ny;
		domain >> Nz;
		domain >> ndiscs;
		domain >> Lx;
		domain >> Ly;
		domain >> Lz;
		//.......................................................................
	}
	// **************************************************************
	// Broadcast simulation parameters from rank 0 to all other procs
	MPI_Barrier(comm);
	//.................................................
	// Computational domain
	MPI_Bcast(&Nx,1,MPI_INT,0,comm);
	MPI_Bcast(&Ny,1,MPI_INT,0,comm);
	MPI_Bcast(&Nz,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	MPI_Bcast(&ndiscs,1,MPI_INT,0,comm);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
	//.................................................
	MPI_Barrier(comm);

	// **************************************************************
	if (argc > 1)	depth=atoi(argv[1]);
	else			depth=Nz;
	
	if (nprocs != nprocx*nprocy*nprocz){
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
		INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
	}

	 InitializeRanks( rank, nprocx, nprocy, nprocz, iproc, jproc, kproc, 
			 	 	 rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z,
			 	 	 rank_xy, rank_XY, rank_xY, rank_Xy, rank_xz, rank_XZ, rank_xZ, rank_Xz,
			 	 	 rank_yz, rank_YZ, rank_yZ, rank_Yz );
	 
	 MPI_Barrier(comm);

	Nx += 2;
	Ny += 2;
	Nz += 2;
	//Nx = Ny = Nz;	// Cubic domain

	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

	if (rank==0){
		printf("Process grid = %ix%ix%i \n", nprocx,nprocy,nprocz);
		printf("Sub-domain size = %ix%ix%i \n", Nx,Ny,Nz);
		printf("Physical domain size = %fx%fx%f \n",Lx,Ly,Lz);
		printf("Micromodel depth = %i voxels",depth);
	}


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
	if (rank==1) printf("ndiscs =%i \n",ndiscs);
	//.......................................................................
	double *cx,*cy,*cz,*rad;
	cx = new double[ndiscs];
	cy = new double[ndiscs];
	rad = new double[ndiscs];
	//.......................................................................
	if (rank == 0)	printf("Reading the disc packing \n");
	if (rank == 0)	ReadDiscPacking(ndiscs,cx,cy,rad);
	MPI_Barrier(comm);
	// Broadcast the sphere packing to all processes
	MPI_Bcast(cx,ndiscs,MPI_DOUBLE,0,comm);
	MPI_Bcast(cy,ndiscs,MPI_DOUBLE,0,comm);
	MPI_Bcast(rad,ndiscs,MPI_DOUBLE,0,comm);
	//...........................................................................
	MPI_Barrier(comm);
	if (rank == 0){
		cout << "Domain set." << endl;
		printf("************ \n");
		printf("Discs are: \n");
		for (int disc=0; disc<ndiscs; disc++){
			printf("%f,%f,%f\n",cx[disc],cy[disc],rad[disc]);
		}
		printf("************ \n");
	}

	if (nprocx > 1 && rank==0) printf("Disc packs are 2D -- are you sure you want nprocx > 1? \n");
	//.......................................................................
	SignedDistanceDiscPack(SignDist.data(),ndiscs,cx,cy,rad,Lx,Ly,Lz,Nx,Ny,Nz,
					   iproc,jproc,kproc,nprocx,nprocy,nprocz);
	//.......................................................................
	// Assign walls in the signed distance functions (x,y boundaries)
	double dst;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				// Assign the micromodel depth
				dst = 1.0*(i-2);
				if (1.0*(depth+2-i) < dst) dst = 1.0*(depth+2-i);
				if (dst < SignDist(i,j,k)) 		SignDist(i,j,k) = dst;
//				dst = (iproc*(Nx-2)+i-1)*1.0;
//				if ((Nx-2)*nprocx-2-iproc*(Nx-2)-i+1 < dst) 		dst = 1.0*((Nx-2)*nprocx-2-iproc*(Nx-2)-i+1);
				// Add walls at the boundary
				if ( (jproc*(Ny-2)+j-2)*1.0 < dst) 				dst = (jproc*(Ny-2)+j-2)*1.0;
				if ((Ny-2)*nprocy-3-(jproc*(Ny-2)+j-1)*1.0 < dst) 	dst = ((Ny-2)*nprocy-3-(jproc*(Ny-2)+j-1))*1.0;
				// Assign the Signed Distance where valid
				if (dst < SignDist(i,j,k)) 			SignDist(i,j,k) = dst;

				n = k*Nx*Ny+j*Nx+i;
				id[n] = 0;
			}
		}
	}

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


        sprintf(LocalRankFilename,"ID.%05i",rank);
        FILE *ID = fopen(LocalRankFilename,"wb");
        fwrite(id,1,N,ID);
        fclose(ID);


	//.......................................................................
	sprintf(LocalRankString,"%05d",rank);
	sprintf(LocalRankFilename,"%s%s","SignDist.",LocalRankString);
	WriteLocalSolidDistance(LocalRankFilename, SignDist.data(), N);
	//......................................................................

	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
