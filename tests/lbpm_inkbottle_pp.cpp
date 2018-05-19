#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "analysis/TwoPhase.h"
#include "common/MPI_Helpers.h"

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
	{
	// parallel domain size (# of sub-domains)
	int nprocx,nprocy,nprocz;
	int iproc,jproc,kproc;
	//*****************************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************

	double UpperTubeRadius =15.0;
	double LowerTubeRadius =15.0;
	int BC;
	int BubbleTop,BubbleBottom;
	double BulbRadius;
	LowerTubeRadius=strtod(argv[1],NULL);
	LowerTubeRadius=strtod(argv[2],NULL);
	BulbRadius=strtod(argv[3],NULL);
	BubbleBottom = atoi(argv[4]);
	BubbleTop = atoi(argv[5]);

	BC=0;

	if (rank == 0){
		printf("********************************************************\n");
		printf("Generate ink bottle geometry");
		printf(	" lower tube radius = %f, upper tube radius=%f, bulb radius %f voxels \n",LowerTubeRadius,UpperTubeRadius,BulbRadius);
		printf("********************************************************\n");
	}

	// Variables that specify the computational domain  
	string FILENAME;
	int Nx,Ny,Nz;		// local sub-domain size
	int nspheres;		// number of spheres in the packing
	double Lx,Ly,Lz;	// Domain length
	int i,j,k,n;

	// pmmc threshold values

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
		domain >> nspheres;
		domain >> Lx;
		domain >> Ly;
		domain >> Lz;
		//.......................................................................
	}
	// **************************************************************
	// Broadcast simulation parameters from rank 0 to all other procs
	MPI_Barrier(comm);
	// Computational domain
	MPI_Bcast(&Nx,1,MPI_INT,0,comm);
	MPI_Bcast(&Ny,1,MPI_INT,0,comm);
	MPI_Bcast(&Nz,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocx,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocy,1,MPI_INT,0,comm);
	MPI_Bcast(&nprocz,1,MPI_INT,0,comm);
	MPI_Bcast(&nspheres,1,MPI_INT,0,comm);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,comm);
	MPI_Bcast(&Lz,1,MPI_DOUBLE,0,comm);
	//.................................................
	MPI_Barrier(comm);
	
	// **************************************************************
	if (nprocs != nprocx*nprocy*nprocz){
		printf("nprocx =  %i \n",nprocx);
		printf("nprocy =  %i \n",nprocy);
		printf("nprocz =  %i \n",nprocz);
		INSIST(nprocs == nprocx*nprocy*nprocz,"Fatal error in processor count!");
	}

	if (rank==0){
		printf("********************************************************\n");
		printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
		printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
		printf("********************************************************\n");
	}

	// Initialized domain and averaging framework for Two-Phase Flow
        std::shared_ptr<Domain> Dm (new Domain(Nx,Ny,Nz,rank,nprocx,nprocy,nprocz,Lx,Ly,Lz,BC));
        Dm->CommInit();
        std::shared_ptr<TwoPhase> Averages( new TwoPhase(Dm) );

	InitializeRanks( rank, nprocx, nprocy, nprocz, iproc, jproc, kproc,
			 	 	 rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z,
			 	 	 rank_xy, rank_XY, rank_xY, rank_Xy, rank_xz, rank_XZ, rank_xZ, rank_Xz,
			 	 	 rank_yz, rank_YZ, rank_yZ, rank_Yz );
	 
	MPI_Barrier(comm);

	Nz += 2;
	Nx = Ny = Nz;	// Cubic domain

	int N = Nx*Ny*Nz;
	int dist_mem_size = N*sizeof(double);

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
	//if (pBC) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
	double pore_vol;

	sum=0;
	double BulbDist;
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				// Cylindrical capillary tube aligned with the z direction
				if (k<Nz/2)
					Averages->SDs(i,j,k) = LowerTubeRadius-sqrt(1.0*((i-Nx/2)*(i-Nx/2)
									+ (j-Ny/2)*(j-Ny/2)));
				else
					Averages->SDs(i,j,k) = UpperTubeRadius-sqrt(1.0*((i-Nx/2)*(i-Nx/2)
									+ (j-Ny/2)*(j-Ny/2)));

				BulbDist = BulbRadius-sqrt(1.0*((i-Nx/2)*(i-Nx/2)+ (j-Ny/2)*(j-Ny/2) + (k-Nz/2)*(k-Nz/2)));

				if (BulbDist > Averages->SDs(i,j,k)) Averages->SDs(i,j,k) = BulbDist;

				// Initialize phase positions
				if (Averages->SDs(i,j,k) < 0.0){
					id[n] = 0;
				}
				else if (Dm->kproc()*Nz+k<BubbleBottom){
					id[n] = 2;
					sum++;
				}
				else if (Dm->kproc()*Nz+k<BubbleTop){
					id[n] = 1;
					sum++;
				}
				else{
					id[n] = 2;
					sum++;
				}
			}
		}
	}
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

    sprintf(LocalRankFilename,"SignDist.%05i",rank);
    FILE *DIST = fopen(LocalRankFilename,"wb");
    fwrite(Averages->SDs.data(),8,Averages->SDs.length(),DIST);
    fclose(DIST);

	sprintf(LocalRankFilename,"ID.%05i",rank);
	FILE *ID = fopen(LocalRankFilename,"wb");
	fwrite(id,1,N,ID);
	fclose(ID);

	}
	// ****************************************************
	MPI_Barrier(comm);
	MPI_Finalize();
	// ****************************************************
}
