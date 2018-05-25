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


std::shared_ptr<Database> loadInputs( )
{
    auto db = std::make_shared<Database>( "Domain.in" );
    const int dim = 50;
    db->putScalar<int>( "BC", 0 );
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
	{
	//*****************************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************

	double TubeRadius =15.0;
	int BC;
	int BubbleTop,BubbleBottom;
	TubeRadius=strtod(argv[1],NULL);
	BC=atoi(argv[2]);
	BubbleBottom = atoi(argv[3]);
	BubbleTop = atoi(argv[4]);

	if (rank == 0){
		printf("********************************************************\n");
		printf("Generate 3D cylindrical capillary tube geometry with radius = %f voxels \n",TubeRadius);
		printf("********************************************************\n");
	}

	// Variables that specify the computational domain  
	string FILENAME;
	int i,j,k,n;

	// pmmc threshold values

	
	// **************************************************************
        // Load inputs
        auto db = loadInputs( );
        int Nx = db->getVector<int>( "n" )[0];
        int Ny = db->getVector<int>( "n" )[1];
        int Nz = db->getVector<int>( "n" )[2];
        int nprocx = db->getVector<int>( "nproc" )[0];
        int nprocy = db->getVector<int>( "nproc" )[1];
        int nprocz = db->getVector<int>( "nproc" )[2];
		int kproc = rank/(nprocx*nprocy);
		int jproc = (rank-nprocx*nprocy*kproc)/nprocx;
		int iproc = rank-nprocx*nprocy*kproc-nprocx*jproc;
		
	if (rank==0){
		printf("********************************************************\n");
		printf("Sub-domain size = %i x %i x %i\n",Nz,Nz,Nz);
		printf("Parallel domain size = %i x %i x %i\n",nprocx,nprocy,nprocz);
		printf("********************************************************\n");
	}

	double Lx=1.f;
	double Ly=1.f;
	double Lz=1.f;
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
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				n = k*Nx*Ny + j*Nz + i;
				// Cylindrical capillary tube aligned with the z direction
				Averages->SDs(i,j,k) = TubeRadius-sqrt(1.0*((i-Nx/2)*(i-Nx/2)
									+ (j-Ny/2)*(j-Ny/2)));

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
