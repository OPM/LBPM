// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "analysis/pmmc.h"
#ifdef PROFILE
 #include "ProfilerApp.h"
#endif

//#include "Domain.h"

using namespace std;


inline void ReadBinaryFile(char *FILENAME, double *Data, int N)
{
    int n;
    double value;
    ifstream File(FILENAME,ios::binary);
    for (n=0; n<N; n++){
        // Write the two density values
        File.read((char*) &value, sizeof(value));
        Data[n] = value;

    }
    File.close();
}


void readRankData( int proc, int nx, int ny, int nz, DoubleArray& Phase, DoubleArray& SignDist )
{
    Phase.resize(nx,ny,nz);
    SignDist.resize(nx,ny,nz);
    char file1[40], file2[40];
    sprintf(file1,"SignDist.%05d",proc);
    sprintf(file2,"Phase.%05d",proc);
    ReadBinaryFile(file1, Phase.data(), nx*ny*nz);
    ReadBinaryFile(file2, SignDist.data(), nx*ny*nz);
}


int main(int argc, char **argv)
{
	// Initialize MPI
	int rank, nprocs;
	MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nprocs);
#ifdef PROFILE
	PROFILE_ENABLE(0);
    PROFILE_DISABLE_TRACE();
    PROFILE_SYNCHRONIZE();
    PROFILE_START("main");
#endif

    if ( rank==0 ) {
        printf("-----------------------------------------------------------\n");
        printf("Labeling Blobs from Two-Phase Lattice Boltzmann Simulation \n");
        printf("-----------------------------------------------------------\n");
    }

    //.......................................................................
    // Reading the domain information file
    //.......................................................................
    int nprocx, nprocy, nprocz, nx, ny, nz, nspheres;
    double Lx, Ly, Lz;
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

    // Check that the number of processors >= the number of ranks
    if ( rank==0 ) {
        printf("Number of MPI ranks required: %i \n", nprocx*nprocy*nprocz);
        printf("Number of MPI ranks used: %i \n", nprocs);
        printf("Full domain size: %i x %i x %i  \n",nx*nprocx,ny*nprocy,nz*nprocz);
    }
    if ( nprocs < nprocx*nprocy*nprocz )
        ERROR("Insufficient number of processors");

    // Get the rank info
    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

    // Read the local file
    DoubleArray Phase;
    DoubleArray SignDist;
    readRankData( rank, nx+2, ny+2, nz+2, Phase, SignDist );

    // Communication the halos
    fillHalo<double> fillData(comm,rank_info,{nx,ny,nz},{1,1,1},0,1);
    fillData.fill(Phase);
    fillData.fill(SignDist);

    // Find blob domains
    if ( rank==0 ) { printf("Finding blob domains\n"); }
    double vF=0.0;
    double vS=0.0;
    IntArray GlobalBlobID;
    int nblobs = ComputeGlobalBlobIDs(nx,ny,nz,rank_info,
        Phase,SignDist,vF,vS,GlobalBlobID,comm);
    if ( rank==0 ) { printf("Identified %i blobs\n",nblobs); }

    // Write the local blob ids
    char LocalRankFilename[100];
    sprintf(LocalRankFilename,"BlobLabel.%05i",rank);
    FILE *BLOBLOCAL = fopen(LocalRankFilename,"wb");
    fwrite(GlobalBlobID.data(),4,GlobalBlobID.length(),BLOBLOCAL);
    fclose(BLOBLOCAL);
    printf("Wrote BlobLabel.%05i \n",rank);


    /*FILE *BLOBS = fopen("Blobs.dat","wb");
    fwrite(GlobalBlobID.data(),4,Nx*Ny*Nz,BLOBS);
    fclose(BLOBS);*/
#ifdef PROFILE
    PROFILE_STOP("main");
    PROFILE_SAVE("BlobIdentifyParallel",false);
#endif
    MPI_Barrier(comm);
	MPI_Finalize();
    return 0;  
}

