#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "common/Array.h"

using namespace std;



//***************************************************************************************

int test_communication( MPI_Comm comm, int nprocx, int nprocy, int nprocz )
{
    int rank,nprocs;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    int iproc,jproc,kproc;
    int sendtag,recvtag;
    if (rank==0)    printf("\nRunning test %i %i %i\n",nprocx,nprocy,nprocz);
    //*****************************************
    // MPI ranks for all 18 neighbors
    //**********************************
    int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
    int rank_xy,rank_XY,rank_xY,rank_Xy;
    int rank_xz,rank_XZ,rank_xZ,rank_Xz;
    int rank_yz,rank_YZ,rank_yZ,rank_Yz;
    //**********************************
    InitializeRanks( rank, nprocx, nprocy, nprocz,
	    iproc, jproc, kproc, 
        rank_x, rank_y, rank_z, 
	    rank_X, rank_Y, rank_Z,
	    rank_xy, rank_XY, rank_xY, rank_Xy,
	    rank_xz, rank_XZ, rank_xZ, rank_Xz,
	    rank_yz, rank_YZ, rank_yZ, rank_Yz );
    MPI_Barrier(comm);

    //**********************************

    // Set up MPI communication structurese
    if (rank==0)    printf ("Setting up communication control structures \n");
    //......................................................................................
    // Get the actual D3Q19 communication counts (based on location of solid phase)
    // Discrete velocity set symmetry implies the sendcount = recvcount
    int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
    int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
    int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
    sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
    sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
    sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;

    //......................................................................................
    int *sendList_x, *sendList_y, *sendList_z, *sendList_X, *sendList_Y, *sendList_Z;
    int *sendList_xy, *sendList_yz, *sendList_xz, *sendList_Xy, *sendList_Yz, *sendList_xZ;
    int *sendList_xY, *sendList_yZ, *sendList_Xz, *sendList_XY, *sendList_YZ, *sendList_XZ;
    //......................................................................................
    // send buffers
    sendList_x = new int [sendCount_x];
    sendList_y = new int [sendCount_y];
    sendList_z = new int [sendCount_z];
    sendList_X = new int [sendCount_X];
    sendList_Y = new int [sendCount_Y];
    sendList_Z = new int [sendCount_Z];
    sendList_xy = new int [sendCount_xy];
    sendList_yz = new int [sendCount_yz];
    sendList_xz = new int [sendCount_xz];
    sendList_Xy = new int [sendCount_Xy];
    sendList_Yz = new int [sendCount_Yz];
    sendList_xZ = new int [sendCount_xZ];
    sendList_xY = new int [sendCount_xY];
    sendList_yZ = new int [sendCount_yZ];
    sendList_Xz = new int [sendCount_Xz];
    sendList_XY = new int [sendCount_XY];
    sendList_YZ = new int [sendCount_YZ];
    sendList_XZ = new int [sendCount_XZ];
    if (rank==0)    printf ("Preparing the sendlists \n");
    //......................................................................................
    // Populate the send list
    sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
    sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
    sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;

    MPI_Barrier(comm);
    if (rank==0)    printf ("SendLists are ready on host\n");
    //......................................................................................
    // Use MPI to fill in the recvCounts form the associated processes
    int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
    int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
    int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
    //......................................................................................
    //**********************************************************************************
    // Fill in the recieve counts using MPI
    sendtag = recvtag = 3;
    CommunicateSendRecvCounts( comm, sendtag, recvtag, 
        rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z,
        rank_xy, rank_XY, rank_xY, rank_Xy,
        rank_xz, rank_XZ, rank_xZ, rank_Xz,
        rank_yz, rank_YZ, rank_yZ, rank_Yz,
        sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z,
        sendCount_xy, sendCount_XY, sendCount_xY, sendCount_Xy,
        sendCount_xz, sendCount_XZ, sendCount_xZ, sendCount_Xz,
        sendCount_yz, sendCount_YZ, sendCount_yZ, sendCount_Yz,
        recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z,
        recvCount_xy, recvCount_XY, recvCount_xY, recvCount_Xy,
        recvCount_xz, recvCount_XZ, recvCount_xZ, recvCount_Xz,
        recvCount_yz, recvCount_YZ, recvCount_yZ, recvCount_Yz );

    //......................................................................................
    int *recvList_x, *recvList_y, *recvList_z, *recvList_X, *recvList_Y, *recvList_Z;
    int *recvList_xy, *recvList_yz, *recvList_xz, *recvList_Xy, *recvList_Yz, *recvList_xZ;
    int *recvList_xY, *recvList_yZ, *recvList_Xz, *recvList_XY, *recvList_YZ, *recvList_XZ;
    //......................................................................................
    // recv buffers
    recvList_x = new int [recvCount_x];
    recvList_y = new int [recvCount_y];
    recvList_z = new int [recvCount_z];
    recvList_X = new int [recvCount_X];
    recvList_Y = new int [recvCount_Y];
    recvList_Z = new int [recvCount_Z];
    recvList_xy = new int [recvCount_xy];
    recvList_yz = new int [recvCount_yz];
    recvList_xz = new int [recvCount_xz];
    recvList_Xy = new int [recvCount_Xy];
    recvList_Yz = new int [recvCount_Yz];
    recvList_xZ = new int [recvCount_xZ];
    recvList_xY = new int [recvCount_xY];
    recvList_yZ = new int [recvCount_yZ];
    recvList_Xz = new int [recvCount_Xz];
    recvList_XY = new int [recvCount_XY];
    recvList_YZ = new int [recvCount_YZ];
    recvList_XZ = new int [recvCount_XZ];
    //......................................................................................
    //......................................................................................
    // Use MPI to fill in the appropriate values for recvList
    // Fill in the recieve lists using MPI
    sendtag = recvtag = 4;
    CommunicateRecvLists( comm, sendtag, recvtag, 
        sendList_x, sendList_y, sendList_z, sendList_X, sendList_Y, sendList_Z,
        sendList_xy, sendList_XY, sendList_xY, sendList_Xy,
        sendList_xz, sendList_XZ, sendList_xZ, sendList_Xz,
        sendList_yz, sendList_YZ, sendList_yZ, sendList_Yz,
        sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z,
        sendCount_xy, sendCount_XY, sendCount_xY, sendCount_Xy,
        sendCount_xz, sendCount_XZ, sendCount_xZ, sendCount_Xz,
        sendCount_yz, sendCount_YZ, sendCount_yZ, sendCount_Yz,
        recvList_x, recvList_y, recvList_z, recvList_X, recvList_Y, recvList_Z,
        recvList_xy, recvList_XY, recvList_xY, recvList_Xy,
        recvList_xz, recvList_XZ, recvList_xZ, recvList_Xz,
        recvList_yz, recvList_YZ, recvList_yZ, recvList_Yz,
        recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z,
        recvCount_xy, recvCount_XY, recvCount_xY, recvCount_Xy,
        recvCount_xz, recvCount_XZ, recvCount_xZ, recvCount_Xz,
        recvCount_yz, recvCount_YZ, recvCount_yZ, recvCount_Yz,
        rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z, rank_xy, rank_XY, rank_xY,
        rank_Xy, rank_xz, rank_XZ, rank_xZ, rank_Xz, rank_yz, rank_YZ, rank_yZ, rank_Yz );
    MPI_Barrier(comm);
    if (rank==0)    printf ("RecvLists finished\n");
    
    // Free memory
    delete [] sendList_x,  delete [] sendList_y,  delete [] sendList_z;
    delete [] sendList_X,  delete [] sendList_Y,  delete [] sendList_Z;
    delete [] sendList_xy, delete [] sendList_xz, delete [] sendList_yz;
    delete [] sendList_xY, delete [] sendList_xZ, delete [] sendList_yZ;
    delete [] sendList_Xy, delete [] sendList_Xz, delete [] sendList_Yz;
    delete [] sendList_XY, delete [] sendList_XZ, delete [] sendList_YZ;
    delete [] recvList_x,  delete [] recvList_y,  delete [] recvList_z;
    delete [] recvList_X,  delete [] recvList_Y,  delete [] recvList_Z;
    delete [] recvList_xy, delete [] recvList_xz, delete [] recvList_yz;
    delete [] recvList_xY, delete [] recvList_xZ, delete [] recvList_yZ;
    delete [] recvList_Xy, delete [] recvList_Xz, delete [] recvList_Yz;
    delete [] recvList_XY, delete [] recvList_XZ, delete [] recvList_YZ;

    // Finished with no errors
    return 0;
}


template<class TYPE>
int testHalo( MPI_Comm comm, int nprocx, int nprocy, int nprocz, int depth )
{
    int rank,nprocs;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    if ( rank==0 )
        printf("\nRunning Halo test %i %i %i %i\n",nprocx,nprocy,nprocz,depth);

    const RankInfoStruct rank_info(rank,nprocx,nprocy,nprocz);

    int nx = 10;
    int ny = 11;
    int nz = 7;
    std::vector<size_t> size(4);
    size[0] = nx+2;
    size[1] = ny+2;
    size[2] = nz+2;
    size[3] = depth;
    Array<TYPE> array(size);
    array.fill(-1);
    
    // Fill the local array
    int Nx = nx*nprocx;
    int Ny = ny*nprocy;
    int Nz = nz*nprocz;
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int k=0; k<nz; k++) {
                for (int d=0; d<depth; d++) {
                    int iglobal = i + rank_info.ix*nx;
                    int jglobal = j + rank_info.jy*ny;
                    int kglobal = k + rank_info.kz*nz;
                    int ijk = iglobal + jglobal*Nx + kglobal*Nx*Ny + d*Nx*Ny*Nz;
                    array(i+1,j+1,k+1,d) = ijk;
                }
            }
        }
    }

    // Communicate the halo
    fillHalo<TYPE> fillData(comm,rank_info,{nx,ny,nz},{1,1,1},0,depth);
    fillData.fill(array);

    // Check the results
    bool pass = true;
    for (int i=-1; i<nx+1; i++) {
        for (int j=-1; j<ny+1; j++) {
            for (int k=-1; k<nz+1; k++) {
                for (int d=0; d<depth; d++) {
                    int iglobal = i + rank_info.ix*nx;
                    int jglobal = j + rank_info.jy*ny;
                    int kglobal = k + rank_info.kz*nz;
                    iglobal = (iglobal+Nx)%Nx;
                    jglobal = (jglobal+Ny)%Ny;
                    kglobal = (kglobal+Nz)%Nz;
                    int ijk = iglobal + jglobal*Nx + kglobal*Nx*Ny + d*Nx*Ny*Nz;
                    if ( array(i+1,j+1,k+1,d) != ijk )
                        pass = false;
                }
            }
        }
    }
    int N_errors = 0;
    if ( !pass ) {
        std::cout << "Failed halo test\n";
        N_errors++;
    }
    return N_errors;
}


int main(int argc, char **argv)
{
    // Initialize MPI
    int rank,nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);

    // Run the test with different domains
    int N_errors = 0;
    N_errors += test_communication( comm, nprocs, 1, 1 );
    N_errors += test_communication( comm, 1, nprocs, 1 );
    N_errors += test_communication( comm, 1, 1, nprocs );
    if ( nprocs==4 ) {
        N_errors += test_communication( comm, 2, 2, 1 );
        N_errors += test_communication( comm, 2, 1, 2 );
        N_errors += test_communication( comm, 1, 2, 2 );
    }

    // Run the halo tests with different domains
    N_errors += testHalo<int>( comm, nprocs, 1, 1, 1 );
    N_errors += testHalo<int>( comm, 1, nprocs, 1, 1 );
    N_errors += testHalo<int>( comm, 1, 1, nprocs, 1 );
    N_errors += testHalo<double>( comm, nprocs, 1, 1, 3 );
    N_errors += testHalo<double>( comm, 1, nprocs, 1, 3 );
    N_errors += testHalo<double>( comm, 1, 1, nprocs, 3 );
    if ( nprocs==4 ) {
        N_errors += testHalo<int>( comm, 2, 2, 1, 1 );
        N_errors += testHalo<int>( comm, 2, 1, 2, 1 );
        N_errors += testHalo<int>( comm, 1, 2, 2, 1 );
    }
    if ( nprocs==8 ) {
        N_errors += testHalo<int>( comm, 2, 2, 2, 1 );
    }

    // Finished
    MPI_Barrier(comm);
    int N_errors_global=0;
    MPI_Allreduce( &N_errors, &N_errors_global, 1, MPI_INT, MPI_SUM, comm );
    MPI_Barrier(comm);
    MPI_Finalize();
    if ( rank==0 ) {
        if ( N_errors_global==0 )
            std::cout << "All tests passed\n";
        else
            std::cout << "Some tests failed\n";
    }
    return N_errors_global;
}
