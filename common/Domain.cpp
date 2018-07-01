// Created by James McClure
// Copyright 2008-2013
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <exception>      // std::exception
#include <stdexcept>

#include "common/Domain.h"
#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"

// Inline function to read line without a return argument
static inline void fgetl( char * str, int num, FILE * stream )
{
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}

/********************************************************
 * Constructors/Destructor                               *
 ********************************************************/
Domain::Domain( int nx, int ny, int nz, int rnk, int npx, int npy, int npz, 
               double lx, double ly, double lz, int BC):
	Nx(0), Ny(0), Nz(0), 
	Lx(0), Ly(0), Lz(0), Volume(0), BoundaryCondition(0),
	Comm(MPI_COMM_WORLD),
	sendCount_x(0), sendCount_y(0), sendCount_z(0), sendCount_X(0), sendCount_Y(0), sendCount_Z(0),
	sendCount_xy(0), sendCount_yz(0), sendCount_xz(0), sendCount_Xy(0), sendCount_Yz(0), sendCount_xZ(0),
	sendCount_xY(0), sendCount_yZ(0), sendCount_Xz(0), sendCount_XY(0), sendCount_YZ(0), sendCount_XZ(0),
	sendList_x(NULL), sendList_y(NULL), sendList_z(NULL), sendList_X(NULL), sendList_Y(NULL), sendList_Z(NULL),
	sendList_xy(NULL), sendList_yz(NULL), sendList_xz(NULL), sendList_Xy(NULL), sendList_Yz(NULL), sendList_xZ(NULL),
	sendList_xY(NULL), sendList_yZ(NULL), sendList_Xz(NULL), sendList_XY(NULL), sendList_YZ(NULL), sendList_XZ(NULL),
	sendBuf_x(NULL), sendBuf_y(NULL), sendBuf_z(NULL), sendBuf_X(NULL), sendBuf_Y(NULL), sendBuf_Z(NULL),
	sendBuf_xy(NULL), sendBuf_yz(NULL), sendBuf_xz(NULL), sendBuf_Xy(NULL), sendBuf_Yz(NULL), sendBuf_xZ(NULL),
	sendBuf_xY(NULL), sendBuf_yZ(NULL), sendBuf_Xz(NULL), sendBuf_XY(NULL), sendBuf_YZ(NULL), sendBuf_XZ(NULL),
	recvCount_x(0), recvCount_y(0), recvCount_z(0), recvCount_X(0), recvCount_Y(0), recvCount_Z(0),
	recvCount_xy(0), recvCount_yz(0), recvCount_xz(0), recvCount_Xy(0), recvCount_Yz(0), recvCount_xZ(0),
	recvCount_xY(0), recvCount_yZ(0), recvCount_Xz(0), recvCount_XY(0), recvCount_YZ(0), recvCount_XZ(0),
	recvList_x(NULL), recvList_y(NULL), recvList_z(NULL), recvList_X(NULL), recvList_Y(NULL), recvList_Z(NULL),
	recvList_xy(NULL), recvList_yz(NULL), recvList_xz(NULL), recvList_Xy(NULL), recvList_Yz(NULL), recvList_xZ(NULL),
	recvList_xY(NULL), recvList_yZ(NULL), recvList_Xz(NULL), recvList_XY(NULL), recvList_YZ(NULL), recvList_XZ(NULL),
	recvBuf_x(NULL), recvBuf_y(NULL), recvBuf_z(NULL), recvBuf_X(NULL), recvBuf_Y(NULL), recvBuf_Z(NULL),
	recvBuf_xy(NULL), recvBuf_yz(NULL), recvBuf_xz(NULL), recvBuf_Xy(NULL), recvBuf_Yz(NULL), recvBuf_xZ(NULL),
	recvBuf_xY(NULL), recvBuf_yZ(NULL), recvBuf_Xz(NULL), recvBuf_XY(NULL), recvBuf_YZ(NULL), recvBuf_XZ(NULL),
	sendData_x(NULL), sendData_y(NULL), sendData_z(NULL), sendData_X(NULL), sendData_Y(NULL), sendData_Z(NULL),
	sendData_xy(NULL), sendData_yz(NULL), sendData_xz(NULL), sendData_Xy(NULL), sendData_Yz(NULL), sendData_xZ(NULL),
	sendData_xY(NULL), sendData_yZ(NULL), sendData_Xz(NULL), sendData_XY(NULL), sendData_YZ(NULL), sendData_XZ(NULL),
	recvData_x(NULL), recvData_y(NULL), recvData_z(NULL), recvData_X(NULL), recvData_Y(NULL), recvData_Z(NULL),
	recvData_xy(NULL), recvData_yz(NULL), recvData_xz(NULL), recvData_Xy(NULL), recvData_Yz(NULL), recvData_xZ(NULL),
	recvData_xY(NULL), recvData_yZ(NULL), recvData_Xz(NULL), recvData_XY(NULL), recvData_YZ(NULL), recvData_XZ(NULL),
	id(NULL)
{	
	// set up the neighbor ranks
    int myrank;
    MPI_Comm_rank( Comm, &myrank );
	rank_info = RankInfoStruct( myrank, rank_info.nx, rank_info.ny, rank_info.nz );
	
	MPI_Barrier(Comm);
	
    auto db = std::make_shared<Database>( );
    db->putScalar<int>( "BC", BC );
    db->putVector<int>( "nproc", { npx, npx, npx } );
    db->putVector<int>( "n", { nx, ny, nz } );
    db->putScalar<int>( "nspheres", 0 );
    db->putVector<double>( "L", { lx, ly, lz } );
    initialize( db );
}
Domain::Domain( std::shared_ptr<Database> db, MPI_Comm Communicator):
	Nx(0), Ny(0), Nz(0), 
	Lx(0), Ly(0), Lz(0), Volume(0), BoundaryCondition(0),
	Comm(MPI_COMM_NULL),
	sendCount_x(0), sendCount_y(0), sendCount_z(0), sendCount_X(0), sendCount_Y(0), sendCount_Z(0),
	sendCount_xy(0), sendCount_yz(0), sendCount_xz(0), sendCount_Xy(0), sendCount_Yz(0), sendCount_xZ(0),
	sendCount_xY(0), sendCount_yZ(0), sendCount_Xz(0), sendCount_XY(0), sendCount_YZ(0), sendCount_XZ(0),
	sendList_x(NULL), sendList_y(NULL), sendList_z(NULL), sendList_X(NULL), sendList_Y(NULL), sendList_Z(NULL),
	sendList_xy(NULL), sendList_yz(NULL), sendList_xz(NULL), sendList_Xy(NULL), sendList_Yz(NULL), sendList_xZ(NULL),
	sendList_xY(NULL), sendList_yZ(NULL), sendList_Xz(NULL), sendList_XY(NULL), sendList_YZ(NULL), sendList_XZ(NULL),
	sendBuf_x(NULL), sendBuf_y(NULL), sendBuf_z(NULL), sendBuf_X(NULL), sendBuf_Y(NULL), sendBuf_Z(NULL),
	sendBuf_xy(NULL), sendBuf_yz(NULL), sendBuf_xz(NULL), sendBuf_Xy(NULL), sendBuf_Yz(NULL), sendBuf_xZ(NULL),
	sendBuf_xY(NULL), sendBuf_yZ(NULL), sendBuf_Xz(NULL), sendBuf_XY(NULL), sendBuf_YZ(NULL), sendBuf_XZ(NULL),
	recvCount_x(0), recvCount_y(0), recvCount_z(0), recvCount_X(0), recvCount_Y(0), recvCount_Z(0),
	recvCount_xy(0), recvCount_yz(0), recvCount_xz(0), recvCount_Xy(0), recvCount_Yz(0), recvCount_xZ(0),
	recvCount_xY(0), recvCount_yZ(0), recvCount_Xz(0), recvCount_XY(0), recvCount_YZ(0), recvCount_XZ(0),
	recvList_x(NULL), recvList_y(NULL), recvList_z(NULL), recvList_X(NULL), recvList_Y(NULL), recvList_Z(NULL),
	recvList_xy(NULL), recvList_yz(NULL), recvList_xz(NULL), recvList_Xy(NULL), recvList_Yz(NULL), recvList_xZ(NULL),
	recvList_xY(NULL), recvList_yZ(NULL), recvList_Xz(NULL), recvList_XY(NULL), recvList_YZ(NULL), recvList_XZ(NULL),
	recvBuf_x(NULL), recvBuf_y(NULL), recvBuf_z(NULL), recvBuf_X(NULL), recvBuf_Y(NULL), recvBuf_Z(NULL),
	recvBuf_xy(NULL), recvBuf_yz(NULL), recvBuf_xz(NULL), recvBuf_Xy(NULL), recvBuf_Yz(NULL), recvBuf_xZ(NULL),
	recvBuf_xY(NULL), recvBuf_yZ(NULL), recvBuf_Xz(NULL), recvBuf_XY(NULL), recvBuf_YZ(NULL), recvBuf_XZ(NULL),
	sendData_x(NULL), sendData_y(NULL), sendData_z(NULL), sendData_X(NULL), sendData_Y(NULL), sendData_Z(NULL),
	sendData_xy(NULL), sendData_yz(NULL), sendData_xz(NULL), sendData_Xy(NULL), sendData_Yz(NULL), sendData_xZ(NULL),
	sendData_xY(NULL), sendData_yZ(NULL), sendData_Xz(NULL), sendData_XY(NULL), sendData_YZ(NULL), sendData_XZ(NULL),
	recvData_x(NULL), recvData_y(NULL), recvData_z(NULL), recvData_X(NULL), recvData_Y(NULL), recvData_Z(NULL),
	recvData_xy(NULL), recvData_yz(NULL), recvData_xz(NULL), recvData_Xy(NULL), recvData_Yz(NULL), recvData_xZ(NULL),
	recvData_xY(NULL), recvData_yZ(NULL), recvData_Xz(NULL), recvData_XY(NULL), recvData_YZ(NULL), recvData_XZ(NULL),
	id(NULL)
{
    MPI_Comm_dup(Communicator,&Comm);

	// set up the neighbor ranks
    int myrank;
    MPI_Comm_rank( Comm, &myrank );
    initialize( db );
	rank_info = RankInfoStruct( myrank, rank_info.nx, rank_info.ny, rank_info.nz );
	MPI_Barrier(Comm);
}
void Domain::initialize( std::shared_ptr<Database> db )
{	
    d_db = db;
    auto nproc = d_db->getVector<int>("nproc");
    auto n = d_db->getVector<int>("n");
    auto L = d_db->getVector<double>("L");
    //nspheres = d_db->getScalar<int>("nspheres");
    ASSERT( n.size() == 3u );
    ASSERT( L.size() == 3u );
    ASSERT( nproc.size() == 3u );
    int nx = n[0];
    int ny = n[1];
    int nz = n[2];
    Lx = L[0];
    Ly = L[1];
    Lz = L[2];
    Nx = nx+2;
    Ny = ny+2;
    Nz = nz+2;
    // Initialize ranks
    int myrank;
    MPI_Comm_rank( Comm, &myrank );
	rank_info = RankInfoStruct(myrank,nproc[0],nproc[1],nproc[2]);
    // Fill remaining variables
	N = Nx*Ny*Nz;
	Volume = nx*ny*nx*nproc[0]*nproc[1]*nproc[2]*1.0;
	id = new char[N];
	memset(id,0,N);
	BoundaryCondition = d_db->getScalar<int>("BC");
    int nprocs;
    MPI_Comm_size( Comm, &nprocs );
	INSIST(nprocs == nproc[0]*nproc[1]*nproc[2],"Fatal error in processor count!");
}
Domain::~Domain()
{
	// Free sendList
	delete [] sendList_x;   delete [] sendList_y;   delete [] sendList_z;
	delete [] sendList_X;   delete [] sendList_Y;   delete [] sendList_Z;
	delete [] sendList_xy;  delete [] sendList_yz;  delete [] sendList_xz;
	delete [] sendList_Xy;  delete [] sendList_Yz;  delete [] sendList_xZ;
	delete [] sendList_xY;  delete [] sendList_yZ;  delete [] sendList_Xz;
	delete [] sendList_XY;  delete [] sendList_YZ;  delete [] sendList_XZ;
	// Free sendBuf
	delete [] sendBuf_x;    delete [] sendBuf_y;    delete [] sendBuf_z;
	delete [] sendBuf_X;    delete [] sendBuf_Y;    delete [] sendBuf_Z;
	delete [] sendBuf_xy;   delete [] sendBuf_yz;   delete [] sendBuf_xz;
	delete [] sendBuf_Xy;   delete [] sendBuf_Yz;   delete [] sendBuf_xZ;
	delete [] sendBuf_xY;   delete [] sendBuf_yZ;   delete [] sendBuf_Xz;
	delete [] sendBuf_XY;   delete [] sendBuf_YZ;   delete [] sendBuf_XZ;
	// Free recvList
	delete [] recvList_x;   delete [] recvList_y;   delete [] recvList_z;
	delete [] recvList_X;   delete [] recvList_Y;   delete [] recvList_Z;
	delete [] recvList_xy;  delete [] recvList_yz;  delete [] recvList_xz;
	delete [] recvList_Xy;  delete [] recvList_Yz;  delete [] recvList_xZ;
	delete [] recvList_xY;  delete [] recvList_yZ;  delete [] recvList_Xz;
	delete [] recvList_XY;  delete [] recvList_YZ;  delete [] recvList_XZ;
	// Free recvBuf
	delete [] recvBuf_x;    delete [] recvBuf_y;    delete [] recvBuf_z;
	delete [] recvBuf_X;    delete [] recvBuf_Y;    delete [] recvBuf_Z;
	delete [] recvBuf_xy;   delete [] recvBuf_yz;   delete [] recvBuf_xz;
	delete [] recvBuf_Xy;   delete [] recvBuf_Yz;   delete [] recvBuf_xZ;
	delete [] recvBuf_xY;   delete [] recvBuf_yZ;   delete [] recvBuf_Xz;
	delete [] recvBuf_XY;   delete [] recvBuf_YZ;   delete [] recvBuf_XZ;
	// Free sendData
	delete [] sendData_x;   delete [] sendData_y;   delete [] sendData_z;
	delete [] sendData_X;   delete [] sendData_Y;   delete [] sendData_Z;
	delete [] sendData_xy;  delete [] sendData_xY;  delete [] sendData_Xy;
	delete [] sendData_XY;  delete [] sendData_xz;  delete [] sendData_xZ;
	delete [] sendData_Xz;  delete [] sendData_XZ;  delete [] sendData_yz;
	delete [] sendData_yZ;  delete [] sendData_Yz;  delete [] sendData_YZ;
	// Free recvData
	delete [] recvData_x;   delete [] recvData_y;   delete [] recvData_z;
	delete [] recvData_X;   delete [] recvData_Y;   delete [] recvData_Z;
	delete [] recvData_xy;  delete [] recvData_xY;  delete [] recvData_Xy;
	delete [] recvData_XY;  delete [] recvData_xz;  delete [] recvData_xZ;
	delete [] recvData_Xz;  delete [] recvData_XZ;  delete [] recvData_yz;
	delete [] recvData_yZ;  delete [] recvData_Yz;  delete [] recvData_YZ;
	// Free id
	delete [] id;
	// Free the communicator
	if ( Comm != MPI_COMM_WORLD && Comm != MPI_COMM_NULL ) {
		MPI_Comm_free(&Comm);
	}
}


/********************************************************
 * Initialize communication                              *
 ********************************************************/
void Domain::CommInit()
{
	int i,j,k,n;
	int sendtag = 21;
	int recvtag = 21;
	//......................................................................................
	sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
	sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
	sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
	//......................................................................................
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				// Check the phase ID
				if (id[k*Nx*Ny+j*Nx+i] != 0){
					// Counts for the six faces
					if (i==1)    sendCount_x++;
					if (j==1)    sendCount_y++;
					if (k==1)    sendCount_z++;
					if (i==Nx-2)    sendCount_X++;
					if (j==Ny-2)    sendCount_Y++;
					if (k==Nz-2)    sendCount_Z++;
					// Counts for the twelve edges
					if (i==1 && j==1)    sendCount_xy++;
					if (i==1 && j==Ny-2)    sendCount_xY++;
					if (i==Nx-2 && j==1)    sendCount_Xy++;
					if (i==Nx-2 && j==Ny-2)    sendCount_XY++;

					if (i==1 && k==1)    sendCount_xz++;
					if (i==1 && k==Nz-2)    sendCount_xZ++;
					if (i==Nx-2 && k==1)    sendCount_Xz++;
					if (i==Nx-2 && k==Nz-2)    sendCount_XZ++;

					if (j==1 && k==1)    sendCount_yz++;
					if (j==1 && k==Nz-2)    sendCount_yZ++;
					if (j==Ny-2 && k==1)    sendCount_Yz++;
					if (j==Ny-2 && k==Nz-2)    sendCount_YZ++;
				}
			}
		}
	}

	// allocate send lists
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
	// Populate the send list
	sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
	sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
	sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				// Local value to send
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] != 0){
					// Counts for the six faces
					if (i==1)        sendList_x[sendCount_x++]=n;
					if (j==1)        sendList_y[sendCount_y++]=n;
					if (k==1)        sendList_z[sendCount_z++]=n;
					if (i==Nx-2)    sendList_X[sendCount_X++]=n;
					if (j==Ny-2)    sendList_Y[sendCount_Y++]=n;
					if (k==Nz-2)    sendList_Z[sendCount_Z++]=n;
					// Counts for the twelve edges
					if (i==1 && j==1)        sendList_xy[sendCount_xy++]=n;
					if (i==1 && j==Ny-2)    sendList_xY[sendCount_xY++]=n;
					if (i==Nx-2 && j==1)    sendList_Xy[sendCount_Xy++]=n;
					if (i==Nx-2 && j==Ny-2)    sendList_XY[sendCount_XY++]=n;

					if (i==1 && k==1)        sendList_xz[sendCount_xz++]=n;
					if (i==1 && k==Nz-2)    sendList_xZ[sendCount_xZ++]=n;
					if (i==Nx-2 && k==1)    sendList_Xz[sendCount_Xz++]=n;
					if (i==Nx-2 && k==Nz-2)    sendList_XZ[sendCount_XZ++]=n;

					if (j==1 && k==1)        sendList_yz[sendCount_yz++]=n;
					if (j==1 && k==Nz-2)    sendList_yZ[sendCount_yZ++]=n;
					if (j==Ny-2 && k==1)    sendList_Yz[sendCount_Yz++]=n;
					if (j==Ny-2 && k==Nz-2)    sendList_YZ[sendCount_YZ++]=n;
				}
			}
		}
	}

	// allocate send buffers
	sendBuf_x = new int [sendCount_x];
	sendBuf_y = new int [sendCount_y];
	sendBuf_z = new int [sendCount_z];
	sendBuf_X = new int [sendCount_X];
	sendBuf_Y = new int [sendCount_Y];
	sendBuf_Z = new int [sendCount_Z];
	sendBuf_xy = new int [sendCount_xy];
	sendBuf_yz = new int [sendCount_yz];
	sendBuf_xz = new int [sendCount_xz];
	sendBuf_Xy = new int [sendCount_Xy];
	sendBuf_Yz = new int [sendCount_Yz];
	sendBuf_xZ = new int [sendCount_xZ];
	sendBuf_xY = new int [sendCount_xY];
	sendBuf_yZ = new int [sendCount_yZ];
	sendBuf_Xz = new int [sendCount_Xz];
	sendBuf_XY = new int [sendCount_XY];
	sendBuf_YZ = new int [sendCount_YZ];
	sendBuf_XZ = new int [sendCount_XZ];
	//......................................................................................
	MPI_Isend(&sendCount_x, 1,MPI_INT,rank_x(),sendtag+0,Comm,&req1[0]);
	MPI_Irecv(&recvCount_X, 1,MPI_INT,rank_X(),recvtag+0,Comm,&req2[0]);
	MPI_Isend(&sendCount_X, 1,MPI_INT,rank_X(),sendtag+1,Comm,&req1[1]);
	MPI_Irecv(&recvCount_x, 1,MPI_INT,rank_x(),recvtag+1,Comm,&req2[1]);
	MPI_Isend(&sendCount_y, 1,MPI_INT,rank_y(),sendtag+2,Comm,&req1[2]);
	MPI_Irecv(&recvCount_Y, 1,MPI_INT,rank_Y(),recvtag+2,Comm,&req2[2]);
	MPI_Isend(&sendCount_Y, 1,MPI_INT,rank_Y(),sendtag+3,Comm,&req1[3]);
	MPI_Irecv(&recvCount_y, 1,MPI_INT,rank_y(),recvtag+3,Comm,&req2[3]);
	MPI_Isend(&sendCount_z, 1,MPI_INT,rank_z(),sendtag+4,Comm,&req1[4]);
	MPI_Irecv(&recvCount_Z, 1,MPI_INT,rank_Z(),recvtag+4,Comm,&req2[4]);
	MPI_Isend(&sendCount_Z, 1,MPI_INT,rank_Z(),sendtag+5,Comm,&req1[5]);
	MPI_Irecv(&recvCount_z, 1,MPI_INT,rank_z(),recvtag+5,Comm,&req2[5]);
	MPI_Isend(&sendCount_xy, 1,MPI_INT,rank_xy(),sendtag+6,Comm,&req1[6]);
	MPI_Irecv(&recvCount_XY, 1,MPI_INT,rank_XY(),recvtag+6,Comm,&req2[6]);
	MPI_Isend(&sendCount_XY, 1,MPI_INT,rank_XY(),sendtag+7,Comm,&req1[7]);
	MPI_Irecv(&recvCount_xy, 1,MPI_INT,rank_xy(),recvtag+7,Comm,&req2[7]);
	MPI_Isend(&sendCount_Xy, 1,MPI_INT,rank_Xy(),sendtag+8,Comm,&req1[8]);
	MPI_Irecv(&recvCount_xY, 1,MPI_INT,rank_xY(),recvtag+8,Comm,&req2[8]);
	MPI_Isend(&sendCount_xY, 1,MPI_INT,rank_xY(),sendtag+9,Comm,&req1[9]);
	MPI_Irecv(&recvCount_Xy, 1,MPI_INT,rank_Xy(),recvtag+9,Comm,&req2[9]);
	MPI_Isend(&sendCount_xz, 1,MPI_INT,rank_xz(),sendtag+10,Comm,&req1[10]);
	MPI_Irecv(&recvCount_XZ, 1,MPI_INT,rank_XZ(),recvtag+10,Comm,&req2[10]);
	MPI_Isend(&sendCount_XZ, 1,MPI_INT,rank_XZ(),sendtag+11,Comm,&req1[11]);
	MPI_Irecv(&recvCount_xz, 1,MPI_INT,rank_xz(),recvtag+11,Comm,&req2[11]);
	MPI_Isend(&sendCount_Xz, 1,MPI_INT,rank_Xz(),sendtag+12,Comm,&req1[12]);
	MPI_Irecv(&recvCount_xZ, 1,MPI_INT,rank_xZ(),recvtag+12,Comm,&req2[12]);
	MPI_Isend(&sendCount_xZ, 1,MPI_INT,rank_xZ(),sendtag+13,Comm,&req1[13]);
	MPI_Irecv(&recvCount_Xz, 1,MPI_INT,rank_Xz(),recvtag+13,Comm,&req2[13]);
	MPI_Isend(&sendCount_yz, 1,MPI_INT,rank_yz(),sendtag+14,Comm,&req1[14]);
	MPI_Irecv(&recvCount_YZ, 1,MPI_INT,rank_YZ(),recvtag+14,Comm,&req2[14]);
	MPI_Isend(&sendCount_YZ, 1,MPI_INT,rank_YZ(),sendtag+15,Comm,&req1[15]);
	MPI_Irecv(&recvCount_yz, 1,MPI_INT,rank_yz(),recvtag+15,Comm,&req2[15]);
	MPI_Isend(&sendCount_Yz, 1,MPI_INT,rank_Yz(),sendtag+16,Comm,&req1[16]);
	MPI_Irecv(&recvCount_yZ, 1,MPI_INT,rank_yZ(),recvtag+16,Comm,&req2[16]);
	MPI_Isend(&sendCount_yZ, 1,MPI_INT,rank_yZ(),sendtag+17,Comm,&req1[17]);
	MPI_Irecv(&recvCount_Yz, 1,MPI_INT,rank_Yz(),recvtag+17,Comm,&req2[17]);
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	MPI_Barrier(Comm);
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
	MPI_Isend(sendList_x, sendCount_x,MPI_INT,rank_x(),sendtag,Comm,&req1[0]);
	MPI_Irecv(recvList_X, recvCount_X,MPI_INT,rank_X(),recvtag,Comm,&req2[0]);
	MPI_Isend(sendList_X, sendCount_X,MPI_INT,rank_X(),sendtag,Comm,&req1[1]);
	MPI_Irecv(recvList_x, recvCount_x,MPI_INT,rank_x(),recvtag,Comm,&req2[1]);
	MPI_Isend(sendList_y, sendCount_y,MPI_INT,rank_y(),sendtag,Comm,&req1[2]);
	MPI_Irecv(recvList_Y, recvCount_Y,MPI_INT,rank_Y(),recvtag,Comm,&req2[2]);
	MPI_Isend(sendList_Y, sendCount_Y,MPI_INT,rank_Y(),sendtag,Comm,&req1[3]);
	MPI_Irecv(recvList_y, recvCount_y,MPI_INT,rank_y(),recvtag,Comm,&req2[3]);
	MPI_Isend(sendList_z, sendCount_z,MPI_INT,rank_z(),sendtag,Comm,&req1[4]);
	MPI_Irecv(recvList_Z, recvCount_Z,MPI_INT,rank_Z(),recvtag,Comm,&req2[4]);
	MPI_Isend(sendList_Z, sendCount_Z,MPI_INT,rank_Z(),sendtag,Comm,&req1[5]);
	MPI_Irecv(recvList_z, recvCount_z,MPI_INT,rank_z(),recvtag,Comm,&req2[5]);
	MPI_Isend(sendList_xy, sendCount_xy,MPI_INT,rank_xy(),sendtag,Comm,&req1[6]);
	MPI_Irecv(recvList_XY, recvCount_XY,MPI_INT,rank_XY(),recvtag,Comm,&req2[6]);
	MPI_Isend(sendList_XY, sendCount_XY,MPI_INT,rank_XY(),sendtag,Comm,&req1[7]);
	MPI_Irecv(recvList_xy, recvCount_xy,MPI_INT,rank_xy(),recvtag,Comm,&req2[7]);
	MPI_Isend(sendList_Xy, sendCount_Xy,MPI_INT,rank_Xy(),sendtag,Comm,&req1[8]);
	MPI_Irecv(recvList_xY, recvCount_xY,MPI_INT,rank_xY(),recvtag,Comm,&req2[8]);
	MPI_Isend(sendList_xY, sendCount_xY,MPI_INT,rank_xY(),sendtag,Comm,&req1[9]);
	MPI_Irecv(recvList_Xy, recvCount_Xy,MPI_INT,rank_Xy(),recvtag,Comm,&req2[9]);
	MPI_Isend(sendList_xz, sendCount_xz,MPI_INT,rank_xz(),sendtag,Comm,&req1[10]);
	MPI_Irecv(recvList_XZ, recvCount_XZ,MPI_INT,rank_XZ(),recvtag,Comm,&req2[10]);
	MPI_Isend(sendList_XZ, sendCount_XZ,MPI_INT,rank_XZ(),sendtag,Comm,&req1[11]);
	MPI_Irecv(recvList_xz, recvCount_xz,MPI_INT,rank_xz(),recvtag,Comm,&req2[11]);
	MPI_Isend(sendList_Xz, sendCount_Xz,MPI_INT,rank_Xz(),sendtag,Comm,&req1[12]);
	MPI_Irecv(recvList_xZ, recvCount_xZ,MPI_INT,rank_xZ(),recvtag,Comm,&req2[12]);
	MPI_Isend(sendList_xZ, sendCount_xZ,MPI_INT,rank_xZ(),sendtag,Comm,&req1[13]);
	MPI_Irecv(recvList_Xz, recvCount_Xz,MPI_INT,rank_Xz(),recvtag,Comm,&req2[13]);
	MPI_Isend(sendList_yz, sendCount_yz,MPI_INT,rank_yz(),sendtag,Comm,&req1[14]);
	MPI_Irecv(recvList_YZ, recvCount_YZ,MPI_INT,rank_YZ(),recvtag,Comm,&req2[14]);
	MPI_Isend(sendList_YZ, sendCount_YZ,MPI_INT,rank_YZ(),sendtag,Comm,&req1[15]);
	MPI_Irecv(recvList_yz, recvCount_yz,MPI_INT,rank_yz(),recvtag,Comm,&req2[15]);
	MPI_Isend(sendList_Yz, sendCount_Yz,MPI_INT,rank_Yz(),sendtag,Comm,&req1[16]);
	MPI_Irecv(recvList_yZ, recvCount_yZ,MPI_INT,rank_yZ(),recvtag,Comm,&req2[16]);
	MPI_Isend(sendList_yZ, sendCount_yZ,MPI_INT,rank_yZ(),sendtag,Comm,&req1[17]);
	MPI_Irecv(recvList_Yz, recvCount_Yz,MPI_INT,rank_Yz(),recvtag,Comm,&req2[17]);
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	//......................................................................................
	for (int idx=0; idx<recvCount_x; idx++)    recvList_x[idx] -= (Nx-2);
	for (int idx=0; idx<recvCount_X; idx++)    recvList_X[idx] += (Nx-2);
	for (int idx=0; idx<recvCount_y; idx++)    recvList_y[idx] -= (Ny-2)*Nx;
	for (int idx=0; idx<recvCount_Y; idx++)    recvList_Y[idx] += (Ny-2)*Nx;
	for (int idx=0; idx<recvCount_z; idx++)    recvList_z[idx] -= (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Z; idx++)    recvList_Z[idx] += (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_xy; idx++)    recvList_xy[idx] -= (Nx-2)+(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_XY; idx++)    recvList_XY[idx] += (Nx-2)+(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_xY; idx++)    recvList_xY[idx] -= (Nx-2)-(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_Xy; idx++)    recvList_Xy[idx] += (Nx-2)-(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_xz; idx++)    recvList_xz[idx] -= (Nx-2)+(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_XZ; idx++)    recvList_XZ[idx] += (Nx-2)+(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_xZ; idx++)    recvList_xZ[idx] -= (Nx-2)-(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Xz; idx++)    recvList_Xz[idx] += (Nx-2)-(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_yz; idx++)    recvList_yz[idx] -= (Ny-2)*Nx + (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_YZ; idx++)    recvList_YZ[idx] += (Ny-2)*Nx + (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_yZ; idx++)    recvList_yZ[idx] -= (Ny-2)*Nx - (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Yz; idx++)    recvList_Yz[idx] += (Ny-2)*Nx - (Nz-2)*Nx*Ny;
	//......................................................................................
	// allocate recv buffers
	recvBuf_x = new int [recvCount_x];
	recvBuf_y = new int [recvCount_y];
	recvBuf_z = new int [recvCount_z];
	recvBuf_X = new int [recvCount_X];
	recvBuf_Y = new int [recvCount_Y];
	recvBuf_Z = new int [recvCount_Z];
	recvBuf_xy = new int [recvCount_xy];
	recvBuf_yz = new int [recvCount_yz];
	recvBuf_xz = new int [recvCount_xz];
	recvBuf_Xy = new int [recvCount_Xy];
	recvBuf_Yz = new int [recvCount_Yz];
	recvBuf_xZ = new int [recvCount_xZ];
	recvBuf_xY = new int [recvCount_xY];
	recvBuf_yZ = new int [recvCount_yZ];
	recvBuf_Xz = new int [recvCount_Xz];
	recvBuf_XY = new int [recvCount_XY];
	recvBuf_YZ = new int [recvCount_YZ];
	recvBuf_XZ = new int [recvCount_XZ];
	//......................................................................................
	// send buffers
	sendData_x = new double [sendCount_x];
	sendData_y = new double [sendCount_y];
	sendData_z = new double [sendCount_z];
	sendData_X = new double [sendCount_X];
	sendData_Y = new double [sendCount_Y];
	sendData_Z = new double [sendCount_Z];
	sendData_xy = new double [sendCount_xy];
	sendData_yz = new double [sendCount_yz];
	sendData_xz = new double [sendCount_xz];
	sendData_Xy = new double [sendCount_Xy];
	sendData_Yz = new double [sendCount_Yz];
	sendData_xZ = new double [sendCount_xZ];
	sendData_xY = new double [sendCount_xY];
	sendData_yZ = new double [sendCount_yZ];
	sendData_Xz = new double [sendCount_Xz];
	sendData_XY = new double [sendCount_XY];
	sendData_YZ = new double [sendCount_YZ];
	sendData_XZ = new double [sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvData_x = new double [recvCount_x];
	recvData_y = new double [recvCount_y];
	recvData_z = new double [recvCount_z];
	recvData_X = new double [recvCount_X];
	recvData_Y = new double [recvCount_Y];
	recvData_Z = new double [recvCount_Z];
	recvData_xy = new double [recvCount_xy];
	recvData_yz = new double [recvCount_yz];
	recvData_xz = new double [recvCount_xz];
	recvData_Xy = new double [recvCount_Xy];
	recvData_xZ = new double [recvCount_xZ];
	recvData_xY = new double [recvCount_xY];
	recvData_yZ = new double [recvCount_yZ];
	recvData_Yz = new double [recvCount_Yz];
	recvData_Xz = new double [recvCount_Xz];
	recvData_XY = new double [recvCount_XY];
	recvData_YZ = new double [recvCount_YZ];
	recvData_XZ = new double [recvCount_XZ];
	//......................................................................................
	// Guarantee consistency of ID labels in halo
	CommHaloIDs();
}

void Domain::ReadIDs(){
	// Read the IDs from input file
  int nprocs=nprocx()*nprocy()*nprocz();
    size_t readID;
    char LocalRankString[8];
    char LocalRankFilename[40];
    //.......................................................................
    if (rank() == 0)    printf("Read input media... \n");
    //.......................................................................
    sprintf(LocalRankString,"%05d",rank());
    sprintf(LocalRankFilename,"%s%s","ID.",LocalRankString);
    // .......... READ THE INPUT FILE .......................................
    if (rank()==0) printf("Initialize from segmented data: solid=0, NWP=1, WP=2 \n");
    sprintf(LocalRankFilename,"ID.%05i",rank());
    FILE *IDFILE = fopen(LocalRankFilename,"rb");
    if (IDFILE==NULL) ERROR("Domain::ReadIDs --  Error opening file: ID.xxxxx");
    readID=fread(id,1,N,IDFILE);
    if (readID != size_t(N)) printf("Domain::ReadIDs -- Error reading ID (rank=%i) \n",rank());
    fclose(IDFILE);
    // Compute the porosity
    double sum;
    double porosity;
    double sum_local=0.0;
    double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
    if (BoundaryCondition > 0) iVol_global = 1.0/(1.0*(Nx-2)*nprocx()*(Ny-2)*nprocy()*((Nz-2)*nprocz()-6));
 	//.........................................................
 	// If external boundary conditions are applied remove solid
    if (BoundaryCondition >  0  && kproc() == 0){
 		for (int k=0; k<3; k++){
 			for (int j=0;j<Ny;j++){
 				for (int i=0;i<Nx;i++){
 					int n = k*Nx*Ny+j*Nx+i;
 					id[n] = 1;
 				}                    
 			}
 		}
 	}
    if (BoundaryCondition >  0  && kproc() == nprocz()-1){
 		for (int k=Nz-3; k<Nz; k++){
 			for (int j=0;j<Ny;j++){
 				for (int i=0;i<Nx;i++){
 					int n = k*Nx*Ny+j*Nx+i;
 					id[n] = 2;
 				}                    
 			}
 		}
 	}
    for (int k=1;k<Nz-1;k++){
        for (int j=1;j<Ny-1;j++){
            for (int i=1;i<Nx-1;i++){
                int n = k*Nx*Ny+j*Nx+i;
                if (id[n] > 0){
                    sum_local+=1.0;
                }
            }
        }
    }
    MPI_Allreduce(&sum_local,&sum,1,MPI_DOUBLE,MPI_SUM,Comm);
    porosity = sum*iVol_global;
    if (rank()==0) printf("Media porosity = %f \n",porosity);
 	//.........................................................
}
int Domain::PoreCount(){
	/*
	 * count the number of nodes occupied by mobile phases
	 */
    int Npore=0;  // number of local pore nodes
    for (int k=1;k<Nz-1;k++){
        for (int j=1;j<Ny-1;j++){
            for (int i=1;i<Nx-1;i++){
                int n = k*Nx*Ny+j*Nx+i;
                if (id[n] > 0){
                    Npore++;
                }
            }
        }
    }
    return Npore;
}

void Domain::CommunicateMeshHalo(DoubleArray &Mesh)
{
	int sendtag, recvtag;
	sendtag = recvtag = 7;
	double *MeshData = Mesh.data();
	PackMeshData(sendList_x, sendCount_x ,sendData_x, MeshData);
	PackMeshData(sendList_X, sendCount_X ,sendData_X, MeshData);
	PackMeshData(sendList_y, sendCount_y ,sendData_y, MeshData);
	PackMeshData(sendList_Y, sendCount_Y ,sendData_Y, MeshData);
	PackMeshData(sendList_z, sendCount_z ,sendData_z, MeshData);
	PackMeshData(sendList_Z, sendCount_Z ,sendData_Z, MeshData);
	PackMeshData(sendList_xy, sendCount_xy ,sendData_xy, MeshData);
	PackMeshData(sendList_Xy, sendCount_Xy ,sendData_Xy, MeshData);
	PackMeshData(sendList_xY, sendCount_xY ,sendData_xY, MeshData);
	PackMeshData(sendList_XY, sendCount_XY ,sendData_XY, MeshData);
	PackMeshData(sendList_xz, sendCount_xz ,sendData_xz, MeshData);
	PackMeshData(sendList_Xz, sendCount_Xz ,sendData_Xz, MeshData);
	PackMeshData(sendList_xZ, sendCount_xZ ,sendData_xZ, MeshData);
	PackMeshData(sendList_XZ, sendCount_XZ ,sendData_XZ, MeshData);
	PackMeshData(sendList_yz, sendCount_yz ,sendData_yz, MeshData);
	PackMeshData(sendList_Yz, sendCount_Yz ,sendData_Yz, MeshData);
	PackMeshData(sendList_yZ, sendCount_yZ ,sendData_yZ, MeshData);
	PackMeshData(sendList_YZ, sendCount_YZ ,sendData_YZ, MeshData);
	//......................................................................................
	MPI_Sendrecv(sendData_x,sendCount_x,MPI_DOUBLE,rank_x(),sendtag,
			recvData_X,recvCount_X,MPI_DOUBLE,rank_X(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_X,sendCount_X,MPI_DOUBLE,rank_X(),sendtag,
			recvData_x,recvCount_x,MPI_DOUBLE,rank_x(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_y,sendCount_y,MPI_DOUBLE,rank_y(),sendtag,
			recvData_Y,recvCount_Y,MPI_DOUBLE,rank_Y(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Y,sendCount_Y,MPI_DOUBLE,rank_Y(),sendtag,
			recvData_y,recvCount_y,MPI_DOUBLE,rank_y(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_z,sendCount_z,MPI_DOUBLE,rank_z(),sendtag,
			recvData_Z,recvCount_Z,MPI_DOUBLE,rank_Z(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Z,sendCount_Z,MPI_DOUBLE,rank_Z(),sendtag,
			recvData_z,recvCount_z,MPI_DOUBLE,rank_z(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xy,sendCount_xy,MPI_DOUBLE,rank_xy(),sendtag,
			recvData_XY,recvCount_XY,MPI_DOUBLE,rank_XY(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_XY,sendCount_XY,MPI_DOUBLE,rank_XY(),sendtag,
			recvData_xy,recvCount_xy,MPI_DOUBLE,rank_xy(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Xy,sendCount_Xy,MPI_DOUBLE,rank_Xy(),sendtag,
			recvData_xY,recvCount_xY,MPI_DOUBLE,rank_xY(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xY,sendCount_xY,MPI_DOUBLE,rank_xY(),sendtag,
			recvData_Xy,recvCount_Xy,MPI_DOUBLE,rank_Xy(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xz,sendCount_xz,MPI_DOUBLE,rank_xz(),sendtag,
			recvData_XZ,recvCount_XZ,MPI_DOUBLE,rank_XZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_XZ,sendCount_XZ,MPI_DOUBLE,rank_XZ(),sendtag,
			recvData_xz,recvCount_xz,MPI_DOUBLE,rank_xz(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Xz,sendCount_Xz,MPI_DOUBLE,rank_Xz(),sendtag,
			recvData_xZ,recvCount_xZ,MPI_DOUBLE,rank_xZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xZ,sendCount_xZ,MPI_DOUBLE,rank_xZ(),sendtag,
			recvData_Xz,recvCount_Xz,MPI_DOUBLE,rank_Xz(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_yz,sendCount_yz,MPI_DOUBLE,rank_yz(),sendtag,
			recvData_YZ,recvCount_YZ,MPI_DOUBLE,rank_YZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_YZ,sendCount_YZ,MPI_DOUBLE,rank_YZ(),sendtag,
			recvData_yz,recvCount_yz,MPI_DOUBLE,rank_yz(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Yz,sendCount_Yz,MPI_DOUBLE,rank_Yz(),sendtag,
			recvData_yZ,recvCount_yZ,MPI_DOUBLE,rank_yZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_yZ,sendCount_yZ,MPI_DOUBLE,rank_yZ(),sendtag,
			recvData_Yz,recvCount_Yz,MPI_DOUBLE,rank_Yz(),recvtag,Comm,MPI_STATUS_IGNORE);
	//........................................................................................
	UnpackMeshData(recvList_x, recvCount_x ,recvData_x, MeshData);
	UnpackMeshData(recvList_X, recvCount_X ,recvData_X, MeshData);
	UnpackMeshData(recvList_y, recvCount_y ,recvData_y, MeshData);
	UnpackMeshData(recvList_Y, recvCount_Y ,recvData_Y, MeshData);
	UnpackMeshData(recvList_z, recvCount_z ,recvData_z, MeshData);
	UnpackMeshData(recvList_Z, recvCount_Z ,recvData_Z, MeshData);
	UnpackMeshData(recvList_xy, recvCount_xy ,recvData_xy, MeshData);
	UnpackMeshData(recvList_Xy, recvCount_Xy ,recvData_Xy, MeshData);
	UnpackMeshData(recvList_xY, recvCount_xY ,recvData_xY, MeshData);
	UnpackMeshData(recvList_XY, recvCount_XY ,recvData_XY, MeshData);
	UnpackMeshData(recvList_xz, recvCount_xz ,recvData_xz, MeshData);
	UnpackMeshData(recvList_Xz, recvCount_Xz ,recvData_Xz, MeshData);
	UnpackMeshData(recvList_xZ, recvCount_xZ ,recvData_xZ, MeshData);
	UnpackMeshData(recvList_XZ, recvCount_XZ ,recvData_XZ, MeshData);
	UnpackMeshData(recvList_yz, recvCount_yz ,recvData_yz, MeshData);
	UnpackMeshData(recvList_Yz, recvCount_Yz ,recvData_Yz, MeshData);
	UnpackMeshData(recvList_yZ, recvCount_yZ ,recvData_yZ, MeshData);
	UnpackMeshData(recvList_YZ, recvCount_YZ ,recvData_YZ, MeshData);
}

void Domain::PackID(int *list, int count, char *sendbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

void Domain::UnpackID(int *list, int count, char *recvbuf, char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}

void Domain::CommHaloIDs(){
	// Ensure consistency for halo IDs
	char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
    //......................................................................................
    int sendCountID_x, sendCountID_y, sendCountID_z, sendCountID_X, sendCountID_Y, sendCountID_Z;
    int sendCountID_xy, sendCountID_yz, sendCountID_xz, sendCountID_Xy, sendCountID_Yz, sendCountID_xZ;
    int sendCountID_xY, sendCountID_yZ, sendCountID_Xz, sendCountID_XY, sendCountID_YZ, sendCountID_XZ;
    //......................................................................................
    int recvCountID_x, recvCountID_y, recvCountID_z, recvCountID_X, recvCountID_Y, recvCountID_Z;
    int recvCountID_xy, recvCountID_yz, recvCountID_xz, recvCountID_Xy, recvCountID_Yz, recvCountID_xZ;
    int recvCountID_xY, recvCountID_yZ, recvCountID_Xz, recvCountID_XY, recvCountID_YZ, recvCountID_XZ;
    //......................................................................................
	//......................................................................................
	sendCountID_x = sendCountID_y = sendCountID_z = sendCountID_X = sendCountID_Y = sendCountID_Z = 0;
	sendCountID_xy = sendCountID_yz = sendCountID_xz = sendCountID_Xy = sendCountID_Yz = sendCountID_xZ = 0;
	sendCountID_xY = sendCountID_yZ = sendCountID_Xz = sendCountID_XY = sendCountID_YZ = sendCountID_XZ = 0;
	//......................................................................................
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				// Check the phase ID
				// CountIDs for the six faces
				if (i==1)    sendCountID_x++;
				if (j==1)    sendCountID_y++;
				if (k==1)    sendCountID_z++;
				if (i==Nx-2)    sendCountID_X++;
				if (j==Ny-2)    sendCountID_Y++;
				if (k==Nz-2)    sendCountID_Z++;
				// CountIDs for the twelve edges
				if (i==1 && j==1)    sendCountID_xy++;
				if (i==1 && j==Ny-2)    sendCountID_xY++;
				if (i==Nx-2 && j==1)    sendCountID_Xy++;
				if (i==Nx-2 && j==Ny-2)    sendCountID_XY++;

				if (i==1 && k==1)    sendCountID_xz++;
				if (i==1 && k==Nz-2)    sendCountID_xZ++;
				if (i==Nx-2 && k==1)    sendCountID_Xz++;
				if (i==Nx-2 && k==Nz-2)    sendCountID_XZ++;

				if (j==1 && k==1)    sendCountID_yz++;
				if (j==1 && k==Nz-2)    sendCountID_yZ++;
				if (j==Ny-2 && k==1)    sendCountID_Yz++;
				if (j==Ny-2 && k==Nz-2)    sendCountID_YZ++;
			}
		}
	}
	// send buffers
	sendID_x = new char [sendCountID_x];
	sendID_y = new char [sendCountID_y];
	sendID_z = new char [sendCountID_z];
	sendID_X = new char [sendCountID_X];
	sendID_Y = new char [sendCountID_Y];
	sendID_Z = new char [sendCountID_Z];
	sendID_xy = new char [sendCountID_xy];
	sendID_yz = new char [sendCountID_yz];
	sendID_xz = new char [sendCountID_xz];
	sendID_Xy = new char [sendCountID_Xy];
	sendID_Yz = new char [sendCountID_Yz];
	sendID_xZ = new char [sendCountID_xZ];
	sendID_xY = new char [sendCountID_xY];
	sendID_yZ = new char [sendCountID_yZ];
	sendID_Xz = new char [sendCountID_Xz];
	sendID_XY = new char [sendCountID_XY];
	sendID_YZ = new char [sendCountID_YZ];
	sendID_XZ = new char [sendCountID_XZ];
	//......................................................................................
	// recv buffers
	recvID_x = new char [recvCountID_x];
	recvID_y = new char [recvCountID_y];
	recvID_z = new char [recvCountID_z];
	recvID_X = new char [recvCountID_X];
	recvID_Y = new char [recvCountID_Y];
	recvID_Z = new char [recvCountID_Z];
	recvID_xy = new char [recvCountID_xy];
	recvID_yz = new char [recvCountID_yz];
	recvID_xz = new char [recvCountID_xz];
	recvID_Xy = new char [recvCountID_Xy];
	recvID_xZ = new char [recvCountID_xZ];
	recvID_xY = new char [recvCountID_xY];
	recvID_yZ = new char [recvCountID_yZ];
	recvID_Yz = new char [recvCountID_Yz];
	recvID_Xz = new char [recvCountID_Xz];
	recvID_XY = new char [recvCountID_XY];
	recvID_YZ = new char [recvCountID_YZ];
	recvID_XZ = new char [recvCountID_XZ];
	//......................................................................................
	int sendtag,recvtag;
	sendtag = recvtag = 7;

	// Pack and send the updated ID values
	PackID(sendList_x, sendCountID_x ,sendID_x, id);
	PackID(sendList_X, sendCountID_X ,sendID_X, id);
	PackID(sendList_y, sendCountID_y ,sendID_y, id);
	PackID(sendList_Y, sendCountID_Y ,sendID_Y, id);
	PackID(sendList_z, sendCountID_z ,sendID_z, id);
	PackID(sendList_Z, sendCountID_Z ,sendID_Z, id);
	PackID(sendList_xy, sendCountID_xy ,sendID_xy, id);
	PackID(sendList_Xy, sendCountID_Xy ,sendID_Xy, id);
	PackID(sendList_xY, sendCountID_xY ,sendID_xY, id);
	PackID(sendList_XY, sendCountID_XY ,sendID_XY, id);
	PackID(sendList_xz, sendCountID_xz ,sendID_xz, id);
	PackID(sendList_Xz, sendCountID_Xz ,sendID_Xz, id);
	PackID(sendList_xZ, sendCountID_xZ ,sendID_xZ, id);
	PackID(sendList_XZ, sendCountID_XZ ,sendID_XZ, id);
	PackID(sendList_yz, sendCountID_yz ,sendID_yz, id);
	PackID(sendList_Yz, sendCountID_Yz ,sendID_Yz, id);
	PackID(sendList_yZ, sendCountID_yZ ,sendID_yZ, id);
	PackID(sendList_YZ, sendCountID_YZ ,sendID_YZ, id);
	//......................................................................................
	MPI_Sendrecv(sendID_x,sendCountID_x,MPI_CHAR,rank_x(),sendtag,
			recvID_X,recvCountID_X,MPI_CHAR,rank_X(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_X,sendCountID_X,MPI_CHAR,rank_X(),sendtag,
			recvID_x,recvCountID_x,MPI_CHAR,rank_x(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_y,sendCountID_y,MPI_CHAR,rank_y(),sendtag,
			recvID_Y,recvCountID_Y,MPI_CHAR,rank_Y(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Y,sendCountID_Y,MPI_CHAR,rank_Y(),sendtag,
			recvID_y,recvCountID_y,MPI_CHAR,rank_y(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_z,sendCountID_z,MPI_CHAR,rank_z(),sendtag,
			recvID_Z,recvCountID_Z,MPI_CHAR,rank_Z(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Z,sendCountID_Z,MPI_CHAR,rank_Z(),sendtag,
			recvID_z,recvCountID_z,MPI_CHAR,rank_z(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xy,sendCountID_xy,MPI_CHAR,rank_xy(),sendtag,
			recvID_XY,recvCountID_XY,MPI_CHAR,rank_XY(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XY,sendCountID_XY,MPI_CHAR,rank_XY(),sendtag,
			recvID_xy,recvCountID_xy,MPI_CHAR,rank_xy(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xy,sendCountID_Xy,MPI_CHAR,rank_Xy(),sendtag,
			recvID_xY,recvCountID_xY,MPI_CHAR,rank_xY(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xY,sendCountID_xY,MPI_CHAR,rank_xY(),sendtag,
			recvID_Xy,recvCountID_Xy,MPI_CHAR,rank_Xy(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xz,sendCountID_xz,MPI_CHAR,rank_xz(),sendtag,
			recvID_XZ,recvCountID_XZ,MPI_CHAR,rank_XZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_XZ,sendCountID_XZ,MPI_CHAR,rank_XZ(),sendtag,
			recvID_xz,recvCountID_xz,MPI_CHAR,rank_xz(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Xz,sendCountID_Xz,MPI_CHAR,rank_Xz(),sendtag,
			recvID_xZ,recvCountID_xZ,MPI_CHAR,rank_xZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_xZ,sendCountID_xZ,MPI_CHAR,rank_xZ(),sendtag,
			recvID_Xz,recvCountID_Xz,MPI_CHAR,rank_Xz(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yz,sendCountID_yz,MPI_CHAR,rank_yz(),sendtag,
			recvID_YZ,recvCountID_YZ,MPI_CHAR,rank_YZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_YZ,sendCountID_YZ,MPI_CHAR,rank_YZ(),sendtag,
			recvID_yz,recvCountID_yz,MPI_CHAR,rank_yz(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_Yz,sendCountID_Yz,MPI_CHAR,rank_Yz(),sendtag,
			recvID_yZ,recvCountID_yZ,MPI_CHAR,rank_yZ(),recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendID_yZ,sendCountID_yZ,MPI_CHAR,rank_yZ(),sendtag,
			recvID_Yz,recvCountID_Yz,MPI_CHAR,rank_Yz(),recvtag,Comm,MPI_STATUS_IGNORE);
	//......................................................................................
	UnpackID(recvList_x, recvCountID_x ,recvID_x, id);
	UnpackID(recvList_X, recvCountID_X ,recvID_X, id);
	UnpackID(recvList_y, recvCountID_y ,recvID_y, id);
	UnpackID(recvList_Y, recvCountID_Y ,recvID_Y, id);
	UnpackID(recvList_z, recvCountID_z ,recvID_z, id);
	UnpackID(recvList_Z, recvCountID_Z ,recvID_Z, id);
	UnpackID(recvList_xy, recvCountID_xy ,recvID_xy, id);
	UnpackID(recvList_Xy, recvCountID_Xy ,recvID_Xy, id);
	UnpackID(recvList_xY, recvCountID_xY ,recvID_xY, id);
	UnpackID(recvList_XY, recvCountID_XY ,recvID_XY, id);
	UnpackID(recvList_xz, recvCountID_xz ,recvID_xz, id);
	UnpackID(recvList_Xz, recvCountID_Xz ,recvID_Xz, id);
	UnpackID(recvList_xZ, recvCountID_xZ ,recvID_xZ, id);
	UnpackID(recvList_XZ, recvCountID_XZ ,recvID_XZ, id);
	UnpackID(recvList_yz, recvCountID_yz ,recvID_yz, id);
	UnpackID(recvList_Yz, recvCountID_Yz ,recvID_Yz, id);
	UnpackID(recvList_yZ, recvCountID_yZ ,recvID_yZ, id);
	UnpackID(recvList_YZ, recvCountID_YZ ,recvID_YZ, id);
	//......................................................................................
}

// Ideally stuff below here should be moved somewhere else -- doesn't really belong here
void WriteCheckpoint(const char *FILENAME, const double *cDen, const double *cfq, int Np)
{
    int q,n;
    double value;
    ofstream File(FILENAME,ios::binary);
    for (n=0; n<Np; n++){
        // Write the two density values
        value = cDen[n];
        File.write((char*) &value, sizeof(value));
        value = cDen[Np+n];
        File.write((char*) &value, sizeof(value));
        // Write the even distributions
        for (q=0; q<19; q++){
            value = cfq[q*Np+n];
            File.write((char*) &value, sizeof(value));
        }
    }
    File.close();

}

void ReadCheckpoint(char *FILENAME, double *cPhi, double *cfq, int Np)
{
    int q=0, n=0;
    double value=0;
    ifstream File(FILENAME,ios::binary);
    for (n=0; n<Np; n++){
        File.read((char*) &value, sizeof(value));
        cPhi[n] = value;
        // Read the distributions
        for (q=0; q<19; q++){
            File.read((char*) &value, sizeof(value));
            cfq[q*Np+n] = value;
        }
    }
    File.close();
}

void ReadBinaryFile(char *FILENAME, double *Data, int N)
{
  int n;
  double value;
  ifstream File(FILENAME,ios::binary);
  if (File.good()){
    for (n=0; n<N; n++){
      // Write the two density values                                                                                
      File.read((char*) &value, sizeof(value));
      Data[n] = value;

    }
  }
  else {
    for (n=0; n<N; n++) Data[n] = 1.2e-34;
  }
  File.close();
}

