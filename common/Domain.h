#ifndef Domain_INC
#define Domain_INC
// Created by James McClure
// Copyright 2008-2013
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <Array.h>
#include <math.h>
#include <time.h>
#include <exception>      // std::exception
#include <stdexcept>

#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"

int MAX_BLOB_COUNT=50;

using namespace std;

struct Domain{

	Domain(int nx, int ny, int nz, int rnk, int npx, int npy, int npz, 
			double lx, double ly, double lz, int BC){
	        Volume = nx*ny*nx*npx*npy*npz*1.0;
		Nx = nx+2; Ny = ny+2; Nz = nz+2; 
		Lx = lx, Ly = ly, Lz = lz;
		rank = rnk;
		nprocx=npx; nprocy=npy; nprocz=npz;
		N = Nx*Ny*Nz;
		id = new char [N];
		BlobLabel.resize(Nx,Ny,Nz);
		BlobGraph.resize(18,MAX_BLOB_COUNT,MAX_BLOB_COUNT);
		BoundaryCondition = BC;
		rank_info=RankInfoStruct(rank,nprocx,nprocy,nprocz);
	}
	~Domain();

	// Basic domain information
	int Nx,Ny,Nz,N;
	int iproc,jproc,kproc;
 	int nprocx,nprocy,nprocz;
    double Lx,Ly,Lz,Volume;
	int rank;
	int BoundaryCondition;
	RankInfoStruct rank_info;
	MPI_Group Group;	// Group of processors associated with this domain
	MPI_Comm Comm;		// MPI Communicator for this domain

	//**********************************
	// MPI ranks for all 18 neighbors
	//**********************************
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//**********************************
	//......................................................................................
	// Get the actual D3Q19 communication counts (based on location of solid phase)
	// Discrete velocity set symmetry implies the sendcount = recvcount
	//......................................................................................
	int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
	int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
	int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
	//......................................................................................
	int *sendList_x, *sendList_y, *sendList_z, *sendList_X, *sendList_Y, *sendList_Z;
	int *sendList_xy, *sendList_yz, *sendList_xz, *sendList_Xy, *sendList_Yz, *sendList_xZ;
	int *sendList_xY, *sendList_yZ, *sendList_Xz, *sendList_XY, *sendList_YZ, *sendList_XZ;
	//......................................................................................
	int *sendBuf_x, *sendBuf_y, *sendBuf_z, *sendBuf_X, *sendBuf_Y, *sendBuf_Z;
	int *sendBuf_xy, *sendBuf_yz, *sendBuf_xz, *sendBuf_Xy, *sendBuf_Yz, *sendBuf_xZ;
	int *sendBuf_xY, *sendBuf_yZ, *sendBuf_Xz, *sendBuf_XY, *sendBuf_YZ, *sendBuf_XZ;
	//......................................................................................
	int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
	int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
	int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
	//......................................................................................
	int *recvList_x, *recvList_y, *recvList_z, *recvList_X, *recvList_Y, *recvList_Z;
	int *recvList_xy, *recvList_yz, *recvList_xz, *recvList_Xy, *recvList_Yz, *recvList_xZ;
	int *recvList_xY, *recvList_yZ, *recvList_Xz, *recvList_XY, *recvList_YZ, *recvList_XZ;
	//......................................................................................
	int *recvBuf_x, *recvBuf_y, *recvBuf_z, *recvBuf_X, *recvBuf_Y, *recvBuf_Z;
	int *recvBuf_xy, *recvBuf_yz, *recvBuf_xz, *recvBuf_Xy, *recvBuf_Yz, *recvBuf_xZ;
	int *recvBuf_xY, *recvBuf_yZ, *recvBuf_Xz, *recvBuf_XY, *recvBuf_YZ, *recvBuf_XZ;
	//......................................................................................	
	double *sendData_x, *sendData_y, *sendData_z, *sendData_X, *sendData_Y, *sendData_Z;
	double *sendData_xy, *sendData_yz, *sendData_xz, *sendData_Xy, *sendData_Yz, *sendData_xZ;
	double *sendData_xY, *sendData_yZ, *sendData_Xz, *sendData_XY, *sendData_YZ, *sendData_XZ;
	double *recvData_x, *recvData_y, *recvData_z, *recvData_X, *recvData_Y, *recvData_Z;
	double *recvData_xy, *recvData_yz, *recvData_xz, *recvData_Xy, *recvData_Yz, *recvData_xZ;
	double *recvData_xY, *recvData_yZ, *recvData_Xz, *recvData_XY, *recvData_YZ, *recvData_XZ;

	// Solid indicator function
	char *id;
	// Blob information
	IntArray BlobLabel;
	IntArray BlobGraph;

	void InitializeRanks();
	void CommInit(MPI_Comm comm);
	void CommunicateMeshHalo(DoubleArray &Mesh);
	void BlobComm(MPI_Comm comm);

	void AssignBlobConnections(){
		getBlobConnections(recvList_x, recvCount_x, rank_x, 0);
		getBlobConnections(recvList_y, recvCount_y, rank_y, 1);
		getBlobConnections(recvList_z, recvCount_z, rank_z, 2);
		getBlobConnections(recvList_X, recvCount_X, rank_X, 3);
		getBlobConnections(recvList_Y, recvCount_y, rank_Y, 4);
		getBlobConnections(recvList_Z, recvCount_Z, rank_Z, 5);
		getBlobConnections(recvList_xy, recvCount_xy, rank_xy, 6);
		getBlobConnections(recvList_xY, recvCount_xY, rank_xY, 7);
		getBlobConnections(recvList_Xy, recvCount_Xy, rank_Xy, 8);
		getBlobConnections(recvList_XY, recvCount_XY, rank_XY, 9);
		getBlobConnections(recvList_xz, recvCount_xz, rank_xz, 10);
		getBlobConnections(recvList_xZ, recvCount_xZ, rank_xZ, 11);
		getBlobConnections(recvList_Xz, recvCount_Xz, rank_Xz, 12);
		getBlobConnections(recvList_XZ, recvCount_XZ, rank_XZ, 13);
		getBlobConnections(recvList_yz, recvCount_yz, rank_yz, 14);
		getBlobConnections(recvList_yZ, recvCount_yZ, rank_yZ, 15);
		getBlobConnections(recvList_Yz, recvCount_Yz, rank_Yz, 16);
		getBlobConnections(recvList_YZ, recvCount_YZ, rank_YZ, 17);
	}

private:

	int getRankForBlock( int i, int j, int k )
	{
		int i2 = (i+nprocx)%nprocx;
		int j2 = (j+nprocy)%nprocy;
		int k2 = (k+nprocz)%nprocz;
		return i2 + j2*nprocx + k2*nprocx*nprocy;
	}
	void PackBlobData(int *list, int count, int *sendbuf, int *data){
		// Fill in the phase ID values from neighboring processors
		// This packs up the values that need to be sent from one processor to another
		int idx,n;
		for (idx=0; idx<count; idx++){
			n = list[idx];
			sendbuf[idx] = data[n];
		}
	}
	void UnpackBlobData(int *list, int count, int *recvbuf, int *data){
		// Fill in the phase ID values from neighboring processors
		// This unpacks the values once they have been recieved from neighbors
		int idx,n;
		for (idx=0; idx<count; idx++){
			n = list[idx];
			data[n] = recvbuf[idx];
		}
	}
	
	int VoxelConnection(int n){
		int d[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
		{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
		{1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
		{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
		{-1,-1,1},{-1,-1,-1}};   // directions to neighbors
		
		int returnVal = -1;
		int x,y,z;
		// Get the 3-D indices
		x = n%Nx;
		y = (n/Nx)%Ny;
		z = n/(Nx*Ny);
		int nodx,nody,nodz;
		for (int p=0;p<26;p++){
			nodx=x+d[p][0];
			// Get the neighbor and guarantee it is in the domain 
			if (nodx < 0 ){ nodx = 0; }		
			if (nodx > Nx-1 ){ nodx = Nx-1; }
			nody=y+d[p][1];
			if (nody < 0 ){ nody = 0; }		
			if (nody > Ny-1 ){ nody = Ny-1; }
			nodz=z+d[p][2];
			if (nodz < 0 ){ nodz = 0; }		
			if (nodz > Nz-1 ){ nodz = Nz-1; }
			
			if (BlobLabel(nodx,nody,nodz) > returnVal )	returnVal = BlobLabel(nodx,nody,nodz);
		}
		return returnVal;
	}
	
	void getBlobConnections(int * List, int count, int neighbor, int direction){

		int idx,n,localValue,neighborValue;
		int x,y,z;
		for (idx=0; idx<count; idx++){
			n = List[idx];
			// Get the 3-D indices
			x = n%Nx;
			y = (n/Nx)%Ny;
			z = n/(Nx*Ny);
			neighborValue = BlobLabel(x,y,z);
			if (neighborValue > -1){
				localValue = VoxelConnection(n);
				printf("Blob (%i,%i) connects to neighbor blob (%i,%i)", localValue, rank, neighborValue, neighbor);
				BlobGraph(direction,localValue,neighbor) = 1; // Set the BlobGraph to TRUE for this pair		
			}
		}
	}
};

// Inline function to read line without a return argument
static inline void fgetl( char * str, int num, FILE * stream )
{
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}

Domain::~Domain(){
	delete sendData_x;
	delete sendData_y;
	delete sendData_z;
	delete sendData_X;
	delete sendData_Y;
	delete sendData_Z;
	delete sendData_xy;
	delete sendData_xY;
	delete sendData_Xy;
	delete sendData_XY;
	delete sendData_xz;
	delete sendData_xZ;
	delete sendData_Xz;
	delete sendData_XZ;
	delete sendData_yz;
	delete sendData_yZ;
	delete sendData_Yz;
	delete sendData_YZ;
	delete recvData_x;
	delete recvData_y;
	delete recvData_z;
	delete recvData_X;
	delete recvData_Y;
	delete recvData_Z;
	delete recvData_xy;
	delete recvData_xY;
	delete recvData_Xy;
	delete recvData_XY;
	delete recvData_xz;
	delete recvData_xZ;
	delete recvData_Xz;
	delete recvData_XZ;
	delete recvData_yz;
	delete recvData_yZ;
	delete recvData_Yz;
	delete recvData_YZ;
}

void Domain::InitializeRanks()
{
	// map the rank to the block index
	iproc = rank%nprocx;
	jproc = (rank/nprocx)%nprocy;
	kproc = rank/(nprocx*nprocy);

	// set up the neighbor ranks
    int i = iproc;
    int j = jproc;
    int k = kproc;
	rank_X = getRankForBlock(i+1,j,k);
	rank_x = getRankForBlock(i-1,j,k);
	rank_Y = getRankForBlock(i,j+1,k);
	rank_y = getRankForBlock(i,j-1,k);
	rank_Z = getRankForBlock(i,j,k+1);
	rank_z = getRankForBlock(i,j,k-1);
	rank_XY = getRankForBlock(i+1,j+1,k);
	rank_xy = getRankForBlock(i-1,j-1,k);
	rank_Xy = getRankForBlock(i+1,j-1,k);
	rank_xY = getRankForBlock(i-1,j+1,k);
	rank_XZ = getRankForBlock(i+1,j,k+1);
	rank_xz = getRankForBlock(i-1,j,k-1);
	rank_Xz = getRankForBlock(i+1,j,k-1);
	rank_xZ = getRankForBlock(i-1,j,k+1);
	rank_YZ = getRankForBlock(i,j+1,k+1);
	rank_yz = getRankForBlock(i,j-1,k-1);
	rank_Yz = getRankForBlock(i,j+1,k-1);
	rank_yZ = getRankForBlock(i,j-1,k+1);
}


void Domain::CommInit(MPI_Comm Communicator){
	int i,j,k,n;
	int sendtag = 21;
	int recvtag = 21;

	//......................................................................................
	//Get the ranks of each process and it's neighbors
	// map the rank to the block index
	iproc = rank%nprocx;
	jproc = (rank/nprocx)%nprocy;
	kproc = rank/(nprocx*nprocy);
	// set up the neighbor ranks
    i = iproc;
    j = jproc;
    k = kproc;
	rank_X = getRankForBlock(i+1,j,k);
	rank_x = getRankForBlock(i-1,j,k);
	rank_Y = getRankForBlock(i,j+1,k);
	rank_y = getRankForBlock(i,j-1,k);
	rank_Z = getRankForBlock(i,j,k+1);
	rank_z = getRankForBlock(i,j,k-1);
	rank_XY = getRankForBlock(i+1,j+1,k);
	rank_xy = getRankForBlock(i-1,j-1,k);
	rank_Xy = getRankForBlock(i+1,j-1,k);
	rank_xY = getRankForBlock(i-1,j+1,k);
	rank_XZ = getRankForBlock(i+1,j,k+1);
	rank_xz = getRankForBlock(i-1,j,k-1);
	rank_Xz = getRankForBlock(i+1,j,k-1);
	rank_xZ = getRankForBlock(i-1,j,k+1);
	rank_YZ = getRankForBlock(i,j+1,k+1);
	rank_yz = getRankForBlock(i,j-1,k-1);
	rank_Yz = getRankForBlock(i,j+1,k-1);
	rank_yZ = getRankForBlock(i,j-1,k+1);
	//......................................................................................

	MPI_Comm_group(Communicator,&Group);
	MPI_Comm_create(Communicator,Group,&Comm);

	//......................................................................................
	MPI_Request req1[18], req2[18];
	MPI_Status stat1[18],stat2[18];
	//......................................................................................
	sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
	sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
	sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
	//......................................................................................
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				// Check the phase ID
				if (id[k*Nx*Ny+j*Nx+i] != 0){
					// Counts for the six faces
					if (i==1)	sendCount_x++;
					if (j==1)	sendCount_y++;
					if (k==1)	sendCount_z++;
					if (i==Nx-2)	sendCount_X++;
					if (j==Ny-2)	sendCount_Y++;
					if (k==Nz-2)	sendCount_Z++;
					// Counts for the twelve edges
					if (i==1 && j==1)	sendCount_xy++;
					if (i==1 && j==Ny-2)	sendCount_xY++;
					if (i==Nx-2 && j==1)	sendCount_Xy++;
					if (i==Nx-2 && j==Ny-2)	sendCount_XY++;

					if (i==1 && k==1)	sendCount_xz++;
					if (i==1 && k==Nz-2)	sendCount_xZ++;
					if (i==Nx-2 && k==1)	sendCount_Xz++;
					if (i==Nx-2 && k==Nz-2)	sendCount_XZ++;

					if (j==1 && k==1)	sendCount_yz++;
					if (j==1 && k==Nz-2)	sendCount_yZ++;
					if (j==Ny-2 && k==1)	sendCount_Yz++;
					if (j==Ny-2 && k==Nz-2)	sendCount_YZ++;
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
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				// Local value to send
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] != 0){
					// Counts for the six faces
					if (i==1)		sendList_x[sendCount_x++]=n;
					if (j==1)		sendList_y[sendCount_y++]=n;
					if (k==1)		sendList_z[sendCount_z++]=n;
					if (i==Nx-2)	sendList_X[sendCount_X++]=n;
					if (j==Ny-2)	sendList_Y[sendCount_Y++]=n;
					if (k==Nz-2)	sendList_Z[sendCount_Z++]=n;
					// Counts for the twelve edges
					if (i==1 && j==1)		sendList_xy[sendCount_xy++]=n;
					if (i==1 && j==Ny-2)	sendList_xY[sendCount_xY++]=n;
					if (i==Nx-2 && j==1)	sendList_Xy[sendCount_Xy++]=n;
					if (i==Nx-2 && j==Ny-2)	sendList_XY[sendCount_XY++]=n;

					if (i==1 && k==1)		sendList_xz[sendCount_xz++]=n;
					if (i==1 && k==Nz-2)	sendList_xZ[sendCount_xZ++]=n;
					if (i==Nx-2 && k==1)	sendList_Xz[sendCount_Xz++]=n;
					if (i==Nx-2 && k==Nz-2)	sendList_XZ[sendCount_XZ++]=n;

					if (j==1 && k==1)		sendList_yz[sendCount_yz++]=n;
					if (j==1 && k==Nz-2)	sendList_yZ[sendCount_yZ++]=n;
					if (j==Ny-2 && k==1)	sendList_Yz[sendCount_Yz++]=n;
					if (j==Ny-2 && k==Nz-2)	sendList_YZ[sendCount_YZ++]=n;
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
	MPI_Isend(&sendCount_x, 1,MPI_INT,rank_x,sendtag+0,Communicator,&req1[0]);
	MPI_Irecv(&recvCount_X, 1,MPI_INT,rank_X,recvtag+0,Communicator,&req2[0]);
	MPI_Isend(&sendCount_X, 1,MPI_INT,rank_X,sendtag+1,Communicator,&req1[1]);
	MPI_Irecv(&recvCount_x, 1,MPI_INT,rank_x,recvtag+1,Communicator,&req2[1]);
	MPI_Isend(&sendCount_y, 1,MPI_INT,rank_y,sendtag+2,Communicator,&req1[2]);
	MPI_Irecv(&recvCount_Y, 1,MPI_INT,rank_Y,recvtag+2,Communicator,&req2[2]);
	MPI_Isend(&sendCount_Y, 1,MPI_INT,rank_Y,sendtag+3,Communicator,&req1[3]);
	MPI_Irecv(&recvCount_y, 1,MPI_INT,rank_y,recvtag+3,Communicator,&req2[3]);
	MPI_Isend(&sendCount_z, 1,MPI_INT,rank_z,sendtag+4,Communicator,&req1[4]);
	MPI_Irecv(&recvCount_Z, 1,MPI_INT,rank_Z,recvtag+4,Communicator,&req2[4]);
	MPI_Isend(&sendCount_Z, 1,MPI_INT,rank_Z,sendtag+5,Communicator,&req1[5]);
	MPI_Irecv(&recvCount_z, 1,MPI_INT,rank_z,recvtag+5,Communicator,&req2[5]);
	MPI_Isend(&sendCount_xy, 1,MPI_INT,rank_xy,sendtag+6,Communicator,&req1[6]);
	MPI_Irecv(&recvCount_XY, 1,MPI_INT,rank_XY,recvtag+6,Communicator,&req2[6]);
	MPI_Isend(&sendCount_XY, 1,MPI_INT,rank_XY,sendtag+7,Communicator,&req1[7]);
	MPI_Irecv(&recvCount_xy, 1,MPI_INT,rank_xy,recvtag+7,Communicator,&req2[7]);
	MPI_Isend(&sendCount_Xy, 1,MPI_INT,rank_Xy,sendtag+8,Communicator,&req1[8]);
	MPI_Irecv(&recvCount_xY, 1,MPI_INT,rank_xY,recvtag+8,Communicator,&req2[8]);
	MPI_Isend(&sendCount_xY, 1,MPI_INT,rank_xY,sendtag+9,Communicator,&req1[9]);
	MPI_Irecv(&recvCount_Xy, 1,MPI_INT,rank_Xy,recvtag+9,Communicator,&req2[9]);
	MPI_Isend(&sendCount_xz, 1,MPI_INT,rank_xz,sendtag+10,Communicator,&req1[10]);
	MPI_Irecv(&recvCount_XZ, 1,MPI_INT,rank_XZ,recvtag+10,Communicator,&req2[10]);
	MPI_Isend(&sendCount_XZ, 1,MPI_INT,rank_XZ,sendtag+11,Communicator,&req1[11]);
	MPI_Irecv(&recvCount_xz, 1,MPI_INT,rank_xz,recvtag+11,Communicator,&req2[11]);
	MPI_Isend(&sendCount_Xz, 1,MPI_INT,rank_Xz,sendtag+12,Communicator,&req1[12]);
	MPI_Irecv(&recvCount_xZ, 1,MPI_INT,rank_xZ,recvtag+12,Communicator,&req2[12]);
	MPI_Isend(&sendCount_xZ, 1,MPI_INT,rank_xZ,sendtag+13,Communicator,&req1[13]);
	MPI_Irecv(&recvCount_Xz, 1,MPI_INT,rank_Xz,recvtag+13,Communicator,&req2[13]);
	MPI_Isend(&sendCount_yz, 1,MPI_INT,rank_yz,sendtag+14,Communicator,&req1[14]);
	MPI_Irecv(&recvCount_YZ, 1,MPI_INT,rank_YZ,recvtag+14,Communicator,&req2[14]);
	MPI_Isend(&sendCount_YZ, 1,MPI_INT,rank_YZ,sendtag+15,Communicator,&req1[15]);
	MPI_Irecv(&recvCount_yz, 1,MPI_INT,rank_yz,recvtag+15,Communicator,&req2[15]);
	MPI_Isend(&sendCount_Yz, 1,MPI_INT,rank_Yz,sendtag+16,Communicator,&req1[16]);
	MPI_Irecv(&recvCount_yZ, 1,MPI_INT,rank_yZ,recvtag+16,Communicator,&req2[16]);
	MPI_Isend(&sendCount_yZ, 1,MPI_INT,rank_yZ,sendtag+17,Communicator,&req1[17]);
	MPI_Irecv(&recvCount_Yz, 1,MPI_INT,rank_Yz,recvtag+17,Communicator,&req2[17]);
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	MPI_Barrier(Communicator);
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
	MPI_Isend(sendList_x, sendCount_x,MPI_INT,rank_x,sendtag,Communicator,&req1[0]);
	MPI_Irecv(recvList_X, recvCount_X,MPI_INT,rank_X,recvtag,Communicator,&req2[0]);
	MPI_Isend(sendList_X, sendCount_X,MPI_INT,rank_X,sendtag,Communicator,&req1[1]);
	MPI_Irecv(recvList_x, recvCount_x,MPI_INT,rank_x,recvtag,Communicator,&req2[1]);
	MPI_Isend(sendList_y, sendCount_y,MPI_INT,rank_y,sendtag,Communicator,&req1[2]);
	MPI_Irecv(recvList_Y, recvCount_Y,MPI_INT,rank_Y,recvtag,Communicator,&req2[2]);
	MPI_Isend(sendList_Y, sendCount_Y,MPI_INT,rank_Y,sendtag,Communicator,&req1[3]);
	MPI_Irecv(recvList_y, recvCount_y,MPI_INT,rank_y,recvtag,Communicator,&req2[3]);
	MPI_Isend(sendList_z, sendCount_z,MPI_INT,rank_z,sendtag,Communicator,&req1[4]);
	MPI_Irecv(recvList_Z, recvCount_Z,MPI_INT,rank_Z,recvtag,Communicator,&req2[4]);
	MPI_Isend(sendList_Z, sendCount_Z,MPI_INT,rank_Z,sendtag,Communicator,&req1[5]);
	MPI_Irecv(recvList_z, recvCount_z,MPI_INT,rank_z,recvtag,Communicator,&req2[5]);
	MPI_Isend(sendList_xy, sendCount_xy,MPI_INT,rank_xy,sendtag,Communicator,&req1[6]);
	MPI_Irecv(recvList_XY, recvCount_XY,MPI_INT,rank_XY,recvtag,Communicator,&req2[6]);
	MPI_Isend(sendList_XY, sendCount_XY,MPI_INT,rank_XY,sendtag,Communicator,&req1[7]);
	MPI_Irecv(recvList_xy, recvCount_xy,MPI_INT,rank_xy,recvtag,Communicator,&req2[7]);
	MPI_Isend(sendList_Xy, sendCount_Xy,MPI_INT,rank_Xy,sendtag,Communicator,&req1[8]);
	MPI_Irecv(recvList_xY, recvCount_xY,MPI_INT,rank_xY,recvtag,Communicator,&req2[8]);
	MPI_Isend(sendList_xY, sendCount_xY,MPI_INT,rank_xY,sendtag,Communicator,&req1[9]);
	MPI_Irecv(recvList_Xy, recvCount_Xy,MPI_INT,rank_Xy,recvtag,Communicator,&req2[9]);
	MPI_Isend(sendList_xz, sendCount_xz,MPI_INT,rank_xz,sendtag,Communicator,&req1[10]);
	MPI_Irecv(recvList_XZ, recvCount_XZ,MPI_INT,rank_XZ,recvtag,Communicator,&req2[10]);
	MPI_Isend(sendList_XZ, sendCount_XZ,MPI_INT,rank_XZ,sendtag,Communicator,&req1[11]);
	MPI_Irecv(recvList_xz, recvCount_xz,MPI_INT,rank_xz,recvtag,Communicator,&req2[11]);
	MPI_Isend(sendList_Xz, sendCount_Xz,MPI_INT,rank_Xz,sendtag,Communicator,&req1[12]);
	MPI_Irecv(recvList_xZ, recvCount_xZ,MPI_INT,rank_xZ,recvtag,Communicator,&req2[12]);
	MPI_Isend(sendList_xZ, sendCount_xZ,MPI_INT,rank_xZ,sendtag,Communicator,&req1[13]);
	MPI_Irecv(recvList_Xz, recvCount_Xz,MPI_INT,rank_Xz,recvtag,Communicator,&req2[13]);
	MPI_Isend(sendList_yz, sendCount_yz,MPI_INT,rank_yz,sendtag,Communicator,&req1[14]);
	MPI_Irecv(recvList_YZ, recvCount_YZ,MPI_INT,rank_YZ,recvtag,Communicator,&req2[14]);
	MPI_Isend(sendList_YZ, sendCount_YZ,MPI_INT,rank_YZ,sendtag,Communicator,&req1[15]);
	MPI_Irecv(recvList_yz, recvCount_yz,MPI_INT,rank_yz,recvtag,Communicator,&req2[15]);
	MPI_Isend(sendList_Yz, sendCount_Yz,MPI_INT,rank_Yz,sendtag,Communicator,&req1[16]);
	MPI_Irecv(recvList_yZ, recvCount_yZ,MPI_INT,rank_yZ,recvtag,Communicator,&req2[16]);
	MPI_Isend(sendList_yZ, sendCount_yZ,MPI_INT,rank_yZ,sendtag,Communicator,&req1[17]);
	MPI_Irecv(recvList_Yz, recvCount_Yz,MPI_INT,rank_Yz,recvtag,Communicator,&req2[17]);
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	//......................................................................................
	for (int idx=0; idx<recvCount_x; idx++)	recvList_x[idx] -= (Nx-2);
	for (int idx=0; idx<recvCount_X; idx++)	recvList_X[idx] += (Nx-2);
	for (int idx=0; idx<recvCount_y; idx++)	recvList_y[idx] -= (Ny-2)*Nx;
	for (int idx=0; idx<recvCount_Y; idx++)	recvList_Y[idx] += (Ny-2)*Nx;
	for (int idx=0; idx<recvCount_z; idx++)	recvList_z[idx] -= (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Z; idx++)	recvList_Z[idx] += (Nz-2)*Nx*Ny;	
	for (int idx=0; idx<recvCount_xy; idx++)	recvList_xy[idx] -= (Nx-2)+(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_XY; idx++)	recvList_XY[idx] += (Nx-2)+(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_xY; idx++)	recvList_xY[idx] -= (Nx-2)-(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_Xy; idx++)	recvList_Xy[idx] += (Nx-2)-(Ny-2)*Nx;
	for (int idx=0; idx<recvCount_xz; idx++)	recvList_xz[idx] -= (Nx-2)+(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_XZ; idx++)	recvList_XZ[idx] += (Nx-2)+(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_xZ; idx++)	recvList_xZ[idx] -= (Nx-2)-(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Xz; idx++)	recvList_Xz[idx] += (Nx-2)-(Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_yz; idx++)	recvList_yz[idx] -= (Ny-2)*Nx + (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_YZ; idx++)	recvList_YZ[idx] += (Ny-2)*Nx + (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_yZ; idx++)	recvList_yZ[idx] -= (Ny-2)*Nx - (Nz-2)*Nx*Ny;
	for (int idx=0; idx<recvCount_Yz; idx++)	recvList_Yz[idx] += (Ny-2)*Nx - (Nz-2)*Nx*Ny;
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

}

inline void Domain::CommunicateMeshHalo(DoubleArray &Mesh)
{
	int sendtag, recvtag;
	sendtag = recvtag = 7;
    double *MeshData = Mesh.get();
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
	MPI_Sendrecv(sendData_x,sendCount_x,MPI_DOUBLE,rank_x,sendtag,
			recvData_X,recvCount_X,MPI_DOUBLE,rank_X,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_X,sendCount_X,MPI_DOUBLE,rank_X,sendtag,
			recvData_x,recvCount_x,MPI_DOUBLE,rank_x,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_y,sendCount_y,MPI_DOUBLE,rank_y,sendtag,
			recvData_Y,recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Y,sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,
			recvData_y,recvCount_y,MPI_DOUBLE,rank_y,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_z,sendCount_z,MPI_DOUBLE,rank_z,sendtag,
			recvData_Z,recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Z,sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,
			recvData_z,recvCount_z,MPI_DOUBLE,rank_z,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xy,sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,
			recvData_XY,recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_XY,sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,
			recvData_xy,recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Xy,sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,
			recvData_xY,recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xY,sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,
			recvData_Xy,recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xz,sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,
			recvData_XZ,recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_XZ,sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,
			recvData_xz,recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Xz,sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,
			recvData_xZ,recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_xZ,sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,
			recvData_Xz,recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_yz,sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,
			recvData_YZ,recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_YZ,sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,
			recvData_yz,recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_Yz,sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,
			recvData_yZ,recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendData_yZ,sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,
			recvData_Yz,recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,Comm,MPI_STATUS_IGNORE);
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

void Domain::BlobComm(MPI_Comm Communicator){
	//......................................................................................
	int sendtag, recvtag;
	sendtag = recvtag = 51;
	//......................................................................................
    int *BlobLabelData = BlobLabel.get();
	PackBlobData(sendList_x, sendCount_x ,sendBuf_x, BlobLabelData);
	PackBlobData(sendList_X, sendCount_X ,sendBuf_X, BlobLabelData);
	PackBlobData(sendList_y, sendCount_y ,sendBuf_y, BlobLabelData);
	PackBlobData(sendList_Y, sendCount_Y ,sendBuf_Y, BlobLabelData);
	PackBlobData(sendList_z, sendCount_z ,sendBuf_z, BlobLabelData);
	PackBlobData(sendList_Z, sendCount_Z ,sendBuf_Z, BlobLabelData);
	PackBlobData(sendList_xy, sendCount_xy ,sendBuf_xy, BlobLabelData);
	PackBlobData(sendList_Xy, sendCount_Xy ,sendBuf_Xy, BlobLabelData);
	PackBlobData(sendList_xY, sendCount_xY ,sendBuf_xY, BlobLabelData);
	PackBlobData(sendList_XY, sendCount_XY ,sendBuf_XY, BlobLabelData);
	PackBlobData(sendList_xz, sendCount_xz ,sendBuf_xz, BlobLabelData);
	PackBlobData(sendList_Xz, sendCount_Xz ,sendBuf_Xz, BlobLabelData);
	PackBlobData(sendList_xZ, sendCount_xZ ,sendBuf_xZ, BlobLabelData);
	PackBlobData(sendList_XZ, sendCount_XZ ,sendBuf_XZ, BlobLabelData);
	PackBlobData(sendList_yz, sendCount_yz ,sendBuf_yz, BlobLabelData);
	PackBlobData(sendList_Yz, sendCount_Yz ,sendBuf_Yz, BlobLabelData);
	PackBlobData(sendList_yZ, sendCount_yZ ,sendBuf_yZ, BlobLabelData);
	PackBlobData(sendList_YZ, sendCount_YZ ,sendBuf_YZ, BlobLabelData);
	//......................................................................................
	MPI_Sendrecv(sendBuf_x,sendCount_x,MPI_INT,rank_x,sendtag,
			recvBuf_X,recvCount_X,MPI_INT,rank_X,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_X,sendCount_X,MPI_INT,rank_X,sendtag,
			recvBuf_x,recvCount_x,MPI_INT,rank_x,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_y,sendCount_y,MPI_INT,rank_y,sendtag,
			recvBuf_Y,recvCount_Y,MPI_INT,rank_Y,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Y,sendCount_Y,MPI_INT,rank_Y,sendtag,
			recvBuf_y,recvCount_y,MPI_INT,rank_y,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_z,sendCount_z,MPI_INT,rank_z,sendtag,
			recvBuf_Z,recvCount_Z,MPI_INT,rank_Z,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Z,sendCount_Z,MPI_INT,rank_Z,sendtag,
			recvBuf_z,recvCount_z,MPI_INT,rank_z,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xy,sendCount_xy,MPI_INT,rank_xy,sendtag,
			recvBuf_XY,recvCount_XY,MPI_INT,rank_XY,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_XY,sendCount_XY,MPI_INT,rank_XY,sendtag,
			recvBuf_xy,recvCount_xy,MPI_INT,rank_xy,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Xy,sendCount_Xy,MPI_INT,rank_Xy,sendtag,
			recvBuf_xY,recvCount_xY,MPI_INT,rank_xY,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xY,sendCount_xY,MPI_INT,rank_xY,sendtag,
			recvBuf_Xy,recvCount_Xy,MPI_INT,rank_Xy,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xz,sendCount_xz,MPI_INT,rank_xz,sendtag,
			recvBuf_XZ,recvCount_XZ,MPI_INT,rank_XZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_XZ,sendCount_XZ,MPI_INT,rank_XZ,sendtag,
			recvBuf_xz,recvCount_xz,MPI_INT,rank_xz,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Xz,sendCount_Xz,MPI_INT,rank_Xz,sendtag,
			recvBuf_xZ,recvCount_xZ,MPI_INT,rank_xZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xZ,sendCount_xZ,MPI_INT,rank_xZ,sendtag,
			recvBuf_Xz,recvCount_Xz,MPI_INT,rank_Xz,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_yz,sendCount_yz,MPI_INT,rank_yz,sendtag,
			recvBuf_YZ,recvCount_YZ,MPI_INT,rank_YZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_YZ,sendCount_YZ,MPI_INT,rank_YZ,sendtag,
			recvBuf_yz,recvCount_yz,MPI_INT,rank_yz,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Yz,sendCount_Yz,MPI_INT,rank_Yz,sendtag,
			recvBuf_yZ,recvCount_yZ,MPI_INT,rank_yZ,recvtag,Comm,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_yZ,sendCount_yZ,MPI_INT,rank_yZ,sendtag,
			recvBuf_Yz,recvCount_Yz,MPI_INT,rank_Yz,recvtag,Comm,MPI_STATUS_IGNORE);
	//........................................................................................
	UnpackBlobData(recvList_x, recvCount_x ,recvBuf_x, BlobLabelData);
	UnpackBlobData(recvList_X, recvCount_X ,recvBuf_X, BlobLabelData);
	UnpackBlobData(recvList_y, recvCount_y ,recvBuf_y, BlobLabelData);
	UnpackBlobData(recvList_Y, recvCount_Y ,recvBuf_Y, BlobLabelData);
	UnpackBlobData(recvList_z, recvCount_z ,recvBuf_z, BlobLabelData);
	UnpackBlobData(recvList_Z, recvCount_Z ,recvBuf_Z, BlobLabelData);
	UnpackBlobData(recvList_xy, recvCount_xy ,recvBuf_xy, BlobLabelData);
	UnpackBlobData(recvList_Xy, recvCount_Xy ,recvBuf_Xy, BlobLabelData);
	UnpackBlobData(recvList_xY, recvCount_xY ,recvBuf_xY, BlobLabelData);
	UnpackBlobData(recvList_XY, recvCount_XY ,recvBuf_XY, BlobLabelData);
	UnpackBlobData(recvList_xz, recvCount_xz ,recvBuf_xz, BlobLabelData);
	UnpackBlobData(recvList_Xz, recvCount_Xz ,recvBuf_Xz, BlobLabelData);
	UnpackBlobData(recvList_xZ, recvCount_xZ ,recvBuf_xZ, BlobLabelData);
	UnpackBlobData(recvList_XZ, recvCount_XZ ,recvBuf_XZ, BlobLabelData);
	UnpackBlobData(recvList_yz, recvCount_yz ,recvBuf_yz, BlobLabelData);
	UnpackBlobData(recvList_Yz, recvCount_Yz ,recvBuf_Yz, BlobLabelData);
	UnpackBlobData(recvList_yZ, recvCount_yZ ,recvBuf_yZ, BlobLabelData);
	UnpackBlobData(recvList_YZ, recvCount_YZ ,recvBuf_YZ, BlobLabelData);
	//......................................................................................
}
inline void SSO(DoubleArray &Distance, char *ID, Domain &Dm, int timesteps){
	/*
	 * This routine converts the data in the Distance array to a signed distance
	 * by solving the equation df/dt = sign(1-|grad f|), where Distance provides
	 * the values of f on the mesh associated with domain Dm
	 * It has been tested with segmented data initialized to values [-1,1]
	 * and will converge toward the signed distance to the surface bounding the associated phases
	 */

	int Q=26;
    int q,i,j,k,n;
    double dt=0.1;
    int in,jn,kn,nn;
    double Dqx,Dqy,Dqz,Dx,Dy,Dz,W;
    double nx,ny,nz,Cqx,Cqy,Cqz,sign,norm;

    const static int D3Q27[26][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
        {1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
        {0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},
        {-1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1}};

    double weights[26];
    // Compute the weights from the finite differences
    for (q=0; q<Q; q++){
        weights[q] = sqrt(1.0*(D3Q27[q][0]*D3Q27[q][0]) + 1.0*(D3Q27[q][1]*D3Q27[q][1]) + 1.0*(D3Q27[q][2]*D3Q27[q][2]));
    }

    int xdim,ydim,zdim;
    xdim=Dm.Nx-2;
    ydim=Dm.Ny-2;
    zdim=Dm.Nz-2;
    fillHalo<double> fillData(Dm.rank_info,xdim,ydim,zdim,1,1,1,0,1);

    int count = 0;
    while (count < timesteps){

        // Communicate the halo of values
        fillData.fill(Distance);

        // Execute the next timestep
        for (k=1;k<Dm.Nz-1;k++){
            for (j=1;j<Dm.Ny-1;j++){
                for (i=1;i<Dm.Nx-1;i++){

                	n = k*Dm.Nx*Dm.Ny+j*Dm.Nx+i;
                    sign = Distance(i,j,k) / fabs(Distance(i,j,k));

              /*
                    if (!(i+1<Nx)) 	nx=0.5*Distance(i,j,k);
                    else 			nx=0.5*Distance(i+1,j,k);;
                    if (!(j+1<Ny)) 	ny=0.5*Distance(i,j,k);
                    else 			ny=0.5*Distance(i,j+1,k);
                    if (!(k+1<Nz)) 	nz=0.5*Distance(i,j,k);
                    else 			nz=0.5*Distance(i,j,k+1);
                    if (i<1)  		nx-=0.5*Distance(i,j,k);
                    else 			nx-=0.5*Distance(i-1,j,k);
                    if (j<1) 		ny-=0.5*Distance(i,j,k);
                    else 			ny-=0.5*Distance(i,j-1,k);
                    if (k<1)  		nz-=0.5*Distance(i,j,k);
                    else 			nz-=0.5*Distance(i,j,k-1);
                    */

                    //............Compute the Gradient...................................
                    nx = 0.5*(Distance(i+1,j,k) - Distance(i-1,j,k));
                    ny = 0.5*(Distance(i,j+1,k) - Distance(i,j-1,k));
                    nz = 0.5*(Distance(i,j,k+1) - Distance(i,j,k-1));

                    W = 0.0;	Dx = Dy = Dz = 0.0;
                    // Ignore any values that have distances less than a lattice unit
                    // since sometimes the positions may be guessed more accurately from
                    // another source (such as a simulation)
                    // also ignore places where the gradient is zero since this will not
                    // result in any local change to Distance
                    if (nx*nx+ny*ny+nz*nz > 0.0 && !(Distance(i,j,k)*Distance(i,j,k) < 1.0) ){
                        for (q=0; q<26; q++){
                            Cqx = 1.0*D3Q27[q][0];
                            Cqy = 1.0*D3Q27[q][1];
                            Cqz = 1.0*D3Q27[q][2];
                            // get the associated neighbor
                            in = i + D3Q27[q][0];
                            jn = j + D3Q27[q][1];
                            kn = k + D3Q27[q][2];

                            // make sure the neighbor is in the domain (periodic BC)
                            /*				if (in < 0 ) in +=Nx;
                             * 	don't need this in parallel since MPI handles the halos
                             if (jn < 0 ) jn +=Ny;
                             if (kn < 0 ) kn +=Nz;
                             if (!(in < Nx) ) in -=Nx;
                             if (!(jn < Ny) ) jn -=Ny;
                             if (!(kn < Nz) ) kn -=Nz;
                             			// symmetric boundary
                            if (in < 0 ) in = i;
                            if (jn < 0 ) jn = j;
                            if (kn < 0 ) kn = k;
                            if (!(in < Nx) ) in = i;
                            if (!(jn < Ny) ) jn = k;
                            if (!(kn < Nz) ) kn = k;
                            */

                            // Compute the gradient using upwind finite differences
                            Dqx = weights[q]*(Distance(i,j,k) - Distance(in,jn,kn))*Cqx;
                            Dqy = weights[q]*(Distance(i,j,k) - Distance(in,jn,kn))*Cqy;
                            Dqz = weights[q]*(Distance(i,j,k) - Distance(in,jn,kn))*Cqz;

                            // Only include upwind derivatives
                            if (sign*(nx*Cqx + ny*Cqy + nz*Cqz) < 0.0 ){

                                Dx += Dqx;
                                Dy += Dqy;
                                Dz += Dqz;
                                W += weights[q];
                            }
                        }
                        // Normalize by the weight to get the approximation to the gradient
                        Dx /= W;
                        Dy /= W;
                        Dz /= W;

                        norm = sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
                    }
                    else{
                        norm = 0.0;
                    }
                    Distance(i,j,k) += dt*sign*(1.0 - norm);

                    // Disallow any change in phase
                    if (Distance(i,j,k)*2.0*(ID[n]-1.0) < 0) Distance(i,j,k) = -Distance(i,j,k);
                }
            }
        }

        count++;
    }
}


inline void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad)
{
	// Read in the full sphere pack
	//...... READ IN THE SPHERES...................................
	cout << "Reading the packing file..." << endl;
	FILE *fid = fopen("pack.out","rb");
	INSIST(fid!=NULL,"Error opening pack.out");
	//.........Trash the header lines..........
	char * line = new char[100];
	fgetl(line, 100, fid);
	fgetl(line, 100, fid);
	fgetl(line, 100, fid);
	fgetl(line, 100, fid);
	fgetl(line, 100, fid);
	//........read the spheres..................
    // We will read until a blank like or end-of-file is reached
	int count = 0;
	while ( !feof(fid) && fgets(line,100,fid)>0 ) {
		char* line2 = line;
		List_cx[count] = strtod(line2,&line2);
		List_cy[count] = strtod(line2,&line2);
		List_cz[count] = strtod(line2,&line2);
		List_rad[count] = strtod(line2,&line2);
		count++;
	}
	cout << "Number of spheres extracted is: " << count << endl;
    INSIST( count==nspheres, "Specified number of spheres is probably incorrect!" );
	// .............................................................
}

inline void AssignLocalSolidID(char *ID, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
						  double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, 
						  int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz)
{
	// Use sphere lists to determine which nodes are in porespace
	// Write out binary file for nodes
	char value;
	int N = Nx*Ny*Nz; 	// Domain size, including the halo
	double hx,hy,hz;
	double x,y,z;
	double cx,cy,cz,r;
	int imin,imax,jmin,jmax,kmin,kmax;
	int p,i,j,k,n;
	//............................................
	double min_x,min_y,min_z;
//	double max_x,max_y,max_z;
	//............................................
	// Lattice spacing for the entire domain
	// It should generally be true that hx=hy=hz
	// Otherwise, you will end up with ellipsoids
	hx = Lx/(Nx*nprocx-1);
	hy = Ly/(Ny*nprocy-1);
	hz = Lz/(Nz*nprocz-1);	
	//............................................
	// Get maximum and minimum for this domain
	// Halo is included !
	min_x = double(iproc*Nx-1)*hx;
	min_y = double(jproc*Ny-1)*hy;
	min_z = double(kproc*Nz-1)*hz;
//	max_x = ((iproc+1)*Nx+1)*hx;
//	max_y = ((jproc+1)*Ny+1)*hy;
//	max_z = ((kproc+1)*Nz+1)*hz;
	//............................................

	//............................................
		// Pre-initialize local ID 
	for (n=0;n<N;n++){
		ID[n]=1;
	}
	//............................................

	//............................................
	// .........Loop over the spheres.............
	for (p=0;p<nspheres;p++){
		// Get the sphere from the list, map to local min
		cx = List_cx[p] - min_x;
		cy = List_cy[p] - min_y;
		cz = List_cz[p] - min_z;
		r = List_rad[p];
		// Check if
		// Range for this sphere in global indexing
		imin = int ((cx-r)/hx)-1;
		imax = int ((cx+r)/hx)+1;
		jmin = int ((cy-r)/hy)-1;
		jmax = int ((cy+r)/hy)+1;
		kmin = int ((cz-r)/hz)-1;
		kmax = int ((cz+r)/hz)+1;
		// Obviously we have to do something at the edges
		if (imin<0)		imin = 0;
		if (imin>Nx)	imin = Nx;
		if (imax<0)		imax = 0;
		if (imax>Nx)	imax = Nx;
		if (jmin<0)		jmin = 0;
		if (jmin>Ny)	jmin = Ny;
		if (jmax<0)		jmax = 0;
		if (jmax>Ny)	jmax = Ny;
		if (kmin<0)		kmin = 0;
		if (kmin>Nz)	kmin = Nz;
		if (kmax<0)		kmax = 0;
		if (kmax>Nz)	kmax = Nz;
		// Loop over the domain for this sphere (may be null)
		for (i=imin;i<imax;i++){
			for (j=jmin;j<jmax;j++){
				for (k=kmin;k<kmax;k++){
					// Initialize ID value to 'fluid (=1)'
					x = i*hx;
					y = j*hy;
					z = k*hz;
					value = 1;
					// if inside sphere, set to zero
					if ( (cx-x)*(cx-x)+(cy-y)*(cy-y)+(cz-z)*(cz-z) < r*r){
						value=0;
					}
					// get the position in the list
					n = k*Nx*Ny+j*Nx+i;
					if ( ID[n] != 0 ){
						ID[n] = value;
					}
				}
			}
		}
	}
}

inline void SignedDistance(double *Distance, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
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
	//............................................

	//............................................
		// Pre-initialize Distance 
	for (n=0;n<N;n++){
		Distance[n]=100.0;
	}
	//............................................

	//............................................
	// .........Loop over the spheres.............
	for (p=0;p<nspheres;p++){
		// Get the sphere from the list, map to local min
		cx = List_cx[p] - min_x;
		cy = List_cy[p] - min_y;
		cz = List_cz[p] - min_z;
		r = List_rad[p];
		// Check if
		// Range for this sphere in global indexing
		imin = int ((cx-2*r)/hx);
		imax = int ((cx+2*r)/hx)+2;
		jmin = int ((cy-2*r)/hy);
		jmax = int ((cy+2*r)/hy)+2;
		kmin = int ((cz-2*r)/hz);
		kmax = int ((cz+2*r)/hz)+2;
		// Obviously we have to do something at the edges
		if (imin<0)		imin = 0;
		if (imin>Nx)	imin = Nx;
		if (imax<0)		imax = 0;
		if (imax>Nx)	imax = Nx;
		if (jmin<0)		jmin = 0;
		if (jmin>Ny)	jmin = Ny;
		if (jmax<0)		jmax = 0;
		if (jmax>Ny)	jmax = Ny;
		if (kmin<0)		kmin = 0;
		if (kmin>Nz)	kmin = Nz;
		if (kmax<0)		kmax = 0;
		if (kmax>Nz)	kmax = Nz;
		// Loop over the domain for this sphere (may be null)
		for (i=imin;i<imax;i++){
			for (j=jmin;j<jmax;j++){
				for (k=kmin;k<kmax;k++){
					// x,y,z is distance in physical units
					x = i*hx;
					y = j*hy;
					z = k*hz;
					// if inside sphere, set to zero
					// get the position in the list
					n = k*Nx*Ny+j*Nx+i;			
					// Compute the distance
					distance = sqrt((cx-x)*(cx-x)+(cy-y)*(cy-y)+(cz-z)*(cz-z)) - r;
					// Assign the minimum distance
					if (distance < Distance[n])		Distance[n] = distance;
				
				}
			}
		}
	}
	
	// Map the distance to lattice units
	for (n=0; n<N; n++)	Distance[n] = Distance[n]/hx;
}

inline void GenerateResidual(char *ID, int Nx, int Ny, int Nz, double Saturation)
{
	//.......................................................................
	int i,j,k,n,Number,N;
	int x,y,z,ii,jj,kk;
	int sizeX,sizeY,sizeZ;
	int *SizeX, *SizeY, *SizeZ;

#ifdef NORANDOM
	srand(10009);
#else
	srand(time(NULL));
#endif
//	float bin;
	//.......................................................................
	N = Nx*Ny*Nz;
	
	int bin, binCount;	
	ifstream Dist("BlobSize.in");
	Dist >> binCount;
//	printf("Number of blob sizes: %i \n",binCount);
	SizeX = new int [binCount];
	SizeY = new int [binCount];
	SizeZ = new int [binCount];
	for (bin=0; bin<binCount; bin++){
		Dist >> SizeX[bin];
		Dist >> SizeY[bin];
		Dist >> SizeZ[bin];
	//	printf("Blob %i dimension: %i x %i x %i \n",bin, SizeX[bin], SizeY[bin], SizeZ[bin]);
	}
	Dist.close();
	//.......................................................................
//	cout << "Generating blocks... " << endl;	
	// Count for the total number of oil nodes 
	int count = 0;
	// Count the total number of non-solid nodes
	int total = 0;
	for (i=0;i<N;i++){
		if (ID[i] != 0) total++;
	}
	
	float sat = 0.f;
	Number = 0;		// number of features
	while (sat < Saturation){
		Number++;
		// Randomly generate a point in the domain
		x = Nx*float(rand())/float(RAND_MAX);
		y = Ny*float(rand())/float(RAND_MAX);
		z = Nz*float(rand())/float(RAND_MAX);
		
		bin = int(floor(binCount*float(rand())/float(RAND_MAX)));
		sizeX = SizeX[bin];
		sizeY = SizeY[bin];
		sizeZ = SizeZ[bin];
		
//		cout << "Sampling from bin no. " << floor(bin) << endl; 
//		cout << "Feature size is: " << sizeX << "x" << sizeY << "x" << sizeZ << endl; 
		
		for (k=z;k<z+sizeZ;k++){
			for (j=y;j<y+sizeY;j++){
				for (i=x;i<x+sizeX;i++){
					// Identify nodes in the domain (periodic BC)
					ii = i;
					jj = j;
					kk = k;					
					if (ii < 1)			ii+=(Nx-2);
					if (jj < 1)			jj+=(Ny-2);
					if (kk < 1)			kk+=(Nz-2);
					if (!(ii < Nx-1))		ii-=(Nx-2);
					if (!(jj < Ny-1))		jj-=(Ny-2);
					if (!(kk < Nz-1))		kk-=(Nz-2);
					
					n = kk*Nx*Ny+jj*Nx+ii;
					
					if (ID[n] == 2){
						ID[n] = 1;
						count++;
					}
				}
			}
		}
		sat = float(count)/total;
	}
	//.......................................................................
}



inline void FlipID(char *ID, int N)
{
	for (int n=0; n<N; n++){
		if  	 (ID[n] == 1)	ID[n] = 2;
		else if  (ID[n] == 2)	ID[n] = 1;
	}
}

inline void WriteLocalSolidID(char *FILENAME, char *ID, int N)
{
	char value;
	ofstream File(FILENAME,ios::binary);
	for (int n=0; n<N; n++){
		value = ID[n];
		File.write((char*) &value, sizeof(value));
	}
	File.close();
}

inline void WriteLocalSolidDistance(char *FILENAME, double *Distance, int N)
{
	double value;
	ofstream File(FILENAME,ios::binary);
	for (int n=0; n<N; n++){
		value = Distance[n];
		File.write((char*) &value, sizeof(value));
	}
	File.close();
}


inline void WriteCheckpoint(char *FILENAME, double *cDen, double *cDistEven, double *cDistOdd, int N)
{
	int q,n;
	double value;
	ofstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		value = cDen[2*n];
		File.write((char*) &value, sizeof(value));
		value = cDen[2*n+1];
		File.write((char*) &value, sizeof(value));
		// Write the even distributions
		for (q=0; q<10; q++){
			value = cDistEven[q*N+n];
			File.write((char*) &value, sizeof(value));
		}
		// Write the odd distributions
		for (q=0; q<9; q++){
			value = cDistOdd[q*N+n];
			File.write((char*) &value, sizeof(value));
		}
	}
	File.close();

}

inline void ReadCheckpoint(char *FILENAME, double *cDen, double *cDistEven, double *cDistOdd, int N)
{
	int q=0, n=0;
	double value=0;
	ifstream File(FILENAME,ios::binary);
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		cDen[2*n] = value;
	//	if (n== 66276)	printf("Density a  = %f \n",value);
		File.read((char*) &value, sizeof(value));
		cDen[2*n+1] = value;
	//	if (n== 66276)	printf("Density b  = %f \n",value);
		// Read the even distributions
		for (q=0; q<10; q++){
			File.read((char*) &value, sizeof(value));
			cDistEven[q*N+n] = value;
	//		if (n== 66276)	printf("dist even %i  = %f \n",q,value);
		}
		// Read the odd distributions
		for (q=0; q<9; q++){
			File.read((char*) &value, sizeof(value));
			cDistOdd[q*N+n] = value;
	//		if (n== 66276)	printf("dist even %i  = %f \n",q,value);
		}
	}
	File.close();
}

inline void ReadBinaryFile(char *FILENAME, double *Data, int N)
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

#endif
