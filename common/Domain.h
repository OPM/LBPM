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
#include <mpi.h>
#include "common/Utilities.h"


using namespace std;

struct Domain{
	
	Domain(int nx, int ny, int nz, int rnk, int npx, int npy, int npz, 
			double lx, double ly, double lz){
		Nx = nx+2; Ny = ny+2; Nz = nz+2; 
		Lx = lx, Ly = ly, Lz = lz;
		rank = rnk;
		nprocx=npx; nprocy=npy; nprocz=npz;
		N = Nx*Ny*Nz;
		id = new char [N];
		Blobs.New(Nx,Ny,Nz);
	}
		
	// Basic domain information
	int Nx,Ny,Nz,N;
	int iproc,jproc,kproc;
	int nprocx,nprocy,nprocz;
	double Lx,Ly,Lz;
	int rank;
	MPI_Comm Communicator;

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
	
	// Solid indicator function
	char *id;
	// Blob information
	IntArray Blobs;
	
	void InitializeRanks();
	void CommInit(MPI_Comm comm);
	void BlobComm(MPI_Comm comm);
	
private:
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
	int getRankForBlock( int i, int j, int k )
	{
		int i2 = (i+nprocx)%nprocx;
		int j2 = (j+nprocy)%nprocy;
		int k2 = (k+nprocz)%nprocz;
		return i2 + j2*nprocx + k2*nprocx*nprocy;
	}
	
};

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

}

void Domain::BlobComm(MPI_Comm Communicator){
	//......................................................................................
	int sendtag, recvtag;
	sendtag = recvtag = 51;
	//......................................................................................
	PackBlobData(sendList_x, sendCount_x ,sendBuf_x, Blobs.data);
	PackBlobData(sendList_X, sendCount_X ,sendBuf_X, Blobs.data);
	PackBlobData(sendList_y, sendCount_y ,sendBuf_y, Blobs.data);
	PackBlobData(sendList_Y, sendCount_Y ,sendBuf_Y, Blobs.data);
	PackBlobData(sendList_z, sendCount_z ,sendBuf_z, Blobs.data);
	PackBlobData(sendList_Z, sendCount_Z ,sendBuf_Z, Blobs.data);
	PackBlobData(sendList_xy, sendCount_xy ,sendBuf_xy, Blobs.data);
	PackBlobData(sendList_Xy, sendCount_Xy ,sendBuf_Xy, Blobs.data);
	PackBlobData(sendList_xY, sendCount_xY ,sendBuf_xY, Blobs.data);
	PackBlobData(sendList_XY, sendCount_XY ,sendBuf_XY, Blobs.data);
	PackBlobData(sendList_xz, sendCount_xz ,sendBuf_xz, Blobs.data);
	PackBlobData(sendList_Xz, sendCount_Xz ,sendBuf_Xz, Blobs.data);
	PackBlobData(sendList_xZ, sendCount_xZ ,sendBuf_xZ, Blobs.data);
	PackBlobData(sendList_XZ, sendCount_XZ ,sendBuf_XZ, Blobs.data);
	PackBlobData(sendList_yz, sendCount_yz ,sendBuf_yz, Blobs.data);
	PackBlobData(sendList_Yz, sendCount_Yz ,sendBuf_Yz, Blobs.data);
	PackBlobData(sendList_yZ, sendCount_yZ ,sendBuf_yZ, Blobs.data);
	PackBlobData(sendList_YZ, sendCount_YZ ,sendBuf_YZ, Blobs.data);
	//......................................................................................
	MPI_Sendrecv(sendBuf_x,sendCount_x,MPI_INT,rank_x,sendtag,
			recvBuf_X,recvCount_X,MPI_INT,rank_X,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_X,sendCount_X,MPI_INT,rank_X,sendtag,
			recvBuf_x,recvCount_x,MPI_INT,rank_x,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_y,sendCount_y,MPI_INT,rank_y,sendtag,
			recvBuf_Y,recvCount_Y,MPI_INT,rank_Y,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Y,sendCount_Y,MPI_INT,rank_Y,sendtag,
			recvBuf_y,recvCount_y,MPI_INT,rank_y,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_z,sendCount_z,MPI_INT,rank_z,sendtag,
			recvBuf_Z,recvCount_Z,MPI_INT,rank_Z,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Z,sendCount_Z,MPI_INT,rank_Z,sendtag,
			recvBuf_z,recvCount_z,MPI_INT,rank_z,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xy,sendCount_xy,MPI_INT,rank_xy,sendtag,
			recvBuf_XY,recvCount_XY,MPI_INT,rank_XY,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_XY,sendCount_XY,MPI_INT,rank_XY,sendtag,
			recvBuf_xy,recvCount_xy,MPI_INT,rank_xy,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Xy,sendCount_Xy,MPI_INT,rank_Xy,sendtag,
			recvBuf_xY,recvCount_xY,MPI_INT,rank_xY,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xY,sendCount_xY,MPI_INT,rank_xY,sendtag,
			recvBuf_Xy,recvCount_Xy,MPI_INT,rank_Xy,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xz,sendCount_xz,MPI_INT,rank_xz,sendtag,
			recvBuf_XZ,recvCount_XZ,MPI_INT,rank_XZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_XZ,sendCount_XZ,MPI_INT,rank_XZ,sendtag,
			recvBuf_xz,recvCount_xz,MPI_INT,rank_xz,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Xz,sendCount_Xz,MPI_INT,rank_Xz,sendtag,
			recvBuf_xZ,recvCount_xZ,MPI_INT,rank_xZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_xZ,sendCount_xZ,MPI_INT,rank_xZ,sendtag,
			recvBuf_Xz,recvCount_Xz,MPI_INT,rank_Xz,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_yz,sendCount_yz,MPI_INT,rank_yz,sendtag,
			recvBuf_YZ,recvCount_YZ,MPI_INT,rank_YZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_YZ,sendCount_YZ,MPI_INT,rank_YZ,sendtag,
			recvBuf_yz,recvCount_yz,MPI_INT,rank_yz,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_Yz,sendCount_Yz,MPI_INT,rank_Yz,sendtag,
			recvBuf_yZ,recvCount_yZ,MPI_INT,rank_yZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendBuf_yZ,sendCount_yZ,MPI_INT,rank_yZ,sendtag,
			recvBuf_Yz,recvCount_Yz,MPI_INT,rank_Yz,recvtag,Communicator,MPI_STATUS_IGNORE);
	//........................................................................................
	UnpackBlobData(recvList_x, recvCount_x ,recvBuf_x, Blobs.data);
	UnpackBlobData(recvList_X, recvCount_X ,recvBuf_X, Blobs.data);
	UnpackBlobData(recvList_y, recvCount_y ,recvBuf_y, Blobs.data);
	UnpackBlobData(recvList_Y, recvCount_Y ,recvBuf_Y, Blobs.data);
	UnpackBlobData(recvList_z, recvCount_z ,recvBuf_z, Blobs.data);
	UnpackBlobData(recvList_Z, recvCount_Z ,recvBuf_Z, Blobs.data);
	UnpackBlobData(recvList_xy, recvCount_xy ,recvBuf_xy, Blobs.data);
	UnpackBlobData(recvList_Xy, recvCount_Xy ,recvBuf_Xy, Blobs.data);
	UnpackBlobData(recvList_xY, recvCount_xY ,recvBuf_xY, Blobs.data);
	UnpackBlobData(recvList_XY, recvCount_XY ,recvBuf_XY, Blobs.data);
	UnpackBlobData(recvList_xz, recvCount_xz ,recvBuf_xz, Blobs.data);
	UnpackBlobData(recvList_Xz, recvCount_Xz ,recvBuf_Xz, Blobs.data);
	UnpackBlobData(recvList_xZ, recvCount_xZ ,recvBuf_xZ, Blobs.data);
	UnpackBlobData(recvList_XZ, recvCount_XZ ,recvBuf_XZ, Blobs.data);
	UnpackBlobData(recvList_yz, recvCount_yz ,recvBuf_yz, Blobs.data);
	UnpackBlobData(recvList_Yz, recvCount_Yz ,recvBuf_Yz, Blobs.data);
	UnpackBlobData(recvList_yZ, recvCount_yZ ,recvBuf_yZ, Blobs.data);
	UnpackBlobData(recvList_YZ, recvCount_YZ ,recvBuf_YZ, Blobs.data);
	//......................................................................................
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
	fgets(line, 100, fid);
	fgets(line, 100, fid);
	fgets(line, 100, fid);
	fgets(line, 100, fid);
	fgets(line, 100, fid);
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
		imin = int ((cx-r)/hx);
		imax = int ((cx+r)/hx)+2;
		jmin = int ((cy-r)/hy);
		jmax = int ((cy+r)/hy)+2;
		kmin = int ((cz-r)/hz);
		kmax = int ((cz+r)/hz)+2;
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
	int q,n;
	double value;
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
	for (n=0; n<N; n++){
		// Write the two density values
		File.read((char*) &value, sizeof(value));
		Data[n] = value;

	}
	File.close();
}

