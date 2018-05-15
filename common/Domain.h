#ifndef Domain_INC
#define Domain_INC
// Created by James McClure
// Copyright 2008-2013

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <exception>
#include <stdexcept>

#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "common/Database.h"



//! Read the domain information file
std::shared_ptr<Database> read_domain( );


//! Class to hold domain info
class Domain{
public:
    //! Default constructor
	Domain( std::shared_ptr<Database> db );

    //! Obsolete constructor
    Domain( int nx, int ny, int nz, int rnk, int npx, int npy, int npz, 
                   double lx, double ly, double lz, int BC);

    //! Empty constructor
    Domain() = delete;

    //! Copy constructor
    Domain( const Domain& ) = delete;

    //! Assignment operator
    Domain& operator=( const Domain& ) = delete;

    //! Destructor
	~Domain();
    
    //! Get the database
    inline std::shared_ptr<const Database> getDatabase() const { return d_db; }

private:

    void initialize( std::shared_ptr<Database> db );

    std::shared_ptr<Database> d_db;

public:

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

	void InitializeRanks();
	void CommInit(MPI_Comm comm);
	void CommunicateMeshHalo(DoubleArray &Mesh);
	void AssignComponentLabels(double *phase);

	void TestCommInit(MPI_Comm comm);

    //void MemoryOptimizedLayout(IntArray &Map, int *neighborList, int Np);

private:

	inline int getRankForBlock( int i, int j, int k )
	{
		int i2 = (i+nprocx)%nprocx;
		int j2 = (j+nprocy)%nprocy;
		int k2 = (k+nprocz)%nprocz;
		return i2 + j2*nprocx + k2*nprocx*nprocy;
	}
};


double SSO(DoubleArray &Distance, char *ID, Domain &Dm, int timesteps);

void ReadSpherePacking(int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad);

void AssignLocalSolidID(char *ID, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
                        double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, 
                        int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz);

void SignedDistance(double *Distance, int nspheres, double *List_cx, double *List_cy, double *List_cz, double *List_rad,
                    double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, 
                    int iproc, int jproc, int kproc, int nprocx, int nprocy, int nprocz);

void WriteLocalSolidID(char *FILENAME, char *ID, int N);

void WriteLocalSolidDistance(char *FILENAME, double *Distance, int N);

void WriteCheckpoint(const char *FILENAME, const double *cDen, const double *cfq, int Np);

void ReadCheckpoint(char *FILENAME, double *cDen, double *cfq, int Np);

void ReadBinaryFile(char *FILENAME, double *Data, int N);


#endif
