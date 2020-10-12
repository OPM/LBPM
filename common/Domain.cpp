/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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
	database(NULL), Nx(0), Ny(0), Nz(0), 
	Lx(0), Ly(0), Lz(0), Volume(0), BoundaryCondition(0), voxel_length(1),
	Comm(MPI_COMM_WORLD),
	inlet_layers_x(0), inlet_layers_y(0), inlet_layers_z(0),
    inlet_layers_phase(1),outlet_layers_phase(2),
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
    NULL_USE( rnk );
    NULL_USE( npy );
    NULL_USE( npz );
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
	database(db), Nx(0), Ny(0), Nz(0), 
	Lx(0), Ly(0), Lz(0), Volume(0), BoundaryCondition(0),
	Comm(MPI_COMM_NULL),
	inlet_layers_x(0), inlet_layers_y(0), inlet_layers_z(0),
	outlet_layers_x(0), outlet_layers_y(0), outlet_layers_z(0),
    inlet_layers_phase(1),outlet_layers_phase(2),
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

void Domain::initialize( std::shared_ptr<Database> db )
{	
    d_db = db;
    auto nproc = d_db->getVector<int>("nproc");
    auto n = d_db->getVector<int>("n");

    ASSERT( n.size() == 3u );
    ASSERT( nproc.size() == 3u );
    int nx = n[0];
    int ny = n[1];
    int nz = n[2];
    
	if (d_db->keyExists( "InletLayers" )){
		auto InletCount = d_db->getVector<int>( "InletLayers" );
		inlet_layers_x = InletCount[0];
		inlet_layers_y = InletCount[1];
		inlet_layers_z = InletCount[2];
	}
	if (d_db->keyExists( "OutletLayers" )){
		auto OutletCount = d_db->getVector<int>( "OutletLayers" );
		outlet_layers_x = OutletCount[0];
		outlet_layers_y = OutletCount[1];
		outlet_layers_z = OutletCount[2];
	}
	if (d_db->keyExists( "InletLayersPhase" )){
		inlet_layers_phase = d_db->getScalar<int>( "InletLayersPhase" );
	}
	if (d_db->keyExists( "OutletLayersPhase" )){
		outlet_layers_phase = d_db->getScalar<int>( "OutletLayersPhase" );
	}
    voxel_length = 1.0;
    if (d_db->keyExists( "voxel_length" )){
    	voxel_length = d_db->getScalar<double>("voxel_length");
    }
    else if (d_db->keyExists( "L" )){
    	auto Length = d_db->getVector<double>("L");
    	Lx = Length[0];
    	Ly = Length[1];
    	Lz = Length[2];
    	voxel_length = Lx/(nx*nproc[0]);
    }
    Lx = nx*nproc[0]*voxel_length;
    Ly = ny*nproc[1]*voxel_length;
    Lz = nz*nproc[2]*voxel_length;
    Nx = nx+2;
    Ny = ny+2;
    Nz = nz+2;
    // Initialize ranks
    int myrank;
    MPI_Comm_rank( Comm, &myrank );
	rank_info = RankInfoStruct(myrank,nproc[0],nproc[1],nproc[2]);
	// inlet layers only apply to lower part of domain
	if (rank_info.ix > 0) inlet_layers_x = 0;
	if (rank_info.jy > 0) inlet_layers_y = 0;
	if (rank_info.kz > 0) inlet_layers_z = 0;
	// outlet layers only apply to top part of domain
	if (rank_info.ix < nproc[0]-1 ) outlet_layers_x = 0;
	if (rank_info.jy < nproc[1]-1) outlet_layers_y = 0;
	if (rank_info.kz < nproc[2]-1) outlet_layers_z = 0;
    // Fill remaining variables
	N = Nx*Ny*Nz;
	Volume = nx*ny*nz*nproc[0]*nproc[1]*nproc[2]*1.0;

	if (myrank==0) printf("voxel length = %f micron \n", voxel_length);

	id = new signed char[N];
	memset(id,0,N);
	BoundaryCondition = d_db->getScalar<int>("BC");
    int nprocs;
    MPI_Comm_size( Comm, &nprocs );
	INSIST(nprocs == nproc[0]*nproc[1]*nproc[2],"Fatal error in processor count!");
}

void Domain::Decomp( const std::string& Filename )
{
	//.......................................................................
	// Reading the domain information file
	//.......................................................................
	int rank_offset = 0;
	int RANK = rank();
	int nprocs, nprocx, nprocy, nprocz, nx, ny, nz;
	int64_t global_Nx,global_Ny,global_Nz;
	int64_t i,j,k,n;
	int64_t xStart,yStart,zStart;
	int checkerSize;
	bool USE_CHECKER = false;
	//int inlet_layers_x, inlet_layers_y, inlet_layers_z;
	//int outlet_layers_x, outlet_layers_y, outlet_layers_z;
	xStart=yStart=zStart=0;
	inlet_layers_x = 0;
	inlet_layers_y = 0;
	inlet_layers_z = 0;
	outlet_layers_x = 0;
	outlet_layers_y = 0;
	outlet_layers_z = 0;
	inlet_layers_phase=1;
	outlet_layers_phase=2;
	checkerSize = 32;

	// Read domain parameters
	//auto Filename = database->getScalar<std::string>( "Filename" );
	//auto L = database->getVector<double>( "L" );
	auto size = database->getVector<int>( "n" );
	auto SIZE = database->getVector<int>( "N" );
	auto nproc = database->getVector<int>( "nproc" );
	if (database->keyExists( "offset" )){
		auto offset = database->getVector<int>( "offset" );
		xStart = offset[0];
		yStart = offset[1];
		zStart = offset[2];
	}
	if (database->keyExists( "InletLayers" )){
		auto InletCount = database->getVector<int>( "InletLayers" );
		inlet_layers_x = InletCount[0];
		inlet_layers_y = InletCount[1];
		inlet_layers_z = InletCount[2];
	}
	if (database->keyExists( "OutletLayers" )){
		auto OutletCount = database->getVector<int>( "OutletLayers" );
		outlet_layers_x = OutletCount[0];
		outlet_layers_y = OutletCount[1];
		outlet_layers_z = OutletCount[2];
	}
	if (database->keyExists( "checkerSize" )){
		checkerSize = database->getScalar<int>( "checkerSize" );
		USE_CHECKER = true;
	}
	else {
		checkerSize = SIZE[0];
	}
	if (database->keyExists( "InletLayersPhase" )){
		inlet_layers_phase = database->getScalar<int>( "InletLayersPhase" );
	}
	if (database->keyExists( "OutletLayersPhase" )){
		outlet_layers_phase = database->getScalar<int>( "OutletLayersPhase" );
	}
	auto ReadValues = database->getVector<int>( "ReadValues" );
	auto WriteValues = database->getVector<int>( "WriteValues" );
	auto ReadType = database->getScalar<std::string>( "ReadType" );
	
	if (ReadType == "8bit"){
	}
	else if (ReadType == "16bit"){
	}
	else{
		//printf("INPUT ERROR: Valid ReadType are 8bit, 16bit \n");
		ReadType = "8bit";
	}

	nx = size[0];
	ny = size[1];
	nz = size[2];
	nprocx = nproc[0];
	nprocy = nproc[1];
	nprocz = nproc[2];
	global_Nx = SIZE[0];
	global_Ny = SIZE[1];
	global_Nz = SIZE[2];
	nprocs=nprocx*nprocy*nprocz;
	char *SegData = NULL;
	
	if (RANK==0){
		printf("Input media: %s\n",Filename.c_str());
		printf("Relabeling %lu values\n",ReadValues.size());
		for (size_t idx=0; idx<ReadValues.size(); idx++){
			int oldvalue=ReadValues[idx];
			int newvalue=WriteValues[idx];
			printf("oldvalue=%d, newvalue =%d \n",oldvalue,newvalue);
		}

		// Rank=0 reads the entire segmented data and distributes to worker processes
		printf("Dimensions of segmented image: %ld x %ld x %ld \n",global_Nx,global_Ny,global_Nz);
		int64_t SIZE = global_Nx*global_Ny*global_Nz;
		SegData = new char[SIZE];
		if (ReadType == "8bit"){
			printf("Reading 8-bit input data \n");
			FILE *SEGDAT = fopen(Filename.c_str(),"rb");
			if (SEGDAT==NULL) ERROR("Domain.cpp: Error reading segmented data");
			size_t ReadSeg;
			ReadSeg=fread(SegData,1,SIZE,SEGDAT);
			if (ReadSeg != size_t(SIZE)) printf("Domain.cpp: Error reading segmented data \n");
			fclose(SEGDAT);
		}
		else if (ReadType == "16bit"){
			printf("Reading 16-bit input data \n");
			short int *InputData;
			InputData = new short int[SIZE];
			FILE *SEGDAT = fopen(Filename.c_str(),"rb");
			if (SEGDAT==NULL) ERROR("Domain.cpp: Error reading segmented data");
			size_t ReadSeg;
			ReadSeg=fread(InputData,2,SIZE,SEGDAT);
			if (ReadSeg != size_t(SIZE)) printf("Domain.cpp: Error reading segmented data \n");
			fclose(SEGDAT);
			for (int n=0; n<SIZE; n++){
				SegData[n] = char(InputData[n]);
			}
		}
		printf("Read segmented data from %s \n",Filename.c_str());

		// relabel the data
		std::vector<long int> LabelCount(ReadValues.size(),0);
		for (int k = 0; k<global_Nz; k++){
			for (int j = 0; j<global_Ny; j++){
				for (int i = 0; i<global_Nx; i++){
					n = k*global_Nx*global_Ny+j*global_Nx+i;
					//char locval = loc_id[n];
					char locval = SegData[n];
					for (size_t idx=0; idx<ReadValues.size(); idx++){
						signed char oldvalue=ReadValues[idx];
						signed char newvalue=WriteValues[idx];
						if (locval == oldvalue){
							SegData[n] = newvalue;
							LabelCount[idx]++;
							idx = ReadValues.size();
						}
					}
				}
			}
		}
		for (size_t idx=0; idx<ReadValues.size(); idx++){
			long int label=ReadValues[idx];
			long int count=LabelCount[idx];
			printf("Label=%ld, Count=%ld \n",label,count);
		}
		if (USE_CHECKER) {
			if (inlet_layers_x > 0){
				// use checkerboard pattern
				printf("Checkerboard pattern at x inlet for %i layers \n",inlet_layers_x);
				for (int k = 0; k<global_Nz; k++){
					for (int j = 0; j<global_Ny; j++){
						for (int i = xStart; i < xStart+inlet_layers_x; i++){
							if ( (j/checkerSize + k/checkerSize)%2 == 0){
								// void checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
							}
							else{
								// solid checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 0;
							}
						}
					}
				}
			}

			if (inlet_layers_y > 0){
				printf("Checkerboard pattern at y inlet for %i layers \n",inlet_layers_y);
				// use checkerboard pattern
				for (int k = 0; k<global_Nz; k++){
					for (int j = yStart; j < yStart+inlet_layers_y; j++){
						for (int i = 0; i<global_Nx; i++){
							if ( (i/checkerSize + k/checkerSize)%2 == 0){
								// void checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
							}
							else{
								// solid checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 0;
							}
						}
					}
				}
			}

			if (inlet_layers_z > 0){
				printf("Checkerboard pattern at z inlet for %i layers, saturated with phase label=%i \n",inlet_layers_z,inlet_layers_phase);
				// use checkerboard pattern
				for (int k = zStart; k < zStart+inlet_layers_z; k++){
					for (int j = 0; j<global_Ny; j++){
						for (int i = 0; i<global_Nx; i++){
							if ( (i/checkerSize+j/checkerSize)%2 == 0){
								// void checkers
								//SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = inlet_layers_phase;
							}
							else{
								// solid checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 0;
							}
						}
					}
				}
			}

			if (outlet_layers_x > 0){
				// use checkerboard pattern
				printf("Checkerboard pattern at x outlet for %i layers \n",outlet_layers_x);
				for (int k = 0; k<global_Nz; k++){
					for (int j = 0; j<global_Ny; j++){
						for (int i = xStart + nx*nprocx - outlet_layers_x; i <  xStart + nx*nprocx; i++){
							if ( (j/checkerSize + k/checkerSize)%2 == 0){
								// void checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
							}
							else{
								// solid checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 0;
							}
						}
					}
				}
			}

			if (outlet_layers_y > 0){
				printf("Checkerboard pattern at y outlet for %i layers \n",outlet_layers_y);
				// use checkerboard pattern
				for (int k = 0; k<global_Nz; k++){
					for (int j = yStart + ny*nprocy - outlet_layers_y; j < yStart + ny*nprocy; j++){
						for (int i = 0; i<global_Nx; i++){
							if ( (i/checkerSize + k/checkerSize)%2 == 0){
								// void checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
							}
							else{
								// solid checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 0;

							}
						}
					}
				}
			}

			if (outlet_layers_z > 0){
				printf("Checkerboard pattern at z outlet for %i layers, saturated with phase label=%i \n",outlet_layers_z,outlet_layers_phase);
				// use checkerboard pattern
				for (int k = zStart + nz*nprocz - outlet_layers_z; k < zStart + nz*nprocz; k++){
					for (int j = 0; j<global_Ny; j++){
						for (int i = 0; i<global_Nx; i++){
							if ( (i/checkerSize+j/checkerSize)%2 == 0){
								// void checkers
								//SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 2;
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = outlet_layers_phase;
							}
							else{
								// solid checkers
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = 0;
							}
						}
					}
				}
			}
		}
		else {
			if (inlet_layers_z > 0){
				printf("Mixed reflection pattern at z inlet for %i layers, saturated with phase label=%i \n",inlet_layers_z,inlet_layers_phase);
				for (int k = zStart; k < zStart+inlet_layers_z; k++){
					for (int j = 0; j<global_Ny; j++){
						for (int i = 0; i<global_Nx; i++){
							signed char local_id = SegData[k*global_Nx*global_Ny+j*global_Nx+i];
							signed char reflection_id = SegData[(zStart + nz*nprocz - 1)*global_Nx*global_Ny+j*global_Nx+i];
							if ( local_id < 1 && reflection_id > 0){
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = reflection_id;
							}
						}
					}
				}
			}
			if (outlet_layers_z > 0){
				printf("Mixed reflection pattern at z outlet for %i layers, saturated with phase label=%i \n",outlet_layers_z,outlet_layers_phase);
				for (int k = zStart + nz*nprocz - outlet_layers_z; k < zStart + nz*nprocz; k++){
					for (int j = 0; j<global_Ny; j++){
						for (int i = 0; i<global_Nx; i++){
							signed char local_id = SegData[k*global_Nx*global_Ny+j*global_Nx+i];
							signed char reflection_id = SegData[zStart*global_Nx*global_Ny+j*global_Nx+i];
							if ( local_id < 1 && reflection_id > 0){
								SegData[k*global_Nx*global_Ny+j*global_Nx+i] = reflection_id;
							}
						}
					}
				}
			}
		}
	}
	
	// Get the rank info
	int64_t N = (nx+2)*(ny+2)*(nz+2);

	// number of sites to use for periodic boundary condition transition zone
	int64_t z_transition_size = (nprocz*nz - (global_Nz - zStart))/2;
	if (z_transition_size < 0) z_transition_size=0;

	char LocalRankFilename[40];
	char *loc_id;
	loc_id = new char [(nx+2)*(ny+2)*(nz+2)];

	// Set up the sub-domains
	if (RANK==0){
		printf("Distributing subdomains across %i processors \n",nprocs);
		printf("Process grid: %i x %i x %i \n",nprocx,nprocy,nprocz);
		printf("Subdomain size: %i x %i x %i \n",nx,ny,nz);
		printf("Size of transition region: %ld \n", z_transition_size);

		for (int kp=0; kp<nprocz; kp++){
			for (int jp=0; jp<nprocy; jp++){
				for (int ip=0; ip<nprocx; ip++){
					// rank of the process that gets this subdomain
					int rnk = kp*nprocx*nprocy + jp*nprocx + ip;
					// Pack and send the subdomain for rnk
					for (k=0;k<nz+2;k++){
						for (j=0;j<ny+2;j++){
							for (i=0;i<nx+2;i++){
								int64_t x = xStart + ip*nx + i-1;
								int64_t y = yStart + jp*ny + j-1;
								// int64_t z = zStart + kp*nz + k-1;
								int64_t z = zStart + kp*nz + k-1 - z_transition_size;
								if (x<xStart) 	x=xStart;
								if (!(x<global_Nx))	x=global_Nx-1;
								if (y<yStart) 	y=yStart;
								if (!(y<global_Ny))	y=global_Ny-1;
								if (z<zStart) 	z=zStart;
								if (!(z<global_Nz))	z=global_Nz-1;
								int64_t nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
								int64_t nglobal = z*global_Nx*global_Ny+y*global_Nx+x;
								loc_id[nlocal] = SegData[nglobal];
							}
						}
					}
					if (rnk==0){
						for (k=0;k<nz+2;k++){
							for (j=0;j<ny+2;j++){
								for (i=0;i<nx+2;i++){
									int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
									id[nlocal] = loc_id[nlocal];
								}
							}
						}
					}
					else{
						//printf("Sending data to process %i \n", rnk);
						MPI_Send(loc_id,N,MPI_CHAR,rnk,15,Comm);
					}
					// Write the data for this rank data 
					sprintf(LocalRankFilename,"ID.%05i",rnk+rank_offset);
					FILE *ID = fopen(LocalRankFilename,"wb");
					fwrite(loc_id,1,(nx+2)*(ny+2)*(nz+2),ID);
					fclose(ID);
				}
			}
		}

	}
	else{
		// Recieve the subdomain from rank = 0
		//printf("Ready to recieve data %i at process %i \n", N,rank);
		MPI_Recv(id,N,MPI_CHAR,0,15,Comm,MPI_STATUS_IGNORE);
	}
	//Comm.barrier();
	MPI_Barrier(Comm);
	// Compute the porosity
	double sum;
	double sum_local=0.0;
	double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
	if (BoundaryCondition > 0 && BoundaryCondition !=5) iVol_global = 1.0/(1.0*(Nx-2)*nprocx*(Ny-2)*nprocy*((Nz-2)*nprocz-6));
	//.........................................................
	// If external boundary conditions are applied remove solid
//	if (BoundaryCondition >  0 && BoundaryCondition !=5 && kproc() == 0){
//    	if (inlet_layers_z < 4){
//            inlet_layers_z=4;
//            if(RANK==0){
//                printf("NOTE:Non-periodic BC is applied, but the number of Z-inlet layers is not specified (or is smaller than 3 voxels) \n     the number of Z-inlet layer is reset to %i voxels, saturated with phase label=%i \n",inlet_layers_z-1,inlet_layers_phase);
//            } 
//        }	
//		for (int k=0; k<inlet_layers_z; k++){
//			for (int j=0;j<Ny;j++){
//				for (int i=0;i<Nx;i++){
//					int n = k*Nx*Ny+j*Nx+i;
//					id[n] = inlet_layers_phase;
//				}                    
//			}
// 		}
// 	}
//    if (BoundaryCondition >  0 && BoundaryCondition !=5 && kproc() == nprocz-1){
//    	if (outlet_layers_z < 4){
//            outlet_layers_z=4;
//            if(RANK==nprocs-1){
//                printf("NOTE:Non-periodic BC is applied, but the number of Z-outlet layers is not specified (or is smaller than 3 voxels) \n     the number of Z-outlet layer is reset to %i voxels, saturated with phase label=%i \n",outlet_layers_z-1,outlet_layers_phase);
//            } 
//        }	
// 		for (int k=Nz-outlet_layers_z; k<Nz; k++){
// 			for (int j=0;j<Ny;j++){
// 				for (int i=0;i<Nx;i++){
// 					int n = k*Nx*Ny+j*Nx+i;
// 					id[n] = outlet_layers_phase;
// 				}                    
// 			}
// 		}
// 	}
    for (int k=inlet_layers_z+1; k<Nz-outlet_layers_z-1;k++){
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
    //sum = Comm.sumReduce(sum_local);
    porosity = sum*iVol_global;
    if (rank()==0) printf("Media porosity = %f \n",porosity);
 	//.........................................................
}

void Domain::AggregateLabels( const std::string& filename ){
	
	int nx = Nx;
	int ny = Ny;
	int nz = Nz;
	
	int npx = nprocx();
	int npy = nprocy();
	int npz = nprocz();
	
	int ipx = iproc();
	int ipy = jproc();
	int ipz = kproc();		
	
	int nprocs = nprocx()*nprocy()*nprocz();
		
	int full_nx = npx*(nx-2);
	int full_ny = npy*(ny-2);
	int full_nz = npz*(nz-2);
	int local_size = (nx-2)*(ny-2)*(nz-2);
	unsigned long int full_size = long(full_nx)*long(full_ny)*long(full_nz);
	
	signed char *LocalID;
	LocalID = new signed char [local_size];
		
	//printf("aggregate labels: local size=%i, global size = %i",local_size, full_size);
	// assign the ID for the local sub-region
	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				int n = k*nx*ny+j*nx+i;
				signed char local_id_val = id[n]; 
				LocalID[(k-1)*(nx-2)*(ny-2) + (j-1)*(nx-2) + i-1] = local_id_val;
			}
		}
	}
	MPI_Barrier(Comm);

	// populate the FullID 
	if (rank() == 0){
		signed char *FullID;
		FullID = new signed char [full_size];
		// first handle local ID for rank 0
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					int x = i-1;
					int y = j-1;
					int z = k-1;
					int n_local = (k-1)*(nx-2)*(ny-2) + (j-1)*(nx-2) + i-1;
					unsigned long int n_full = z*long(full_nx)*long(full_ny) + y*long(full_nx) + x;
					FullID[n_full] = LocalID[n_local];
				}
			}
		}
		// next get the local ID from the other ranks
		for (int rnk = 1; rnk<nprocs; rnk++){
			ipz = rnk / (npx*npy);
			ipy = (rnk - ipz*npx*npy) / npx;
			ipx = (rnk - ipz*npx*npy - ipy*npx); 
			//printf("ipx=%i ipy=%i ipz=%i\n", ipx, ipy, ipz);
			int tag = 15+rnk;
			MPI_Recv(LocalID,local_size,MPI_CHAR,rnk,tag,Comm,MPI_STATUS_IGNORE);
			for (int k=1; k<nz-1; k++){
				for (int j=1; j<ny-1; j++){
					for (int i=1; i<nx-1; i++){
						int x = i-1 + ipx*(nx-2);
						int y = j-1 + ipy*(ny-2);
						int z = k-1 + ipz*(nz-2);
						int n_local = (k-1)*(nx-2)*(ny-2) + (j-1)*(nx-2) + i-1;
						unsigned long int n_full = z*long(full_nx)*long(full_ny) + y*long(full_nx) + x;
						FullID[n_full] = LocalID[n_local];
					}
				}
			}
		}
		// write the output
		FILE *OUTFILE = fopen(filename.c_str(),"wb");
		fwrite(FullID,1,full_size,OUTFILE);
		fclose(OUTFILE);
	}
	else{
		// send LocalID to rank=0
		int tag = 15+ rank();
		int dstrank = 0;
		MPI_Send(LocalID,local_size,MPI_CHAR,dstrank,tag,Comm);
	}
	MPI_Barrier(Comm);

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
				if (id[k*Nx*Ny+j*Nx+i] > 0){
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
				if (id[n] > 0){
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
	double sum_local=0.0;
	double iVol_global = 1.0/(1.0*(Nx-2)*(Ny-2)*(Nz-2)*nprocs);
	if (BoundaryCondition > 0 && BoundaryCondition !=5) iVol_global = 1.0/(1.0*(Nx-2)*nprocx()*(Ny-2)*nprocy()*((Nz-2)*nprocz()-6));
	//.........................................................
	// If external boundary conditions are applied remove solid
	if (BoundaryCondition >  0 && BoundaryCondition !=5 && kproc() == 0){
    	if (inlet_layers_z < 4)	inlet_layers_z=4;
		for (int k=0; k<inlet_layers_z; k++){
			for (int j=0;j<Ny;j++){
				for (int i=0;i<Nx;i++){
					int n = k*Nx*Ny+j*Nx+i;
					id[n] = 1;
				}                    
			}
 		}
 	}
    if (BoundaryCondition >  0 && BoundaryCondition !=5 && kproc() == nprocz()-1){
    	if (outlet_layers_z < 4)	outlet_layers_z=4;
 		for (int k=Nz-outlet_layers_z; k<Nz; k++){
 			for (int j=0;j<Ny;j++){
 				for (int i=0;i<Nx;i++){
 					int n = k*Nx*Ny+j*Nx+i;
 					id[n] = 2;
 				}                    
 			}
 		}
 	}
    for (int k=inlet_layers_z+1; k<Nz-outlet_layers_z-1;k++){
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

// Ideally stuff below here should be moved somewhere else -- doesn't really belong here
void WriteCheckpoint(const char *FILENAME, const double *cDen, const double *cfq, size_t Np)
{
    double value;
    ofstream File(FILENAME,ios::binary);
    for (size_t n=0; n<Np; n++){
        // Write the two density values
        value = cDen[n];
        File.write((char*) &value, sizeof(value));
        value = cDen[Np+n];
        File.write((char*) &value, sizeof(value));
        // Write the even distributions
        for (size_t q=0; q<19; q++){
            value = cfq[q*Np+n];
            File.write((char*) &value, sizeof(value));
        }
    }
    File.close();

}

void ReadCheckpoint(char *FILENAME, double *cPhi, double *cfq, size_t Np)
{
    double value=0;
    ifstream File(FILENAME,ios::binary);
    for (size_t n=0; n<Np; n++){
        File.read((char*) &value, sizeof(value));
        cPhi[n] = value;
        // Read the distributions
        for (size_t q=0; q<19; q++){
            File.read((char*) &value, sizeof(value));
            cfq[q*Np+n] = value;
        }
    }
    File.close();
}

void ReadBinaryFile(char *FILENAME, double *Data, size_t N)
{
  double value;
  ifstream File(FILENAME,ios::binary);
  if (File.good()){
    for (size_t n=0; n<N; n++){
      // Write the two density values                                                                                
      File.read((char*) &value, sizeof(value));
      Data[n] = value;

    }
  }
  else {
    for (size_t n=0; n<N; n++) Data[n] = 1.2e-34;
  }
  File.close();
}

void Domain::ReadFromFile(const std::string& Filename,const std::string& Datatype, double *UserData)
{
	//........................................................................................
	// Reading the user-defined input file
    // NOTE: so far it only supports BC=0 (periodic) and BC=5 (mixed reflection)
    //       because if checkerboard or inlet/outlet buffer layers are added, the
    //       value of the void space is undefined. 
    // NOTE: if BC=5 is used, where the inlet and outlet layers of the domain are modified, 
    //       user needs to modify the input file accordingly before LBPM simulator read
    //       the input file.
	//........................................................................................
	int rank_offset = 0;
	int RANK = rank();
	int nprocs, nprocx, nprocy, nprocz, nx, ny, nz;
	int64_t global_Nx,global_Ny,global_Nz;
	int64_t i,j,k,n;
    //TODO These offset we may still need them
	int64_t xStart,yStart,zStart;
	xStart=yStart=zStart=0;

	// Read domain parameters
    // TODO currently the size of the data is still read from Domain{}; 
    //      but user may have a user-specified size
	auto size = database->getVector<int>( "n" );
	auto SIZE = database->getVector<int>( "N" );
	auto nproc = database->getVector<int>( "nproc" );
    //TODO currently the funcationality "offset" is disabled as the user-defined input data may have a different size from that of the input domain 
	if (database->keyExists( "offset" )){
		auto offset = database->getVector<int>( "offset" );
		xStart = offset[0];
		yStart = offset[1];
		zStart = offset[2];
	}
	
	nx = size[0];
	ny = size[1];
	nz = size[2];
	nprocx = nproc[0];
	nprocy = nproc[1];
	nprocz = nproc[2];
	global_Nx = SIZE[0];
	global_Ny = SIZE[1];
	global_Nz = SIZE[2];
	nprocs=nprocx*nprocy*nprocz;

	double *SegData = NULL;
	if (RANK==0){
		printf("User-defined input file: %s (data type: %s)\n",Filename.c_str(),Datatype.c_str());
        printf("NOTE: currently only BC=0 or 5 supports user-defined input file!\n");
		// Rank=0 reads the entire segmented data and distributes to worker processes
		printf("Dimensions of the user-defined input file: %ld x %ld x %ld \n",global_Nx,global_Ny,global_Nz);
		int64_t SIZE = global_Nx*global_Ny*global_Nz;

		if (Datatype == "double"){
			printf("Reading input data as double precision floating number\n");
		    SegData = new double[SIZE];
			FILE *SEGDAT = fopen(Filename.c_str(),"rb");
			if (SEGDAT==NULL) ERROR("Domain.cpp: Error reading user-defined file!\n");
			size_t ReadSeg;
			ReadSeg=fread(SegData,8,SIZE,SEGDAT);
			if (ReadSeg != size_t(SIZE)) printf("Domain.cpp: Error reading file: %s\n",Filename.c_str());
			fclose(SEGDAT);
		}
        else{
            ERROR("Error: User-defined input file only supports double-precision floating number!\n"); 
        }
		printf("Read file successfully from %s \n",Filename.c_str());
	}
	
	// Get the rank info
	int64_t N = (nx+2)*(ny+2)*(nz+2);

	// number of sites to use for periodic boundary condition transition zone
	//int64_t z_transition_size = (nprocz*nz - (global_Nz - zStart))/2;
	//if (z_transition_size < 0) z_transition_size=0;
	int64_t z_transition_size = 0;

	//char LocalRankFilename[1000];//just for debug
	double *loc_id;
	loc_id = new double [(nx+2)*(ny+2)*(nz+2)];

	// Set up the sub-domains
	if (RANK==0){
        printf("Decomposing user-defined input file\n");
		printf("Distributing subdomains across %i processors \n",nprocs);
		printf("Process grid: %i x %i x %i \n",nprocx,nprocy,nprocz);
		printf("Subdomain size: %i x %i x %i \n",nx,ny,nz);
		printf("Size of transition region: %ld \n", z_transition_size);

		for (int kp=0; kp<nprocz; kp++){
			for (int jp=0; jp<nprocy; jp++){
				for (int ip=0; ip<nprocx; ip++){
					// rank of the process that gets this subdomain
					int rnk = kp*nprocx*nprocy + jp*nprocx + ip;
					// Pack and send the subdomain for rnk
					for (k=0;k<nz+2;k++){
						for (j=0;j<ny+2;j++){
							for (i=0;i<nx+2;i++){
								int64_t x = xStart + ip*nx + i-1;
								int64_t y = yStart + jp*ny + j-1;
								// int64_t z = zStart + kp*nz + k-1;
								int64_t z = zStart + kp*nz + k-1 - z_transition_size;
								if (x<xStart) 	x=xStart;
								if (!(x<global_Nx))	x=global_Nx-1;
								if (y<yStart) 	y=yStart;
								if (!(y<global_Ny))	y=global_Ny-1;
								if (z<zStart) 	z=zStart;
								if (!(z<global_Nz))	z=global_Nz-1;
								int64_t nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
								int64_t nglobal = z*global_Nx*global_Ny+y*global_Nx+x;
								loc_id[nlocal] = SegData[nglobal];
							}
						}
					}
					if (rnk==0){
						for (k=0;k<nz+2;k++){
							for (j=0;j<ny+2;j++){
								for (i=0;i<nx+2;i++){
									int nlocal = k*(nx+2)*(ny+2) + j*(nx+2) + i;
									UserData[nlocal] = loc_id[nlocal];
								}
							}
						}
					}
					else{
						//printf("Sending data to process %i \n", rnk);
						MPI_Send(loc_id,N,MPI_DOUBLE,rnk,15,Comm);
					}
					// Write the data for this rank data 
                    // NOTE just for debug
					//sprintf(LocalRankFilename,"%s.%05i",Filename.c_str(),rnk+rank_offset);
					//FILE *ID = fopen(LocalRankFilename,"wb");
					//fwrite(loc_id,8,(nx+2)*(ny+2)*(nz+2),ID);
					//fclose(ID);
				}
			}
		}

	}
	else{
		// Recieve the subdomain from rank = 0
		//printf("Ready to recieve data %i at process %i \n", N,rank);
		MPI_Recv(UserData,N,MPI_DOUBLE,0,15,Comm,MPI_STATUS_IGNORE);
	}
	//Comm.barrier();
	MPI_Barrier(Comm);
}
