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
#include "common/ScaLBL.h"

ScaLBL_Communicator::ScaLBL_Communicator(std::shared_ptr <Domain> Dm){
	//......................................................................................
	Lock=false; // unlock the communicator
	//......................................................................................
	// Create a separate copy of the communicator for the device
	//MPI_Comm_group(Dm->Comm,&Group);
	//MPI_Comm_create(Dm->Comm,Group,&MPI_COMM_SCALBL);
	MPI_Comm_dup(Dm->Comm,&MPI_COMM_SCALBL);
	//......................................................................................
	// Copy the domain size and communication information directly from Dm
	Nx = Dm->Nx;
	Ny = Dm->Ny;
	Nz = Dm->Nz;
	N = Nx*Ny*Nz;
	next=0;
	rank=Dm->rank();
	rank_x=Dm->rank_x();
	rank_y=Dm->rank_y();
	rank_z=Dm->rank_z();
	rank_X=Dm->rank_X();
	rank_Y=Dm->rank_Y();
	rank_Z=Dm->rank_Z();
	rank_xy=Dm->rank_xy();
	rank_XY=Dm->rank_XY();
	rank_xY=Dm->rank_xY();
	rank_Xy=Dm->rank_Xy();
	rank_xz=Dm->rank_xz();
	rank_XZ=Dm->rank_XZ();
	rank_xZ=Dm->rank_xZ();
	rank_Xz=Dm->rank_Xz();
	rank_yz=Dm->rank_yz();
	rank_YZ=Dm->rank_YZ();
	rank_yZ=Dm->rank_yZ();
	rank_Yz=Dm->rank_Yz();
	sendCount_x=Dm->sendCount_x;
	sendCount_y=Dm->sendCount_y;
	sendCount_z=Dm->sendCount_z;
	sendCount_X=Dm->sendCount_X;
	sendCount_Y=Dm->sendCount_Y;
	sendCount_Z=Dm->sendCount_Z;
	sendCount_xy=Dm->sendCount_xy;
	sendCount_yz=Dm->sendCount_yz;
	sendCount_xz=Dm->sendCount_xz;
	sendCount_Xy=Dm->sendCount_Xy;
	sendCount_Yz=Dm->sendCount_Yz;
	sendCount_xZ=Dm->sendCount_xZ;
	sendCount_xY=Dm->sendCount_xY;
	sendCount_yZ=Dm->sendCount_yZ;
	sendCount_Xz=Dm->sendCount_Xz;
	sendCount_XY=Dm->sendCount_XY;
	sendCount_YZ=Dm->sendCount_YZ;
	sendCount_XZ=Dm->sendCount_XZ;
	recvCount_x=Dm->recvCount_x;
	recvCount_y=Dm->recvCount_y;
	recvCount_z=Dm->recvCount_z;
	recvCount_X=Dm->recvCount_X;
	recvCount_Y=Dm->recvCount_Y;
	recvCount_Z=Dm->recvCount_Z;
	recvCount_xy=Dm->recvCount_xy;
	recvCount_yz=Dm->recvCount_yz;
	recvCount_xz=Dm->recvCount_xz;
	recvCount_Xy=Dm->recvCount_Xy;
	recvCount_Yz=Dm->recvCount_Yz;
	recvCount_xZ=Dm->recvCount_xZ;
	recvCount_xY=Dm->recvCount_xY;
	recvCount_yZ=Dm->recvCount_yZ;
	recvCount_Xz=Dm->recvCount_Xz;
	recvCount_XY=Dm->recvCount_XY;
	recvCount_YZ=Dm->recvCount_YZ;
	recvCount_XZ=Dm->recvCount_XZ;
	
	iproc = Dm->iproc();
	jproc = Dm->jproc();
	kproc = Dm->kproc();
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
	BoundaryCondition = Dm->BoundaryCondition;
	//......................................................................................

	ScaLBL_AllocateZeroCopy((void **) &sendbuf_x, 5*sendCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_X, 5*sendCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_y, 5*sendCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Y, 5*sendCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_z, 5*sendCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Z, 5*sendCount_Z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xy, sendCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xY, sendCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Xy, sendCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_XY, sendCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xz, sendCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xZ, sendCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Xz, sendCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_XZ, sendCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_yz, sendCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_yZ, sendCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Yz, sendCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_YZ, sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_x, 5*recvCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_X, 5*recvCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_y, 5*recvCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Y, 5*recvCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_z, 5*recvCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Z, 5*recvCount_Z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xy, recvCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xY, recvCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Xy, recvCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_XY, recvCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xz, recvCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xZ, recvCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Xz, recvCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_XZ, recvCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_yz, recvCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_yZ, recvCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Yz, recvCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_YZ, recvCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_x, sendCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_X, sendCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_y, sendCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_Y, sendCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_z, sendCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_Z, sendCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xy, sendCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xY, sendCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_Xy, sendCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_XY, sendCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xz, sendCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xZ, sendCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_Xz, sendCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_XZ, sendCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_yz, sendCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_yZ, sendCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_Yz, sendCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_YZ, sendCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_x, recvCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_X, recvCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_y, recvCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_Y, recvCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_z, recvCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_Z, recvCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_x, 5*recvCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_X, 5*recvCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_y, 5*recvCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_Y, 5*recvCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_z, 5*recvCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_Z, 5*recvCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvDist_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
	//......................................................................................

	ScaLBL_CopyToZeroCopy(dvcSendList_x,Dm->sendList_x,sendCount_x*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_X,Dm->sendList_X,sendCount_X*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_y,Dm->sendList_y,sendCount_y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Y,Dm->sendList_Y,sendCount_Y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_z,Dm->sendList_z,sendCount_z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Z,Dm->sendList_Z,sendCount_Z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xy,Dm->sendList_xy,sendCount_xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_XY,Dm->sendList_XY,sendCount_XY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xY,Dm->sendList_xY,sendCount_xY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Xy,Dm->sendList_Xy,sendCount_Xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xz,Dm->sendList_xz,sendCount_xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_XZ,Dm->sendList_XZ,sendCount_XZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xZ,Dm->sendList_xZ,sendCount_xZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Xz,Dm->sendList_Xz,sendCount_Xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_yz,Dm->sendList_yz,sendCount_yz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_YZ,Dm->sendList_YZ,sendCount_YZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_yZ,Dm->sendList_yZ,sendCount_yZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Yz,Dm->sendList_Yz,sendCount_Yz*sizeof(int));
	//......................................................................................
	ScaLBL_CopyToZeroCopy(dvcRecvList_x,Dm->recvList_x,recvCount_x*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_X,Dm->recvList_X,recvCount_X*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_y,Dm->recvList_y,recvCount_y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Y,Dm->recvList_Y,recvCount_Y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_z,Dm->recvList_z,recvCount_z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Z,Dm->recvList_Z,recvCount_Z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xy,Dm->recvList_xy,recvCount_xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_XY,Dm->recvList_XY,recvCount_XY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xY,Dm->recvList_xY,recvCount_xY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Xy,Dm->recvList_Xy,recvCount_Xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xz,Dm->recvList_xz,recvCount_xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_XZ,Dm->recvList_XZ,recvCount_XZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xZ,Dm->recvList_xZ,recvCount_xZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Xz,Dm->recvList_Xz,recvCount_Xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_yz,Dm->recvList_yz,recvCount_yz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_YZ,Dm->recvList_YZ,recvCount_YZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_yZ,Dm->recvList_yZ,recvCount_yZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Yz,Dm->recvList_Yz,recvCount_Yz*sizeof(int));
	//......................................................................................

	MPI_Barrier(MPI_COMM_SCALBL);

	//...................................................................................
	// Set up the recieve distribution lists
	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,14 .................................
	D3Q19_MapRecv(-1,0,0,Dm->recvList_X,0,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,-1,0,Dm->recvList_X,recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,1,0,Dm->recvList_X,2*recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,0,-1,Dm->recvList_X,3*recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,0,1,Dm->recvList_X,4*recvCount_X,recvCount_X,dvcRecvDist_X);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	D3Q19_MapRecv(1,0,0,Dm->recvList_x,0,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,1,0,Dm->recvList_x,recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,-1,0,Dm->recvList_x,2*recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,0,1,Dm->recvList_x,3*recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,0,-1,Dm->recvList_x,4*recvCount_x,recvCount_x,dvcRecvDist_x);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
	D3Q19_MapRecv(0,-1,0,Dm->recvList_Y,0,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(-1,-1,0,Dm->recvList_Y,recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(1,-1,0,Dm->recvList_Y,2*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(0,-1,-1,Dm->recvList_Y,3*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(0,-1,1,Dm->recvList_Y,4*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
	D3Q19_MapRecv(0,1,0,Dm->recvList_y,0,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(1,1,0,Dm->recvList_y,recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(-1,1,0,Dm->recvList_y,2*recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(0,1,1,Dm->recvList_y,3*recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(0,1,-1,Dm->recvList_y,4*recvCount_y,recvCount_y,dvcRecvDist_y);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
	D3Q19_MapRecv(0,0,-1,Dm->recvList_Z,0,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(-1,0,-1,Dm->recvList_Z,recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(1,0,-1,Dm->recvList_Z,2*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(0,-1,-1,Dm->recvList_Z,3*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(0,1,-1,Dm->recvList_Z,4*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
	D3Q19_MapRecv(0,0,1,Dm->recvList_z,0,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(1,0,1,Dm->recvList_z,recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(-1,0,1,Dm->recvList_z,2*recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(0,1,1,Dm->recvList_z,3*recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(0,-1,1,Dm->recvList_z,4*recvCount_z,recvCount_z,dvcRecvDist_z);
	//..................................................................................
	//...Map recieve list for the xy edge <<<8)................................
	D3Q19_MapRecv(-1,-1,0,Dm->recvList_XY,0,recvCount_XY,dvcRecvDist_XY);
	//...Map recieve list for the Xy edge <<<9)................................
	D3Q19_MapRecv(1,-1,0,Dm->recvList_xY,0,recvCount_xY,dvcRecvDist_xY);
	//...Map recieve list for the xY edge <<<10)................................
	D3Q19_MapRecv(-1,1,0,Dm->recvList_Xy,0,recvCount_Xy,dvcRecvDist_Xy);
	//...Map recieve list for the XY edge <<<7)................................
	D3Q19_MapRecv(1,1,0,Dm->recvList_xy,0,recvCount_xy,dvcRecvDist_xy);
	//...Map recieve list for the xz edge <<<12)................................
	D3Q19_MapRecv(-1,0,-1,Dm->recvList_XZ,0,recvCount_XZ,dvcRecvDist_XZ);
	//...Map recieve list for the xZ edge <<<14)................................
	D3Q19_MapRecv(-1,0,1,Dm->recvList_Xz,0,recvCount_Xz,dvcRecvDist_Xz);
	//...Map recieve list for the Xz edge <<<13)................................
	D3Q19_MapRecv(1,0,-1,Dm->recvList_xZ,0,recvCount_xZ,dvcRecvDist_xZ);
	//...Map recieve list for the XZ edge <<<11)................................
	D3Q19_MapRecv(1,0,1,Dm->recvList_xz,0,recvCount_xz,dvcRecvDist_xz);
	//...Map recieve list for the yz edge <<<16)................................
	D3Q19_MapRecv(0,-1,-1,Dm->recvList_YZ,0,recvCount_YZ,dvcRecvDist_YZ);
	//...Map recieve list for the yZ edge <<<18)................................
	D3Q19_MapRecv(0,-1,1,Dm->recvList_Yz,0,recvCount_Yz,dvcRecvDist_Yz);
	//...Map recieve list for the Yz edge <<<17)................................
	D3Q19_MapRecv(0,1,-1,Dm->recvList_yZ,0,recvCount_yZ,dvcRecvDist_yZ);
	//...Map recieve list for the YZ edge <<<15)................................
	D3Q19_MapRecv(0,1,1,Dm->recvList_yz,0,recvCount_yz,dvcRecvDist_yz);
	//...................................................................................

	//......................................................................................
	MPI_Barrier(MPI_COMM_SCALBL);
	ScaLBL_DeviceBarrier();
	//......................................................................................
	SendCount = sendCount_x+sendCount_X+sendCount_y+sendCount_Y+sendCount_z+sendCount_Z+
			sendCount_xy+sendCount_Xy+sendCount_xY+sendCount_XY+
			sendCount_xZ+sendCount_Xz+sendCount_xZ+sendCount_XZ+
			sendCount_yz+sendCount_Yz+sendCount_yZ+sendCount_YZ;

	RecvCount = recvCount_x+recvCount_X+recvCount_y+recvCount_Y+recvCount_z+recvCount_Z+
			recvCount_xy+recvCount_Xy+recvCount_xY+recvCount_XY+
			recvCount_xZ+recvCount_Xz+recvCount_xZ+recvCount_XZ+
			recvCount_yz+recvCount_Yz+recvCount_yZ+recvCount_YZ;

	CommunicationCount = SendCount+RecvCount;
	//......................................................................................

}


ScaLBL_Communicator::~ScaLBL_Communicator(){
	// destrutor does nothing (bad idea)
	// -- note that there needs to be a way to free memory allocated on the device!!!
}
int ScaLBL_Communicator::LastExterior(){
	return next;
}
int ScaLBL_Communicator::FirstInterior(){
	return first_interior;
}
int ScaLBL_Communicator::LastInterior(){
	return last_interior;
}

void ScaLBL_Communicator::D3Q19_MapRecv(int Cqx, int Cqy, int Cqz, int *list,  int start, int count,
		int *d3q19_recvlist){
	int i,j,k,n,nn,idx;
	int * ReturnDist;
	ReturnDist=new int [count];
	
	for (idx=0; idx<count; idx++){

		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx]; // if (rank == 0) printf("@ rank:%d n=%d\n",rank,n);
		// Get the 3-D indices from the send process
		k = n/(Nx*Ny); j = (n-Nx*Ny*k)/Nx; i = n-Nx*Ny*k-Nx*j;
		// if (rank ==0) printf("@ Get 3D indices from the send process: i=%d, j=%d, k=%d\n",i,j,k);

		// Streaming for the non-local distribution
		i += Cqx; j += Cqy; k += Cqz;
		// if (rank == 0) printf("@ Streaming for the non-local distribution: i=%d, j=%d, k=%d\n",i,j,k);

		// Compute 1D index for the neighbor and save
		nn = k*Nx*Ny+j*Nx+i;
		// if (rank == 0) printf("@ rank:%d: neighbor=%d\n",rank,nn);
		ReturnDist[idx] = nn;
	}

	// Return updated version to the device
	ScaLBL_CopyToDevice(&d3q19_recvlist[start], ReturnDist, count*sizeof(int));

	// clean up the work arrays
	delete [] ReturnDist;
}

int ScaLBL_Communicator::MemoryOptimizedLayoutAA(IntArray &Map, int *neighborList, char *id, int Np){
	/*
	 * Generate a memory optimized layout
	 *   id[n] == 0 implies that site n should be ignored (treat as a mask)
	 *   Map(i,j,k) = idx  <- this is the index for the memory optimized layout
	 *   neighborList(idx) <-stores the neighbors for the D3Q19 model
	 *   note that the number of communications remains the same
	 *   the index in the Send and Recv lists is also updated
	 *   this means that the commuincations are no longer valid for regular data structures
	 */
	int idx,i,j,k,n;

	// Check that Map has size matching sub-domain
	if (Map.size(0) != Nx)
		ERROR("ScaLBL_Communicator::MemoryOptimizedLayout: Map array dimensions do not match! \n");

	// Initialize Map
	for (k=0;k<Nz;k++){
		for (j=0;j<Ny;j++){
			for (i=0;i<Nx;i++){
				Map(i,j,k) = -2;
			}
		}
	}

	//printf("Exterior... \n");

	// ********* Exterior **********
	// Step 1/2: Index the outer walls of the grid only
	idx=0;	next=0;
	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				// domain interior
				Map(i,j,k) = -1;
				// Local index
				n = k*Nx*Ny+j*Nx+i;
				if (id[n] > 0){
					// Counts for the six faces
					if (i==1)       Map(n)=idx++;
					else if (j==1)  Map(n)=idx++;
					else if (k==1)  Map(n)=idx++;
					else if (i==Nx-2)  Map(n)=idx++;
					else if (j==Ny-2)  Map(n)=idx++;
					else if (k==Nz-2)  Map(n)=idx++;
				}
			}
		}
	}
	next=idx;
	
	//printf("Interior... \n");
	
	// ********* Interior **********
	// align the next read
	first_interior=(next/16 + 1)*16;
	idx = first_interior;
	// Step 2/2: Next loop over the domain interior in block-cyclic fashion
	for (k=2; k<Nz-2; k++){
		for (j=2; j<Ny-2; j++){
			for (i=2; i<Nx-2; i++){
				// Local index (regular layout)
				n = k*Nx*Ny + j*Nx + i;
				if (id[n] > 0 ){
					Map(n) = idx++;
					//neighborList[idx++] = n; // index of self in regular layout
				}
			}
		}
	}
	last_interior=idx;
	
	Np = (last_interior/16 + 1)*16;
	//printf("    Np=%i \n",Np);
		
	// Now use Map to determine the neighbors for each lattice direction
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				n=k*Nx*Ny+j*Nx+i;
				idx=Map(i,j,k);
				if (idx > Np) printf("ScaLBL_Communicator::MemoryOptimizedLayout: Map(%i,%i,%i) = %i > %i \n",i,j,k,Map(i,j,k),Np);
				else if (!(idx<0)){
					// store the idx associated with each neighbor
					// store idx for self if neighbor is in solid or out of domain
					//D3Q19 = {{1,0,0},{-1,0,0}
					//         {0,1,0},{0,-1,0}
					//         {0,0,1},{0,0,-1},
					//	       {1,1,0},{-1,-1,0},
					//         {1,-1,0},{-1,1,0},
					//         {1,0,1},{-1,0,-1},
					//         {1,0,-1},{-1,0,1},
					//	       {0,1,1},{0,-1,-1},
					//         {0,1,-1},{0,-1,1}};
					int neighbor;    // cycle through the neighbors of lattice site idx
					neighbor=Map(i-1,j,k);
					if (neighbor<0)	   neighborList[idx]=idx + 2*Np;
					else		       neighborList[idx]=neighbor + 1*Np;

					neighbor=Map(i+1,j,k);
					if (neighbor<0)	   neighborList[Np+idx] = idx + 1*Np;
					else			   neighborList[Np+idx]= neighbor + 2*Np;

					neighbor=Map(i,j-1,k);
					if (neighbor<0)	   neighborList[2*Np+idx]=idx + 4*Np;
					else	       	   neighborList[2*Np+idx]=neighbor + 3*Np;

					neighbor=Map(i,j+1,k);
					if (neighbor<0)	   neighborList[3*Np+idx]=idx + 3*Np;
					else	      	   neighborList[3*Np+idx]=neighbor + 4*Np;

					neighbor=Map(i,j,k-1);
					if (neighbor<0)	   neighborList[4*Np+idx]=idx + 6*Np;
					else	       	   neighborList[4*Np+idx]=neighbor + 5*Np;

					neighbor=Map(i,j,k+1);
					if (neighbor<0)	   neighborList[5*Np+idx]=idx + 5*Np;
					else			   neighborList[5*Np+idx]=neighbor + 6*Np;

					neighbor=Map(i-1,j-1,k);
					if (neighbor<0)	   neighborList[6*Np+idx]=idx + 8*Np;
					else	      	   neighborList[6*Np+idx]=neighbor + 7*Np;

					neighbor=Map(i+1,j+1,k);
					if (neighbor<0)	   neighborList[7*Np+idx]=idx + 7*Np;
					else		       neighborList[7*Np+idx]=neighbor+8*Np;

					neighbor=Map(i-1,j+1,k);
					if (neighbor<0)	   neighborList[8*Np+idx]=idx + 10*Np;
					else		       neighborList[8*Np+idx]=neighbor + 9*Np;

					neighbor=Map(i+1,j-1,k);
					if (neighbor<0)	   neighborList[9*Np+idx]=idx + 9*Np;
					else		       neighborList[9*Np+idx]=neighbor + 10*Np;

					neighbor=Map(i-1,j,k-1);
					if (neighbor<0)    neighborList[10*Np+idx]=idx + 12*Np;
					else	    	   neighborList[10*Np+idx]=neighbor + 11*Np;

					neighbor=Map(i+1,j,k+1);
					if (neighbor<0)	   neighborList[11*Np+idx]=idx + 11*Np;
					else		       neighborList[11*Np+idx]=neighbor + 12*Np;

					neighbor=Map(i-1,j,k+1);
					if (neighbor<0)    neighborList[12*Np+idx]=idx + 14*Np;
					else	  	       neighborList[12*Np+idx]=neighbor + 13*Np;

					neighbor=Map(i+1,j,k-1);
					if (neighbor<0)    neighborList[13*Np+idx]=idx + 13*Np;
					else		       neighborList[13*Np+idx]=neighbor + 14*Np;

					neighbor=Map(i,j-1,k-1);
					if (neighbor<0)	   neighborList[14*Np+idx]=idx + 16*Np;
					else		       neighborList[14*Np+idx]=neighbor + 15*Np;

					neighbor=Map(i,j+1,k+1);
					if (neighbor<0)	   neighborList[15*Np+idx]=idx + 15*Np;
					else		       neighborList[15*Np+idx]=neighbor + 16*Np;

					neighbor=Map(i,j-1,k+1);
					if (neighbor<0)	   neighborList[16*Np+idx]=idx + 18*Np;
					else		       neighborList[16*Np+idx]=neighbor + 17*Np;

					neighbor=Map(i,j+1,k-1);
					if (neighbor<0)	   neighborList[17*Np+idx]=idx + 17*Np;
					else		       neighborList[17*Np+idx]=neighbor + 18*Np;
				}
			}
		}
	}

	//for (idx=0; idx<Np; idx++)	printf("%i: %i %i\n", idx, neighborList[Np],  neighborList[Np+idx]);
	//.......................................................................
	// Now map through  SendList and RecvList to update indices
	// First loop over the send lists

	int *TempBuffer;
	TempBuffer = new int [5*RecvCount];

	//.......................................................................
	// Re-index the send lists
	ScaLBL_CopyToHost(TempBuffer,dvcSendList_x,sendCount_x*sizeof(int));
	for (i=0; i<sendCount_x; i++){
		n = TempBuffer[i];
		//if (rank==0) printf("s: n=%d ",n);
		idx=Map(n);
		//if (rank == 0) printf("s: mapped n=%d\n",idx);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_x,TempBuffer,sendCount_x*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_y,sendCount_y*sizeof(int));
	for (i=0; i<sendCount_y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_y,TempBuffer,sendCount_y*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_z,sendCount_z*sizeof(int));
	for (i=0; i<sendCount_z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_z,TempBuffer,sendCount_z*sizeof(int));


	ScaLBL_CopyToHost(TempBuffer,dvcSendList_X,sendCount_X*sizeof(int));
	for (i=0; i<sendCount_X; i++){
		n = TempBuffer[i];
		//if (rank==0) printf("r: n=%d ",n);
		idx=Map(n);
		//if (rank == 0) printf("r: mapped n=%d\n",idx);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_X,TempBuffer,sendCount_X*sizeof(int));


	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Y,sendCount_Y*sizeof(int));
	for (i=0; i<sendCount_Y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Y,TempBuffer,sendCount_Y*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Z,sendCount_Z*sizeof(int));
	for (i=0; i<sendCount_Z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Z,TempBuffer,sendCount_Z*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xy,sendCount_xy*sizeof(int));
	for (i=0; i<sendCount_xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xy,TempBuffer,sendCount_xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xY,sendCount_xY*sizeof(int));
	for (i=0; i<sendCount_xY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xY,TempBuffer,sendCount_xY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Xy,sendCount_Xy*sizeof(int));
	for (i=0; i<sendCount_Xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Xy,TempBuffer,sendCount_Xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_XY,sendCount_XY*sizeof(int));
	for (i=0; i<sendCount_XY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_XY,TempBuffer,sendCount_XY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xz,sendCount_xz*sizeof(int));
	for (i=0; i<sendCount_xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xz,TempBuffer,sendCount_xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_xZ,sendCount_xZ*sizeof(int));
	for (i=0; i<sendCount_xZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_xZ,TempBuffer,sendCount_xZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Xz,sendCount_Xz*sizeof(int));
	for (i=0; i<sendCount_Xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Xz,TempBuffer,sendCount_Xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_XZ,sendCount_XZ*sizeof(int));
	for (i=0; i<sendCount_XZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_XZ,TempBuffer,sendCount_XZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_yz,sendCount_yz*sizeof(int));
	for (i=0; i<sendCount_yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_yz,TempBuffer,sendCount_yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_Yz,sendCount_Yz*sizeof(int));
	for (i=0; i<sendCount_Yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_Yz,TempBuffer,sendCount_Yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_yZ,sendCount_yZ*sizeof(int));
	for (i=0; i<sendCount_yZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_yZ,TempBuffer,sendCount_yZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcSendList_YZ,sendCount_YZ*sizeof(int));
	for (i=0; i<sendCount_YZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcSendList_YZ,TempBuffer,sendCount_YZ*sizeof(int));
	//.......................................................................
	// Re-index the recieve lists for the D3Q19 distributions
	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_x,5*recvCount_x*sizeof(int));
	for (i=0; i<5*recvCount_x; i++){
		n = TempBuffer[i];
		//if (rank==0) printf("r: n=%d ",n);
		idx=Map(n);
		//if (rank == 0) printf("r: mapped n=%d\n",idx);
		TempBuffer[i]=idx;


	}
	ScaLBL_CopyToDevice(dvcRecvDist_x,TempBuffer,5*recvCount_x*sizeof(int));


	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_y,5*recvCount_y*sizeof(int));
	for (i=0; i<5*recvCount_y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}

	ScaLBL_CopyToDevice(dvcRecvDist_y,TempBuffer,5*recvCount_y*sizeof(int));


	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_z,5*recvCount_z*sizeof(int));
	for (i=0; i<5*recvCount_z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_z,TempBuffer,5*recvCount_z*sizeof(int));


	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_X,5*recvCount_X*sizeof(int));
	for (i=0; i<5*recvCount_X; i++){
		n = TempBuffer[i];
		//if (rank==0) printf("r: n=%d ",n);
		idx=Map(n);
		//if (rank == 0) printf("r: mapped n=%d\n",idx);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_X,TempBuffer,5*recvCount_X*sizeof(int));


	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Y,5*recvCount_Y*sizeof(int));
	for (i=0; i<5*recvCount_Y; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Y,TempBuffer,5*recvCount_Y*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Z,5*recvCount_Z*sizeof(int));
	for (i=0; i<5*recvCount_Z; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Z,TempBuffer,5*recvCount_Z*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xy,recvCount_xy*sizeof(int));
	for (i=0; i<recvCount_xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xy,TempBuffer,recvCount_xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xY,recvCount_xY*sizeof(int));
	for (i=0; i<recvCount_xY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xY,TempBuffer,recvCount_xY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Xy,recvCount_Xy*sizeof(int));
	for (i=0; i<recvCount_Xy; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Xy,TempBuffer,recvCount_Xy*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_XY,recvCount_XY*sizeof(int));
	for (i=0; i<recvCount_XY; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_XY,TempBuffer,recvCount_XY*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xz,recvCount_xz*sizeof(int));
	for (i=0; i<recvCount_xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xz,TempBuffer,recvCount_xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_xZ,recvCount_xZ*sizeof(int));
	for (i=0; i<recvCount_xZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_xZ,TempBuffer,recvCount_xZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Xz,recvCount_Xz*sizeof(int));
	for (i=0; i<recvCount_Xz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Xz,TempBuffer,recvCount_Xz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_XZ,recvCount_XZ*sizeof(int));
	for (i=0; i<recvCount_XZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_XZ,TempBuffer,recvCount_XZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_yz,recvCount_yz*sizeof(int));
	for (i=0; i<recvCount_yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_yz,TempBuffer,recvCount_yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_Yz,recvCount_Yz*sizeof(int));
	for (i=0; i<recvCount_Yz; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_Yz,TempBuffer,recvCount_Yz*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_yZ,recvCount_yZ*sizeof(int));
	for (i=0; i<recvCount_yZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_yZ,TempBuffer,recvCount_yZ*sizeof(int));

	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_YZ,recvCount_YZ*sizeof(int));
	for (i=0; i<recvCount_YZ; i++){
		n = TempBuffer[i];
		idx=Map(n);
		TempBuffer[i]=idx;
	}
	ScaLBL_CopyToDevice(dvcRecvDist_YZ,TempBuffer,recvCount_YZ*sizeof(int));
	//.......................................................................


	// Reset the value of N to match the dense structure
	N = Np;
	
	// Clean up
	delete [] TempBuffer;
	return(Np);
}

void ScaLBL_Communicator::SendD3Q19AA(double *dist){

	// NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
	if (Lock==true){
		ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	// assign tag of 19 to D3Q19 communication
	sendtag = recvtag = 19;
	ScaLBL_DeviceBarrier();
	// Pack the distributions
	//...Packing for x face(2,8,10,12,14)................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_x,0,sendCount_x,sendbuf_x,dist,N);
	ScaLBL_D3Q19_Pack(8,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,dist,N);
	ScaLBL_D3Q19_Pack(10,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,dist,N);
	ScaLBL_D3Q19_Pack(12,dvcSendList_x,3*sendCount_x,sendCount_x,sendbuf_x,dist,N);
	ScaLBL_D3Q19_Pack(14,dvcSendList_x,4*sendCount_x,sendCount_x,sendbuf_x,dist,N);
	
	MPI_Isend(sendbuf_x, 5*sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_SCALBL,&req1[0]);
	MPI_Irecv(recvbuf_X, 5*recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_SCALBL,&req2[0]);
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,0,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(9,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(11,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(13,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,dist,N);
	
	MPI_Isend(sendbuf_X, 5*sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_SCALBL,&req1[1]);
	MPI_Irecv(recvbuf_x, 5*recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_SCALBL,&req2[1]);
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,0,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(8,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(9,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(16,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(18,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,dist,N);
	
	MPI_Isend(sendbuf_y, 5*sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_SCALBL,&req1[2]);
	MPI_Irecv(recvbuf_Y, 5*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_SCALBL,&req2[2]);
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(10,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(15,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(17,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	
	MPI_Isend(sendbuf_Y, 5*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_SCALBL,&req1[3]);
	MPI_Irecv(recvbuf_y, 5*recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_SCALBL,&req2[3]);
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,0,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(12,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(13,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(16,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(17,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,dist,N);
	
	MPI_Isend(sendbuf_z, 5*sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_SCALBL,&req1[4]);
	MPI_Irecv(recvbuf_Z, 5*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_SCALBL,&req2[4]);
	
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(11,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(14,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(15,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(18,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	
	MPI_Isend(sendbuf_Z, 5*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_SCALBL,&req1[5]);
	MPI_Irecv(recvbuf_z, 5*recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_SCALBL,&req2[5]);
	
	//...Pack the xy edge (8)................................
	ScaLBL_D3Q19_Pack(8,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,dist,N);
	MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,MPI_COMM_SCALBL,&req1[6]);
	MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,MPI_COMM_SCALBL,&req2[6]);
	//...Pack the Xy edge (9)................................
	ScaLBL_D3Q19_Pack(9,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,dist,N);
	MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,MPI_COMM_SCALBL,&req1[8]);
	MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,MPI_COMM_SCALBL,&req2[8]);
	//...Pack the xY edge (10)................................
	ScaLBL_D3Q19_Pack(10,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,dist,N);
	MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,MPI_COMM_SCALBL,&req1[9]);
	MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,MPI_COMM_SCALBL,&req2[9]);
	//...Pack the XY edge (7)................................
	ScaLBL_D3Q19_Pack(7,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,dist,N);
	MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,MPI_COMM_SCALBL,&req1[7]);
	MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,MPI_COMM_SCALBL,&req2[7]);
	//...Pack the xz edge (12)................................
	ScaLBL_D3Q19_Pack(12,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,dist,N);
	MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,MPI_COMM_SCALBL,&req1[10]);
	MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,MPI_COMM_SCALBL,&req2[10]);
	//...Pack the xZ edge (14)................................
	ScaLBL_D3Q19_Pack(14,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,dist,N);
	MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,MPI_COMM_SCALBL,&req1[13]);
	MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,MPI_COMM_SCALBL,&req2[13]);
	//...Pack the Xz edge (13)................................
	ScaLBL_D3Q19_Pack(13,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,dist,N);
	MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,MPI_COMM_SCALBL,&req1[12]);
	MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,MPI_COMM_SCALBL,&req2[12]);
	//...Pack the XZ edge (11)................................
	ScaLBL_D3Q19_Pack(11,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,dist,N);
	MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,MPI_COMM_SCALBL,&req1[11]);
	MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,MPI_COMM_SCALBL,&req2[11]);
	//...Pack the yz edge (16)................................
	ScaLBL_D3Q19_Pack(16,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,dist,N);
	MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,MPI_COMM_SCALBL,&req1[14]);
	MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,MPI_COMM_SCALBL,&req2[14]);
	//...Pack the yZ edge (18)................................
	ScaLBL_D3Q19_Pack(18,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,dist,N);
	MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,MPI_COMM_SCALBL,&req1[17]);
	MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,MPI_COMM_SCALBL,&req2[17]);
	//...Pack the Yz edge (17)................................
	ScaLBL_D3Q19_Pack(17,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,dist,N);
	MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,MPI_COMM_SCALBL,&req1[16]);
	MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,MPI_COMM_SCALBL,&req2[16]);
	//...Pack the YZ edge (15)................................
	ScaLBL_D3Q19_Pack(15,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,dist,N);
	MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,MPI_COMM_SCALBL,&req1[15]);
	MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,MPI_COMM_SCALBL,&req2[15]);
	//...................................................................................

}

void ScaLBL_Communicator::RecvD3Q19AA(double *dist){

	// NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	ScaLBL_DeviceBarrier();

	//...................................................................................
	// NOTE: AA Routine writes to opposite 
	// Unpack the distributions on the device
	//...................................................................................
	//...Unpacking for x face(2,8,10,12,14)................................
	ScaLBL_D3Q19_Unpack(2,dvcRecvDist_x,0,recvCount_x,recvbuf_x,dist,N);
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_x,recvCount_x,recvCount_x,recvbuf_x,dist,N);
	ScaLBL_D3Q19_Unpack(10,dvcRecvDist_x,2*recvCount_x,recvCount_x,recvbuf_x,dist,N);
	ScaLBL_D3Q19_Unpack(12,dvcRecvDist_x,3*recvCount_x,recvCount_x,recvbuf_x,dist,N);
	ScaLBL_D3Q19_Unpack(14,dvcRecvDist_x,4*recvCount_x,recvCount_x,recvbuf_x,dist,N);
	//...................................................................................
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q19_Unpack(1,dvcRecvDist_X,0,recvCount_X,recvbuf_X,dist,N);
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_X,recvCount_X,recvCount_X,recvbuf_X,dist,N);
	ScaLBL_D3Q19_Unpack(9,dvcRecvDist_X,2*recvCount_X,recvCount_X,recvbuf_X,dist,N);
	ScaLBL_D3Q19_Unpack(11,dvcRecvDist_X,3*recvCount_X,recvCount_X,recvbuf_X,dist,N);
	ScaLBL_D3Q19_Unpack(13,dvcRecvDist_X,4*recvCount_X,recvCount_X,recvbuf_X,dist,N);
	//...................................................................................
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q19_Unpack(4,dvcRecvDist_y,0,recvCount_y,recvbuf_y,dist,N);
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_y,recvCount_y,recvCount_y,recvbuf_y,dist,N);
	ScaLBL_D3Q19_Unpack(9,dvcRecvDist_y,2*recvCount_y,recvCount_y,recvbuf_y,dist,N);
	ScaLBL_D3Q19_Unpack(16,dvcRecvDist_y,3*recvCount_y,recvCount_y,recvbuf_y,dist,N);
	ScaLBL_D3Q19_Unpack(18,dvcRecvDist_y,4*recvCount_y,recvCount_y,recvbuf_y,dist,N);
	//...................................................................................
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q19_Unpack(3,dvcRecvDist_Y,0,recvCount_Y,recvbuf_Y,dist,N);
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_Y,recvCount_Y,recvCount_Y,recvbuf_Y,dist,N);
	ScaLBL_D3Q19_Unpack(10,dvcRecvDist_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,dist,N);
	ScaLBL_D3Q19_Unpack(15,dvcRecvDist_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,dist,N);
	ScaLBL_D3Q19_Unpack(17,dvcRecvDist_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,dist,N);
	//...................................................................................
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_D3Q19_Unpack(6,dvcRecvDist_z,0,recvCount_z,recvbuf_z,dist,N);
	ScaLBL_D3Q19_Unpack(12,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,dist,N);
	ScaLBL_D3Q19_Unpack(13,dvcRecvDist_z,2*recvCount_z,recvCount_z,recvbuf_z,dist,N);
	ScaLBL_D3Q19_Unpack(16,dvcRecvDist_z,3*recvCount_z,recvCount_z,recvbuf_z,dist,N);
	ScaLBL_D3Q19_Unpack(17,dvcRecvDist_z,4*recvCount_z,recvCount_z,recvbuf_z,dist,N);
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_D3Q19_Unpack(5,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,dist,N);
	ScaLBL_D3Q19_Unpack(11,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,dist,N);
	ScaLBL_D3Q19_Unpack(14,dvcRecvDist_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,dist,N);
	ScaLBL_D3Q19_Unpack(15,dvcRecvDist_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,dist,N);
	ScaLBL_D3Q19_Unpack(18,dvcRecvDist_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,dist,N);
	//..................................................................................
	//...Pack the xy edge (8)................................
	ScaLBL_D3Q19_Unpack(8,dvcRecvDist_xy,0,recvCount_xy,recvbuf_xy,dist,N);
	//...Pack the Xy edge (9)................................
	ScaLBL_D3Q19_Unpack(9,dvcRecvDist_Xy,0,recvCount_Xy,recvbuf_Xy,dist,N);
	//...Pack the xY edge (10)................................
	ScaLBL_D3Q19_Unpack(10,dvcRecvDist_xY,0,recvCount_xY,recvbuf_xY,dist,N);
	//...Pack the XY edge (7)................................
	ScaLBL_D3Q19_Unpack(7,dvcRecvDist_XY,0,recvCount_XY,recvbuf_XY,dist,N);
	//...Pack the xz edge (12)................................
	ScaLBL_D3Q19_Unpack(12,dvcRecvDist_xz,0,recvCount_xz,recvbuf_xz,dist,N);
	//...Pack the xZ edge (14)................................
	ScaLBL_D3Q19_Unpack(14,dvcRecvDist_xZ,0,recvCount_xZ,recvbuf_xZ,dist,N);
	//...Pack the Xz edge (13)................................
	ScaLBL_D3Q19_Unpack(13,dvcRecvDist_Xz,0,recvCount_Xz,recvbuf_Xz,dist,N);
	//...Pack the XZ edge (11)................................
	ScaLBL_D3Q19_Unpack(11,dvcRecvDist_XZ,0,recvCount_XZ,recvbuf_XZ,dist,N);
	//...Pack the yz edge (16)................................
	ScaLBL_D3Q19_Unpack(16,dvcRecvDist_yz,0,recvCount_yz,recvbuf_yz,dist,N);
	//...Pack the yZ edge (18)................................
	ScaLBL_D3Q19_Unpack(18,dvcRecvDist_yZ,0,recvCount_yZ,recvbuf_yZ,dist,N);
	//...Pack the Yz edge (17)................................
	ScaLBL_D3Q19_Unpack(17,dvcRecvDist_Yz,0,recvCount_Yz,recvbuf_Yz,dist,N);
	//...Pack the YZ edge (15)................................
	ScaLBL_D3Q19_Unpack(15,dvcRecvDist_YZ,0,recvCount_YZ,recvbuf_YZ,dist,N);
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................

}

void ScaLBL_Communicator::RecvGrad(double *phi, double *grad){

	// Recieves halo and incorporates into D3Q19 based stencil gradient computation
	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	ScaLBL_DeviceBarrier();

	//...................................................................................
	// Unpack the gradributions on the device
	//...................................................................................
	//...Unpacking for x face(2,8,10,12,14)................................
	ScaLBL_Gradient_Unpack(1.0,-1,0,0,dvcRecvDist_x,0,recvCount_x,recvbuf_x,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,-1,-1,0,dvcRecvDist_x,recvCount_x,recvCount_x,recvbuf_x,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,-1,1,0,dvcRecvDist_x,2*recvCount_x,recvCount_x,recvbuf_x,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,-1,0,1,dvcRecvDist_x,4*recvCount_x,recvCount_x,recvbuf_x,phi,grad,N);
	//...................................................................................
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_Gradient_Unpack(1.0,1,0,0,dvcRecvDist_X,0,recvCount_X,recvbuf_X,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,1,0,dvcRecvDist_X,recvCount_X,recvCount_X,recvbuf_X,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,-1,0,dvcRecvDist_X,2*recvCount_X,recvCount_X,recvbuf_X,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,0,1,dvcRecvDist_X,3*recvCount_X,recvCount_X,recvbuf_X,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,0,-1,dvcRecvDist_X,4*recvCount_X,recvCount_X,recvbuf_X,phi,grad,N);
	//...................................................................................
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_Gradient_Unpack(1.0,0,-1,0,dvcRecvDist_y,0,recvCount_y,recvbuf_y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,-1,-1,0,dvcRecvDist_y,recvCount_y,recvCount_y,recvbuf_y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,-1,0,dvcRecvDist_y,2*recvCount_y,recvCount_y,recvbuf_y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,-1,-1,dvcRecvDist_y,3*recvCount_y,recvCount_y,recvbuf_y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,-1,1,dvcRecvDist_y,4*recvCount_y,recvCount_y,recvbuf_y,phi,grad,N);
	//...................................................................................
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_Gradient_Unpack(1.0,0,1,0,dvcRecvDist_Y,0,recvCount_Y,recvbuf_Y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,1,0,dvcRecvDist_Y,recvCount_Y,recvCount_Y,recvbuf_Y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,-1,1,0,dvcRecvDist_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,1,1,dvcRecvDist_Y,3*recvCount_Y,recvCount_Y,recvbuf_Y,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,1,-1,dvcRecvDist_Y,4*recvCount_Y,recvCount_Y,recvbuf_Y,phi,grad,N);
	//...................................................................................
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_Gradient_Unpack(1.0,0,0,-1,dvcRecvDist_z,0,recvCount_z,recvbuf_z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,-1,0,-1,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,0,-1,dvcRecvDist_z,2*recvCount_z,recvCount_z,recvbuf_z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,-1,-1,dvcRecvDist_z,3*recvCount_z,recvCount_z,recvbuf_z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,1,-1,dvcRecvDist_z,4*recvCount_z,recvCount_z,recvbuf_z,phi,grad,N);
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_Gradient_Unpack(1.0,0,0,1,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,1,0,1,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,-1,0,1,dvcRecvDist_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,1,1,dvcRecvDist_Z,3*recvCount_Z,recvCount_Z,recvbuf_Z,phi,grad,N);
	ScaLBL_Gradient_Unpack(0.5,0,-1,1,dvcRecvDist_Z,4*recvCount_Z,recvCount_Z,recvbuf_Z,phi,grad,N);
	//..................................................................................
	//...Pack the xy edge (8)................................
	ScaLBL_Gradient_Unpack(0.5,-1,-1,0,dvcRecvDist_xy,0,recvCount_xy,recvbuf_xy,phi,grad,N);
	//...Pack the Xy edge (9)................................
	ScaLBL_Gradient_Unpack(0.5,1,-1,0,dvcRecvDist_Xy,0,recvCount_Xy,recvbuf_Xy,phi,grad,N);
	//...Pack the xY edge (10)................................
	ScaLBL_Gradient_Unpack(0.5,-1,1,0,dvcRecvDist_xY,0,recvCount_xY,recvbuf_xY,phi,grad,N);
	//...Pack the XY edge (7)................................
	ScaLBL_Gradient_Unpack(0.5,1,1,0,dvcRecvDist_XY,0,recvCount_XY,recvbuf_XY,phi,grad,N);
	//...Pack the xz edge (12)................................
	ScaLBL_Gradient_Unpack(0.5,-1,0,-1,dvcRecvDist_xz,0,recvCount_xz,recvbuf_xz,phi,grad,N);
	//...Pack the xZ edge (14)................................
	ScaLBL_Gradient_Unpack(0.5,-1,0,1,dvcRecvDist_xZ,0,recvCount_xZ,recvbuf_xZ,phi,grad,N);
	//...Pack the Xz edge (13)................................
	ScaLBL_Gradient_Unpack(0.5,1,0,-1,dvcRecvDist_Xz,0,recvCount_Xz,recvbuf_Xz,phi,grad,N);
	//...Pack the XZ edge (11)................................
	ScaLBL_Gradient_Unpack(0.5,1,0,1,dvcRecvDist_XZ,0,recvCount_XZ,recvbuf_XZ,phi,grad,N);
	//...Pack the yz edge (16)................................
	ScaLBL_Gradient_Unpack(0.5,0,-1,-1,dvcRecvDist_yz,0,recvCount_yz,recvbuf_yz,phi,grad,N);
	//...Pack the yZ edge (18)................................
	ScaLBL_Gradient_Unpack(0.5,0,-1,1,dvcRecvDist_yZ,0,recvCount_yZ,recvbuf_yZ,phi,grad,N);
	//...Pack the Yz edge (17)................................
	ScaLBL_Gradient_Unpack(0.5,0,1,-1,dvcRecvDist_Yz,0,recvCount_Yz,recvbuf_Yz,phi,grad,N);
	//...Pack the YZ edge (15)................................
	ScaLBL_Gradient_Unpack(0.5,0,1,1,dvcRecvDist_YZ,0,recvCount_YZ,recvbuf_YZ,phi,grad,N);
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................

}

void ScaLBL_Communicator::BiSendD3Q7AA(double *Aq, double *Bq){

	// NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
	if (Lock==true){
		ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	// assign tag of 19 to D3Q19 communication
	sendtag = recvtag = 14;
	ScaLBL_DeviceBarrier();
	// Pack the distributions
	//...Packing for x face(2,8,10,12,14)................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_x,0,sendCount_x,sendbuf_x,Aq,N);
	ScaLBL_D3Q19_Pack(2,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,Bq,N);

	MPI_Isend(sendbuf_x, 2*sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_SCALBL,&req1[0]);
	MPI_Irecv(recvbuf_X, 2*recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_SCALBL,&req2[0]);
	
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,0,sendCount_X,sendbuf_X,Aq,N);
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,Bq,N);
	
	MPI_Isend(sendbuf_X, 2*sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_SCALBL,&req1[1]);
	MPI_Irecv(recvbuf_x, 2*recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_SCALBL,&req2[1]);

	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,0,sendCount_y,sendbuf_y,Aq,N);
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,Bq,N);

	MPI_Isend(sendbuf_y, 2*sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_SCALBL,&req1[2]);
	MPI_Irecv(recvbuf_Y, 2*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_SCALBL,&req2[2]);
	
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,Aq,N);
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,Bq,N);

	MPI_Isend(sendbuf_Y, 2*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_SCALBL,&req1[3]);
	MPI_Irecv(recvbuf_y, 2*recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_SCALBL,&req2[3]);
	
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,0,sendCount_z,sendbuf_z,Aq,N);
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,Bq,N);
	
	MPI_Isend(sendbuf_z, 2*sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_SCALBL,&req1[4]);
	MPI_Irecv(recvbuf_Z, 2*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_SCALBL,&req2[4]);
	
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,Aq,N);
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,Bq,N);

	//...................................................................................
	// Send all the distributions
	MPI_Isend(sendbuf_Z, 2*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_SCALBL,&req1[5]);
	MPI_Irecv(recvbuf_z, 2*recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_SCALBL,&req2[5]);

}


void ScaLBL_Communicator::BiRecvD3Q7AA(double *Aq, double *Bq){

	// NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_Waitall(6,req1,stat1);
	MPI_Waitall(6,req2,stat2);
	ScaLBL_DeviceBarrier();

	//...................................................................................
	// NOTE: AA Routine writes to opposite
	// Unpack the distributions on the device
	//...................................................................................
	//...Unpacking for x face(2,8,10,12,14)................................
	ScaLBL_D3Q7_Unpack(2,dvcRecvDist_x,0,recvCount_x,recvbuf_x,Aq,N);
	ScaLBL_D3Q7_Unpack(2,dvcRecvDist_x,recvCount_x,recvCount_x,recvbuf_x,Bq,N);
	//...................................................................................
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q7_Unpack(1,dvcRecvDist_X,0,recvCount_X,recvbuf_X,Aq,N);
	ScaLBL_D3Q7_Unpack(1,dvcRecvDist_X,recvCount_X,recvCount_X,recvbuf_X,Bq,N);
	//...................................................................................
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q7_Unpack(4,dvcRecvDist_y,0,recvCount_y,recvbuf_y,Aq,N);
	ScaLBL_D3Q7_Unpack(4,dvcRecvDist_y,recvCount_y,recvCount_y,recvbuf_y,Bq,N);
	//...................................................................................
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q7_Unpack(3,dvcRecvDist_Y,0,recvCount_Y,recvbuf_Y,Aq,N);
	ScaLBL_D3Q7_Unpack(3,dvcRecvDist_Y,recvCount_Y,recvCount_Y,recvbuf_Y,Bq,N);
	//...................................................................................
	
	if (BoundaryCondition > 0 && kproc == 0){
		// don't unpack little z
		//...Packing for Z face(5,11,14,15,18)................................
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,Aq,N);
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,Bq,N);
	}
	else if (BoundaryCondition > 0 && kproc == nprocz-1){
		// don't unpack big z
		//...Packing for z face(6,12,13,16,17)................................
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,0,recvCount_z,recvbuf_z,Aq,N);
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,Bq,N);
	}
	else {
		//...Packing for z face(6,12,13,16,17)................................
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,0,recvCount_z,recvbuf_z,Aq,N);
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,Bq,N);
		//...Packing for Z face(5,11,14,15,18)................................
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,Aq,N);
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,Bq,N);
	}
	
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................

}

void ScaLBL_Communicator::TriSendD3Q7AA(double *Aq, double *Bq, double *Cq){

	// NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
	if (Lock==true){
		ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	// assign tag of 19 to D3Q19 communication
	sendtag = recvtag = 15;
	ScaLBL_DeviceBarrier();
	// Pack the distributions
	//...Packing for x face(2,8,10,12,14)................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_x,0,sendCount_x,sendbuf_x,Aq,N);
	ScaLBL_D3Q19_Pack(2,dvcSendList_x,sendCount_x,sendCount_x,sendbuf_x,Bq,N);
	ScaLBL_D3Q19_Pack(2,dvcSendList_x,2*sendCount_x,sendCount_x,sendbuf_x,Cq,N);
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,0,sendCount_X,sendbuf_X,Aq,N);
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,Bq,N);
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,Cq,N);
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,0,sendCount_y,sendbuf_y,Aq,N);
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,Bq,N);
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,Cq,N);
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,Aq,N);
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,Bq,N);
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,Cq,N);
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,0,sendCount_z,sendbuf_z,Aq,N);
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,Bq,N);
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,Cq,N);
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,Aq,N);
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,Bq,N);
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,Cq,N);

	//...................................................................................
	// Send all the distributions
	MPI_Isend(sendbuf_x, 3*sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_SCALBL,&req1[0]);
	MPI_Irecv(recvbuf_X, 3*recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_SCALBL,&req2[0]);
	MPI_Isend(sendbuf_X, 3*sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_SCALBL,&req1[1]);
	MPI_Irecv(recvbuf_x, 3*recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_SCALBL,&req2[1]);
	MPI_Isend(sendbuf_y, 3*sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_SCALBL,&req1[2]);
	MPI_Irecv(recvbuf_Y, 3*recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_SCALBL,&req2[2]);
	MPI_Isend(sendbuf_Y, 3*sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_SCALBL,&req1[3]);
	MPI_Irecv(recvbuf_y, 3*recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_SCALBL,&req2[3]);
	MPI_Isend(sendbuf_z, 3*sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_SCALBL,&req1[4]);
	MPI_Irecv(recvbuf_Z, 3*recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_SCALBL,&req2[4]);
	MPI_Isend(sendbuf_Z, 3*sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_SCALBL,&req1[5]);
	MPI_Irecv(recvbuf_z, 3*recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_SCALBL,&req2[5]);

}


void ScaLBL_Communicator::TriRecvD3Q7AA(double *Aq, double *Bq, double *Cq){

	// NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_Waitall(6,req1,stat1);
	MPI_Waitall(6,req2,stat2);
	ScaLBL_DeviceBarrier();

	//...................................................................................
	// NOTE: AA Routine writes to opposite
	// Unpack the distributions on the device
	//...................................................................................
	//...Unpacking for x face(2,8,10,12,14)................................
	ScaLBL_D3Q7_Unpack(2,dvcRecvDist_x,0,recvCount_x,recvbuf_x,Aq,N);
	ScaLBL_D3Q7_Unpack(2,dvcRecvDist_x,recvCount_x,recvCount_x,recvbuf_x,Bq,N);
	ScaLBL_D3Q7_Unpack(2,dvcRecvDist_x,2*recvCount_x,recvCount_x,recvbuf_x,Cq,N);
	//...................................................................................
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q7_Unpack(1,dvcRecvDist_X,0,recvCount_X,recvbuf_X,Aq,N);
	ScaLBL_D3Q7_Unpack(1,dvcRecvDist_X,recvCount_X,recvCount_X,recvbuf_X,Bq,N);
	ScaLBL_D3Q7_Unpack(1,dvcRecvDist_X,2*recvCount_X,recvCount_X,recvbuf_X,Cq,N);
	//...................................................................................
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q7_Unpack(4,dvcRecvDist_y,0,recvCount_y,recvbuf_y,Aq,N);
	ScaLBL_D3Q7_Unpack(4,dvcRecvDist_y,recvCount_y,recvCount_y,recvbuf_y,Bq,N);
	ScaLBL_D3Q7_Unpack(4,dvcRecvDist_y,2*recvCount_y,recvCount_y,recvbuf_y,Cq,N);
	//...................................................................................
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q7_Unpack(3,dvcRecvDist_Y,0,recvCount_Y,recvbuf_Y,Aq,N);
	ScaLBL_D3Q7_Unpack(3,dvcRecvDist_Y,recvCount_Y,recvCount_Y,recvbuf_Y,Bq,N);
	ScaLBL_D3Q7_Unpack(3,dvcRecvDist_Y,2*recvCount_Y,recvCount_Y,recvbuf_Y,Cq,N);
	//...................................................................................
	
	if (BoundaryCondition > 0 && kproc == 0){
		// don't unpack little z
		//...Packing for Z face(5,11,14,15,18)................................
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,Aq,N);
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,Bq,N);
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,Cq,N);
	}
	else if (BoundaryCondition > 0 && kproc == nprocz-1){
		// don't unpack big z
		//...Packing for z face(6,12,13,16,17)................................
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,0,recvCount_z,recvbuf_z,Aq,N);
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,Bq,N);
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,2*recvCount_z,recvCount_z,recvbuf_z,Cq,N);
	}
	else {
		//...Packing for z face(6,12,13,16,17)................................
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,0,recvCount_z,recvbuf_z,Aq,N);
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,recvCount_z,recvCount_z,recvbuf_z,Bq,N);
		ScaLBL_D3Q7_Unpack(6,dvcRecvDist_z,2*recvCount_z,recvCount_z,recvbuf_z,Cq,N);
		//...Packing for Z face(5,11,14,15,18)................................
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,0,recvCount_Z,recvbuf_Z,Aq,N);
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,recvCount_Z,recvCount_Z,recvbuf_Z,Bq,N);
		ScaLBL_D3Q7_Unpack(5,dvcRecvDist_Z,2*recvCount_Z,recvCount_Z,recvbuf_Z,Cq,N);
	}
	
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................

}


void ScaLBL_Communicator::SendHalo(double *data){
	//...................................................................................
	if (Lock==true){
		ERROR("ScaLBL Error (SendHalo): ScaLBL_Communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	ScaLBL_DeviceBarrier();
	//...................................................................................
	sendtag = recvtag = 1;
	//...................................................................................
	ScaLBL_Scalar_Pack(dvcSendList_x, sendCount_x,sendbuf_x, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_y, sendCount_y,sendbuf_y, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_z, sendCount_z,sendbuf_z, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_X, sendCount_X,sendbuf_X, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Y, sendCount_Y,sendbuf_Y, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Z, sendCount_Z,sendbuf_Z, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xy, sendCount_xy,sendbuf_xy, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xY, sendCount_xY,sendbuf_xY, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Xy, sendCount_Xy,sendbuf_Xy, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_XY, sendCount_XY,sendbuf_XY, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xz, sendCount_xz,sendbuf_xz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_xZ, sendCount_xZ,sendbuf_xZ, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Xz, sendCount_Xz,sendbuf_Xz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_XZ, sendCount_XZ,sendbuf_XZ, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_yz, sendCount_yz,sendbuf_yz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_yZ, sendCount_yZ,sendbuf_yZ, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_Yz, sendCount_Yz,sendbuf_Yz, data, N);
	ScaLBL_Scalar_Pack(dvcSendList_YZ, sendCount_YZ,sendbuf_YZ, data, N);
	//...................................................................................
	// Send / Recv all the phase indcator field values
	//...................................................................................

	MPI_Isend(sendbuf_x, sendCount_x,MPI_DOUBLE,rank_x,sendtag,MPI_COMM_SCALBL,&req1[0]);
	MPI_Irecv(recvbuf_X, recvCount_X,MPI_DOUBLE,rank_X,recvtag,MPI_COMM_SCALBL,&req2[0]);
	MPI_Isend(sendbuf_X, sendCount_X,MPI_DOUBLE,rank_X,sendtag,MPI_COMM_SCALBL,&req1[1]);
	MPI_Irecv(recvbuf_x, recvCount_x,MPI_DOUBLE,rank_x,recvtag,MPI_COMM_SCALBL,&req2[1]);
	MPI_Isend(sendbuf_y, sendCount_y,MPI_DOUBLE,rank_y,sendtag,MPI_COMM_SCALBL,&req1[2]);
	MPI_Irecv(recvbuf_Y, recvCount_Y,MPI_DOUBLE,rank_Y,recvtag,MPI_COMM_SCALBL,&req2[2]);
	MPI_Isend(sendbuf_Y, sendCount_Y,MPI_DOUBLE,rank_Y,sendtag,MPI_COMM_SCALBL,&req1[3]);
	MPI_Irecv(recvbuf_y, recvCount_y,MPI_DOUBLE,rank_y,recvtag,MPI_COMM_SCALBL,&req2[3]);
	MPI_Isend(sendbuf_z, sendCount_z,MPI_DOUBLE,rank_z,sendtag,MPI_COMM_SCALBL,&req1[4]);
	MPI_Irecv(recvbuf_Z, recvCount_Z,MPI_DOUBLE,rank_Z,recvtag,MPI_COMM_SCALBL,&req2[4]);
	MPI_Isend(sendbuf_Z, sendCount_Z,MPI_DOUBLE,rank_Z,sendtag,MPI_COMM_SCALBL,&req1[5]);
	MPI_Irecv(recvbuf_z, recvCount_z,MPI_DOUBLE,rank_z,recvtag,MPI_COMM_SCALBL,&req2[5]);
	MPI_Isend(sendbuf_xy, sendCount_xy,MPI_DOUBLE,rank_xy,sendtag,MPI_COMM_SCALBL,&req1[6]);
	MPI_Irecv(recvbuf_XY, recvCount_XY,MPI_DOUBLE,rank_XY,recvtag,MPI_COMM_SCALBL,&req2[6]);
	MPI_Isend(sendbuf_XY, sendCount_XY,MPI_DOUBLE,rank_XY,sendtag,MPI_COMM_SCALBL,&req1[7]);
	MPI_Irecv(recvbuf_xy, recvCount_xy,MPI_DOUBLE,rank_xy,recvtag,MPI_COMM_SCALBL,&req2[7]);
	MPI_Isend(sendbuf_Xy, sendCount_Xy,MPI_DOUBLE,rank_Xy,sendtag,MPI_COMM_SCALBL,&req1[8]);
	MPI_Irecv(recvbuf_xY, recvCount_xY,MPI_DOUBLE,rank_xY,recvtag,MPI_COMM_SCALBL,&req2[8]);
	MPI_Isend(sendbuf_xY, sendCount_xY,MPI_DOUBLE,rank_xY,sendtag,MPI_COMM_SCALBL,&req1[9]);
	MPI_Irecv(recvbuf_Xy, recvCount_Xy,MPI_DOUBLE,rank_Xy,recvtag,MPI_COMM_SCALBL,&req2[9]);
	MPI_Isend(sendbuf_xz, sendCount_xz,MPI_DOUBLE,rank_xz,sendtag,MPI_COMM_SCALBL,&req1[10]);
	MPI_Irecv(recvbuf_XZ, recvCount_XZ,MPI_DOUBLE,rank_XZ,recvtag,MPI_COMM_SCALBL,&req2[10]);
	MPI_Isend(sendbuf_XZ, sendCount_XZ,MPI_DOUBLE,rank_XZ,sendtag,MPI_COMM_SCALBL,&req1[11]);
	MPI_Irecv(recvbuf_xz, recvCount_xz,MPI_DOUBLE,rank_xz,recvtag,MPI_COMM_SCALBL,&req2[11]);
	MPI_Isend(sendbuf_Xz, sendCount_Xz,MPI_DOUBLE,rank_Xz,sendtag,MPI_COMM_SCALBL,&req1[12]);
	MPI_Irecv(recvbuf_xZ, recvCount_xZ,MPI_DOUBLE,rank_xZ,recvtag,MPI_COMM_SCALBL,&req2[12]);
	MPI_Isend(sendbuf_xZ, sendCount_xZ,MPI_DOUBLE,rank_xZ,sendtag,MPI_COMM_SCALBL,&req1[13]);
	MPI_Irecv(recvbuf_Xz, recvCount_Xz,MPI_DOUBLE,rank_Xz,recvtag,MPI_COMM_SCALBL,&req2[13]);
	MPI_Isend(sendbuf_yz, sendCount_yz,MPI_DOUBLE,rank_yz,sendtag,MPI_COMM_SCALBL,&req1[14]);
	MPI_Irecv(recvbuf_YZ, recvCount_YZ,MPI_DOUBLE,rank_YZ,recvtag,MPI_COMM_SCALBL,&req2[14]);
	MPI_Isend(sendbuf_YZ, sendCount_YZ,MPI_DOUBLE,rank_YZ,sendtag,MPI_COMM_SCALBL,&req1[15]);
	MPI_Irecv(recvbuf_yz, recvCount_yz,MPI_DOUBLE,rank_yz,recvtag,MPI_COMM_SCALBL,&req2[15]);
	MPI_Isend(sendbuf_Yz, sendCount_Yz,MPI_DOUBLE,rank_Yz,sendtag,MPI_COMM_SCALBL,&req1[16]);
	MPI_Irecv(recvbuf_yZ, recvCount_yZ,MPI_DOUBLE,rank_yZ,recvtag,MPI_COMM_SCALBL,&req2[16]);
	MPI_Isend(sendbuf_yZ, sendCount_yZ,MPI_DOUBLE,rank_yZ,sendtag,MPI_COMM_SCALBL,&req1[17]);
	MPI_Irecv(recvbuf_Yz, recvCount_Yz,MPI_DOUBLE,rank_Yz,recvtag,MPI_COMM_SCALBL,&req2[17]);
	//...................................................................................
}
void ScaLBL_Communicator::RecvHalo(double *data){

	//...................................................................................
	MPI_Waitall(18,req1,stat1);
	MPI_Waitall(18,req2,stat2);
	ScaLBL_DeviceBarrier();
	//...................................................................................
	//...................................................................................
	ScaLBL_Scalar_Unpack(dvcRecvList_x, recvCount_x,recvbuf_x, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_y, recvCount_y,recvbuf_y, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_z, recvCount_z,recvbuf_z, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_X, recvCount_X,recvbuf_X, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Y, recvCount_Y,recvbuf_Y, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Z, recvCount_Z,recvbuf_Z, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xy, recvCount_xy,recvbuf_xy, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xY, recvCount_xY,recvbuf_xY, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XY, recvCount_XY,recvbuf_XY, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xz, recvCount_xz,recvbuf_xz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yz, recvCount_yz,recvbuf_yz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, data, N);
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................
}

void ScaLBL_Communicator::RegularLayout(IntArray map, const double *data, DoubleArray &regdata){
	// Gets data from the device and stores in regular layout
	int i,j,k,n,idx;
	int Nx = map.size(0);
	int Ny = map.size(1);
	int Nz = map.size(2);

	// initialize the array
	regdata.fill(0.f);
	
	double *TmpDat;
	double value;
	TmpDat = new double [N];
	ScaLBL_CopyToHost(&TmpDat[0],&data[0], N*sizeof(double));
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n=k*Nx*Ny+j*Nx+i;
				idx=map(i,j,k);
				if (!(idx<0)){
					value=TmpDat[idx];
					regdata(i,j,k)=value;
				}
			}
		}
	}
	//printf("r=%i, value=%f   ",rank,value);
	
	delete [] TmpDat;
}


void ScaLBL_Communicator::Color_BC_z(int *Map, double *Phi, double *Den, double vA, double vB){
	double Value=(vA-vB)/(vA+vB);
	if (kproc == 0) {
		// Set the phase indicator field and density on the z inlet
		ScaLBL_Color_BC_z(dvcSendList_z, Map, Phi, Den, vA, vB, sendCount_z, N);
		//ScaLBL_SetSlice_z(Phi,Value,Nx,Ny,Nz,0);
	}
}

void ScaLBL_Communicator::Color_BC_Z(int *Map, double *Phi, double *Den, double vA, double vB){
	double Value=(vA-vB)/(vA+vB);
	if (kproc == nprocz-1){
		// Set the phase indicator field and density on the Z outlet
		ScaLBL_Color_BC_Z(dvcSendList_Z, Map, Phi, Den, vA, vB, sendCount_Z, N);
		//ScaLBL_SetSlice_z(Phi,Value,Nx,Ny,Nz,Nz-1);
	}
}

void ScaLBL_Communicator::D3Q19_Pressure_BC_z(int *neighborList, double *fq, double din, int time){
    //ScaLBL_D3Q19_Pressure_BC_z(int *LIST,fq,din,Nx,Ny,Nz);
	if (kproc == 0) {
		if (time%2==0){
			ScaLBL_D3Q19_AAeven_Pressure_BC_z(dvcSendList_z, fq, din, sendCount_z, N);
		}
		else{
			ScaLBL_D3Q19_AAodd_Pressure_BC_z(neighborList, dvcSendList_z, fq, din, sendCount_z, N);
		}
	}
}

void ScaLBL_Communicator::D3Q19_Pressure_BC_Z(int *neighborList, double *fq, double dout, int time){
    //ScaLBL_D3Q19_Pressure_BC_Z(int *LIST,fq,dout,Nx,Ny,Nz);
	if (kproc == nprocz-1){
		if (time%2==0){
			ScaLBL_D3Q19_AAeven_Pressure_BC_Z(dvcSendList_Z, fq, dout, sendCount_Z, N);
		}
		else{
			ScaLBL_D3Q19_AAodd_Pressure_BC_Z(neighborList, dvcSendList_Z, fq, dout, sendCount_Z, N);
		}
	}
}

double ScaLBL_Communicator::D3Q19_Flux_BC_z(int *neighborList, double *fq, double flux, int time){
	double sum, locsum, din;
	double LocInletArea, InletArea;
	
	// Note that flux = rho_0 * Q
	
	// Compute the inlet area
	if (kproc == 0) 
		LocInletArea = double(sendCount_z);
	else LocInletArea = 0.f;
	
	MPI_Allreduce(&LocInletArea,&InletArea,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_SCALBL);
	//printf("Inlet area = %f \n", InletArea);

	// Set the flux BC
	locsum = 0.f;
	if (time%2==0){
		if (kproc == 0) 
			locsum = ScaLBL_D3Q19_AAeven_Flux_BC_z(dvcSendList_z, fq, flux, InletArea, sendCount_z, N);
		
		MPI_Allreduce(&locsum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_SCALBL);
		din = flux/InletArea + sum;
		//if (rank==0) printf("computed din (even) =%f \n",din);
		if (kproc == 0)
			ScaLBL_D3Q19_AAeven_Pressure_BC_z(dvcSendList_z, fq, din, sendCount_z, N);
	}
	else{
		if (kproc == 0) 
			locsum = ScaLBL_D3Q19_AAodd_Flux_BC_z(neighborList, dvcSendList_z, fq, flux, InletArea, sendCount_z, N);

		MPI_Allreduce(&locsum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_SCALBL);
		din = flux/InletArea + sum;
		
		//if (rank==0) printf("computed din (odd)=%f \n",din);
		if (kproc == 0) 
			ScaLBL_D3Q19_AAodd_Pressure_BC_z(neighborList, dvcSendList_z, fq, din, sendCount_z, N);
	}
	//printf("Inlet pressure = %f \n", din);
	return din;
}

void ScaLBL_Communicator::PrintD3Q19(){
	printf("Printing D3Q19 communication buffer contents \n");

	int i,n;
	double f;
	int *TempBuffer;
	TempBuffer = new int [5*recvCount_x];

	//.......................................................................
	// Re-index the send lists
	ScaLBL_CopyToHost(TempBuffer,dvcRecvDist_x,5*recvCount_x*sizeof(int));
	for (i=0; i<5*recvCount_x; i++){
		n = TempBuffer[i];
		f = recvbuf_x[i];
		printf("Receive %f to %i from buffer index %i \n", f, n, i);
	}

	delete [] TempBuffer;
}

