/* Flow adaptor class for multiphase flow methods */

#include "analysis/Membrane.h"
#include "analysis/distance.h"

Membrane::Membrane(std::shared_ptr <Domain> Dm, int *initialNeighborList, int Nsites) {

	Np = Nsites;
	neighborList = new int[18*Np];
	/* Copy neighborList */
	for (int idx=0; idx<Np; idx++){
		for (int q = 0; q<18; q++){
			neighborList[q*Np+idx] = initialNeighborList[q*18+idx];
		}
	}
	/* Copy communication lists */
	//......................................................................................
	//Lock=false; // unlock the communicator
	//......................................................................................
	// Create a separate copy of the communicator for the device
    //MPI_COMM_SCALBL = Dm->Comm.dup();
	//......................................................................................
	// Copy the domain size and communication information directly from Dm
	Nx = Dm->Nx;
	Ny = Dm->Ny;
	Nz = Dm->Nz;
	N = Nx*Ny*Nz;
	//next=0;
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
	sendCount_x=Dm->sendCount("x");
	sendCount_y=Dm->sendCount("y");
	sendCount_z=Dm->sendCount("z");
	sendCount_X=Dm->sendCount("X");
	sendCount_Y=Dm->sendCount("Y");
	sendCount_Z=Dm->sendCount("Z");
	sendCount_xy=Dm->sendCount("xy");
	sendCount_yz=Dm->sendCount("yz");
	sendCount_xz=Dm->sendCount("xz");
	sendCount_Xy=Dm->sendCount("Xy");
	sendCount_Yz=Dm->sendCount("Yz");
	sendCount_xZ=Dm->sendCount("xZ");
	sendCount_xY=Dm->sendCount("xY");
	sendCount_yZ=Dm->sendCount("yZ");
	sendCount_Xz=Dm->sendCount("Xz");
	sendCount_XY=Dm->sendCount("XY");
	sendCount_YZ=Dm->sendCount("YZ");
	sendCount_XZ=Dm->sendCount("XZ");
	recvCount_x=Dm->recvCount("x");
	recvCount_y=Dm->recvCount("y");
	recvCount_z=Dm->recvCount("z");
	recvCount_X=Dm->recvCount("X");
	recvCount_Y=Dm->recvCount("Y");
	recvCount_Z=Dm->recvCount("Z");
	recvCount_xy=Dm->recvCount("xy");
	recvCount_yz=Dm->recvCount("yz");
	recvCount_xz=Dm->recvCount("xz");
	recvCount_Xy=Dm->recvCount("Xy");
	recvCount_Yz=Dm->recvCount("Yz");
	recvCount_xZ=Dm->recvCount("xZ");
	recvCount_xY=Dm->recvCount("xY");
	recvCount_yZ=Dm->recvCount("yZ");
	recvCount_Xz=Dm->recvCount("Xz");
	recvCount_XY=Dm->recvCount("XY");
	recvCount_YZ=Dm->recvCount("YZ");
	recvCount_XZ=Dm->recvCount("XZ");\
	
	/* check symmetry for send / recv counts */
	if (sendCount_x != recvCount_X)   printf("WARNING: rank %i send/recv mismatch (x/X)! \n",rank);
	if (sendCount_y != recvCount_Y)   printf("WARNING: rank %i send/recv mismatch (y/Y)! \n",rank);
	if (sendCount_z != recvCount_Z)   printf("WARNING: rank %i send/recv mismatch (z/Z)! \n",rank);
	if (sendCount_X != recvCount_x)   printf("WARNING: rank %i send/recv mismatch (X/x)! \n",rank);
	if (sendCount_Y != recvCount_y)   printf("WARNING: rank %i send/recv mismatch (Y/y)! \n",rank);
	if (sendCount_x != recvCount_z)   printf("WARNING: rank %i send/recv mismatch (Z/z)! \n",rank);	
	if (sendCount_xy != recvCount_XY) printf("WARNING: rank %i send/recv mismatch (xy/XY)! \n",rank);
	if (sendCount_Xy != recvCount_xY) printf("WARNING: rank %i send/recv mismatch (Xy/xY)! \n",rank);
	if (sendCount_xY != recvCount_Xy) printf("WARNING: rank %i send/recv mismatch (xY/Xy)! \n",rank);
	if (sendCount_XY != recvCount_xy) printf("WARNING: rank %i send/recv mismatch (XY/xy)! \n",rank);
	if (sendCount_xz != recvCount_XZ) printf("WARNING: rank %i send/recv mismatch (xz/XZ)! \n",rank);
	if (sendCount_Xz != recvCount_xZ) printf("WARNING: rank %i send/recv mismatch (Xz/xZ)! \n",rank);
	if (sendCount_xZ != recvCount_Xz) printf("WARNING: rank %i send/recv mismatch (xZ/Xz)! \n",rank);
	if (sendCount_XZ != recvCount_xz) printf("WARNING: rank %i send/recv mismatch (XZ/xz)! \n",rank);
	if (sendCount_yz != recvCount_YZ) printf("WARNING: rank %i send/recv mismatch (yz/YZ)! \n",rank);
	if (sendCount_Yz != recvCount_yZ) printf("WARNING: rank %i send/recv mismatch (Yz/yZ)! \n",rank);
	if (sendCount_yZ != recvCount_Yz) printf("WARNING: rank %i send/recv mismatch (yZ/Yz)! \n",rank);
	if (sendCount_YZ != recvCount_yz) printf("WARNING: rank %i send/recv mismatch (YZ/yz)! \n",rank);

	iproc = Dm->iproc();
	jproc = Dm->jproc();
	kproc = Dm->kproc();
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
	//BoundaryCondition = Dm->BoundaryCondition;
	//......................................................................................

	ScaLBL_AllocateZeroCopy((void **) &sendbuf_x, 2*5*sendCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_X, 2*5*sendCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_y, 2*5*sendCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Y, 2*5*sendCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_z, 2*5*sendCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Z, 2*5*sendCount_Z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xy, 2*sendCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xY, 2*sendCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Xy, 2*sendCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_XY, 2*sendCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xz, 2*sendCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xZ, 2*sendCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Xz, 2*sendCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_XZ, 2*sendCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_yz, 2*sendCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_yZ, 2*sendCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Yz, 2*sendCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_YZ, 2*sendCount_YZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_x, 2*5*recvCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_X, 2*5*recvCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_y, 2*5*recvCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Y, 2*5*recvCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_z, 2*5*recvCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Z, 2*5*recvCount_Z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xy, 2*recvCount_xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xY, 2*recvCount_xY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Xy, 2*recvCount_Xy*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_XY, 2*recvCount_XY*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xz, 2*recvCount_xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xZ, 2*recvCount_xZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Xz, 2*recvCount_Xz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_XZ, 2*recvCount_XZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_yz, 2*recvCount_yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_yZ, 2*recvCount_yZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Yz, 2*recvCount_Yz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_YZ, 2*recvCount_YZ*sizeof(double));	// Allocate device memory
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

	//......................................................................................
	ScaLBL_CopyToZeroCopy(dvcSendList_x,Dm->sendList("x"),sendCount_x*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_X,Dm->sendList("X"),sendCount_X*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_y,Dm->sendList("y"),sendCount_y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Y,Dm->sendList("Y"),sendCount_Y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_z,Dm->sendList("z"),sendCount_z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Z,Dm->sendList("Z"),sendCount_Z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xy,Dm->sendList("xy"),sendCount_xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_XY,Dm->sendList("XY"),sendCount_XY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xY,Dm->sendList("xY"),sendCount_xY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Xy,Dm->sendList("Xy"),sendCount_Xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xz,Dm->sendList("xz"),sendCount_xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_XZ,Dm->sendList("XZ"),sendCount_XZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_xZ,Dm->sendList("xZ"),sendCount_xZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Xz,Dm->sendList("Xz"),sendCount_Xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_yz,Dm->sendList("yz"),sendCount_yz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_YZ,Dm->sendList("YZ"),sendCount_YZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_yZ,Dm->sendList("yZ"),sendCount_yZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcSendList_Yz,Dm->sendList("Yz"),sendCount_Yz*sizeof(int));
	//......................................................................................
	ScaLBL_CopyToZeroCopy(dvcRecvList_x,Dm->recvList("x"),recvCount_x*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_X,Dm->recvList("X"),recvCount_X*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_y,Dm->recvList("y"),recvCount_y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Y,Dm->recvList("Y"),recvCount_Y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_z,Dm->recvList("z"),recvCount_z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Z,Dm->recvList("Z"),recvCount_Z*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xy,Dm->recvList("xy"),recvCount_xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_XY,Dm->recvList("XY"),recvCount_XY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xY,Dm->recvList("xY"),recvCount_xY*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Xy,Dm->recvList("Xy"),recvCount_Xy*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xz,Dm->recvList("xz"),recvCount_xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_XZ,Dm->recvList("XZ"),recvCount_XZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_xZ,Dm->recvList("xZ"),recvCount_xZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Xz,Dm->recvList("Xz"),recvCount_Xz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_yz,Dm->recvList("yz"),recvCount_yz*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_YZ,Dm->recvList("YZ"),recvCount_YZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_yZ,Dm->recvList("yZ"),recvCount_yZ*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Yz,Dm->recvList("Yz"),recvCount_Yz*sizeof(int));
	//......................................................................................
	
}

Membrane::~Membrane() {
	
	ScaLBL_FreeDeviceMemory( sendbuf_x );
	ScaLBL_FreeDeviceMemory( sendbuf_X );
	ScaLBL_FreeDeviceMemory( sendbuf_y );
	ScaLBL_FreeDeviceMemory( sendbuf_Y );
	ScaLBL_FreeDeviceMemory( sendbuf_z );
	ScaLBL_FreeDeviceMemory( sendbuf_Z );
	ScaLBL_FreeDeviceMemory( sendbuf_xy );
	ScaLBL_FreeDeviceMemory( sendbuf_xY );
	ScaLBL_FreeDeviceMemory( sendbuf_Xy );
	ScaLBL_FreeDeviceMemory( sendbuf_XY );
	ScaLBL_FreeDeviceMemory( sendbuf_xz );
	ScaLBL_FreeDeviceMemory( sendbuf_xZ );
	ScaLBL_FreeDeviceMemory( sendbuf_Xz );
	ScaLBL_FreeDeviceMemory( sendbuf_XZ );
	ScaLBL_FreeDeviceMemory( sendbuf_yz );
	ScaLBL_FreeDeviceMemory( sendbuf_yZ );
	ScaLBL_FreeDeviceMemory( sendbuf_Yz );
	ScaLBL_FreeDeviceMemory( sendbuf_YZ );
	ScaLBL_FreeDeviceMemory( recvbuf_x );
	ScaLBL_FreeDeviceMemory( recvbuf_X );
	ScaLBL_FreeDeviceMemory( recvbuf_y );
	ScaLBL_FreeDeviceMemory( recvbuf_Y );
	ScaLBL_FreeDeviceMemory( recvbuf_z );
	ScaLBL_FreeDeviceMemory( recvbuf_Z );
	ScaLBL_FreeDeviceMemory( recvbuf_xy );
	ScaLBL_FreeDeviceMemory( recvbuf_xY );
	ScaLBL_FreeDeviceMemory( recvbuf_Xy );
	ScaLBL_FreeDeviceMemory( recvbuf_XY );
	ScaLBL_FreeDeviceMemory( recvbuf_xz );
	ScaLBL_FreeDeviceMemory( recvbuf_xZ );
	ScaLBL_FreeDeviceMemory( recvbuf_Xz );
	ScaLBL_FreeDeviceMemory( recvbuf_XZ );
	ScaLBL_FreeDeviceMemory( recvbuf_yz );
	ScaLBL_FreeDeviceMemory( recvbuf_yZ );
	ScaLBL_FreeDeviceMemory( recvbuf_Yz );
	ScaLBL_FreeDeviceMemory( recvbuf_YZ );
	ScaLBL_FreeDeviceMemory( dvcSendList_x );
	ScaLBL_FreeDeviceMemory( dvcSendList_X );
	ScaLBL_FreeDeviceMemory( dvcSendList_y );
	ScaLBL_FreeDeviceMemory( dvcSendList_Y );
	ScaLBL_FreeDeviceMemory( dvcSendList_z );
	ScaLBL_FreeDeviceMemory( dvcSendList_Z );
	ScaLBL_FreeDeviceMemory( dvcSendList_xy );
	ScaLBL_FreeDeviceMemory( dvcSendList_xY );
	ScaLBL_FreeDeviceMemory( dvcSendList_Xy );
	ScaLBL_FreeDeviceMemory( dvcSendList_XY );
	ScaLBL_FreeDeviceMemory( dvcSendList_xz );
	ScaLBL_FreeDeviceMemory( dvcSendList_xZ );
	ScaLBL_FreeDeviceMemory( dvcSendList_Xz );
	ScaLBL_FreeDeviceMemory( dvcSendList_XZ );
	ScaLBL_FreeDeviceMemory( dvcSendList_yz );
	ScaLBL_FreeDeviceMemory( dvcSendList_yZ );
	ScaLBL_FreeDeviceMemory( dvcSendList_Yz );
	ScaLBL_FreeDeviceMemory( dvcSendList_YZ );
	ScaLBL_FreeDeviceMemory( dvcRecvList_x );
	ScaLBL_FreeDeviceMemory( dvcRecvList_X );
	ScaLBL_FreeDeviceMemory( dvcRecvList_y );
	ScaLBL_FreeDeviceMemory( dvcRecvList_Y );
	ScaLBL_FreeDeviceMemory( dvcRecvList_z );
	ScaLBL_FreeDeviceMemory( dvcRecvList_Z );
	ScaLBL_FreeDeviceMemory( dvcRecvList_xy );
	ScaLBL_FreeDeviceMemory( dvcRecvList_xY );
	ScaLBL_FreeDeviceMemory( dvcRecvList_Xy );
	ScaLBL_FreeDeviceMemory( dvcRecvList_XY );
	ScaLBL_FreeDeviceMemory( dvcRecvList_xz );
	ScaLBL_FreeDeviceMemory( dvcRecvList_xZ );
	ScaLBL_FreeDeviceMemory( dvcRecvList_Xz );
	ScaLBL_FreeDeviceMemory( dvcRecvList_XZ );
	ScaLBL_FreeDeviceMemory( dvcRecvList_yz );
	ScaLBL_FreeDeviceMemory( dvcRecvList_yZ );
	ScaLBL_FreeDeviceMemory( dvcRecvList_Yz );
	ScaLBL_FreeDeviceMemory( dvcRecvList_YZ );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_x );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_X );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_y );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_Y );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_z );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_Z );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_xy );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_xY );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_Xy );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_XY );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_xz );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_xZ );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_Xz );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_XZ );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_yz );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_yZ );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_Yz );
	ScaLBL_FreeDeviceMemory( dvcRecvDist_YZ );
}

int Membrane::D3Q19_MapSendRecv(const int Cqx, const int Cqy, const int Cqz, 
		const int shift_x, const int shift_y, const int shift_z, 
		const int sendCount, const int recvCount,
		int rank_q, int rank_Q, int startSend, int startRecv,
		std::shared_ptr <Domain> Dm, const int *originalSendList, 
		const DoubleArray &Distance, int *membraneSendList, int *membraneRecvList){

	int sendtag = 2389;
	int recvtag = 2389;
	int i,j,k,n,nn,idx;
	double dist,locdist;
	
	int *SendList;
	int *RecvList;
	int *RecvList_q;
	SendList = new int [sendCount]; // for this process
	RecvList = new int [recvCount];   // filled in by the neighbor process
	RecvList_q = new int [sendCount]; // goes to the neighbor process

	int regularCount = 0;
	int membraneCount = 0;
	for (idx=0; idx<sendCount; idx++){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = originalSendList[idx]; // if (rank == 0) printf("@ rank:%d n=%d\n",rank,n);
		// Get the 3-D indices from the send process
		k = n/(Nx*Ny); j = (n-Nx*Ny*k)/Nx; i = n-Nx*Ny*k-Nx*j;
		// if (rank ==0) printf("@ Get 3D indices from the send process: i=%d, j=%d, k=%d\n",i,j,k);
		
		/* distance to membrane at the send site */
		dist = Distance(i,j,k);

		// Streaming for the non-local distribution
		i += Cqx; j += Cqy; k += Cqz;
		// if (rank == 0) printf("@ Streaming for the non-local distribution: i=%d, j=%d, k=%d\n",i,j,k);
		/* distance to membrane at the recv site */
		locdist = Distance(i,j,k);
		
		// Compute 1D index for the neighbor and save
		i += shift_x;
		j += shift_y;
		k += shift_z;
		nn = k*Nx*Ny+j*Nx+i;
		
		if (dist*locdist < 0.0){
			/* store membrane values at the end */
			membraneCount++;
			RecvList_q[sendCount-membraneCount] = nn;
			SendList[sendCount-membraneCount] = n;
		}
		else {
			RecvList_q[regularCount] = nn;
			SendList[regularCount++] = n;
		}
	}
	// send the recv list to the neighbor process
    req1[0] = Dm->Comm.Isend(RecvList_q, sendCount, rank_q, sendtag);
    req2[0] = Dm->Comm.Irecv(RecvList, recvCount, rank_Q, recvtag);
    Dm->Comm.barrier();
	
	// Return updated version to the device
	ScaLBL_CopyToDevice(&membraneSendList[startSend], SendList, sendCount*sizeof(int));
	ScaLBL_CopyToDevice(&membraneRecvList[startRecv], RecvList, recvCount*sizeof(int));

	// clean up the work arrays
	delete [] SendList;
	delete [] RecvList;
	delete [] RecvList_q;
	return membraneCount;
}

int Membrane::Create(std::shared_ptr <Domain> Dm, DoubleArray &Distance, IntArray &Map){
	int mlink = 0;
	int i,j,k;
	int idx, neighbor;
	double dist, locdist;
	/* go through the neighborlist structure */
	/* count & cut the links */
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				idx=Map(i,j,k);
				locdist=Distance(i,j,k);

				if (!(idx<0)){
					
					neighbor=Map(i-1,j,k);
					dist=Distance(i-1,j,k);
					if (dist*locdist < 0.0){
						neighborList[idx]=idx + 2*Np;
					}

					neighbor=Map(i+1,j,k);
					dist=Distance(i+1,j,k);
					if (dist*locdist < 0.0){
						neighborList[Np+idx] = idx + 1*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k);
					dist=Distance(i,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[2*Np+idx]=idx + 4*Np;
					}

					neighbor=Map(i,j+1,k);
					dist=Distance(i,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[3*Np+idx]=idx + 3*Np;
						mlink++;
					}

					neighbor=Map(i,j,k-1);
					dist=Distance(i,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[4*Np+idx]=idx + 6*Np;
					}

					neighbor=Map(i,j,k+1);
					dist=Distance(i,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[5*Np+idx]=idx + 5*Np;
						mlink++;
					}

					neighbor=Map(i-1,j-1,k);
					dist=Distance(i-1,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[6*Np+idx]=idx + 8*Np;
					}

					neighbor=Map(i+1,j+1,k);
					dist=Distance(i+1,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[7*Np+idx]=idx + 7*Np;
						mlink++;
					}

					neighbor=Map(i-1,j+1,k);
					dist=Distance(i-1,j+1,k);
					if (dist*locdist < 0.0){
						neighborList[8*Np+idx]=idx + 10*Np;
					}

					neighbor=Map(i+1,j-1,k);
					dist=Distance(i+1,j-1,k);
					if (dist*locdist < 0.0){
						neighborList[9*Np+idx]=idx + 9*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k-1);
					dist=Distance(i-1,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[10*Np+idx]=idx + 12*Np;
					}

					neighbor=Map(i+1,j,k+1);
					dist=Distance(i+1,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[11*Np+idx]=idx + 11*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k+1);
					dist=Distance(i-1,j,k+1);
					if (dist*locdist < 0.0){
						neighborList[12*Np+idx]=idx + 14*Np;
					}

					neighbor=Map(i+1,j,k-1);
					dist=Distance(i+1,j,k-1);
					if (dist*locdist < 0.0){
						neighborList[13*Np+idx]=idx + 13*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k-1);
					dist=Distance(i,j-1,k-1);
					if (dist*locdist < 0.0){
						neighborList[14*Np+idx]=idx + 16*Np;
					}

					neighbor=Map(i,j+1,k+1);
					dist=Distance(i,j+1,k+1);
					if (dist*locdist < 0.0){
						neighborList[15*Np+idx]=idx + 15*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k+1);
					dist=Distance(i,j-1,k+1);
					if (dist*locdist < 0.0){
						neighborList[16*Np+idx]=idx + 18*Np;
					}

					neighbor=Map(i,j+1,k-1);
					dist=Distance(i,j+1,k-1);
					if (dist*locdist < 0.0){
						neighborList[17*Np+idx]=idx + 17*Np;
						mlink++;
					}
				}
			}
		}
	}
	
	/* allocate memory */
	membraneLinks = new int [2*mlink];
	membraneDist = new double [2*mlink];
	membraneTag = new int [mlink];

	/* construct the membrane*/
	mlink = 0;
	int insideMem = 0; int outsideMem = 0;
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				idx=Map(i,j,k);
				locdist=Distance(i,j,k);

				if (!(idx<0)){

					neighbor=Map(i+1,j,k);
					dist=Distance(i+1,j,k);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 1*Np;
						membraneLinks[outsideMem] = neighbor + 2*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i,j+1,k);
					dist=Distance(i,j+1,k);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 3*Np;
						membraneLinks[outsideMem] = neighbor + 4*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i,j,k+1);
					dist=Distance(i,j,k+1);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 5*Np;
						membraneLinks[outsideMem] = neighbor + 6*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j+1,k);
					dist=Distance(i+1,j+1,k);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 7*Np;
						membraneLinks[outsideMem] = neighbor+8*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j-1,k);
					dist=Distance(i+1,j-1,k);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 9*Np;
						membraneLinks[outsideMem] = neighbor + 10*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j,k+1);
					dist=Distance(i+1,j,k+1);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 11*Np;
						membraneLinks[outsideMem] = neighbor + 12*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j,k-1);
					dist=Distance(i+1,j,k-1);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 13*Np;
						membraneLinks[outsideMem] = neighbor + 14*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i,j+1,k+1);
					dist=Distance(i,j+1,k+1);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 15*Np;
						membraneLinks[outsideMem] =neighbor + 16*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}

					neighbor=Map(i,j+1,k-1);
					dist=Distance(i,j+1,k-1);
					if (dist*locdist < 0.0){
						if (locdist < 0.0){
							insideMem = 2*mlink;
							outsideMem = 2*mlink+1;
						}
						else{
							insideMem = 2*mlink+1;
							outsideMem = 2*mlink;
						}
						membraneLinks[insideMem] = idx + 17*Np;
						membraneLinks[outsideMem] = neighbor + 18*Np;
						membraneDist[insideMem] = locdist;
						membraneDist[outsideMem] = dist;
						mlink++;
					}
				}
			}
		}
	}
	
	/* Re-organize communication based on membrane structure*/
	//membraneCount_x = D3Q19_MapSendRecv(Cx, Cy, Cz, Dm->sendList("x"), Distance, );
	// MPI_COMM_SCALBL.barrier();
	

	//...................................................................................
	// Set up the recieve distribution lists
	//...................................................................................
	//...Map recieve list for the X face: q=2,8,10,12,14 .................................
/*	D3Q19_MapRecv(-1,0,0, Dm->recvList("X"),0,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,-1,0,Dm->recvList("X"),recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,1,0, Dm->recvList("X"),2*recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,0,-1,Dm->recvList("X"),3*recvCount_X,recvCount_X,dvcRecvDist_X);
	D3Q19_MapRecv(-1,0,1, Dm->recvList("X"),4*recvCount_X,recvCount_X,dvcRecvDist_X);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	D3Q19_MapRecv(1,0,0, Dm->recvList("x"),0,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,1,0, Dm->recvList("x"),recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,-1,0,Dm->recvList("x"),2*recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,0,1, Dm->recvList("x"),3*recvCount_x,recvCount_x,dvcRecvDist_x);
	D3Q19_MapRecv(1,0,-1,Dm->recvList("x"),4*recvCount_x,recvCount_x,dvcRecvDist_x);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
	D3Q19_MapRecv(0,-1,0, Dm->recvList("Y"),0,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(-1,-1,0,Dm->recvList("Y"),recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(1,-1,0, Dm->recvList("Y"),2*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(0,-1,-1,Dm->recvList("Y"),3*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	D3Q19_MapRecv(0,-1,1, Dm->recvList("Y"),4*recvCount_Y,recvCount_Y,dvcRecvDist_Y);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
	D3Q19_MapRecv(0,1,0, Dm->recvList("y"),0,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(1,1,0, Dm->recvList("y"),recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(-1,1,0,Dm->recvList("y"),2*recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(0,1,1, Dm->recvList("y"),3*recvCount_y,recvCount_y,dvcRecvDist_y);
	D3Q19_MapRecv(0,1,-1,Dm->recvList("y"),4*recvCount_y,recvCount_y,dvcRecvDist_y);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
	D3Q19_MapRecv(0,0,-1, Dm->recvList("Z"),0,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(-1,0,-1,Dm->recvList("Z"),recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(1,0,-1, Dm->recvList("Z"),2*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(0,-1,-1,Dm->recvList("Z"),3*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	D3Q19_MapRecv(0,1,-1, Dm->recvList("Z"),4*recvCount_Z,recvCount_Z,dvcRecvDist_Z);
	//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
	D3Q19_MapRecv(0,0,1, Dm->recvList("z"),0,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(1,0,1, Dm->recvList("z"),recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(-1,0,1,Dm->recvList("z"),2*recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(0,1,1, Dm->recvList("z"),3*recvCount_z,recvCount_z,dvcRecvDist_z);
	D3Q19_MapRecv(0,-1,1,Dm->recvList("z"),4*recvCount_z,recvCount_z,dvcRecvDist_z);
	//..................................................................................
	//...Map recieve list for the xy edge <<<8)................................
	D3Q19_MapRecv(-1,-1,0,Dm->recvList("XY"),0,recvCount_XY,dvcRecvDist_XY);
	//...Map recieve list for the Xy edge <<<9)................................
	D3Q19_MapRecv(1,-1,0,Dm->recvList("xY"),0,recvCount_xY,dvcRecvDist_xY);
	//...Map recieve list for the xY edge <<<10)................................
	D3Q19_MapRecv(-1,1,0,Dm->recvList("Xy"),0,recvCount_Xy,dvcRecvDist_Xy);
	//...Map recieve list for the XY edge <<<7)................................
	D3Q19_MapRecv(1,1,0,Dm->recvList("xy"),0,recvCount_xy,dvcRecvDist_xy);
	//...Map recieve list for the xz edge <<<12)................................
	D3Q19_MapRecv(-1,0,-1,Dm->recvList("XZ"),0,recvCount_XZ,dvcRecvDist_XZ);
	//...Map recieve list for the xZ edge <<<14)................................
	D3Q19_MapRecv(-1,0,1,Dm->recvList("Xz"),0,recvCount_Xz,dvcRecvDist_Xz);
	//...Map recieve list for the Xz edge <<<13)................................
	D3Q19_MapRecv(1,0,-1,Dm->recvList("xZ"),0,recvCount_xZ,dvcRecvDist_xZ);
	//...Map recieve list for the XZ edge <<<11)................................
	D3Q19_MapRecv(1,0,1,Dm->recvList("xz"),0,recvCount_xz,dvcRecvDist_xz);
	//...Map recieve list for the yz edge <<<16)................................
	D3Q19_MapRecv(0,-1,-1,Dm->recvList("YZ"),0,recvCount_YZ,dvcRecvDist_YZ);
	//...Map recieve list for the yZ edge <<<18)................................
	D3Q19_MapRecv(0,-1,1,Dm->recvList("Yz"),0,recvCount_Yz,dvcRecvDist_Yz);
	//...Map recieve list for the Yz edge <<<17)................................
	D3Q19_MapRecv(0,1,-1,Dm->recvList("yZ"),0,recvCount_yZ,dvcRecvDist_yZ);
	//...Map recieve list for the YZ edge <<<15)................................
	D3Q19_MapRecv(0,1,1,Dm->recvList("yz"),0,recvCount_yz,dvcRecvDist_yz);
	//...................................................................................

	//......................................................................................
	MPI_COMM_SCALBL.barrier();
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
	*/
	return mlink;
}
