/* Flow adaptor class for multiphase flow methods */

#include "common/Membrane.h"
#include "analysis/distance.h"

Membrane::Membrane(std::shared_ptr <Domain> Dm, int *dvcNeighborList, int Nsites) {

	Np = Nsites;
	initialNeighborList = new int[18*Np];
    ScaLBL_AllocateDeviceMemory((void **)&NeighborList, 18*Np*sizeof(int));    
    Lock=false; // unlock the communicator
	//......................................................................................
	// Create a separate copy of the communicator for the device
    MPI_COMM_SCALBL = Dm->Comm.dup();
    
    ScaLBL_CopyToHost(initialNeighborList, dvcNeighborList, 18*Np*sizeof(int));
    Dm->Comm.barrier();
    ScaLBL_CopyToDevice(NeighborList, initialNeighborList, 18*Np*sizeof(int));
    
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
	recvCount_XZ=Dm->recvCount("XZ");
	
	if (rank == 0){
		printf("**** Creating membrane data structure ****** \n");
	}
	printf("   Number of active lattice sites (rank = %i): %i \n",rank, Np);
	
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
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_x, 5*recvCount_x*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_X, 5*recvCount_X*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_y, 5*recvCount_y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_Y, 5*recvCount_Y*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_z, 5*recvCount_z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_Z, 5*recvCount_Z*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_xy, recvCount_xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_xY, recvCount_xY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_Xy, recvCount_Xy*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_XY, recvCount_XY*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_xz, recvCount_xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_xZ, recvCount_xZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_Xz, recvCount_Xz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_XZ, recvCount_XZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_yz, recvCount_yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_yZ, recvCount_yZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_Yz, recvCount_Yz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvLinks_YZ, recvCount_YZ*sizeof(int));	// Allocate device memory
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
	
	delete [] initialNeighborList;
	delete [] membraneLinks;    
	delete [] membraneTag; 
	delete [] membraneDist;
	
	ScaLBL_FreeDeviceMemory( coefficient_x );
	ScaLBL_FreeDeviceMemory( coefficient_X );
	ScaLBL_FreeDeviceMemory( coefficient_y );
	ScaLBL_FreeDeviceMemory( coefficient_Y );
	ScaLBL_FreeDeviceMemory( coefficient_z );
	ScaLBL_FreeDeviceMemory( coefficient_Z );
	
	ScaLBL_FreeDeviceMemory( NeighborList );
	ScaLBL_FreeDeviceMemory( MembraneLinks );
	ScaLBL_FreeDeviceMemory( MembraneCoef );
	
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
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_x );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_X );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_y );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_Y );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_z );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_Z );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_xy );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_xY );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_Xy );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_XY );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_xz );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_xZ );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_Xz );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_XZ );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_yz );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_yZ );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_Yz );
	ScaLBL_FreeDeviceMemory( dvcRecvLinks_YZ );
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

int Membrane::Create(std::shared_ptr <Domain> Dm, DoubleArray &Distance, IntArray &Map){
	int mlink = 0;
	int i,j,k;
	int idx, neighbor;
	double dist, locdist;
	
	if (rank == 0) printf("   Copy initial neighborlist... \n");
	int * neighborList = new int[18*Np];
	/* Copy neighborList */
	for (int idx=0; idx<Np; idx++){
		for (int q = 0; q<18; q++){
			neighborList[q*Np+idx] = initialNeighborList[q*Np+idx];
		}
	}
	
	/* go through the neighborlist structure */
	/* count & cut the links */
	if (rank == 0) printf("   Cut membrane links... \n");
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				idx=Map(i,j,k);
				locdist=Distance(i,j,k);

				if (!(idx<0)){
					
					neighbor=Map(i-1,j,k);
					dist=Distance(i-1,j,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[idx]=idx + 2*Np;
					}

					neighbor=Map(i+1,j,k);
					dist=Distance(i+1,j,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[Np+idx] = idx + 1*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k);
					dist=Distance(i,j-1,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[2*Np+idx]=idx + 4*Np;
					}

					neighbor=Map(i,j+1,k);
					dist=Distance(i,j+1,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[3*Np+idx]=idx + 3*Np;
						mlink++;
					}

					neighbor=Map(i,j,k-1);
					dist=Distance(i,j,k-1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[4*Np+idx]=idx + 6*Np;
					}

					neighbor=Map(i,j,k+1);
					dist=Distance(i,j,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[5*Np+idx]=idx + 5*Np;
						mlink++;
					}

					neighbor=Map(i-1,j-1,k);
					dist=Distance(i-1,j-1,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
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
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[8*Np+idx]=idx + 10*Np;
					}

					neighbor=Map(i+1,j-1,k);
					dist=Distance(i+1,j-1,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[9*Np+idx]=idx + 9*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k-1);
					dist=Distance(i-1,j,k-1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[10*Np+idx]=idx + 12*Np;
					}

					neighbor=Map(i+1,j,k+1);
					dist=Distance(i+1,j,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[11*Np+idx]=idx + 11*Np;
						mlink++;
					}

					neighbor=Map(i-1,j,k+1);
					dist=Distance(i-1,j,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[12*Np+idx]=idx + 14*Np;
					}

					neighbor=Map(i+1,j,k-1);
					dist=Distance(i+1,j,k-1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[13*Np+idx]=idx + 13*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k-1);
					dist=Distance(i,j-1,k-1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[14*Np+idx]=idx + 16*Np;
					}

					neighbor=Map(i,j+1,k+1);
					dist=Distance(i,j+1,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[15*Np+idx]=idx + 15*Np;
						mlink++;
					}

					neighbor=Map(i,j-1,k+1);
					dist=Distance(i,j-1,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[16*Np+idx]=idx + 18*Np;
					}

					neighbor=Map(i,j+1,k-1);
					dist=Distance(i,j+1,k-1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						neighborList[17*Np+idx]=idx + 17*Np;
						mlink++;
					}
				}
			}
		}
	}

	/* allocate memory */
	membraneTag = new int [mlink];
	membraneLinks = new int [2*mlink];
	membraneDist = new double [2*mlink];
	membraneLinkCount = mlink;

	if (rank == 0) printf("   (cut %i links crossing membrane) \n",mlink);

	/* construct the membrane*/
	/* *
	 *  Sites inside the membrane (negative distance) -- store at 2*mlink
	 *  Sites outside the membrane (positive distance) -- store at 2*mlink+1
	 */
	if (rank == 0) printf("   Construct membrane data structures... \n");
	mlink = 0;
	int localSite = 0; int neighborSite = 0;
	for (k=1;k<Nz-1;k++){
		for (j=1;j<Ny-1;j++){
			for (i=1;i<Nx-1;i++){
				idx=Map(i,j,k);
				locdist=Distance(i,j,k);

				if (!(idx<0)){

					neighbor=Map(i+1,j,k);
					dist=Distance(i+1,j,k);
					if (dist*locdist < 0.0){
						if (locdist < 0.0 && !(neighbor<0)){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 1*Np;
						membraneLinks[neighborSite] = neighbor + 2*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i,j+1,k);
					dist=Distance(i,j+1,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 3*Np;
						membraneLinks[neighborSite] = neighbor + 4*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i,j,k+1);
					dist=Distance(i,j,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 5*Np;
						membraneLinks[neighborSite] = neighbor + 6*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j+1,k);
					dist=Distance(i+1,j+1,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 7*Np;
						membraneLinks[neighborSite] = neighbor+8*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j-1,k);
					dist=Distance(i+1,j-1,k);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 9*Np;
						membraneLinks[neighborSite] = neighbor + 10*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j,k+1);
					dist=Distance(i+1,j,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 11*Np;
						membraneLinks[neighborSite] = neighbor + 12*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i+1,j,k-1);
					dist=Distance(i+1,j,k-1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 13*Np;
						membraneLinks[neighborSite] = neighbor + 14*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i,j+1,k+1);
					dist=Distance(i,j+1,k+1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 15*Np;
						membraneLinks[neighborSite] = neighbor + 16*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					}

					neighbor=Map(i,j+1,k-1);
					dist=Distance(i,j+1,k-1);
					if (dist*locdist < 0.0 && !(neighbor<0)){
						if (locdist < 0.0){
							localSite = 2*mlink;
							neighborSite = 2*mlink+1;
						}
						else{
							localSite = 2*mlink+1;
							neighborSite = 2*mlink;
						}
						membraneLinks[localSite] = idx + 17*Np;
						membraneLinks[neighborSite] = neighbor + 18*Np;
						membraneDist[localSite] = locdist;
						membraneDist[neighborSite] = dist;
						mlink++;
					} 
				}
			}
		}
	}
	
	if (rank == 0) printf("   Create device data structures... \n");

	/* Create device copies of data structures */
    ScaLBL_AllocateDeviceMemory((void **)&MembraneLinks, 2*mlink*sizeof(int));
    ScaLBL_AllocateDeviceMemory((void **)&MembraneCoef, 2*mlink*sizeof(double));
    //ScaLBL_AllocateDeviceMemory((void **)&MembraneDistance, 2*mlink*sizeof(double));
    ScaLBL_AllocateDeviceMemory((void **)&MembraneDistance, Nx*Ny*Nz*sizeof(double));
    
    ScaLBL_CopyToDevice(NeighborList, neighborList, 18*Np*sizeof(int));
    ScaLBL_CopyToDevice(MembraneLinks, membraneLinks, 2*mlink*sizeof(int));
    //ScaLBL_CopyToDevice(MembraneDistance, membraneDist, 2*mlink*sizeof(double));
    ScaLBL_CopyToDevice(MembraneDistance, Distance.data(), Nx*Ny*Nz*sizeof(double));

	
	if (rank == 0) printf("   Construct communication data structures... \n");
	/* Re-organize communication based on membrane structure*/
	//...Map recieve list for the X face: q=2,8,10,12,14 .................................
	linkCount_X[0] = D3Q19_MapRecv(-1,0,0, Dm->recvList("X"),0,recvCount_X,dvcRecvDist_X,dvcRecvLinks_X,Distance,Map);
	linkCount_X[1] = D3Q19_MapRecv(-1,-1,0,Dm->recvList("X"),recvCount_X,recvCount_X,dvcRecvDist_X,dvcRecvLinks_X,Distance,Map);
	linkCount_X[2] = D3Q19_MapRecv(-1,1,0, Dm->recvList("X"),2*recvCount_X,recvCount_X,dvcRecvDist_X,dvcRecvLinks_X,Distance,Map);
	linkCount_X[3] = D3Q19_MapRecv(-1,0,-1,Dm->recvList("X"),3*recvCount_X,recvCount_X,dvcRecvDist_X,dvcRecvLinks_X,Distance,Map);
	linkCount_X[4] = D3Q19_MapRecv(-1,0,1, Dm->recvList("X"),4*recvCount_X,recvCount_X,dvcRecvDist_X,dvcRecvLinks_X,Distance,Map);
	//...................................................................................
	//...Map recieve list for the x face: q=1,7,9,11,13..................................
	linkCount_x[0] = D3Q19_MapRecv(1,0,0, Dm->recvList("x"),0,recvCount_x,dvcRecvDist_x,dvcRecvLinks_x,Distance,Map);
	linkCount_x[1] = D3Q19_MapRecv(1,1,0, Dm->recvList("x"),recvCount_x,recvCount_x,dvcRecvDist_x,dvcRecvLinks_x,Distance,Map);
	linkCount_x[2] = D3Q19_MapRecv(1,-1,0,Dm->recvList("x"),2*recvCount_x,recvCount_x,dvcRecvDist_x,dvcRecvLinks_x,Distance,Map);
	linkCount_x[3] = D3Q19_MapRecv(1,0,1, Dm->recvList("x"),3*recvCount_x,recvCount_x,dvcRecvDist_x,dvcRecvLinks_x,Distance,Map);
	linkCount_x[4] = D3Q19_MapRecv(1,0,-1,Dm->recvList("x"),4*recvCount_x,recvCount_x,dvcRecvDist_x,dvcRecvLinks_x,Distance,Map);
	//...................................................................................
	//...Map recieve list for the y face: q=4,8,9,16,18 ...................................
	linkCount_Y[0] = D3Q19_MapRecv(0,-1,0, Dm->recvList("Y"),0,recvCount_Y,dvcRecvDist_Y,dvcRecvLinks_Y,Distance,Map);
	linkCount_Y[1] = D3Q19_MapRecv(-1,-1,0,Dm->recvList("Y"),recvCount_Y,recvCount_Y,dvcRecvDist_Y,dvcRecvLinks_Y,Distance,Map);
	linkCount_Y[2] = D3Q19_MapRecv(1,-1,0, Dm->recvList("Y"),2*recvCount_Y,recvCount_Y,dvcRecvDist_Y,dvcRecvLinks_Y,Distance,Map);
	linkCount_Y[3] = D3Q19_MapRecv(0,-1,-1,Dm->recvList("Y"),3*recvCount_Y,recvCount_Y,dvcRecvDist_Y,dvcRecvLinks_Y,Distance,Map);
	linkCount_Y[4] = D3Q19_MapRecv(0,-1,1, Dm->recvList("Y"),4*recvCount_Y,recvCount_Y,dvcRecvDist_Y,dvcRecvLinks_Y,Distance,Map);
	//...................................................................................
	//...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
	linkCount_y[0] = D3Q19_MapRecv(0,1,0, Dm->recvList("y"),0,recvCount_y,dvcRecvDist_y,dvcRecvLinks_y,Distance,Map);
	linkCount_y[1] = D3Q19_MapRecv(1,1,0, Dm->recvList("y"),recvCount_y,recvCount_y,dvcRecvDist_y,dvcRecvLinks_y,Distance,Map);
	linkCount_y[2] = D3Q19_MapRecv(-1,1,0,Dm->recvList("y"),2*recvCount_y,recvCount_y,dvcRecvDist_y,dvcRecvLinks_y,Distance,Map);
	linkCount_y[3] = D3Q19_MapRecv(0,1,1, Dm->recvList("y"),3*recvCount_y,recvCount_y,dvcRecvDist_y,dvcRecvLinks_y,Distance,Map);
	linkCount_y[4] = D3Q19_MapRecv(0,1,-1,Dm->recvList("y"),4*recvCount_y,recvCount_y,dvcRecvDist_y,dvcRecvLinks_y,Distance,Map);
	//...................................................................................
	//...Map recieve list for the z face<<<6,12,13,16,17)..............................................
	linkCount_Z[0] = D3Q19_MapRecv(0,0,-1, Dm->recvList("Z"),0,recvCount_Z,dvcRecvDist_Z,dvcRecvLinks_Z,Distance,Map);
	linkCount_Z[1] = D3Q19_MapRecv(-1,0,-1,Dm->recvList("Z"),recvCount_Z,recvCount_Z,dvcRecvDist_Z,dvcRecvLinks_Z,Distance,Map);
	linkCount_Z[2] = D3Q19_MapRecv(1,0,-1, Dm->recvList("Z"),2*recvCount_Z,recvCount_Z,dvcRecvDist_Z,dvcRecvLinks_Z,Distance,Map);
	linkCount_Z[3] = D3Q19_MapRecv(0,-1,-1,Dm->recvList("Z"),3*recvCount_Z,recvCount_Z,dvcRecvDist_Z,dvcRecvLinks_Z,Distance,Map);
	linkCount_Z[4] = D3Q19_MapRecv(0,1,-1, Dm->recvList("Z"),4*recvCount_Z,recvCount_Z,dvcRecvDist_Z,dvcRecvLinks_Z,Distance,Map);
	//...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
	linkCount_z[0] = D3Q19_MapRecv(0,0,1, Dm->recvList("z"),0,recvCount_z,dvcRecvDist_z,dvcRecvLinks_z,Distance,Map);
	linkCount_z[1] = D3Q19_MapRecv(1,0,1, Dm->recvList("z"),recvCount_z,recvCount_z,dvcRecvDist_z,dvcRecvLinks_z,Distance,Map);
	linkCount_z[2] = D3Q19_MapRecv(-1,0,1,Dm->recvList("z"),2*recvCount_z,recvCount_z,dvcRecvDist_z,dvcRecvLinks_z,Distance,Map);
	linkCount_z[3] = D3Q19_MapRecv(0,1,1, Dm->recvList("z"),3*recvCount_z,recvCount_z,dvcRecvDist_z,dvcRecvLinks_z,Distance,Map);
	linkCount_z[4] = D3Q19_MapRecv(0,-1,1,Dm->recvList("z"),4*recvCount_z,recvCount_z,dvcRecvDist_z,dvcRecvLinks_z,Distance,Map);
	//..................................................................................
	//...Map recieve list for the xy edge <<<8)................................
	linkCount_XY = D3Q19_MapRecv(-1,-1,0,Dm->recvList("XY"),0,recvCount_XY,dvcRecvDist_XY,dvcRecvLinks_XY,Distance,Map);
	//...Map recieve list for the Xy edge <<<9)................................
	linkCount_xY = D3Q19_MapRecv(1,-1,0,Dm->recvList("xY"),0,recvCount_xY,dvcRecvDist_xY,dvcRecvLinks_xY,Distance,Map);
	//...Map recieve list for the xY edge <<<10)................................
	linkCount_Xy = D3Q19_MapRecv(-1,1,0,Dm->recvList("Xy"),0,recvCount_Xy,dvcRecvDist_Xy,dvcRecvLinks_Xy,Distance,Map);
	//...Map recieve list for the XY edge <<<7)................................
	linkCount_xy = D3Q19_MapRecv(1,1,0,Dm->recvList("xy"),0,recvCount_xy,dvcRecvDist_xy,dvcRecvLinks_xy,Distance,Map);
	//...Map recieve list for the xz edge <<<12)................................
	linkCount_XZ = D3Q19_MapRecv(-1,0,-1,Dm->recvList("XZ"),0,recvCount_XZ,dvcRecvDist_XZ,dvcRecvLinks_XZ,Distance,Map);
	//...Map recieve list for the xZ edge <<<14)................................
	linkCount_Xz = D3Q19_MapRecv(-1,0,1,Dm->recvList("Xz"),0,recvCount_Xz,dvcRecvDist_Xz,dvcRecvLinks_Xz,Distance,Map);
	//...Map recieve list for the Xz edge <<<13)................................
	linkCount_xZ = D3Q19_MapRecv(1,0,-1,Dm->recvList("xZ"),0,recvCount_xZ,dvcRecvDist_xZ,dvcRecvLinks_xZ,Distance,Map);
	//...Map recieve list for the XZ edge <<<11)................................
	linkCount_xz = D3Q19_MapRecv(1,0,1,Dm->recvList("xz"),0,recvCount_xz,dvcRecvDist_xz,dvcRecvLinks_xz,Distance,Map);
	//...Map recieve list for the yz edge <<<16)................................
	linkCount_YZ = D3Q19_MapRecv(0,-1,-1,Dm->recvList("YZ"),0,recvCount_YZ,dvcRecvDist_YZ,dvcRecvLinks_YZ,Distance,Map);
	//...Map recieve list for the yZ edge <<<18)................................
	linkCount_Yz = D3Q19_MapRecv(0,-1,1,Dm->recvList("Yz"),0,recvCount_Yz,dvcRecvDist_Yz,dvcRecvLinks_Yz,Distance,Map);
	//...Map recieve list for the Yz edge <<<17)................................
	linkCount_yZ = D3Q19_MapRecv(0,1,-1,Dm->recvList("yZ"),0,recvCount_yZ,dvcRecvDist_yZ,dvcRecvLinks_yZ,Distance,Map);
	//...Map recieve list for the YZ edge <<<15)................................
	linkCount_yz = D3Q19_MapRecv(0,1,1,Dm->recvList("yz"),0,recvCount_yz,dvcRecvDist_yz,dvcRecvLinks_yz,Distance,Map);
	//...................................................................................
	if (rank == 0) printf("    x count = %i \n",linkCount_x[0]);
	if (rank == 0) printf("    X count = %i \n",linkCount_X[0]);
	if (rank == 0) printf("    y count = %i \n",linkCount_y[0]);
	if (rank == 0) printf("    Y count = %i \n",linkCount_Y[0]);
	if (rank == 0) printf("    z count = %i \n",linkCount_z[0]);
	if (rank == 0) printf("    Z count = %i \n",linkCount_Z[0]);

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
	// Allocate membrane coefficient buffers (for d3q7 recv)
	ScaLBL_AllocateZeroCopy((void **) &coefficient_x, (recvCount_x - linkCount_x[0])*sizeof(double));
	ScaLBL_AllocateZeroCopy((void **) &coefficient_X, (recvCount_X - linkCount_X[0])*sizeof(double));
	ScaLBL_AllocateZeroCopy((void **) &coefficient_y, (recvCount_y - linkCount_y[0])*sizeof(double));
	ScaLBL_AllocateZeroCopy((void **) &coefficient_Y, (recvCount_Y - linkCount_Y[0])*sizeof(double));
	ScaLBL_AllocateZeroCopy((void **) &coefficient_z, (recvCount_z - linkCount_z[0])*sizeof(double));
	ScaLBL_AllocateZeroCopy((void **) &coefficient_Z, (recvCount_Z - linkCount_Z[0])*sizeof(double));
	//......................................................................................
	
	return mlink;
}

int Membrane::D3Q19_MapRecv(int Cqx, int Cqy, int Cqz, const int *list, int start, int count,
		int *d3q19_recvlist, int *d3q19_linkList, DoubleArray &Distance,  IntArray &Map){
	
	int linkCount = 0;
	int memLinkCount=0;
	int i,j,k,n,nn,idx;
	int * ReturnDist;
	double distanceNonLocal,distanceLocal;
	ReturnDist=new int [count];
	int *LinkList=new int [count];
	int *memLinkList=new int [count];
	
	for (idx=0; idx<count; idx++){

		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx];
		// Get the 3-D indices from the send process
		k = n/(Nx*Ny); j = (n-Nx*Ny*k)/Nx; i = n-Nx*Ny*k-Nx*j;
		// if (rank ==0) printf("@ Get 3D indices from the send process: i=%d, j=%d, k=%d\n",i,j,k);

		distanceLocal = Distance(i,j,k);  // this site should be in the halo

		// Streaming for the non-local distribution
		i += Cqx; j += Cqy; k += Cqz;

		nn = Map(i,j,k);
		distanceNonLocal = Distance(i,j,k);
		
		
		printf("CHECK:  idx=%i, n=%i, (%i, %i, %i) shift {%i, %i, %i}, stored nn=%i \n",idx,n,i,j,k,Cqx,Cqy,Cqz,nn);
		
		ReturnDist[idx] = nn;
		//if (nn < 0){
		//	printf("   Check map for site (%i, %i, %i) based on Cq=(%i, %i, %i) \n", i,j,k,Cqx,Cqy,Cqz);
		//}
		
		/* tag the links to swap out later*/
		if (distanceLocal*distanceNonLocal < 0.0){
			memLinkList[memLinkCount++] = idx;
		}
		else {
			LinkList[linkCount++] = idx;
		}
	}
	
	/* add membrane links at the end */
	for (int link=0; link<memLinkCount; link++){
		idx = memLinkList[link];
		LinkList[linkCount+link] = idx;
	}
	
	/* quick check */
	if (idx != count){
		printf("ERROR forming membrane communication links! \n");
	}
	
	// Return updated version to the device
	ScaLBL_CopyToDevice(&d3q19_recvlist[start], ReturnDist, count*sizeof(int));
	ScaLBL_CopyToDevice(&d3q19_linkList[start], LinkList, count*sizeof(int));

	// clean up the work arrays
	delete [] ReturnDist;
	delete [] memLinkList;
	delete [] LinkList;
	return linkCount;
}

void Membrane::SendD3Q19AA(double *dist){

	if (Lock==true){
		ERROR("Membrane Error (SendD3Q19): Membrane communicator is locked -- did you forget to match Send/Recv calls?");
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
	
	req1[0] = MPI_COMM_SCALBL.Isend(sendbuf_x, 5*sendCount_x,rank_x,sendtag);
	req2[0] = MPI_COMM_SCALBL.Irecv(recvbuf_X, 5*recvCount_X,rank_X,recvtag);
	//...Packing for X face(1,7,9,11,13)................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,0,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_X,sendCount_X,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(9,dvcSendList_X,2*sendCount_X,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(11,dvcSendList_X,3*sendCount_X,sendCount_X,sendbuf_X,dist,N);
	ScaLBL_D3Q19_Pack(13,dvcSendList_X,4*sendCount_X,sendCount_X,sendbuf_X,dist,N);
	
	req1[1] = MPI_COMM_SCALBL.Isend(sendbuf_X, 5*sendCount_X,rank_X,sendtag);
	req2[1] = MPI_COMM_SCALBL.Irecv(recvbuf_x, 5*recvCount_x,rank_x,recvtag);
	
	//...Packing for y face(4,8,9,16,18).................................
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,0,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(8,dvcSendList_y,sendCount_y,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(9,dvcSendList_y,2*sendCount_y,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(16,dvcSendList_y,3*sendCount_y,sendCount_y,sendbuf_y,dist,N);
	ScaLBL_D3Q19_Pack(18,dvcSendList_y,4*sendCount_y,sendCount_y,sendbuf_y,dist,N);
	
	req1[2] = MPI_COMM_SCALBL.Isend(sendbuf_y, 5*sendCount_y,rank_y,sendtag);
	req2[2] = MPI_COMM_SCALBL.Irecv(recvbuf_Y, 5*recvCount_Y,rank_Y,recvtag);
	//...Packing for Y face(3,7,10,15,17).................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(7,dvcSendList_Y,sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(10,dvcSendList_Y,2*sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(15,dvcSendList_Y,3*sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	ScaLBL_D3Q19_Pack(17,dvcSendList_Y,4*sendCount_Y,sendCount_Y,sendbuf_Y,dist,N);
	
	req1[3] = MPI_COMM_SCALBL.Isend(sendbuf_Y, 5*sendCount_Y,rank_Y,sendtag);
	req2[3] = MPI_COMM_SCALBL.Irecv(recvbuf_y, 5*recvCount_y,rank_y,recvtag);
	
	//...Packing for z face(6,12,13,16,17)................................
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,0,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(12,dvcSendList_z,sendCount_z,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(13,dvcSendList_z,2*sendCount_z,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(16,dvcSendList_z,3*sendCount_z,sendCount_z,sendbuf_z,dist,N);
	ScaLBL_D3Q19_Pack(17,dvcSendList_z,4*sendCount_z,sendCount_z,sendbuf_z,dist,N);
	
	req1[4] = MPI_COMM_SCALBL.Isend(sendbuf_z, 5*sendCount_z,rank_z,sendtag);
	req2[4] = MPI_COMM_SCALBL.Irecv(recvbuf_Z, 5*recvCount_Z,rank_Z,recvtag);
	
	//...Packing for Z face(5,11,14,15,18)................................
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(11,dvcSendList_Z,sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(14,dvcSendList_Z,2*sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(15,dvcSendList_Z,3*sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	ScaLBL_D3Q19_Pack(18,dvcSendList_Z,4*sendCount_Z,sendCount_Z,sendbuf_Z,dist,N);
	
	req1[5] = MPI_COMM_SCALBL.Isend(sendbuf_Z, 5*sendCount_Z,rank_Z,sendtag);
	req2[5] = MPI_COMM_SCALBL.Irecv(recvbuf_z, 5*recvCount_z,rank_z,recvtag);
	
	//...Pack the xy edge (8)................................
	ScaLBL_D3Q19_Pack(8,dvcSendList_xy,0,sendCount_xy,sendbuf_xy,dist,N);
	req1[6] = MPI_COMM_SCALBL.Isend(sendbuf_xy, sendCount_xy,rank_xy,sendtag);
	req2[6] = MPI_COMM_SCALBL.Irecv(recvbuf_XY, recvCount_XY,rank_XY,recvtag);
	//...Pack the Xy edge (9)................................
	ScaLBL_D3Q19_Pack(9,dvcSendList_Xy,0,sendCount_Xy,sendbuf_Xy,dist,N);
	req1[8] = MPI_COMM_SCALBL.Isend(sendbuf_Xy, sendCount_Xy,rank_Xy,sendtag);
	req2[8] = MPI_COMM_SCALBL.Irecv(recvbuf_xY, recvCount_xY,rank_xY,recvtag);
	//...Pack the xY edge (10)................................
	ScaLBL_D3Q19_Pack(10,dvcSendList_xY,0,sendCount_xY,sendbuf_xY,dist,N);
	req1[9] = MPI_COMM_SCALBL.Isend(sendbuf_xY, sendCount_xY,rank_xY,sendtag);
	req2[9] = MPI_COMM_SCALBL.Irecv(recvbuf_Xy, recvCount_Xy,rank_Xy,recvtag);
	//...Pack the XY edge (7)................................
	ScaLBL_D3Q19_Pack(7,dvcSendList_XY,0,sendCount_XY,sendbuf_XY,dist,N);
	req1[7] = MPI_COMM_SCALBL.Isend(sendbuf_XY, sendCount_XY,rank_XY,sendtag);
	req2[7] = MPI_COMM_SCALBL.Irecv(recvbuf_xy, recvCount_xy,rank_xy,recvtag);
	//...Pack the xz edge (12)................................
	ScaLBL_D3Q19_Pack(12,dvcSendList_xz,0,sendCount_xz,sendbuf_xz,dist,N);
	req1[10] = MPI_COMM_SCALBL.Isend(sendbuf_xz, sendCount_xz,rank_xz,sendtag);
	req2[10] = MPI_COMM_SCALBL.Irecv(recvbuf_XZ, recvCount_XZ,rank_XZ,recvtag);
	//...Pack the xZ edge (14)................................
	ScaLBL_D3Q19_Pack(14,dvcSendList_xZ,0,sendCount_xZ,sendbuf_xZ,dist,N);
	req1[13] = MPI_COMM_SCALBL.Isend(sendbuf_xZ, sendCount_xZ,rank_xZ,sendtag);
	req2[13] = MPI_COMM_SCALBL.Irecv(recvbuf_Xz, recvCount_Xz,rank_Xz,recvtag);
	//...Pack the Xz edge (13)................................
	ScaLBL_D3Q19_Pack(13,dvcSendList_Xz,0,sendCount_Xz,sendbuf_Xz,dist,N);
	req1[12] = MPI_COMM_SCALBL.Isend(sendbuf_Xz, sendCount_Xz,rank_Xz,sendtag);
	req2[12] = MPI_COMM_SCALBL.Irecv(recvbuf_xZ, recvCount_xZ,rank_xZ,recvtag);
	//...Pack the XZ edge (11)................................
	ScaLBL_D3Q19_Pack(11,dvcSendList_XZ,0,sendCount_XZ,sendbuf_XZ,dist,N);
	req1[11] = MPI_COMM_SCALBL.Isend(sendbuf_XZ, sendCount_XZ,rank_XZ,sendtag);
	req2[11] = MPI_COMM_SCALBL.Irecv(recvbuf_xz, recvCount_xz,rank_xz,recvtag);
	//...Pack the yz edge (16)................................
	ScaLBL_D3Q19_Pack(16,dvcSendList_yz,0,sendCount_yz,sendbuf_yz,dist,N);
	req1[14] = MPI_COMM_SCALBL.Isend(sendbuf_yz, sendCount_yz,rank_yz,sendtag);
	req2[14] = MPI_COMM_SCALBL.Irecv(recvbuf_YZ, recvCount_YZ,rank_YZ,recvtag);
	//...Pack the yZ edge (18)................................
	ScaLBL_D3Q19_Pack(18,dvcSendList_yZ,0,sendCount_yZ,sendbuf_yZ,dist,N);
	req1[17] = MPI_COMM_SCALBL.Isend(sendbuf_yZ, sendCount_yZ,rank_yZ,sendtag);
	req2[17] = MPI_COMM_SCALBL.Irecv(recvbuf_Yz, recvCount_Yz,rank_Yz,recvtag);
	//...Pack the Yz edge (17)................................
	ScaLBL_D3Q19_Pack(17,dvcSendList_Yz,0,sendCount_Yz,sendbuf_Yz,dist,N);
	req1[16] = MPI_COMM_SCALBL.Isend(sendbuf_Yz, sendCount_Yz,rank_Yz,sendtag);
	req2[16] = MPI_COMM_SCALBL.Irecv(recvbuf_yZ, recvCount_yZ,rank_yZ,recvtag);
	//...Pack the YZ edge (15)................................
	ScaLBL_D3Q19_Pack(15,dvcSendList_YZ,0,sendCount_YZ,sendbuf_YZ,dist,N);
	req1[15] = MPI_COMM_SCALBL.Isend(sendbuf_YZ, sendCount_YZ,rank_YZ,sendtag);
	req2[15] = MPI_COMM_SCALBL.Irecv(recvbuf_yz, recvCount_yz,rank_yz,recvtag);
	//...................................................................................

}

void Membrane::RecvD3Q19AA(double *dist){

	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_COMM_SCALBL.waitAll(18,req1);
	MPI_COMM_SCALBL.waitAll(18,req2);
	ScaLBL_DeviceBarrier();

	//...................................................................................
	// NOTE: AA Routine writes to opposite 
	// Unpack the distributions on the device
	//...................................................................................
	//...Unpacking for x face(2,8,10,12,14)................................
	Membrane_D3Q19_Unpack(2,dvcRecvDist_x, dvcRecvLinks_x,0,linkCount_x[0],recvbuf_x,dist,N);
	Membrane_D3Q19_Unpack(8,dvcRecvDist_x, dvcRecvLinks_x,recvCount_x,linkCount_x[1],recvbuf_x,dist,N);
	Membrane_D3Q19_Unpack(10,dvcRecvDist_x, dvcRecvLinks_x,2*recvCount_x,linkCount_x[2],recvbuf_x,dist,N);
	Membrane_D3Q19_Unpack(12,dvcRecvDist_x, dvcRecvLinks_x,3*recvCount_x,linkCount_x[3],recvbuf_x,dist,N);
	Membrane_D3Q19_Unpack(14,dvcRecvDist_x, dvcRecvLinks_x,4*recvCount_x,linkCount_x[4],recvbuf_x,dist,N);
	//...................................................................................
	//...Packing for X face(1,7,9,11,13)................................
	Membrane_D3Q19_Unpack(1,dvcRecvDist_X, dvcRecvLinks_X,0,linkCount_X[0],recvbuf_X,dist,N);
	Membrane_D3Q19_Unpack(7,dvcRecvDist_X, dvcRecvLinks_X,recvCount_X,linkCount_X[1],recvbuf_X,dist,N);
	Membrane_D3Q19_Unpack(9,dvcRecvDist_X, dvcRecvLinks_X,2*recvCount_X,linkCount_X[2],recvbuf_X,dist,N);
	Membrane_D3Q19_Unpack(11,dvcRecvDist_X, dvcRecvLinks_X,3*recvCount_X,linkCount_X[3],recvbuf_X,dist,N);
	Membrane_D3Q19_Unpack(13,dvcRecvDist_X, dvcRecvLinks_X,4*recvCount_X,linkCount_X[4],recvbuf_X,dist,N);
	//...................................................................................
	//...Packing for y face(4,8,9,16,18).................................
	Membrane_D3Q19_Unpack(4,dvcRecvDist_y, dvcRecvLinks_y,0,linkCount_y[0],recvbuf_y,dist,N);
	Membrane_D3Q19_Unpack(8,dvcRecvDist_y, dvcRecvLinks_y,recvCount_y,linkCount_y[1],recvbuf_y,dist,N);
	Membrane_D3Q19_Unpack(9,dvcRecvDist_y, dvcRecvLinks_y,2*recvCount_y,linkCount_y[2],recvbuf_y,dist,N);
	Membrane_D3Q19_Unpack(16,dvcRecvDist_y, dvcRecvLinks_y,3*recvCount_y,linkCount_y[3],recvbuf_y,dist,N);
	Membrane_D3Q19_Unpack(18,dvcRecvDist_y, dvcRecvLinks_y,4*recvCount_y,linkCount_y[4],recvbuf_y,dist,N);
	//...................................................................................
	//...Packing for Y face(3,7,10,15,17).................................
	Membrane_D3Q19_Unpack(3,dvcRecvDist_Y, dvcRecvLinks_Y,0,linkCount_Y[0],recvbuf_Y,dist,N);
	Membrane_D3Q19_Unpack(7,dvcRecvDist_Y, dvcRecvLinks_Y,recvCount_Y,linkCount_Y[1],recvbuf_Y,dist,N);
	Membrane_D3Q19_Unpack(10,dvcRecvDist_Y, dvcRecvLinks_Y,2*recvCount_Y,linkCount_Y[2],recvbuf_Y,dist,N);
	Membrane_D3Q19_Unpack(15,dvcRecvDist_Y, dvcRecvLinks_Y,3*recvCount_Y,linkCount_Y[3],recvbuf_Y,dist,N);
	Membrane_D3Q19_Unpack(17,dvcRecvDist_Y, dvcRecvLinks_Y,4*recvCount_Y,linkCount_Y[4],recvbuf_Y,dist,N);
	//...................................................................................
	//...Packing for z face(6,12,13,16,17)................................
	Membrane_D3Q19_Unpack(6,dvcRecvDist_z, dvcRecvLinks_z,0,linkCount_z[0],recvbuf_z,dist,N);
	Membrane_D3Q19_Unpack(12,dvcRecvDist_z, dvcRecvLinks_z,recvCount_z,linkCount_z[1],recvbuf_z,dist,N);
	Membrane_D3Q19_Unpack(13,dvcRecvDist_z, dvcRecvLinks_z,2*recvCount_z,linkCount_z[2],recvbuf_z,dist,N);
	Membrane_D3Q19_Unpack(16,dvcRecvDist_z, dvcRecvLinks_z,3*recvCount_z,linkCount_z[3],recvbuf_z,dist,N);
	Membrane_D3Q19_Unpack(17,dvcRecvDist_z, dvcRecvLinks_z,4*recvCount_z,linkCount_z[4],recvbuf_z,dist,N);
	//...Packing for Z face(5,11,14,15,18)................................
	Membrane_D3Q19_Unpack(5,dvcRecvDist_Z, dvcRecvLinks_Z,0,linkCount_Z[0],recvbuf_Z,dist,N);
	Membrane_D3Q19_Unpack(11,dvcRecvDist_Z, dvcRecvLinks_Z,recvCount_Z,linkCount_Z[1],recvbuf_Z,dist,N);
	Membrane_D3Q19_Unpack(14,dvcRecvDist_Z, dvcRecvLinks_Z,2*recvCount_Z,linkCount_Z[2],recvbuf_Z,dist,N);
	Membrane_D3Q19_Unpack(15,dvcRecvDist_Z, dvcRecvLinks_Z,3*recvCount_Z,linkCount_Z[3],recvbuf_Z,dist,N);
	Membrane_D3Q19_Unpack(18,dvcRecvDist_Z, dvcRecvLinks_Z,4*recvCount_Z,linkCount_Z[4],recvbuf_Z,dist,N);
	//..................................................................................
	//...Pack the xy edge (8)...............................
	Membrane_D3Q19_Unpack(8,dvcRecvDist_xy, dvcRecvLinks_xy,0,recvCount_xy,recvbuf_xy,dist,N);
	//...Pack the Xy edge (9)................................
	Membrane_D3Q19_Unpack(9,dvcRecvDist_Xy, dvcRecvLinks_Xy,0,recvCount_Xy,recvbuf_Xy,dist,N);
	//...Pack the xY edge (10)................................
	Membrane_D3Q19_Unpack(10,dvcRecvDist_xY, dvcRecvLinks_xY,0,recvCount_xY,recvbuf_xY,dist,N);
	//...Pack the XY edge (7)................................
	Membrane_D3Q19_Unpack(7,dvcRecvDist_XY, dvcRecvLinks_XY,0,recvCount_XY,recvbuf_XY,dist,N);
	//...Pack the xz edge (12)................................
	Membrane_D3Q19_Unpack(12,dvcRecvDist_xz, dvcRecvLinks_xz,0,recvCount_xz,recvbuf_xz,dist,N);
	//...Pack the xZ edge (14)................................
	Membrane_D3Q19_Unpack(14,dvcRecvDist_xZ, dvcRecvLinks_xZ,0,recvCount_xZ,recvbuf_xZ,dist,N);
	//...Pack the Xz edge (13)................................
	Membrane_D3Q19_Unpack(13,dvcRecvDist_Xz, dvcRecvLinks_Xz,0,recvCount_Xz,recvbuf_Xz,dist,N);
	//...Pack the XZ edge (11)................................
	Membrane_D3Q19_Unpack(11,dvcRecvDist_XZ, dvcRecvLinks_XZ,0,recvCount_XZ,recvbuf_XZ,dist,N);
	//...Pack the yz edge (16)................................
	Membrane_D3Q19_Unpack(16,dvcRecvDist_yz, dvcRecvLinks_yz,0,recvCount_yz,recvbuf_yz,dist,N);
	//...Pack the yZ edge (18)................................
	Membrane_D3Q19_Unpack(18,dvcRecvDist_yZ, dvcRecvLinks_yZ,0,recvCount_yZ,recvbuf_yZ,dist,N);
	//...Pack the Yz edge (17)................................
	Membrane_D3Q19_Unpack(17,dvcRecvDist_Yz, dvcRecvLinks_Yz,0,recvCount_Yz,recvbuf_Yz,dist,N);
	//...Pack the YZ edge (15)................................
	Membrane_D3Q19_Unpack(15,dvcRecvDist_YZ, dvcRecvLinks_YZ,0,recvCount_YZ,recvbuf_YZ,dist,N);
	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................
}

void Membrane::SendD3Q7AA(double *dist){

	if (Lock==true){
		ERROR("Membrane Error (SendD3Q7): Membrane communicator is locked -- did you forget to match Send/Recv calls?");
	}
	else{
		Lock=true;
	}
	// assign tag of 37 to D3Q7 communication
	sendtag = recvtag = 37;
	ScaLBL_DeviceBarrier();
	// Pack the distributions
	//...Packing for x face(q=2)................................
	ScaLBL_D3Q19_Pack(2,dvcSendList_x,0,sendCount_x,sendbuf_x,dist,Np);
	req1[0] = MPI_COMM_SCALBL.Isend(sendbuf_x, sendCount_x,rank_x,sendtag);
	req2[0] = MPI_COMM_SCALBL.Irecv(recvbuf_X, recvCount_X,rank_X,recvtag);
	//...Packing for X face(q=1)................................
	ScaLBL_D3Q19_Pack(1,dvcSendList_X,0,sendCount_X,sendbuf_X,dist,Np);
	req1[1] = MPI_COMM_SCALBL.Isend(sendbuf_X, sendCount_X,rank_X,sendtag);
	req2[1] = MPI_COMM_SCALBL.Irecv(recvbuf_x, recvCount_x,rank_x,recvtag);
	//for (int idx=0; idx<sendCount_X; idx++) printf(" SendX(%i)=%e \n",idx,sendbuf_X[idx]);
	//...Packing for y face(q=4).................................
	ScaLBL_D3Q19_Pack(4,dvcSendList_y,0,sendCount_y,sendbuf_y,dist,Np);	
	req1[2] = MPI_COMM_SCALBL.Isend(sendbuf_y, sendCount_y,rank_y,sendtag);
	req2[2] = MPI_COMM_SCALBL.Irecv(recvbuf_Y, recvCount_Y,rank_Y,recvtag);
	//...Packing for Y face(q=3).................................
	ScaLBL_D3Q19_Pack(3,dvcSendList_Y,0,sendCount_Y,sendbuf_Y,dist,Np);
	req1[3] = MPI_COMM_SCALBL.Isend(sendbuf_Y, sendCount_Y,rank_Y,sendtag);
	req2[3] = MPI_COMM_SCALBL.Irecv(recvbuf_y, recvCount_y,rank_y,recvtag);
	//for (int idx=0; idx<sendCount_Y; idx++) printf(" SendY(%i)=%e \n",idx,sendbuf_Y[idx]);
	//...Packing for z face(q=6)................................
	ScaLBL_D3Q19_Pack(6,dvcSendList_z,0,sendCount_z,sendbuf_z,dist,Np);
	req1[4] = MPI_COMM_SCALBL.Isend(sendbuf_z, sendCount_z,rank_z,sendtag);
	req2[4] = MPI_COMM_SCALBL.Irecv(recvbuf_Z, recvCount_Z,rank_Z,recvtag);
	//...Packing for Z face(q=5)................................
	ScaLBL_D3Q19_Pack(5,dvcSendList_Z,0,sendCount_Z,sendbuf_Z,dist,Np);
	req1[5] = MPI_COMM_SCALBL.Isend(sendbuf_Z, sendCount_Z,rank_Z,sendtag);
	req2[5] = MPI_COMM_SCALBL.Irecv(recvbuf_z, recvCount_z,rank_z,recvtag); 

}

void Membrane::RecvD3Q7AA(double *dist){

	//...................................................................................
	// Wait for completion of D3Q19 communication
	MPI_COMM_SCALBL.waitAll(6,req1);
	MPI_COMM_SCALBL.waitAll(6,req2);
	ScaLBL_DeviceBarrier();
	//...................................................................................
	// NOTE: AA Routine writes to opposite 
	// Unpack the distributions on the device
	//...................................................................................
	//...Unpacking for x face(q=2)................................
	ScaLBL_D3Q7_Membrane_Unpack(2,dvcRecvDist_x, dvcRecvLinks_x,0,linkCount_x[0],recvCount_x,recvbuf_x,dist,Np,coefficient_x);
	//...................................................................................
	//...Packing for X face(q=1)................................
	ScaLBL_D3Q7_Membrane_Unpack(1,dvcRecvDist_X, dvcRecvLinks_X,0,linkCount_X[0],recvCount_X,recvbuf_X,dist,Np,coefficient_X);
	//...................................................................................
	//...Packing for y face(q=4).................................
	ScaLBL_D3Q7_Membrane_Unpack(4,dvcRecvDist_y, dvcRecvLinks_y,0,linkCount_y[0],recvCount_y,recvbuf_y,dist,Np,coefficient_y);
	//...................................................................................
	//...Packing for Y face(q=3).................................
	ScaLBL_D3Q7_Membrane_Unpack(3,dvcRecvDist_Y, dvcRecvLinks_Y,0,linkCount_Y[0],recvCount_Y,recvbuf_Y,dist,Np,coefficient_Y);
	//...................................................................................
	//...Packing for z face(q=6)................................
	ScaLBL_D3Q7_Membrane_Unpack(6,dvcRecvDist_z, dvcRecvLinks_z,0,linkCount_z[0],recvCount_z,recvbuf_z,dist,Np,coefficient_z);
	//...Packing for Z face(q=5)................................
	ScaLBL_D3Q7_Membrane_Unpack(5,dvcRecvDist_Z, dvcRecvLinks_Z,0,linkCount_Z[0],recvCount_Z,recvbuf_Z,dist,Np,coefficient_Z);
	//..................................................................................

	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................

}

//	std::shared_ptr<Database> db){
void Membrane::AssignCoefficients(int *Map, double *Psi, string method){
	/* Assign mass transfer coefficients to the membrane data structure */
	
	double Threshold;
	double MassFractionIn,MassFractionOut,ThresholdMassFractionIn,ThresholdMassFractionOut;
	
	Threshold = -55.0;
	MassFractionIn = 0.0;
	MassFractionOut = 0.0;
	ThresholdMassFractionOut = 0.0;				
	ThresholdMassFractionIn = 0.0;
	
	if (method == "Voltage Gated Potassium"){
		MassFractionIn = 0.0;
		MassFractionOut = 0.0;
		ThresholdMassFractionOut = 0.0;				
		ThresholdMassFractionIn = 1.0;
	}
	ScaLBL_D3Q7_Membrane_AssignLinkCoef(MembraneLinks, Map, MembraneDistance, Psi, MembraneCoef,
			Threshold, MassFractionIn, MassFractionOut, ThresholdMassFractionIn, ThresholdMassFractionOut,
			membraneLinkCount, Nx, Ny, Nz, Np);

	ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(-1,0,0,Map,MembraneDistance,Psi,Threshold,
			MassFractionIn,MassFractionOut,ThresholdMassFractionIn,ThresholdMassFractionOut,
			dvcRecvDist_X,dvcRecvLinks_X,coefficient_X,0,linkCount_X[0],recvCount_X,
			Np,Nx,Ny,Nz);

	ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(1,0,0,Map,MembraneDistance,Psi,Threshold,
			MassFractionIn,MassFractionOut,ThresholdMassFractionIn,ThresholdMassFractionOut,
			dvcRecvDist_x,dvcRecvLinks_x,coefficient_x,0,linkCount_x[0],recvCount_x,
			Np,Nx,Ny,Nz);
	
	ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(0,-1,0,Map,MembraneDistance,Psi,Threshold,
			MassFractionIn,MassFractionOut,ThresholdMassFractionIn,ThresholdMassFractionOut,
			dvcRecvDist_Y,dvcRecvLinks_Y,coefficient_Y,0,linkCount_Y[0],recvCount_Y,
			Np,Nx,Ny,Nz);
	
	ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(0,1,0,Map,MembraneDistance,Psi,Threshold,
			MassFractionIn,MassFractionOut,ThresholdMassFractionIn,ThresholdMassFractionOut,
			dvcRecvDist_y,dvcRecvLinks_y,coefficient_y,0,linkCount_y[0],recvCount_y,
			Np,Nx,Ny,Nz);
	
	ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(0,0,-1,Map,MembraneDistance,Psi,Threshold,
			MassFractionIn,MassFractionOut,ThresholdMassFractionIn,ThresholdMassFractionOut,
			dvcRecvDist_Z,dvcRecvLinks_Z,coefficient_Z,0,linkCount_Z[0],recvCount_Z,
			Np,Nx,Ny,Nz);
	
	ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(0,0,1,Map,MembraneDistance,Psi,Threshold,
			MassFractionIn,MassFractionOut,ThresholdMassFractionIn,ThresholdMassFractionOut,
			dvcRecvDist_z,dvcRecvLinks_z,coefficient_z,0,linkCount_z[0],recvCount_z,
			Np,Nx,Ny,Nz);
	
}
