/* 
This class implements support for halo widths larger than 1
 */

ScaLBLWideHalo_Communicator::ScaLBLWideHalo_Communicator(std::shared_ptr <Domain> Dm, int width)
{
	//......................................................................................
	Lock=false; // unlock the communicator
	//......................................................................................
	// Create a separate copy of the communicator for the device
	MPI_Comm_dup(Dm->Comm,&MPI_COMM_SCALBL);
	//......................................................................................
	// Copy the domain size and communication information directly from Dm
	Nx = Dm->Nx;
	Ny = Dm->Ny;
	Nz = Dm->Nz;
	N = Nx*Ny*Nz;
	Nxh = Nx + 2*(width - 1);
	Nyh = Ny + 2*(width - 1);
	Nzh = Nz + 2*(width - 1);
	Nh = Nxh*Nyh*Nzh;
	next=0;
	
	rank=Dm->rank();
	iproc = Dm->iproc();
	jproc = Dm->jproc();
	kproc = Dm->kproc();
	nprocx = Dm->nprocx();
	nprocy = Dm->nprocy();
	nprocz = Dm->nprocz();
	rank_info = RankInfoStruct(myrank,nprocx,nprocy,nprocz);
    rank = rank_info.rank[1][1][1]; 
    rank_X = rank_info.rank[2][1][1]; 
    rank_x = rank_info.rank[0][1][1]; 
    rank_Y = rank_info.rank[1][2][1]; 
    rank_y = rank_info.rank[1][0][1]; 
    rank_Z = rank_info.rank[1][1][2]; 
    rank_z = rank_info.rank[1][1][0]; 
    rank_XY = rank_info.rank[2][2][1]; 
    rank_xy = rank_info.rank[0][0][1]; 
    rank_Xy = rank_info.rank[2][0][1]; 
    rank_xY = rank_info.rank[0][2][1]; 
    rank_XZ = rank_info.rank[2][1][2]; 
    rank_xz = rank_info.rank[0][1][0]; 
    rank_Xz = rank_info.rank[2][1][0]; 
    rank_xZ = rank_info.rank[0][1][2]; 
    rank_YZ = rank_info.rank[1][2][2]; 
    rank_yz = rank_info.rank[1][0][0]; 
    rank_Yz = rank_info.rank[1][2][0]; 
    rank_yZ = rank_info.rank[1][0][2]; 
    rank_XYz = rank_info.rank[2][2][0]; 
    rank_xyz = rank_info.rank[0][0][0]; 
    rank_Xyz = rank_info.rank[2][0][0]; 
    rank_xYz = rank_info.rank[0][2][0]; 
    rank_XYZ = rank_info.rank[2][2][2]; 
    rank_xyZ = rank_info.rank[0][0][2]; 
    rank_XyZ = rank_info.rank[2][0][2]; 
    rank_xYZ = rank_info.rank[0][2][2]; 

	sendCount_x = (Ny-2)*(Nz-2)*width;
	sendCount_y = (Nx-2)*(Nz-2)*width;
	sendCount_z = (Nx-2)*(Ny-2)*width;
	sendCount_X = (Ny-2)*(Nz-2)*width;
	sendCount_Y = (Nx-2)*(Nz-2)*width;
	sendCount_Z = (Nx-2)*(Ny-2)*width;
	sendCount_xy = (Nz-2)*width*width;
	sendCount_yz = (Nx-2)*width*width;
	sendCount_xz = (Ny-2)*width*width;
	sendCount_Xy = (Nz-2)*width*width;
	sendCount_Yz = (Nx-2)*width*width;
	sendCount_xZ = (Ny-2)*width*width;
	sendCount_xY = (Nz-2)*width*width;
	sendCount_yZ = (Nx-2)*width*width;
	sendCount_Xz = (Ny-2)*width*width;
	sendCount_XY = (Nz-2)*width*width;
	sendCount_YZ = (Nx-2)*width*width;
	sendCount_XZ = (Ny-2)*width*width;
	sendCount_xyz = width*width*width;
	sendCount_Xyz = width*width*width;
	sendCount_xYz = width*width*width;
	sendCount_XYz = width*width*width;
	sendCount_xyZ = width*width*width;
	sendCount_XyZ = width*width*width;
	sendCount_xYZ = width*width*width;
	sendCount_XYZ = width*width*width;
	
	RecvCount_x = (Ny-2)*(Nz-2)*width;
	RecvCount_y = (Nx-2)*(Nz-2)*width;
	RecvCount_z = (Nx-2)*(Ny-2)*width;
	RecvCount_X = (Ny-2)*(Nz-2)*width;
	RecvCount_Y = (Nx-2)*(Nz-2)*width;
	RecvCount_Z = (Nx-2)*(Ny-2)*width;
	RecvCount_xy = (Nz-2)*width*width;
	RecvCount_yz = (Nx-2)*width*width;
	RecvCount_xz = (Ny-2)*width*width;
	RecvCount_Xy = (Nz-2)*width*width;
	RecvCount_Yz = (Nx-2)*width*width;
	RecvCount_xZ = (Ny-2)*width*width;
	RecvCount_xY = (Nz-2)*width*width;
	RecvCount_yZ = (Nx-2)*width*width;
	RecvCount_Xz = (Ny-2)*width*width;
	RecvCount_XY = (Nz-2)*width*width;
	RecvCount_YZ = (Nx-2)*width*width;
	RecvCount_XZ = (Ny-2)*width*width;
	RecvCount_xyz = width*width*width;
	RecvCount_Xyz = width*width*width;
	RecvCount_xYz = width*width*width;
	RecvCount_XYz = width*width*width;
	RecvCount_xyZ = width*width*width;
	RecvCount_XyZ = width*width*width;
	RecvCount_xYZ = width*width*width;
	RecvCount_XYZ = width*width*width;
	//......................................................................................
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_x, sendCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_X, sendCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_y, sendCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Y, sendCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_z, sendCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Z, sendCount_Z*sizeof(double));	// Allocate device memory
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
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xyz, sendCount_xyz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xYz, sendCount_xYz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_Xyz, sendCount_Xyz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_XYz, sendCount_XYz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xyZ, sendCount_xyZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_xYZ, sendCount_xYZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_XyZ, sendCount_XyZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &sendbuf_XYZ, sendCount_XYZ*sizeof(double));	// Allocate device memory
	//......................................................................................
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_x, recvCount_x*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_X, recvCount_X*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_y, recvCount_y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Y, recvCount_Y*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_z, recvCount_z*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Z, recvCount_Z*sizeof(double));	// Allocate device memory
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
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xyz, recvCount_xyz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xYz, recvCount_xYz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_Xyz, recvCount_Xyz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_XYz, recvCount_XYz*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xyZ, recvCount_xyZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_xYZ, recvCount_xYZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_XyZ, recvCount_XyZ*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf_XYZ, recvCount_XYZ*sizeof(double));	// Allocate device memory
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
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xyz, sendCount_xyz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xYz, sendCount_xYz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_Xyz, sendCount_Xyz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_XYz, sendCount_XYz*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xyZ, sendCount_xyZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_xYZ, sendCount_xYZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_XyZ, sendCount_XyZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcSendList_XYZ, sendCount_XYZ*sizeof(int));	// Allocate device memory
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
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xyz, recvCount_xyZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xYz, recvCount_xYZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_Xyz, recvCount_XyZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_XYz, recvCount_XYZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xyZ, recvCount_xyZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_xYZ, recvCount_xYZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_XyZ, recvCount_XyZ*sizeof(int));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &dvcRecvList_XYZ, recvCount_XYZ*sizeof(int));	// Allocate device memory
	//......................................................................................

	MPI_Barrier(MPI_COMM_SCALBL);
	
	/*  Fill in communications patterns for the lists */
	int *sendList, *recvList;
	sendList = new int [width*max(Nxh,Nyh)*max(Nxh,Nzh)];
	/* x face */
	int count = 0;
	for (int k=0; k<width; k++){
		for (j=2; j<Nyh-2; j++){
			for (i=2; i<Nxh-2; i++){
				SendList[count] = (k+width)*Nxh*Nyh + j*Nxh + i;
				RecvList[count++] = k*Nxh*Nyh + j*Nxh + i;
			}
		}
	}
	ScaLBL_CopyToZeroCopy(dvcSendList_x,sendList,sendCount_x*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_x,recvList,recvCount_x*sizeof(int));
	/* X face */
	count = 0;
	for (int k=0; k<width; k++){
		for (j=2; j<Nyh-2; j++){
			for (i=2; i<Nxh-2; i++){
				SendList[count] = (Nzh-width-k)*Nxh*Nyh + j*Nxh + i;
				RecvList[count++] = (Nzh-k)*Nxh*Nyh + j*Nxh + i;
			}
		}
	}
	ScaLBL_CopyToZeroCopy(dvcSendList_X,sendList,sendCount_X*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_X,recvList,recvCount_X*sizeof(int));
	/* y face */
	count = 0;
	for (k=2; k<Nzh-2; k++){
		for (int j=0; j<width; j++){
			for (i=2; i<Nxh-2; i++){
				SendList[count] = k*Nxh*Nyh + (j+width)*Nxh + i;
				RecvList[count++] = k*Nxh*Nyh + j*Nxh + i;
			}
		}
	}
	ScaLBL_CopyToZeroCopy(dvcSendList_y,sendList,sendCount_y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_y,recvList,recvCount_y*sizeof(int));
	/* Y face */
	count = 0;
	for (k=2; k<Nzh-2; k++){
		for (int j=0; j<width; j++){
			for (i=2; i<Nxh-2; i++){
				SendList[count] = k*Nxh*Nyh + (Nyh-width-j)*Nxh + i;
				RecvList[count++] = k*Nxh*Nyh + (Nyh-j)*Nxh + i;
			}
		}
	}
	ScaLBL_CopyToZeroCopy(dvcSendList_Y,sendList,sendCount_Y*sizeof(int));
	ScaLBL_CopyToZeroCopy(dvcRecvList_Y,recvList,recvCount_Y*sizeof(int));
	

}

void ScaLBLWideHalo_Communicator::Send(double *data){
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
void ScaLBLWideHalo_Communicator::Recv(double *data){

	//...................................................................................
	MPI_Waitall(26,req1,stat1);
	MPI_Waitall(26,req2,stat2);
	ScaLBL_DeviceBarrier();
	//...................................................................................
	//...................................................................................
	ScaLBL_Scalar_Unpack(dvcRecvList_x, recvCount_x,recvbuf_x, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_y, recvCount_y,recvbuf_y, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_X, recvCount_X,recvbuf_X, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Y, recvCount_Y,recvbuf_Y, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xy, recvCount_xy,recvbuf_xy, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xY, recvCount_xY,recvbuf_xY, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xy, recvCount_Xy,recvbuf_Xy, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XY, recvCount_XY,recvbuf_XY, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_z, recvCount_z,recvbuf_z, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xz, recvCount_xz,recvbuf_xz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Xz, recvCount_Xz,recvbuf_Xz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yz, recvCount_yz,recvbuf_yz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Yz, recvCount_Yz,recvbuf_Yz, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_Z, recvCount_Z,recvbuf_Z, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_xZ, recvCount_xZ,recvbuf_xZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_XZ, recvCount_XZ,recvbuf_XZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_yZ, recvCount_yZ,recvbuf_yZ, data, N);
	ScaLBL_Scalar_Unpack(dvcRecvList_YZ, recvCount_YZ,recvbuf_YZ, data, N);

	//...................................................................................
	Lock=false; // unlock the communicator after communications complete
	//...................................................................................

}

inline int getHaloBlock(){
}
}
