#ifndef COMMUNICATION_H_INC
#define COMMUNICATION_H_INC

// ********** COMMUNICTION **************************************
/*
 //..............................................................
 // Communication helper routines for MPI
 //..............................................................
 */

using namespace std;


//***************************************************************************************
inline void PackMeshData(int *list, int count, double *sendbuf, DoubleArray &Values){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;
	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = Values.data[n];
	}
}
inline void UnpackMeshData(int *list, int count, double *recvbuf, DoubleArray &Values){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		Values.data[n] = recvbuf[idx];
	}
}


//***************************************************************************************
inline void CommunicateMeshHalo(DoubleArray &MeshData, MPI_Comm Communicator,
		double *sendbuf_x,double *sendbuf_y,double *sendbuf_z,double *sendbuf_X,double *sendbuf_Y,double *sendbuf_Z,
		double *sendbuf_xy,double *sendbuf_XY,double *sendbuf_xY,double *sendbuf_Xy,
		double *sendbuf_xz,double *sendbuf_XZ,double *sendbuf_xZ,double *sendbuf_Xz,
		double *sendbuf_yz,double *sendbuf_YZ,double *sendbuf_yZ,double *sendbuf_Yz,
		double *recvbuf_x,double *recvbuf_y,double *recvbuf_z,double *recvbuf_X,double *recvbuf_Y,double *recvbuf_Z,
		double *recvbuf_xy,double *recvbuf_XY,double *recvbuf_xY,double *recvbuf_Xy,
		double *recvbuf_xz,double *recvbuf_XZ,double *recvbuf_xZ,double *recvbuf_Xz,
		double *recvbuf_yz,double *recvbuf_YZ,double *recvbuf_yZ,double *recvbuf_Yz,
		int *sendList_x,int *sendList_y,int *sendList_z,int *sendList_X,int *sendList_Y,int *sendList_Z,
		int *sendList_xy,int *sendList_XY,int *sendList_xY,int *sendList_Xy,
		int *sendList_xz,int *sendList_XZ,int *sendList_xZ,int *sendList_Xz,
		int *sendList_yz,int *sendList_YZ,int *sendList_yZ,int *sendList_Yz,
		int sendCount_x,int sendCount_y,int sendCount_z,int sendCount_X,int sendCount_Y,int sendCount_Z,
		int sendCount_xy,int sendCount_XY,int sendCount_xY,int sendCount_Xy,
		int sendCount_xz,int sendCount_XZ,int sendCount_xZ,int sendCount_Xz,
		int sendCount_yz,int sendCount_YZ,int sendCount_yZ,int sendCount_Yz,
		int *recvList_x,int *recvList_y,int *recvList_z,int *recvList_X,int *recvList_Y,int *recvList_Z,
		int *recvList_xy,int *recvList_XY,int *recvList_xY,int *recvList_Xy,
		int *recvList_xz,int *recvList_XZ,int *recvList_xZ,int *recvList_Xz,
		int *recvList_yz,int *recvList_YZ,int *recvList_yZ,int *recvList_Yz,
		int recvCount_x,int recvCount_y,int recvCount_z,int recvCount_X,int recvCount_Y,int recvCount_Z,
		int recvCount_xy,int recvCount_XY,int recvCount_xY,int recvCount_Xy,
		int recvCount_xz,int recvCount_XZ,int recvCount_xZ,int recvCount_Xz,
		int recvCount_yz,int recvCount_YZ,int recvCount_yZ,int recvCount_Yz,
		int rank_x,int rank_y,int rank_z,int rank_X,int rank_Y,int rank_Z,int rank_xy,int rank_XY,int rank_xY,
		int rank_Xy,int rank_xz,int rank_XZ,int rank_xZ,int rank_Xz,int rank_yz,int rank_YZ,int rank_yZ,int rank_Yz)
{
	int sendtag, recvtag;
	sendtag = recvtag = 7;
	PackMeshData(sendList_x, sendCount_x ,sendbuf_x, MeshData);
	PackMeshData(sendList_X, sendCount_X ,sendbuf_X, MeshData);
	PackMeshData(sendList_y, sendCount_y ,sendbuf_y, MeshData);
	PackMeshData(sendList_Y, sendCount_Y ,sendbuf_Y, MeshData);
	PackMeshData(sendList_z, sendCount_z ,sendbuf_z, MeshData);
	PackMeshData(sendList_Z, sendCount_Z ,sendbuf_Z, MeshData);
	PackMeshData(sendList_xy, sendCount_xy ,sendbuf_xy, MeshData);
	PackMeshData(sendList_Xy, sendCount_Xy ,sendbuf_Xy, MeshData);
	PackMeshData(sendList_xY, sendCount_xY ,sendbuf_xY, MeshData);
	PackMeshData(sendList_XY, sendCount_XY ,sendbuf_XY, MeshData);
	PackMeshData(sendList_xz, sendCount_xz ,sendbuf_xz, MeshData);
	PackMeshData(sendList_Xz, sendCount_Xz ,sendbuf_Xz, MeshData);
	PackMeshData(sendList_xZ, sendCount_xZ ,sendbuf_xZ, MeshData);
	PackMeshData(sendList_XZ, sendCount_XZ ,sendbuf_XZ, MeshData);
	PackMeshData(sendList_yz, sendCount_yz ,sendbuf_yz, MeshData);
	PackMeshData(sendList_Yz, sendCount_Yz ,sendbuf_Yz, MeshData);
	PackMeshData(sendList_yZ, sendCount_yZ ,sendbuf_yZ, MeshData);
	PackMeshData(sendList_YZ, sendCount_YZ ,sendbuf_YZ, MeshData);
	//......................................................................................
	MPI_Sendrecv(sendbuf_x,sendCount_x,MPI_CHAR,rank_x,sendtag,
			recvbuf_X,recvCount_X,MPI_CHAR,rank_X,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_X,sendCount_X,MPI_CHAR,rank_X,sendtag,
			recvbuf_x,recvCount_x,MPI_CHAR,rank_x,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_y,sendCount_y,MPI_CHAR,rank_y,sendtag,
			recvbuf_Y,recvCount_Y,MPI_CHAR,rank_Y,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_Y,sendCount_Y,MPI_CHAR,rank_Y,sendtag,
			recvbuf_y,recvCount_y,MPI_CHAR,rank_y,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_z,sendCount_z,MPI_CHAR,rank_z,sendtag,
			recvbuf_Z,recvCount_Z,MPI_CHAR,rank_Z,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_Z,sendCount_Z,MPI_CHAR,rank_Z,sendtag,
			recvbuf_z,recvCount_z,MPI_CHAR,rank_z,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_xy,sendCount_xy,MPI_CHAR,rank_xy,sendtag,
			recvbuf_XY,recvCount_XY,MPI_CHAR,rank_XY,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_XY,sendCount_XY,MPI_CHAR,rank_XY,sendtag,
			recvbuf_xy,recvCount_xy,MPI_CHAR,rank_xy,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_Xy,sendCount_Xy,MPI_CHAR,rank_Xy,sendtag,
			recvbuf_xY,recvCount_xY,MPI_CHAR,rank_xY,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_xY,sendCount_xY,MPI_CHAR,rank_xY,sendtag,
			recvbuf_Xy,recvCount_Xy,MPI_CHAR,rank_Xy,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_xz,sendCount_xz,MPI_CHAR,rank_xz,sendtag,
			recvbuf_XZ,recvCount_XZ,MPI_CHAR,rank_XZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_XZ,sendCount_XZ,MPI_CHAR,rank_XZ,sendtag,
			recvbuf_xz,recvCount_xz,MPI_CHAR,rank_xz,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_Xz,sendCount_Xz,MPI_CHAR,rank_Xz,sendtag,
			recvbuf_xZ,recvCount_xZ,MPI_CHAR,rank_xZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_xZ,sendCount_xZ,MPI_CHAR,rank_xZ,sendtag,
			recvbuf_Xz,recvCount_Xz,MPI_CHAR,rank_Xz,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_yz,sendCount_yz,MPI_CHAR,rank_yz,sendtag,
			recvbuf_YZ,recvCount_YZ,MPI_CHAR,rank_YZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_YZ,sendCount_YZ,MPI_CHAR,rank_YZ,sendtag,
			recvbuf_yz,recvCount_yz,MPI_CHAR,rank_yz,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_Yz,sendCount_Yz,MPI_CHAR,rank_Yz,sendtag,
			recvbuf_yZ,recvCount_yZ,MPI_CHAR,rank_yZ,recvtag,Communicator,MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbuf_yZ,sendCount_yZ,MPI_CHAR,rank_yZ,sendtag,
			recvbuf_Yz,recvCount_Yz,MPI_CHAR,rank_Yz,recvtag,Communicator,MPI_STATUS_IGNORE);
	//........................................................................................
	UnpackMeshData(recvList_x, recvCount_x ,recvbuf_x, MeshData);
	UnpackMeshData(recvList_X, recvCount_X ,recvbuf_X, MeshData);
	UnpackMeshData(recvList_y, recvCount_y ,recvbuf_y, MeshData);
	UnpackMeshData(recvList_Y, recvCount_Y ,recvbuf_Y, MeshData);
	UnpackMeshData(recvList_z, recvCount_z ,recvbuf_z, MeshData);
	UnpackMeshData(recvList_Z, recvCount_Z ,recvbuf_Z, MeshData);
	UnpackMeshData(recvList_xy, recvCount_xy ,recvbuf_xy, MeshData);
	UnpackMeshData(recvList_Xy, recvCount_Xy ,recvbuf_Xy, MeshData);
	UnpackMeshData(recvList_xY, recvCount_xY ,recvbuf_xY, MeshData);
	UnpackMeshData(recvList_XY, recvCount_XY ,recvbuf_XY, MeshData);
	UnpackMeshData(recvList_xz, recvCount_xz ,recvbuf_xz, MeshData);
	UnpackMeshData(recvList_Xz, recvCount_Xz ,recvbuf_Xz, MeshData);
	UnpackMeshData(recvList_xZ, recvCount_xZ ,recvbuf_xZ, MeshData);
	UnpackMeshData(recvList_XZ, recvCount_XZ ,recvbuf_XZ, MeshData);
	UnpackMeshData(recvList_yz, recvCount_yz ,recvbuf_yz, MeshData);
	UnpackMeshData(recvList_Yz, recvCount_Yz ,recvbuf_Yz, MeshData);
	UnpackMeshData(recvList_yZ, recvCount_yZ ,recvbuf_yZ, MeshData);
	UnpackMeshData(recvList_YZ, recvCount_YZ ,recvbuf_YZ, MeshData);
}



#endif

