/* Flow adaptor class for multiphase flow methods */

#ifndef ScaLBL_Membrane_INC
#define ScaLBL_Membrane_INC
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"

/**
 * \class Membrane
 * @brief 
 * The Membrane class operates on ScaLBL data structures to insert membrane
 * 
 */

class Membrane {
public:
    int Np;
    int Nx,Ny,Nz,N;
    
    int *neighborList;     // modified neighborlist
    int *membraneLinks;    // D3Q19 links that cross membrane
    double *membraneDist;  // distance to membrane for each linked site

    /**
    * \brief Create a flow adaptor to operate on the LB model
    * @param         ScaLBL - originating data structures 
    * @param 		 neighborList - list of neighbors for each site
    */
    Membrane(std::shared_ptr <Domain> Dm, int *initialNeighborList, int Nsites);

    /**
	* \brief Destructor
	*/
    ~Membrane();

    /**
    * \brief  Create membrane 
    * \details  Create membrane structure from signed distance function
    * @param Dm             - domain structure
    * @param Distance       - signed distance to membrane 
    * @param Map            - mapping between regular layout and compact layout
    */
    int Create(std::shared_ptr <Domain> Dm, DoubleArray &Distance, IntArray &Map);
    
    void SendRecv(double *dist);
	//......................................................................................
	// Buffers to store data sent and recieved by this MPI process
	double *sendbuf_x, *sendbuf_y, *sendbuf_z, *sendbuf_X, *sendbuf_Y, *sendbuf_Z;
	double *sendbuf_xy, *sendbuf_yz, *sendbuf_xz, *sendbuf_Xy, *sendbuf_Yz, *sendbuf_xZ;
	double *sendbuf_xY, *sendbuf_yZ, *sendbuf_Xz, *sendbuf_XY, *sendbuf_YZ, *sendbuf_XZ;
	double *recvbuf_x, *recvbuf_y, *recvbuf_z, *recvbuf_X, *recvbuf_Y, *recvbuf_Z;
	double *recvbuf_xy, *recvbuf_yz, *recvbuf_xz, *recvbuf_Xy, *recvbuf_Yz, *recvbuf_xZ;
	double *recvbuf_xY, *recvbuf_yZ, *recvbuf_Xz, *recvbuf_XY, *recvbuf_YZ, *recvbuf_XZ;
	//......................................................................................
    
private:
    int iproc,jproc,kproc;
    int nprocx,nprocy,nprocz;
	// Give the object it's own MPI communicator
	RankInfoStruct rank_info;
	Utilities::MPI MPI_COMM_SCALBL;		// MPI Communicator for this domain
	MPI_Request req1[18],req2[18];
    /**
    * \brief   Set up membrane communication
    * \details associate p2p communication links to membrane where necessary
    *          returns the number of membrane links
    *          regular communications are stored in the first part of the list
    *          membrane communications are stored in the last part of the list
    * @param Cqx            - discrete velocity (x)
    * @param Cqy            - discrete velocity (y)
    * @param Cqz            - discrete velocity (z)
    * @param list           - list of recieved values
    * @param count          - number recieved values
    * @param Distance       - signed distance to membrane 
    * @param recvDistance   - distance values from neighboring processor
    * @param d3q19_recvlist - device array with the saved list
    * */
    int D3Q19_MapSendRecv(const int Cqx, const int Cqy, const int Cqz, 
    		const int shift_x, const int shift_y, const int shift_z, 
    		const int sendCount, const int recvCount,
    		int rank_q, int rank_Q, int startSend, int startRecv,
    		std::shared_ptr <Domain> Dm, const int *originalSendList, 
    		const DoubleArray &Distance, int *membraneSendList, int *membraneRecvList);
	//......................................................................................
	// MPI ranks for all 18 neighbors
	//......................................................................................
	// These variables are all private to prevent external things from modifying them!!
	//......................................................................................
	int rank;
	int rank_x,rank_y,rank_z,rank_X,rank_Y,rank_Z;
	int rank_xy,rank_XY,rank_xY,rank_Xy;
	int rank_xz,rank_XZ,rank_xZ,rank_Xz;
	int rank_yz,rank_YZ,rank_yZ,rank_Yz;
	//......................................................................................
	//......................................................................................
	int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
	int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
	int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
	//......................................................................................
	int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
	int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
	int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
	//......................................................................................
	// Send buffers that reside on the compute device
	int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X, *dvcSendList_Y, *dvcSendList_Z;
	int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy, *dvcSendList_Yz, *dvcSendList_xZ;
	int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY, *dvcSendList_YZ, *dvcSendList_XZ;
	// Recieve buffers that reside on the compute device
	int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X, *dvcRecvList_Y, *dvcRecvList_Z;
	int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy, *dvcRecvList_Yz, *dvcRecvList_xZ;
	int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY, *dvcRecvList_YZ, *dvcRecvList_XZ;
	// Recieve buffers for the distributions
	int *dvcRecvDist_x, *dvcRecvDist_y, *dvcRecvDist_z, *dvcRecvDist_X, *dvcRecvDist_Y, *dvcRecvDist_Z;
	int *dvcRecvDist_xy, *dvcRecvDist_yz, *dvcRecvDist_xz, *dvcRecvDist_Xy, *dvcRecvDist_Yz, *dvcRecvDist_xZ;
	int *dvcRecvDist_xY, *dvcRecvDist_yZ, *dvcRecvDist_Xz, *dvcRecvDist_XY, *dvcRecvDist_YZ, *dvcRecvDist_XZ;
	//......................................................................................
};
#endif