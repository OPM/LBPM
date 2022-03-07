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
* \brief  Unpack D3Q19 distributions after communication using links determined based on membrane location
* @param q  - index for distribution based on D3Q19 discrete velocity structure
* @param list - list of distributions to communicate
* @param links - list of active links based on the membrane location
* @param start -  index to start parsing the list 
* @param count -  number of values to unppack 
* @param recvbuf - memory buffer where recieved values have been stored
* @param dist - memory buffer to hold the distributions
* @param N - size of the distributions (derived from Domain structure)
*/
extern "C" void Membrane_D3Q19_Unpack(int q, int *list, int *links, int start, int count, double *recvbuf, double *dist, int N);


/**
* \brief Set custom link rules for D3Q19 distribution based on membrane location
* @param q  - index for distribution based on D3Q19 discrete velocity structure
* @param list - list of distributions to communicate
* @param links - list of active links based on the membrane location
* @param coef  - coefficient to determine the local mass transport for each membrane link
* @param start -  index to start parsing the list 
* @param offset - offset to start reading membrane links
* @param count -  number of values to unppack 
* @param recvbuf - memory buffer where recieved values have been stored
* @param dist - memory buffer to hold the distributions
* @param N - size of the distributions (derived from Domain structure)
*/
extern "C" void Membrane_D3Q19_Transport(int q, int *list, int *links, double *coef, int start, int offset, 
		int linkCount, double *recvbuf, double *dist, int N);

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
    int *membraneTag;      // label each link in the membrane
    double *membraneDist;  // distance to membrane for each linked site
    double *membraneCoef;  // mass transport coefficient for the membrane 

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
        
	void SendD3Q19AA(double *dist);
	void RecvD3Q19AA(double *dist);
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
	bool Lock; 	// use Lock to make sure only one call at a time to protect data in transit
	int sendtag, recvtag;
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
    * @param d3q19_recvlist - device array with the saved list
    * @param d3q19_linkList - sorted list with regular and membrane links
    * @param Distance       - signed distance to membrane 
    * @param Map      		- data structure used to define mapping between dense and sparse representation
    * */
    int D3Q19_MapRecv(int Cqx, int Cqy, int Cqz, const int *list, int start, int count,
    		int *d3q19_recvlist, int *d3q19_linkList, DoubleArray &Distance,  IntArray &Map);
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
	int SendCount, RecvCount, CommunicationCount;
	//......................................................................................
	int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
	int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
	int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
	//......................................................................................
	int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
	int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
	int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
	//......................................................................................
	int linkCount_x[5], linkCount_y[5], linkCount_z[5], linkCount_X[5], linkCount_Y[5], linkCount_Z[5];
	int linkCount_xy, linkCount_yz, linkCount_xz, linkCount_Xy, linkCount_Yz, linkCount_xZ;
	int linkCount_xY, linkCount_yZ, linkCount_Xz, linkCount_XY, linkCount_YZ, linkCount_XZ;
	//......................................................................................
	// Send buffers that reside on the compute device
	int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X, *dvcSendList_Y, *dvcSendList_Z;
	int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy, *dvcSendList_Yz, *dvcSendList_xZ;
	int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY, *dvcSendList_YZ, *dvcSendList_XZ;
	// Recieve buffers that reside on the compute device
	int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X, *dvcRecvList_Y, *dvcRecvList_Z;
	int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy, *dvcRecvList_Yz, *dvcRecvList_xZ;
	int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY, *dvcRecvList_YZ, *dvcRecvList_XZ;
	// Link lists  that reside on the compute device
	int *dvcRecvLinks_x, *dvcRecvLinks_y, *dvcRecvLinks_z, *dvcRecvLinks_X, *dvcRecvLinks_Y, *dvcRecvLinks_Z;
	int *dvcRecvLinks_xy, *dvcRecvLinks_yz, *dvcRecvLinks_xz, *dvcRecvLinks_Xy, *dvcRecvLinks_Yz, *dvcRecvLinks_xZ;
	int *dvcRecvLinks_xY, *dvcRecvLinks_yZ, *dvcRecvLinks_Xz, *dvcRecvLinks_XY, *dvcRecvLinks_YZ, *dvcRecvLinks_XZ;
	// Recieve buffers for the distributions
	int *dvcRecvDist_x, *dvcRecvDist_y, *dvcRecvDist_z, *dvcRecvDist_X, *dvcRecvDist_Y, *dvcRecvDist_Z;
	int *dvcRecvDist_xy, *dvcRecvDist_yz, *dvcRecvDist_xz, *dvcRecvDist_Xy, *dvcRecvDist_Yz, *dvcRecvDist_xZ;
	int *dvcRecvDist_xY, *dvcRecvDist_yZ, *dvcRecvDist_Xz, *dvcRecvDist_XY, *dvcRecvDist_YZ, *dvcRecvDist_XZ;
	//......................................................................................
};
#endif