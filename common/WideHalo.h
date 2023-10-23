/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

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
/* 
This class implements support for halo widths larger than 1
 */
#ifndef WideHalo_H
#define WideHalo_H
#include "common/ScaLBL.h"
#include "common/MPI.h"

class ScaLBLWideHalo_Communicator {
public:
    //......................................................................................
    ScaLBLWideHalo_Communicator(std::shared_ptr<Domain> Dm, int width);
    ~ScaLBLWideHalo_Communicator();
    //......................................................................................
    //MPI_Comm MPI_COMM_SCALBL;		// MPI Communicator
    Utilities::MPI MPI_COMM_SCALBL;
    unsigned long int CommunicationCount, SendCount, RecvCount;
    int Nx, Ny, Nz, N;     // original domain structure
    int Nxh, Nyh, Nzh, Nh; // with wide halo
    DoubleArray Map;       // map to regular halo
    int first_interior, last_interior;
    //......................................................................................
    //  Set up for D3Q19 distributions -- all 27 neighbors are needed
    //......................................................................................
    // Buffers to store data sent and recieved by this MPI process
    double *sendbuf_x, *sendbuf_y, *sendbuf_z, *sendbuf_X, *sendbuf_Y,
        *sendbuf_Z;
    double *sendbuf_xy, *sendbuf_yz, *sendbuf_xz, *sendbuf_Xy, *sendbuf_Yz,
        *sendbuf_xZ;
    double *sendbuf_xY, *sendbuf_yZ, *sendbuf_Xz, *sendbuf_XY, *sendbuf_YZ,
        *sendbuf_XZ;
    double *sendbuf_xyz, *sendbuf_Xyz, *sendbuf_xYz, *sendbuf_XYz;
    double *sendbuf_xyZ, *sendbuf_XyZ, *sendbuf_xYZ, *sendbuf_XYZ;
    double *recvbuf_x, *recvbuf_y, *recvbuf_z, *recvbuf_X, *recvbuf_Y,
        *recvbuf_Z;
    double *recvbuf_xy, *recvbuf_yz, *recvbuf_xz, *recvbuf_Xy, *recvbuf_Yz,
        *recvbuf_xZ;
    double *recvbuf_xY, *recvbuf_yZ, *recvbuf_Xz, *recvbuf_XY, *recvbuf_YZ,
        *recvbuf_XZ;
    double *recvbuf_xyz, *recvbuf_Xyz, *recvbuf_xYz, *recvbuf_XYz;
    double *recvbuf_xyZ, *recvbuf_XyZ, *recvbuf_xYZ, *recvbuf_XYZ;
    //......................................................................................
    int LastExterior();
    int FirstInterior();
    int LastInterior();

    void Send(double *data);
    void Recv(double *data);

    // Debugging and unit testing functions
    void PrintDebug();

private:
    bool
        Lock; // use Lock to make sure only one call at a time to protect data in transit
    // only one set of Send requests can be active at any time (per instance)
    int i, j, k, n;
    int iproc, jproc, kproc;
    int nprocx, nprocy, nprocz;
    int sendtag, recvtag;
    // Give the object it's own MPI communicator
    RankInfoStruct rank_info;
    MPI_Request req1[26], req2[26];
    //......................................................................................
    // MPI ranks for all 18 neighbors
    //......................................................................................
    // These variables are all private to prevent external things from modifying them!!
    //......................................................................................
    int rank;
    int rank_x, rank_y, rank_z, rank_X, rank_Y, rank_Z;
    int rank_xy, rank_XY, rank_xY, rank_Xy;
    int rank_xz, rank_XZ, rank_xZ, rank_Xz;
    int rank_yz, rank_YZ, rank_yZ, rank_Yz;
    int rank_xyz, rank_Xyz, rank_xYz, rank_XYz;
    int rank_xyZ, rank_XyZ, rank_xYZ, rank_XYZ;
    //......................................................................................
    //......................................................................................
    int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y,
        sendCount_Z;
    int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz,
        sendCount_xZ;
    int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ,
        sendCount_XZ;
    int sendCount_xyz, sendCount_Xyz, sendCount_xYz, sendCount_XYz;
    int sendCount_xyZ, sendCount_XyZ, sendCount_xYZ, sendCount_XYZ;
    //......................................................................................
    int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y,
        recvCount_Z;
    int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz,
        recvCount_xZ;
    int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ,
        recvCount_XZ;
    int recvCount_xyz, recvCount_Xyz, recvCount_xYz, recvCount_XYz;
    int recvCount_xyZ, recvCount_XyZ, recvCount_xYZ, recvCount_XYZ;
    //......................................................................................
    // Send buffers that reside on the compute device
    int *dvcSendList_x, *dvcSendList_y, *dvcSendList_z, *dvcSendList_X,
        *dvcSendList_Y, *dvcSendList_Z;
    int *dvcSendList_xy, *dvcSendList_yz, *dvcSendList_xz, *dvcSendList_Xy,
        *dvcSendList_Yz, *dvcSendList_xZ;
    int *dvcSendList_xY, *dvcSendList_yZ, *dvcSendList_Xz, *dvcSendList_XY,
        *dvcSendList_YZ, *dvcSendList_XZ;
    int *dvcSendList_xyz, *dvcSendList_Xyz, *dvcSendList_xYz, *dvcSendList_XYz;
    int *dvcSendList_xyZ, *dvcSendList_XyZ, *dvcSendList_xYZ, *dvcSendList_XYZ;
    // Recieve buffers that reside on the compute device
    int *dvcRecvList_x, *dvcRecvList_y, *dvcRecvList_z, *dvcRecvList_X,
        *dvcRecvList_Y, *dvcRecvList_Z;
    int *dvcRecvList_xy, *dvcRecvList_yz, *dvcRecvList_xz, *dvcRecvList_Xy,
        *dvcRecvList_Yz, *dvcRecvList_xZ;
    int *dvcRecvList_xY, *dvcRecvList_yZ, *dvcRecvList_Xz, *dvcRecvList_XY,
        *dvcRecvList_YZ, *dvcRecvList_XZ;
    int *dvcRecvList_xyz, *dvcRecvList_Xyz, *dvcRecvList_xYz, *dvcRecvList_XYz;
    int *dvcRecvList_xyZ, *dvcRecvList_XyZ, *dvcRecvList_xYZ, *dvcRecvList_XYZ;
    //......................................................................................

    inline int getHaloBlock(int imin, int imax, int jmin, int jmax, int kmin,
                            int kmax, int *&dvcList) {
        int count = 0;
        int *List;
        List = new int[(imax - imin) * (jmax - jmin) * (kmax - kmin)];
        for (k = kmin; k < kmax; k++) {
            for (j = jmin; j < jmax; j++) {
                for (i = imin; i < imax; i++) {
                    List[count++] = k * Nxh * Nyh + j * Nxh + i;
                }
            }
        }
        size_t numbytes = count * sizeof(int);
        ScaLBL_AllocateZeroCopy((void **)&dvcList,
                                numbytes); // Allocate device memory
        ScaLBL_CopyToZeroCopy(dvcList, List, numbytes);
        return count;
    }
};
#endif
