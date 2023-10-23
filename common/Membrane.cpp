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
/* Membrane class for lattice Boltzmann models */

#include "common/Membrane.h"
#include "analysis/distance.h"

Membrane::Membrane(std::shared_ptr<ScaLBL_Communicator> sComm,
                   int *dvcNeighborList, int Nsites) {

    Np = Nsites;
    initialNeighborList = new int[18 * Np];
    ScaLBL_AllocateDeviceMemory((void **)&NeighborList, 18 * Np * sizeof(int));
    Lock = false; // unlock the communicator
        //......................................................................................
        // Create a separate copy of the communicator for the device
    MPI_COMM_SCALBL = sComm->MPI_COMM_SCALBL.dup();
    int myrank = sComm->MPI_COMM_SCALBL.getRank();
    rank_info =
        RankInfoStruct(myrank, rank_info.nx, rank_info.ny, rank_info.nz);

    ScaLBL_CopyToHost(initialNeighborList, dvcNeighborList,
                      18 * Np * sizeof(int));
    sComm->MPI_COMM_SCALBL.barrier();
    ScaLBL_CopyToDevice(NeighborList, initialNeighborList,
                        18 * Np * sizeof(int));

    /* Copy communication lists */
    //......................................................................................
    //Lock=false; // unlock the communicator
    //......................................................................................
    // Create a separate copy of the communicator for the device
    //MPI_COMM_SCALBL = sComm->Comm.dup();
    //......................................................................................
    // Copy the domain size and communication information directly from sComm
    Nx = sComm->Nx;
    Ny = sComm->Ny;
    Nz = sComm->Nz;
    N = Nx * Ny * Nz;
    //next=0;
    rank = sComm->rank;
    rank_x = sComm->rank_x;
    rank_y = sComm->rank_y;
    rank_z = sComm->rank_z;
    rank_X = sComm->rank_X;
    rank_Y = sComm->rank_Y;
    rank_Z = sComm->rank_Z;

    BoundaryCondition = sComm->BoundaryCondition;

    if (rank == 0) {
        printf("**** Creating membrane data structure ****** \n");
        printf("   Number of active lattice sites (rank = %i): %i \n", rank,
               Np);
    }

    sendCount_x = sComm->sendCount_x;
    sendCount_y = sComm->sendCount_y;
    sendCount_z = sComm->sendCount_z;
    sendCount_X = sComm->sendCount_X;
    sendCount_Y = sComm->sendCount_Y;
    sendCount_Z = sComm->sendCount_Z;

    recvCount_x = sComm->recvCount_x;
    recvCount_y = sComm->recvCount_y;
    recvCount_z = sComm->recvCount_z;
    recvCount_X = sComm->recvCount_X;
    recvCount_Y = sComm->recvCount_Y;
    recvCount_Z = sComm->recvCount_Z;

    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_x, recvCount_x * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_y, recvCount_y * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_z, recvCount_z * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_X, recvCount_X * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_Y, recvCount_Y * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_Z, recvCount_Z * sizeof(int));

    ScaLBL_AllocateZeroCopy((void **)&dvcRecvLinks_x,
                            recvCount_x * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvLinks_y,
                            recvCount_y * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvLinks_z,
                            recvCount_z * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvLinks_X,
                            recvCount_X * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvLinks_Y,
                            recvCount_Y * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvLinks_Z,
                            recvCount_Z * sizeof(int));

    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_x, recvCount_x * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_y, recvCount_y * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_z, recvCount_z * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_X, recvCount_X * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_Y, recvCount_Y * sizeof(int));
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_Z, recvCount_Z * sizeof(int));

    ScaLBL_AllocateZeroCopy((void **)&sendbuf_x, sendCount_x * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_y, sendCount_y * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_z, sendCount_z * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_X, sendCount_X * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_Y, sendCount_Y * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_Z, sendCount_Z * sizeof(double));

    ScaLBL_AllocateZeroCopy((void **)&recvbuf_x, recvCount_x * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_y, recvCount_y * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_z, recvCount_z * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_X, recvCount_X * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_Y, recvCount_Y * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_Z, recvCount_Z * sizeof(double));

    sendCount_x = sComm->copySendList("x", dvcSendList_x);
    sendCount_y = sComm->copySendList("y", dvcSendList_y);
    sendCount_z = sComm->copySendList("z", dvcSendList_z);
    sendCount_X = sComm->copySendList("X", dvcSendList_X);
    sendCount_Y = sComm->copySendList("Y", dvcSendList_Y);
    sendCount_Z = sComm->copySendList("Z", dvcSendList_Z);

    recvCount_x = sComm->copyRecvList("x", dvcRecvDist_x);
    recvCount_y = sComm->copyRecvList("y", dvcRecvDist_y);
    recvCount_z = sComm->copyRecvList("z", dvcRecvDist_z);
    recvCount_X = sComm->copyRecvList("X", dvcRecvDist_X);
    recvCount_Y = sComm->copyRecvList("Y", dvcRecvDist_Y);
    recvCount_Z = sComm->copyRecvList("Z", dvcRecvDist_Z);
}

Membrane::~Membrane() {

    delete[] initialNeighborList;
    delete[] membraneLinks;
    delete[] membraneTag;
    delete[] membraneDist;

    ScaLBL_FreeDeviceMemory(coefficient_x);
    ScaLBL_FreeDeviceMemory(coefficient_X);
    ScaLBL_FreeDeviceMemory(coefficient_y);
    ScaLBL_FreeDeviceMemory(coefficient_Y);
    ScaLBL_FreeDeviceMemory(coefficient_z);
    ScaLBL_FreeDeviceMemory(coefficient_Z);

    ScaLBL_FreeDeviceMemory(NeighborList);
    ScaLBL_FreeDeviceMemory(MembraneLinks);
    ScaLBL_FreeDeviceMemory(MembraneCoef);
    ScaLBL_FreeDeviceMemory(MembraneDistance);

    ScaLBL_FreeDeviceMemory(sendbuf_x);
    ScaLBL_FreeDeviceMemory(sendbuf_X);
    ScaLBL_FreeDeviceMemory(sendbuf_y);
    ScaLBL_FreeDeviceMemory(sendbuf_Y);
    ScaLBL_FreeDeviceMemory(sendbuf_z);
    ScaLBL_FreeDeviceMemory(sendbuf_Z);
    /*	ScaLBL_FreeDeviceMemory( sendbuf_xy );
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
	*/
    ScaLBL_FreeDeviceMemory(recvbuf_x);
    ScaLBL_FreeDeviceMemory(recvbuf_X);
    ScaLBL_FreeDeviceMemory(recvbuf_y);
    ScaLBL_FreeDeviceMemory(recvbuf_Y);
    ScaLBL_FreeDeviceMemory(recvbuf_z);
    ScaLBL_FreeDeviceMemory(recvbuf_Z);
    /*
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
	*/
    ScaLBL_FreeDeviceMemory(dvcSendList_x);
    ScaLBL_FreeDeviceMemory(dvcSendList_X);
    ScaLBL_FreeDeviceMemory(dvcSendList_y);
    ScaLBL_FreeDeviceMemory(dvcSendList_Y);
    ScaLBL_FreeDeviceMemory(dvcSendList_z);
    ScaLBL_FreeDeviceMemory(dvcSendList_Z);
    /*
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
	*/
    ScaLBL_FreeDeviceMemory(dvcRecvLinks_x);
    ScaLBL_FreeDeviceMemory(dvcRecvLinks_X);
    ScaLBL_FreeDeviceMemory(dvcRecvLinks_y);
    ScaLBL_FreeDeviceMemory(dvcRecvLinks_Y);
    ScaLBL_FreeDeviceMemory(dvcRecvLinks_z);
    ScaLBL_FreeDeviceMemory(dvcRecvLinks_Z);
    /*
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
	*/
    ScaLBL_FreeDeviceMemory(dvcRecvDist_x);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_X);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_y);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_Y);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_z);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_Z);
    /*
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
	*/
}

int Membrane::Create(DoubleArray &Distance, IntArray &Map) {
    int mlink = 0;
    int i, j, k;

    int idx, neighbor;
    double dist, locdist;

    if (rank == 0)
        printf("   Copy initial neighborlist... \n");
    int *neighborList = new int[18 * Np];
    /* Copy neighborList */
    for (int idx = 0; idx < Np; idx++) {
        for (int q = 0; q < 18; q++) {
            neighborList[q * Np + idx] = initialNeighborList[q * Np + idx];
        }
    }

    int Q = 7; // for D3Q7 model

    /* go through the neighborlist structure */
    /* count & cut the links */
    if (rank == 0)
        printf("   Cut membrane links... \n");
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                idx = Map(i, j, k);
                locdist = Distance(i, j, k);

                if (!(idx < 0)) {

                    neighbor = Map(i - 1, j, k);
                    dist = Distance(i - 1, j, k);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        neighborList[idx] = idx + 2 * Np;
                    }

                    neighbor = Map(i + 1, j, k);
                    dist = Distance(i + 1, j, k);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        neighborList[Np + idx] = idx + 1 * Np;
                        mlink++;
                    }

                    neighbor = Map(i, j - 1, k);
                    dist = Distance(i, j - 1, k);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        neighborList[2 * Np + idx] = idx + 4 * Np;
                    }

                    neighbor = Map(i, j + 1, k);
                    dist = Distance(i, j + 1, k);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        neighborList[3 * Np + idx] = idx + 3 * Np;
                        mlink++;
                    }

                    neighbor = Map(i, j, k - 1);
                    dist = Distance(i, j, k - 1);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        neighborList[4 * Np + idx] = idx + 6 * Np;
                    }

                    neighbor = Map(i, j, k + 1);
                    dist = Distance(i, j, k + 1);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        neighborList[5 * Np + idx] = idx + 5 * Np;
                        mlink++;
                    }

                    if (Q > 7) {
                        neighbor = Map(i - 1, j - 1, k);
                        dist = Distance(i - 1, j - 1, k);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[6 * Np + idx] = idx + 8 * Np;
                        }

                        neighbor = Map(i + 1, j + 1, k);
                        dist = Distance(i + 1, j + 1, k);
                        if (dist * locdist < 0.0) {
                            neighborList[7 * Np + idx] = idx + 7 * Np;
                            mlink++;
                        }

                        neighbor = Map(i - 1, j + 1, k);
                        dist = Distance(i - 1, j + 1, k);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[8 * Np + idx] = idx + 10 * Np;
                        }

                        neighbor = Map(i + 1, j - 1, k);
                        dist = Distance(i + 1, j - 1, k);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[9 * Np + idx] = idx + 9 * Np;
                            mlink++;
                        }

                        neighbor = Map(i - 1, j, k - 1);
                        dist = Distance(i - 1, j, k - 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[10 * Np + idx] = idx + 12 * Np;
                        }

                        neighbor = Map(i + 1, j, k + 1);
                        dist = Distance(i + 1, j, k + 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[11 * Np + idx] = idx + 11 * Np;
                            mlink++;
                        }

                        neighbor = Map(i - 1, j, k + 1);
                        dist = Distance(i - 1, j, k + 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[12 * Np + idx] = idx + 14 * Np;
                        }

                        neighbor = Map(i + 1, j, k - 1);
                        dist = Distance(i + 1, j, k - 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[13 * Np + idx] = idx + 13 * Np;
                            mlink++;
                        }

                        neighbor = Map(i, j - 1, k - 1);
                        dist = Distance(i, j - 1, k - 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[14 * Np + idx] = idx + 16 * Np;
                        }

                        neighbor = Map(i, j + 1, k + 1);
                        dist = Distance(i, j + 1, k + 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[15 * Np + idx] = idx + 15 * Np;
                            mlink++;
                        }

                        neighbor = Map(i, j - 1, k + 1);
                        dist = Distance(i, j - 1, k + 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[16 * Np + idx] = idx + 18 * Np;
                        }

                        neighbor = Map(i, j + 1, k - 1);
                        dist = Distance(i, j + 1, k - 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            neighborList[17 * Np + idx] = idx + 17 * Np;
                            mlink++;
                        }
                    }
                }
            }
        }
    }

    /* allocate memory */
    membraneTag = new int[mlink];
    membraneLinks = new int[2 * mlink];
    membraneDist = new double[2 * mlink];
    membraneLinkCount = mlink;

    if (rank == 0)
        printf("   (cut %i links crossing membrane) \n", mlink);

    /* construct the membrane*/
    /* *
	 *  Sites inside the membrane (negative distance) -- store at 2*mlink
	 *  Sites outside the membrane (positive distance) -- store at 2*mlink+1
	 */
    if (rank == 0)
        printf("   Construct membrane data structures... \n");
    mlink = 0;
    int localSite = 0;
    int neighborSite = 0;
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                idx = Map(i, j, k);
                locdist = Distance(i, j, k);

                if (!(idx < 0)) {

                    neighbor = Map(i + 1, j, k);
                    dist = Distance(i + 1, j, k);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        if (locdist < 0.0) {
                            localSite = 2 * mlink;
                            neighborSite = 2 * mlink + 1;
                        } else {
                            localSite = 2 * mlink + 1;
                            neighborSite = 2 * mlink;
                        }
                        membraneLinks[localSite] = idx + 1 * Np;
                        membraneLinks[neighborSite] = neighbor + 2 * Np;
                        membraneDist[localSite] = locdist;
                        membraneDist[neighborSite] = dist;
                        mlink++;
                    }

                    neighbor = Map(i, j + 1, k);
                    dist = Distance(i, j + 1, k);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        if (locdist < 0.0) {
                            localSite = 2 * mlink;
                            neighborSite = 2 * mlink + 1;
                        } else {
                            localSite = 2 * mlink + 1;
                            neighborSite = 2 * mlink;
                        }
                        membraneLinks[localSite] = idx + 3 * Np;
                        membraneLinks[neighborSite] = neighbor + 4 * Np;
                        membraneDist[localSite] = locdist;
                        membraneDist[neighborSite] = dist;
                        mlink++;
                    }

                    neighbor = Map(i, j, k + 1);
                    dist = Distance(i, j, k + 1);
                    if (dist * locdist < 0.0 && !(neighbor < 0)) {
                        if (locdist < 0.0) {
                            localSite = 2 * mlink;
                            neighborSite = 2 * mlink + 1;
                        } else {
                            localSite = 2 * mlink + 1;
                            neighborSite = 2 * mlink;
                        }
                        membraneLinks[localSite] = idx + 5 * Np;
                        membraneLinks[neighborSite] = neighbor + 6 * Np;
                        membraneDist[localSite] = locdist;
                        membraneDist[neighborSite] = dist;
                        mlink++;
                    }

                    if (Q > 7) {

                        neighbor = Map(i + 1, j + 1, k);
                        dist = Distance(i + 1, j + 1, k);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            if (locdist < 0.0) {
                                localSite = 2 * mlink;
                                neighborSite = 2 * mlink + 1;
                            } else {
                                localSite = 2 * mlink + 1;
                                neighborSite = 2 * mlink;
                            }
                            membraneLinks[localSite] = idx + 7 * Np;
                            membraneLinks[neighborSite] = neighbor + 8 * Np;
                            membraneDist[localSite] = locdist;
                            membraneDist[neighborSite] = dist;
                            mlink++;
                        }

                        neighbor = Map(i + 1, j - 1, k);
                        dist = Distance(i + 1, j - 1, k);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            if (locdist < 0.0) {
                                localSite = 2 * mlink;
                                neighborSite = 2 * mlink + 1;
                            } else {
                                localSite = 2 * mlink + 1;
                                neighborSite = 2 * mlink;
                            }
                            membraneLinks[localSite] = idx + 9 * Np;
                            membraneLinks[neighborSite] = neighbor + 10 * Np;
                            membraneDist[localSite] = locdist;
                            membraneDist[neighborSite] = dist;
                            mlink++;
                        }

                        neighbor = Map(i + 1, j, k + 1);
                        dist = Distance(i + 1, j, k + 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            if (locdist < 0.0) {
                                localSite = 2 * mlink;
                                neighborSite = 2 * mlink + 1;
                            } else {
                                localSite = 2 * mlink + 1;
                                neighborSite = 2 * mlink;
                            }
                            membraneLinks[localSite] = idx + 11 * Np;
                            membraneLinks[neighborSite] = neighbor + 12 * Np;
                            membraneDist[localSite] = locdist;
                            membraneDist[neighborSite] = dist;
                            mlink++;
                        }

                        neighbor = Map(i + 1, j, k - 1);
                        dist = Distance(i + 1, j, k - 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            if (locdist < 0.0) {
                                localSite = 2 * mlink;
                                neighborSite = 2 * mlink + 1;
                            } else {
                                localSite = 2 * mlink + 1;
                                neighborSite = 2 * mlink;
                            }
                            membraneLinks[localSite] = idx + 13 * Np;
                            membraneLinks[neighborSite] = neighbor + 14 * Np;
                            membraneDist[localSite] = locdist;
                            membraneDist[neighborSite] = dist;
                            mlink++;
                        }

                        neighbor = Map(i, j + 1, k + 1);
                        dist = Distance(i, j + 1, k + 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            if (locdist < 0.0) {
                                localSite = 2 * mlink;
                                neighborSite = 2 * mlink + 1;
                            } else {
                                localSite = 2 * mlink + 1;
                                neighborSite = 2 * mlink;
                            }
                            membraneLinks[localSite] = idx + 15 * Np;
                            membraneLinks[neighborSite] = neighbor + 16 * Np;
                            membraneDist[localSite] = locdist;
                            membraneDist[neighborSite] = dist;
                            mlink++;
                        }

                        neighbor = Map(i, j + 1, k - 1);
                        dist = Distance(i, j + 1, k - 1);
                        if (dist * locdist < 0.0 && !(neighbor < 0)) {
                            if (locdist < 0.0) {
                                localSite = 2 * mlink;
                                neighborSite = 2 * mlink + 1;
                            } else {
                                localSite = 2 * mlink + 1;
                                neighborSite = 2 * mlink;
                            }
                            membraneLinks[localSite] = idx + 17 * Np;
                            membraneLinks[neighborSite] = neighbor + 18 * Np;
                            membraneDist[localSite] = locdist;
                            membraneDist[neighborSite] = dist;
                            mlink++;
                        }
                    }
                }
            }
        }
    }

    if (rank == 0)
        printf("   Create device data structures... \n");
    /* Create device copies of data structures */
    ScaLBL_AllocateDeviceMemory((void **)&MembraneLinks,
                                2 * mlink * sizeof(int));
    ScaLBL_AllocateDeviceMemory((void **)&MembraneCoef,
                                2 * mlink * sizeof(double));
    //ScaLBL_AllocateDeviceMemory((void **)&MembraneDistance, 2*mlink*sizeof(double));
    ScaLBL_AllocateDeviceMemory((void **)&MembraneDistance,
                                Nx * Ny * Nz * sizeof(double));

    ScaLBL_CopyToDevice(NeighborList, neighborList, 18 * Np * sizeof(int));
    ScaLBL_CopyToDevice(MembraneLinks, membraneLinks, 2 * mlink * sizeof(int));
    //ScaLBL_CopyToDevice(MembraneDistance, membraneDist, 2*mlink*sizeof(double));
    ScaLBL_CopyToDevice(MembraneDistance, Distance.data(),
                        Nx * Ny * Nz * sizeof(double));

    int *dvcTmpMap;
    ScaLBL_AllocateDeviceMemory((void **)&dvcTmpMap, sizeof(int) * Np);
    int *TmpMap;
    TmpMap = new int[Np];
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int idx = Map(i, j, k);
                if (!(idx < 0))
                    TmpMap[idx] = k * Nx * Ny + j * Nx + i;
            }
        }
    }
    ScaLBL_CopyToDevice(dvcTmpMap, TmpMap, sizeof(int) * Np);

    //int Membrane::D3Q7_MapRecv(int Cqx, int Cqy, int Cqz, int *d3q19_recvlist,
    // int count, int *membraneRecvLabels, DoubleArray &Distance,  int  *dvcMap){
    if (rank == 0)
        printf("   Construct communication data structures... \n");
    /* Re-organize communication based on membrane structure*/
    //...dvcMap recieve list for the X face: q=2,8,10,12,14 .................................
    linkCount_X[0] = D3Q7_MapRecv(-1, 0, 0, dvcRecvDist_X, recvCount_X,
                                  dvcRecvLinks_X, Distance, dvcTmpMap);
    //...................................................................................
    //...dvcMap recieve list for the x face: q=1,7,9,11,13..................................
    linkCount_x[0] = D3Q7_MapRecv(1, 0, 0, dvcRecvDist_x, recvCount_x,
                                  dvcRecvLinks_x, Distance, dvcTmpMap);
    //...................................................................................
    //...dvcMap recieve list for the y face: q=4,8,9,16,18 ...................................
    linkCount_Y[0] = D3Q7_MapRecv(0, -1, 0, dvcRecvDist_Y, recvCount_Y,
                                  dvcRecvLinks_Y, Distance, dvcTmpMap);
    //...................................................................................
    //...dvcMap recieve list for the Y face: q=3,7,10,15,17 ..................................
    linkCount_y[0] = D3Q7_MapRecv(0, 1, 0, dvcRecvDist_y, recvCount_y,
                                  dvcRecvLinks_y, Distance, dvcTmpMap);
    //...................................................................................
    //...dvcMap recieve list for the z face<<<6,12,13,16,17)..............................................
    linkCount_Z[0] = D3Q7_MapRecv(0, 0, -1, dvcRecvDist_Z, recvCount_Z,
                                  dvcRecvLinks_Z, Distance, dvcTmpMap);

    //...dvcMap recieve list for the Z face<<<5,11,14,15,18)..............................................
    linkCount_z[0] = D3Q7_MapRecv(0, 0, 1, dvcRecvDist_z, recvCount_z,
                                  dvcRecvLinks_z, Distance, dvcTmpMap);
    //..................................................................................

    //......................................................................................
    MPI_COMM_SCALBL.barrier();
    ScaLBL_DeviceBarrier();
    //.......................................................................
    SendCount = sendCount_x + sendCount_X + sendCount_y + sendCount_Y +
                sendCount_z + sendCount_Z;
    RecvCount = recvCount_x + recvCount_X + recvCount_y + recvCount_Y +
                recvCount_z + recvCount_Z;
    CommunicationCount = SendCount + RecvCount;
    //......................................................................................

    //......................................................................................
    // Allocate membrane coefficient buffers (for d3q7 recv)
    ScaLBL_AllocateZeroCopy((void **)&coefficient_x,
                            2 * (recvCount_x) * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&coefficient_X,
                            2 * (recvCount_X) * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&coefficient_y,
                            2 * (recvCount_y) * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&coefficient_Y,
                            2 * (recvCount_Y) * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&coefficient_z,
                            2 * (recvCount_z) * sizeof(double));
    ScaLBL_AllocateZeroCopy((void **)&coefficient_Z,
                            2 * (recvCount_Z) * sizeof(double));
    //......................................................................................

    ScaLBL_FreeDeviceMemory(dvcTmpMap);
    delete[] neighborList;
    delete[] TmpMap;
    return mlink;
}

void Membrane::Write(string filename) {

    int mlink = membraneLinkCount;
    std::ofstream ofs(filename, std::ofstream::out);
    /* Create local copies of membrane data structures */
    double *tmpMembraneCoef; // mass transport coefficient for the membrane
    tmpMembraneCoef = new double[2 * mlink * sizeof(double)];
    ScaLBL_CopyToHost(tmpMembraneCoef, MembraneCoef,
                      2 * mlink * sizeof(double));
    int i, j, k;
    for (int m = 0; m < mlink; m++) {
        double a1 = tmpMembraneCoef[2 * m];
        double a2 = tmpMembraneCoef[2 * m + 1];
        int m1 = membraneLinks[2 * m] % Np;
        int m2 = membraneLinks[2 * m + 1] % Np;
        // map index to global i,j,k
        k = m1 / (Nx * Ny);
        j = (m1 - Nx * Ny * k) / Nx;
        i = m1 - Nx * Ny * k - Nx * j;
        ofs << i << " " << j << " " << k << " " << a1;
        k = m2 / (Nx * Ny);
        j = (m2 - Nx * Ny * k) / Nx;
        i = m2 - Nx * Ny * k - Nx * j;
        ofs << i << " " << j << " " << k << " " << a2 << endl;
    }
    ofs.close();

    /*FILE *VELX_FILE;
    sprintf(LocalRankFilename, "Velocity_X.%05i.raw", rank);
    VELX_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELX_FILE);
    fclose(VELX_FILE);
    */
    delete[] tmpMembraneCoef;
}

void Membrane::Read(string filename) {

    int mlink = membraneLinkCount;
    /* Create local copies of membrane data structures */
    double *tmpMembraneCoef; // mass transport coefficient for the membrane
    tmpMembraneCoef = new double[2 * mlink * sizeof(double)];

    FILE *fid = fopen(filename.c_str(), "r");
    INSIST(fid != NULL, "Error opening membrane file \n");
    //........read the spheres..................
    // We will read until a blank like or end-of-file is reached
    int count = 0;
    int i, j, k;
    int ii, jj, kk;
    double a1, a2;

    while (fscanf(fid, "%i,%i,%i,%lf,%i,%i,%i,%lf,\n", &i, &j, &k, &a1, &ii,
                  &jj, &kk, &a2) == 8) {
        printf("%i, %i, %i, %lf \n", i, j, k, a2);
        count++;
    }
    if (count != mlink) {
        printf("WARNING (Membrane::Read): number of file lines does not match "
               "number of links \n");
    }
    fclose(fid);

    ScaLBL_CopyToDevice(MembraneCoef, tmpMembraneCoef,
                        2 * mlink * sizeof(double));

    /*FILE *VELX_FILE;
    sprintf(LocalRankFilename, "Velocity_X.%05i.raw", rank);
    VELX_FILE = fopen(LocalRankFilename, "wb");
    fwrite(PhaseField.data(), 8, N, VELX_FILE);
    fclose(VELX_FILE);
    */
    delete[] tmpMembraneCoef;
}

int Membrane::D3Q7_MapRecv(int Cqx, int Cqy, int Cqz, int *d3q19_recvlist,
                           int count, int *membraneRecvLabels,
                           DoubleArray &Distance, int *dvcMap) {

    int i, j, k, n, nn, idx;
    double distanceNonLocal, distanceLocal;
    int *ReturnLabels;
    ReturnLabels = new int[count];
    int *list;
    list = new int[count];
    ScaLBL_CopyToHost(list, d3q19_recvlist, count * sizeof(int));

    int *TmpMap;
    TmpMap = new int[Np];
    ScaLBL_CopyToHost(TmpMap, dvcMap, Np * sizeof(int));

    int countMembraneLinks = 0;
    for (idx = 0; idx < count; idx++) {
        //printf("    Read 1 \n");
        // Get the value from the list -- note that n is the index is from the send (non-local) process
        nn = list[idx]; // if (rank == 0) printf("@ rank:%d n=%d\n",rank,n);
        //printf("    Read 2 \n");
        n = TmpMap[nn];
        //printf("    idx= %i(%i), nn=%i, n= %i \n",idx,count,nn,n);

        // Get the 3-D indices from the send process
        k = n / (Nx * Ny);
        j = (n - Nx * Ny * k) / Nx;
        i = n - Nx * Ny * k - Nx * j;
        // if (rank ==0) printf("@ Get 3D indices from the send process: i=%d, j=%d, k=%d\n",i,j,k);

        distanceLocal = Distance(i, j, k); // this site should be in the halo
        //printf("    Local value %i, %i, %i \n",i,j,k);

        // Streaming for the non-local distribution
        i -= Cqx;
        j -= Cqy;
        k -= Cqz;
        distanceNonLocal = Distance(i, j, k);
        //printf("    Nonlocal value %i, %i, %i \n",i,j,k);
        ReturnLabels[idx] = 0;
        if (distanceLocal * distanceNonLocal < 0.0) {
            if (distanceLocal > 0.0)
                ReturnLabels[idx] = 1;
            else
                ReturnLabels[idx] = 2;
            countMembraneLinks++;
        }
    }
    // Return updated version to the device
    ScaLBL_CopyToDevice(membraneRecvLabels, ReturnLabels, count * sizeof(int));

    // clean up the work arrays
    delete[] ReturnLabels;
    delete[] TmpMap;
    delete[] list;
    return countMembraneLinks;
}

void Membrane::SendD3Q7AA(double *dist) {

    if (Lock == true) {
        ERROR("Membrane Error (SendD3Q7): Membrane communicator is locked -- "
              "did you forget to match Send/Recv calls?");
    } else {
        Lock = true;
    }
    // assign tag of 37 to D3Q7 communication
    sendtag = recvtag = 37;
    ScaLBL_DeviceBarrier();
    // Pack the distributions
    //...Packing for x face(q=2)................................
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, 0, sendCount_x, sendbuf_x, dist, Np);
    req2[0] = MPI_COMM_SCALBL.Irecv(recvbuf_X, recvCount_X, rank_X, recvtag);
    req1[0] = MPI_COMM_SCALBL.Isend(sendbuf_x, sendCount_x, rank_x, sendtag);
    //...Packing for X face(q=1)................................
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, 0, sendCount_X, sendbuf_X, dist, Np);
    req2[1] = MPI_COMM_SCALBL.Irecv(recvbuf_x, recvCount_x, rank_x, recvtag);
    req1[1] = MPI_COMM_SCALBL.Isend(sendbuf_X, sendCount_X, rank_X, sendtag);
    //for (int idx=0; idx<sendCount_X; idx++) printf(" SendX(%i)=%e \n",idx,sendbuf_X[idx]);
    //...Packing for y face(q=4).................................
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, 0, sendCount_y, sendbuf_y, dist, Np);
    req2[2] = MPI_COMM_SCALBL.Irecv(recvbuf_Y, recvCount_Y, rank_Y, recvtag);
    req1[2] = MPI_COMM_SCALBL.Isend(sendbuf_y, sendCount_y, rank_y, sendtag);
    //...Packing for Y face(q=3).................................
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, 0, sendCount_Y, sendbuf_Y, dist, Np);
    req2[3] = MPI_COMM_SCALBL.Irecv(recvbuf_y, recvCount_y, rank_y, recvtag);
    req1[3] = MPI_COMM_SCALBL.Isend(sendbuf_Y, sendCount_Y, rank_Y, sendtag);
    //for (int idx=0; idx<sendCount_Y; idx++) printf(" SendY(%i)=%e \n",idx,sendbuf_Y[idx]);
    //...Packing for z face(q=6)................................
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, 0, sendCount_z, sendbuf_z, dist, Np);
    req2[4] = MPI_COMM_SCALBL.Irecv(recvbuf_Z, recvCount_Z, rank_Z, recvtag);
    req1[4] = MPI_COMM_SCALBL.Isend(sendbuf_z, sendCount_z, rank_z, sendtag);
    //...Packing for Z face(q=5)................................
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, 0, sendCount_Z, sendbuf_Z, dist, Np);
    req2[5] = MPI_COMM_SCALBL.Irecv(recvbuf_z, recvCount_z, rank_z, recvtag);
    req1[5] = MPI_COMM_SCALBL.Isend(sendbuf_Z, sendCount_Z, rank_Z, sendtag);
}

void Membrane::RecvD3Q7AA(double *dist, bool BounceBack) {

    //...................................................................................
    // Wait for completion of D3Q19 communication
    MPI_COMM_SCALBL.waitAll(6, req1);
    MPI_COMM_SCALBL.waitAll(6, req2);
    ScaLBL_DeviceBarrier();
    //...................................................................................
    // NOTE: AA Routine writes to opposite
    // Unpack the distributions on the device
    //...................................................................................
    //...Unpacking for x face(q=2)................................
    ScaLBL_D3Q7_Membrane_Unpack(2, dvcRecvDist_x, recvbuf_x, recvCount_x, dist,
                                Np, coefficient_x);
    //...................................................................................
    //...Packing for X face(q=1)................................
    ScaLBL_D3Q7_Membrane_Unpack(1, dvcRecvDist_X, recvbuf_X, recvCount_X, dist,
                                Np, coefficient_X);
    //...................................................................................
    //...Packing for y face(q=4).................................
    ScaLBL_D3Q7_Membrane_Unpack(4, dvcRecvDist_y, recvbuf_y, recvCount_y, dist,
                                Np, coefficient_y);
    //...................................................................................
    //...Packing for Y face(q=3).................................
    ScaLBL_D3Q7_Membrane_Unpack(3, dvcRecvDist_Y, recvbuf_Y, recvCount_Y, dist,
                                Np, coefficient_Y);
    //...................................................................................
    //if (BoundaryCondition > 0 && rank_info.kz == 0)
    if (BounceBack &&
        rank_info.kz == 0) { /* leave the bounce-back distributions in place */
    } else {
        //...Packing for z face(q=6)................................
        ScaLBL_D3Q7_Membrane_Unpack(6, dvcRecvDist_z, recvbuf_z, recvCount_z,
                                    dist, Np, coefficient_z);
    }
    //if (BoundaryCondition > 0 && rank_info.kz == rank_info.nz-1)
    if (BounceBack &&
        rank_info.kz ==
            rank_info.nz -
                1) { /* leave the bounce-back distributions in place */
    } else {
        //...Packing for Z face(q=5)................................
        ScaLBL_D3Q7_Membrane_Unpack(5, dvcRecvDist_Z, recvbuf_Z, recvCount_Z,
                                    dist, Np, coefficient_Z);
        //..................................................................................
    }
    MPI_COMM_SCALBL.barrier();
    //...................................................................................
    Lock = false; // unlock the communicator after communications complete
    //...................................................................................
}

void Membrane::IonTransport(double *dist, double *den) {

    ScaLBL_D3Q7_Membrane_IonTransport(MembraneLinks, MembraneCoef, dist, den,
                                      membraneLinkCount, Np);
}

//	std::shared_ptr<Database> db){
void Membrane::AssignCoefficients(int *Map, double *Psi, double Threshold,
                                  double MassFractionIn, double MassFractionOut,
                                  double ThresholdMassFractionIn,
                                  double ThresholdMassFractionOut) {
    /* Assign mass transfer coefficients to the membrane data structure */

    if (membraneLinkCount > 0)
        ScaLBL_D3Q7_Membrane_AssignLinkCoef(
            MembraneLinks, Map, MembraneDistance, Psi, MembraneCoef, Threshold,
            MassFractionIn, MassFractionOut, ThresholdMassFractionIn,
            ThresholdMassFractionOut, membraneLinkCount, Nx, Ny, Nz, Np);

    if (linkCount_X[0] < recvCount_X)
        ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
            -1, 0, 0, Map, MembraneDistance, Psi, Threshold, MassFractionIn,
            MassFractionOut, ThresholdMassFractionIn, ThresholdMassFractionOut,
            dvcRecvDist_X, dvcRecvLinks_X, coefficient_X, 0, linkCount_X[0],
            recvCount_X, Np, Nx, Ny, Nz);

    if (linkCount_x[0] < recvCount_x)
        ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
            1, 0, 0, Map, MembraneDistance, Psi, Threshold, MassFractionIn,
            MassFractionOut, ThresholdMassFractionIn, ThresholdMassFractionOut,
            dvcRecvDist_x, dvcRecvLinks_x, coefficient_x, 0, linkCount_x[0],
            recvCount_x, Np, Nx, Ny, Nz);

    if (linkCount_Y[0] < recvCount_Y)
        ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
            0, -1, 0, Map, MembraneDistance, Psi, Threshold, MassFractionIn,
            MassFractionOut, ThresholdMassFractionIn, ThresholdMassFractionOut,
            dvcRecvDist_Y, dvcRecvLinks_Y, coefficient_Y, 0, linkCount_Y[0],
            recvCount_Y, Np, Nx, Ny, Nz);

    if (linkCount_y[0] < recvCount_y)
        ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
            0, 1, 0, Map, MembraneDistance, Psi, Threshold, MassFractionIn,
            MassFractionOut, ThresholdMassFractionIn, ThresholdMassFractionOut,
            dvcRecvDist_y, dvcRecvLinks_y, coefficient_y, 0, linkCount_y[0],
            recvCount_y, Np, Nx, Ny, Nz);

    if (linkCount_Z[0] < recvCount_Z)
        ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
            0, 0, -1, Map, MembraneDistance, Psi, Threshold, MassFractionIn,
            MassFractionOut, ThresholdMassFractionIn, ThresholdMassFractionOut,
            dvcRecvDist_Z, dvcRecvLinks_Z, coefficient_Z, 0, linkCount_Z[0],
            recvCount_Z, Np, Nx, Ny, Nz);

    if (linkCount_z[0] < recvCount_z)
        ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
            0, 0, 1, Map, MembraneDistance, Psi, Threshold, MassFractionIn,
            MassFractionOut, ThresholdMassFractionIn, ThresholdMassFractionOut,
            dvcRecvDist_z, dvcRecvLinks_z, coefficient_z, 0, linkCount_z[0],
            recvCount_z, Np, Nx, Ny, Nz);
}
