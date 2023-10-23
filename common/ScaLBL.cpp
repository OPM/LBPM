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
#include "common/ScaLBL.h"
#include <chrono>

ScaLBL_Communicator::ScaLBL_Communicator(std::shared_ptr<Domain> Dm) {
    //......................................................................................
    Lock = false; // unlock the communicator
    //......................................................................................
    // Create a separate copy of the communicator for the device
    MPI_COMM_SCALBL = Dm->Comm.dup();
    int myrank = MPI_COMM_SCALBL.getRank();
    rank_info =
        RankInfoStruct(myrank, rank_info.nx, rank_info.ny, rank_info.nz);

    //......................................................................................
    // Copy the domain size and communication information directly from Dm
    Nx = Dm->Nx;
    Ny = Dm->Ny;
    Nz = Dm->Nz;
    N = Nx * Ny * Nz;
    next = 0;
    rank = Dm->rank();
    rank_x = Dm->rank_x();
    rank_y = Dm->rank_y();
    rank_z = Dm->rank_z();
    rank_X = Dm->rank_X();
    rank_Y = Dm->rank_Y();
    rank_Z = Dm->rank_Z();
    rank_xy = Dm->rank_xy();
    rank_XY = Dm->rank_XY();
    rank_xY = Dm->rank_xY();
    rank_Xy = Dm->rank_Xy();
    rank_xz = Dm->rank_xz();
    rank_XZ = Dm->rank_XZ();
    rank_xZ = Dm->rank_xZ();
    rank_Xz = Dm->rank_Xz();
    rank_yz = Dm->rank_yz();
    rank_YZ = Dm->rank_YZ();
    rank_yZ = Dm->rank_yZ();
    rank_Yz = Dm->rank_Yz();
    sendCount_x = Dm->sendCount("x");
    sendCount_y = Dm->sendCount("y");
    sendCount_z = Dm->sendCount("z");
    sendCount_X = Dm->sendCount("X");
    sendCount_Y = Dm->sendCount("Y");
    sendCount_Z = Dm->sendCount("Z");
    sendCount_xy = Dm->sendCount("xy");
    sendCount_yz = Dm->sendCount("yz");
    sendCount_xz = Dm->sendCount("xz");
    sendCount_Xy = Dm->sendCount("Xy");
    sendCount_Yz = Dm->sendCount("Yz");
    sendCount_xZ = Dm->sendCount("xZ");
    sendCount_xY = Dm->sendCount("xY");
    sendCount_yZ = Dm->sendCount("yZ");
    sendCount_Xz = Dm->sendCount("Xz");
    sendCount_XY = Dm->sendCount("XY");
    sendCount_YZ = Dm->sendCount("YZ");
    sendCount_XZ = Dm->sendCount("XZ");
    recvCount_x = Dm->recvCount("x");
    recvCount_y = Dm->recvCount("y");
    recvCount_z = Dm->recvCount("z");
    recvCount_X = Dm->recvCount("X");
    recvCount_Y = Dm->recvCount("Y");
    recvCount_Z = Dm->recvCount("Z");
    recvCount_xy = Dm->recvCount("xy");
    recvCount_yz = Dm->recvCount("yz");
    recvCount_xz = Dm->recvCount("xz");
    recvCount_Xy = Dm->recvCount("Xy");
    recvCount_Yz = Dm->recvCount("Yz");
    recvCount_xZ = Dm->recvCount("xZ");
    recvCount_xY = Dm->recvCount("xY");
    recvCount_yZ = Dm->recvCount("yZ");
    recvCount_Xz = Dm->recvCount("Xz");
    recvCount_XY = Dm->recvCount("XY");
    recvCount_YZ = Dm->recvCount("YZ");
    recvCount_XZ = Dm->recvCount("XZ");

    iproc = Dm->iproc();
    jproc = Dm->jproc();
    kproc = Dm->kproc();
    nprocx = Dm->nprocx();
    nprocy = Dm->nprocy();
    nprocz = Dm->nprocz();
    BoundaryCondition = Dm->BoundaryCondition;
    //......................................................................................

    ScaLBL_AllocateZeroCopy((void **)&sendbuf_x,
                            2 * 5 * sendCount_x *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_X,
                            2 * 5 * sendCount_X *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_y,
                            2 * 5 * sendCount_y *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_Y,
                            2 * 5 * sendCount_Y *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_z,
                            2 * 5 * sendCount_z *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_Z,
                            2 * 5 * sendCount_Z *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_xy,
                            2 * sendCount_xy *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_xY,
                            2 * sendCount_xY *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_Xy,
                            2 * sendCount_Xy *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_XY,
                            2 * sendCount_XY *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_xz,
                            2 * sendCount_xz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_xZ,
                            2 * sendCount_xZ *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_Xz,
                            2 * sendCount_Xz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_XZ,
                            2 * sendCount_XZ *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_yz,
                            2 * sendCount_yz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_yZ,
                            2 * sendCount_yZ *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_Yz,
                            2 * sendCount_Yz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&sendbuf_YZ,
                            2 * sendCount_YZ *
                                sizeof(double)); // Allocate device memory
    //......................................................................................
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_x,
                            2 * 5 * recvCount_x *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_X,
                            2 * 5 * recvCount_X *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_y,
                            2 * 5 * recvCount_y *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_Y,
                            2 * 5 * recvCount_Y *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_z,
                            2 * 5 * recvCount_z *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_Z,
                            2 * 5 * recvCount_Z *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_xy,
                            2 * recvCount_xy *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_xY,
                            2 * recvCount_xY *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_Xy,
                            2 * recvCount_Xy *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_XY,
                            2 * recvCount_XY *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_xz,
                            2 * recvCount_xz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_xZ,
                            2 * recvCount_xZ *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_Xz,
                            2 * recvCount_Xz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_XZ,
                            2 * recvCount_XZ *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_yz,
                            2 * recvCount_yz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_yZ,
                            2 * recvCount_yZ *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_Yz,
                            2 * recvCount_Yz *
                                sizeof(double)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&recvbuf_YZ,
                            2 * recvCount_YZ *
                                sizeof(double)); // Allocate device memory
    //......................................................................................
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_x,
                            sendCount_x *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_X,
                            sendCount_X *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_y,
                            sendCount_y *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_Y,
                            sendCount_Y *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_z,
                            sendCount_z *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_Z,
                            sendCount_Z *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_xy,
                            sendCount_xy *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_xY,
                            sendCount_xY *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_Xy,
                            sendCount_Xy *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_XY,
                            sendCount_XY *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_xz,
                            sendCount_xz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_xZ,
                            sendCount_xZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_Xz,
                            sendCount_Xz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_XZ,
                            sendCount_XZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_yz,
                            sendCount_yz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_yZ,
                            sendCount_yZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_Yz,
                            sendCount_Yz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcSendList_YZ,
                            sendCount_YZ *
                                sizeof(int)); // Allocate device memory
    //......................................................................................
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_x,
                            recvCount_x *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_X,
                            recvCount_X *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_y,
                            recvCount_y *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_Y,
                            recvCount_Y *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_z,
                            recvCount_z *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_Z,
                            recvCount_Z *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_xy,
                            recvCount_xy *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_xY,
                            recvCount_xY *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_Xy,
                            recvCount_Xy *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_XY,
                            recvCount_XY *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_xz,
                            recvCount_xz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_xZ,
                            recvCount_xZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_Xz,
                            recvCount_Xz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_XZ,
                            recvCount_XZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_yz,
                            recvCount_yz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_yZ,
                            recvCount_yZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_Yz,
                            recvCount_Yz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvList_YZ,
                            recvCount_YZ *
                                sizeof(int)); // Allocate device memory
    //......................................................................................
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_x,
                            5 * recvCount_x *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_X,
                            5 * recvCount_X *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_y,
                            5 * recvCount_y *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_Y,
                            5 * recvCount_Y *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_z,
                            5 * recvCount_z *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_Z,
                            5 * recvCount_Z *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_xy,
                            recvCount_xy *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_xY,
                            recvCount_xY *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_Xy,
                            recvCount_Xy *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_XY,
                            recvCount_XY *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_xz,
                            recvCount_xz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_xZ,
                            recvCount_xZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_Xz,
                            recvCount_Xz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_XZ,
                            recvCount_XZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_yz,
                            recvCount_yz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_yZ,
                            recvCount_yZ *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_Yz,
                            recvCount_Yz *
                                sizeof(int)); // Allocate device memory
    ScaLBL_AllocateZeroCopy((void **)&dvcRecvDist_YZ,
                            recvCount_YZ *
                                sizeof(int)); // Allocate device memory
    //......................................................................................

    //......................................................................................
    ScaLBL_CopyToZeroCopy(dvcSendList_x, Dm->sendList("x"),
                          sendCount_x * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_X, Dm->sendList("X"),
                          sendCount_X * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_y, Dm->sendList("y"),
                          sendCount_y * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_Y, Dm->sendList("Y"),
                          sendCount_Y * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_z, Dm->sendList("z"),
                          sendCount_z * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_Z, Dm->sendList("Z"),
                          sendCount_Z * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_xy, Dm->sendList("xy"),
                          sendCount_xy * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_XY, Dm->sendList("XY"),
                          sendCount_XY * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_xY, Dm->sendList("xY"),
                          sendCount_xY * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_Xy, Dm->sendList("Xy"),
                          sendCount_Xy * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_xz, Dm->sendList("xz"),
                          sendCount_xz * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_XZ, Dm->sendList("XZ"),
                          sendCount_XZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_xZ, Dm->sendList("xZ"),
                          sendCount_xZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_Xz, Dm->sendList("Xz"),
                          sendCount_Xz * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_yz, Dm->sendList("yz"),
                          sendCount_yz * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_YZ, Dm->sendList("YZ"),
                          sendCount_YZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_yZ, Dm->sendList("yZ"),
                          sendCount_yZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcSendList_Yz, Dm->sendList("Yz"),
                          sendCount_Yz * sizeof(int));
    //......................................................................................
    ScaLBL_CopyToZeroCopy(dvcRecvList_x, Dm->recvList("x"),
                          recvCount_x * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_X, Dm->recvList("X"),
                          recvCount_X * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_y, Dm->recvList("y"),
                          recvCount_y * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_Y, Dm->recvList("Y"),
                          recvCount_Y * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_z, Dm->recvList("z"),
                          recvCount_z * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_Z, Dm->recvList("Z"),
                          recvCount_Z * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_xy, Dm->recvList("xy"),
                          recvCount_xy * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_XY, Dm->recvList("XY"),
                          recvCount_XY * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_xY, Dm->recvList("xY"),
                          recvCount_xY * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_Xy, Dm->recvList("Xy"),
                          recvCount_Xy * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_xz, Dm->recvList("xz"),
                          recvCount_xz * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_XZ, Dm->recvList("XZ"),
                          recvCount_XZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_xZ, Dm->recvList("xZ"),
                          recvCount_xZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_Xz, Dm->recvList("Xz"),
                          recvCount_Xz * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_yz, Dm->recvList("yz"),
                          recvCount_yz * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_YZ, Dm->recvList("YZ"),
                          recvCount_YZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_yZ, Dm->recvList("yZ"),
                          recvCount_yZ * sizeof(int));
    ScaLBL_CopyToZeroCopy(dvcRecvList_Yz, Dm->recvList("Yz"),
                          recvCount_Yz * sizeof(int));
    //......................................................................................

    MPI_COMM_SCALBL.barrier();

    //...................................................................................
    // Set up the recieve distribution lists
    //...................................................................................
    //...Map recieve list for the X face: q=2,8,10,12,14 .................................
    D3Q19_MapRecv(-1, 0, 0, Dm->recvList("X"), 0, recvCount_X, dvcRecvDist_X);
    D3Q19_MapRecv(-1, -1, 0, Dm->recvList("X"), recvCount_X, recvCount_X,
                  dvcRecvDist_X);
    D3Q19_MapRecv(-1, 1, 0, Dm->recvList("X"), 2 * recvCount_X, recvCount_X,
                  dvcRecvDist_X);
    D3Q19_MapRecv(-1, 0, -1, Dm->recvList("X"), 3 * recvCount_X, recvCount_X,
                  dvcRecvDist_X);
    D3Q19_MapRecv(-1, 0, 1, Dm->recvList("X"), 4 * recvCount_X, recvCount_X,
                  dvcRecvDist_X);
    //...................................................................................
    //...Map recieve list for the x face: q=1,7,9,11,13..................................
    D3Q19_MapRecv(1, 0, 0, Dm->recvList("x"), 0, recvCount_x, dvcRecvDist_x);
    D3Q19_MapRecv(1, 1, 0, Dm->recvList("x"), recvCount_x, recvCount_x,
                  dvcRecvDist_x);
    D3Q19_MapRecv(1, -1, 0, Dm->recvList("x"), 2 * recvCount_x, recvCount_x,
                  dvcRecvDist_x);
    D3Q19_MapRecv(1, 0, 1, Dm->recvList("x"), 3 * recvCount_x, recvCount_x,
                  dvcRecvDist_x);
    D3Q19_MapRecv(1, 0, -1, Dm->recvList("x"), 4 * recvCount_x, recvCount_x,
                  dvcRecvDist_x);
    //...................................................................................
    //...Map recieve list for the y face: q=4,8,9,16,18 ...................................
    D3Q19_MapRecv(0, -1, 0, Dm->recvList("Y"), 0, recvCount_Y, dvcRecvDist_Y);
    D3Q19_MapRecv(-1, -1, 0, Dm->recvList("Y"), recvCount_Y, recvCount_Y,
                  dvcRecvDist_Y);
    D3Q19_MapRecv(1, -1, 0, Dm->recvList("Y"), 2 * recvCount_Y, recvCount_Y,
                  dvcRecvDist_Y);
    D3Q19_MapRecv(0, -1, -1, Dm->recvList("Y"), 3 * recvCount_Y, recvCount_Y,
                  dvcRecvDist_Y);
    D3Q19_MapRecv(0, -1, 1, Dm->recvList("Y"), 4 * recvCount_Y, recvCount_Y,
                  dvcRecvDist_Y);
    //...................................................................................
    //...Map recieve list for the Y face: q=3,7,10,15,17 ..................................
    D3Q19_MapRecv(0, 1, 0, Dm->recvList("y"), 0, recvCount_y, dvcRecvDist_y);
    D3Q19_MapRecv(1, 1, 0, Dm->recvList("y"), recvCount_y, recvCount_y,
                  dvcRecvDist_y);
    D3Q19_MapRecv(-1, 1, 0, Dm->recvList("y"), 2 * recvCount_y, recvCount_y,
                  dvcRecvDist_y);
    D3Q19_MapRecv(0, 1, 1, Dm->recvList("y"), 3 * recvCount_y, recvCount_y,
                  dvcRecvDist_y);
    D3Q19_MapRecv(0, 1, -1, Dm->recvList("y"), 4 * recvCount_y, recvCount_y,
                  dvcRecvDist_y);
    //...................................................................................
    //...Map recieve list for the z face<<<6,12,13,16,17)..............................................
    D3Q19_MapRecv(0, 0, -1, Dm->recvList("Z"), 0, recvCount_Z, dvcRecvDist_Z);
    D3Q19_MapRecv(-1, 0, -1, Dm->recvList("Z"), recvCount_Z, recvCount_Z,
                  dvcRecvDist_Z);
    D3Q19_MapRecv(1, 0, -1, Dm->recvList("Z"), 2 * recvCount_Z, recvCount_Z,
                  dvcRecvDist_Z);
    D3Q19_MapRecv(0, -1, -1, Dm->recvList("Z"), 3 * recvCount_Z, recvCount_Z,
                  dvcRecvDist_Z);
    D3Q19_MapRecv(0, 1, -1, Dm->recvList("Z"), 4 * recvCount_Z, recvCount_Z,
                  dvcRecvDist_Z);
    //...Map recieve list for the Z face<<<5,11,14,15,18)..............................................
    D3Q19_MapRecv(0, 0, 1, Dm->recvList("z"), 0, recvCount_z, dvcRecvDist_z);
    D3Q19_MapRecv(1, 0, 1, Dm->recvList("z"), recvCount_z, recvCount_z,
                  dvcRecvDist_z);
    D3Q19_MapRecv(-1, 0, 1, Dm->recvList("z"), 2 * recvCount_z, recvCount_z,
                  dvcRecvDist_z);
    D3Q19_MapRecv(0, 1, 1, Dm->recvList("z"), 3 * recvCount_z, recvCount_z,
                  dvcRecvDist_z);
    D3Q19_MapRecv(0, -1, 1, Dm->recvList("z"), 4 * recvCount_z, recvCount_z,
                  dvcRecvDist_z);
    //..................................................................................
    //...Map recieve list for the xy edge <<<8)................................
    D3Q19_MapRecv(-1, -1, 0, Dm->recvList("XY"), 0, recvCount_XY,
                  dvcRecvDist_XY);
    //...Map recieve list for the Xy edge <<<9)................................
    D3Q19_MapRecv(1, -1, 0, Dm->recvList("xY"), 0, recvCount_xY,
                  dvcRecvDist_xY);
    //...Map recieve list for the xY edge <<<10)................................
    D3Q19_MapRecv(-1, 1, 0, Dm->recvList("Xy"), 0, recvCount_Xy,
                  dvcRecvDist_Xy);
    //...Map recieve list for the XY edge <<<7)................................
    D3Q19_MapRecv(1, 1, 0, Dm->recvList("xy"), 0, recvCount_xy, dvcRecvDist_xy);
    //...Map recieve list for the xz edge <<<12)................................
    D3Q19_MapRecv(-1, 0, -1, Dm->recvList("XZ"), 0, recvCount_XZ,
                  dvcRecvDist_XZ);
    //...Map recieve list for the xZ edge <<<14)................................
    D3Q19_MapRecv(-1, 0, 1, Dm->recvList("Xz"), 0, recvCount_Xz,
                  dvcRecvDist_Xz);
    //...Map recieve list for the Xz edge <<<13)................................
    D3Q19_MapRecv(1, 0, -1, Dm->recvList("xZ"), 0, recvCount_xZ,
                  dvcRecvDist_xZ);
    //...Map recieve list for the XZ edge <<<11)................................
    D3Q19_MapRecv(1, 0, 1, Dm->recvList("xz"), 0, recvCount_xz, dvcRecvDist_xz);
    //...Map recieve list for the yz edge <<<16)................................
    D3Q19_MapRecv(0, -1, -1, Dm->recvList("YZ"), 0, recvCount_YZ,
                  dvcRecvDist_YZ);
    //...Map recieve list for the yZ edge <<<18)................................
    D3Q19_MapRecv(0, -1, 1, Dm->recvList("Yz"), 0, recvCount_Yz,
                  dvcRecvDist_Yz);
    //...Map recieve list for the Yz edge <<<17)................................
    D3Q19_MapRecv(0, 1, -1, Dm->recvList("yZ"), 0, recvCount_yZ,
                  dvcRecvDist_yZ);
    //...Map recieve list for the YZ edge <<<15)................................
    D3Q19_MapRecv(0, 1, 1, Dm->recvList("yz"), 0, recvCount_yz, dvcRecvDist_yz);
    //...................................................................................

    //......................................................................................
    MPI_COMM_SCALBL.barrier();
    ScaLBL_DeviceBarrier();
    //......................................................................................
    SendCount = 5 * (sendCount_x + sendCount_X + sendCount_y + sendCount_Y +
                     sendCount_z + sendCount_Z) +
                sendCount_xy + sendCount_Xy + sendCount_xY + sendCount_XY +
                sendCount_xZ + sendCount_Xz + sendCount_xZ + sendCount_XZ +
                sendCount_yz + sendCount_Yz + sendCount_yZ + sendCount_YZ;

    RecvCount = 5 * (recvCount_x + recvCount_X + recvCount_y + recvCount_Y +
                     recvCount_z + recvCount_Z) +
                recvCount_xy + recvCount_Xy + recvCount_xY + recvCount_XY +
                recvCount_xZ + recvCount_Xz + recvCount_xZ + recvCount_XZ +
                recvCount_yz + recvCount_Yz + recvCount_yZ + recvCount_YZ;

    CommunicationCount = SendCount + RecvCount;
    //......................................................................................

    //...................................................................................
    // Set up the persistent communication for D3Q19AA (use tags 130-145)
    //...................................................................................
    req_D3Q19AA.clear();
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_x, 5 * sendCount_x, rank_x, 130));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_X, 5 * recvCount_X, rank_X, 130));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_X, 5 * sendCount_X, rank_X, 131));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_x, 5 * recvCount_x, rank_x, 131));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_y, 5 * sendCount_y, rank_y, 132));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_Y, 5 * recvCount_Y, rank_Y, 132));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_Y, 5 * sendCount_Y, rank_Y, 133));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_y, 5 * recvCount_y, rank_y, 133));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_z, 5 * sendCount_z, rank_z, 134));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_Z, 5 * recvCount_Z, rank_Z, 134));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_Z, 5 * sendCount_Z, rank_Z, 135));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_z, 5 * recvCount_z, rank_z, 135));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_xy, sendCount_xy, rank_xy, 136));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_XY, recvCount_XY, rank_XY, 136));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_XY, sendCount_XY, rank_XY, 137));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_xy, recvCount_xy, rank_xy, 137));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_Xy, sendCount_Xy, rank_Xy, 138));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_xY, recvCount_xY, rank_xY, 138));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_xY, sendCount_xY, rank_xY, 139));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_Xy, recvCount_Xy, rank_Xy, 139));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_xz, sendCount_xz, rank_xz, 140));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_XZ, recvCount_XZ, rank_XZ, 140));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_xZ, sendCount_xZ, rank_xZ, 143));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_Xz, recvCount_Xz, rank_Xz, 143));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_Xz, sendCount_Xz, rank_Xz, 142));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_xZ, recvCount_xZ, rank_xZ, 142));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_XZ, sendCount_XZ, rank_XZ, 141));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_xz, recvCount_xz, rank_xz, 141));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_yz, sendCount_yz, rank_yz, 144));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_YZ, recvCount_YZ, rank_YZ, 144));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_yZ, sendCount_yZ, rank_yZ, 147));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_Yz, recvCount_Yz, rank_Yz, 147));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_Yz, sendCount_Yz, rank_Yz, 146));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_yZ, recvCount_yZ, rank_yZ, 146));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Isend_init(sendbuf_YZ, sendCount_YZ, rank_YZ, 145));
    req_D3Q19AA.push_back(
        MPI_COMM_SCALBL.Irecv_init(recvbuf_yz, recvCount_yz, rank_yz, 145));
}

ScaLBL_Communicator::~ScaLBL_Communicator() {

    ScaLBL_FreeDeviceMemory(sendbuf_x);
    ScaLBL_FreeDeviceMemory(sendbuf_X);
    ScaLBL_FreeDeviceMemory(sendbuf_y);
    ScaLBL_FreeDeviceMemory(sendbuf_Y);
    ScaLBL_FreeDeviceMemory(sendbuf_z);
    ScaLBL_FreeDeviceMemory(sendbuf_Z);
    ScaLBL_FreeDeviceMemory(sendbuf_xy);
    ScaLBL_FreeDeviceMemory(sendbuf_xY);
    ScaLBL_FreeDeviceMemory(sendbuf_Xy);
    ScaLBL_FreeDeviceMemory(sendbuf_XY);
    ScaLBL_FreeDeviceMemory(sendbuf_xz);
    ScaLBL_FreeDeviceMemory(sendbuf_xZ);
    ScaLBL_FreeDeviceMemory(sendbuf_Xz);
    ScaLBL_FreeDeviceMemory(sendbuf_XZ);
    ScaLBL_FreeDeviceMemory(sendbuf_yz);
    ScaLBL_FreeDeviceMemory(sendbuf_yZ);
    ScaLBL_FreeDeviceMemory(sendbuf_Yz);
    ScaLBL_FreeDeviceMemory(sendbuf_YZ);
    ScaLBL_FreeDeviceMemory(recvbuf_x);
    ScaLBL_FreeDeviceMemory(recvbuf_X);
    ScaLBL_FreeDeviceMemory(recvbuf_y);
    ScaLBL_FreeDeviceMemory(recvbuf_Y);
    ScaLBL_FreeDeviceMemory(recvbuf_z);
    ScaLBL_FreeDeviceMemory(recvbuf_Z);
    ScaLBL_FreeDeviceMemory(recvbuf_xy);
    ScaLBL_FreeDeviceMemory(recvbuf_xY);
    ScaLBL_FreeDeviceMemory(recvbuf_Xy);
    ScaLBL_FreeDeviceMemory(recvbuf_XY);
    ScaLBL_FreeDeviceMemory(recvbuf_xz);
    ScaLBL_FreeDeviceMemory(recvbuf_xZ);
    ScaLBL_FreeDeviceMemory(recvbuf_Xz);
    ScaLBL_FreeDeviceMemory(recvbuf_XZ);
    ScaLBL_FreeDeviceMemory(recvbuf_yz);
    ScaLBL_FreeDeviceMemory(recvbuf_yZ);
    ScaLBL_FreeDeviceMemory(recvbuf_Yz);
    ScaLBL_FreeDeviceMemory(recvbuf_YZ);
    ScaLBL_FreeDeviceMemory(dvcSendList_x);
    ScaLBL_FreeDeviceMemory(dvcSendList_X);
    ScaLBL_FreeDeviceMemory(dvcSendList_y);
    ScaLBL_FreeDeviceMemory(dvcSendList_Y);
    ScaLBL_FreeDeviceMemory(dvcSendList_z);
    ScaLBL_FreeDeviceMemory(dvcSendList_Z);
    ScaLBL_FreeDeviceMemory(dvcSendList_xy);
    ScaLBL_FreeDeviceMemory(dvcSendList_xY);
    ScaLBL_FreeDeviceMemory(dvcSendList_Xy);
    ScaLBL_FreeDeviceMemory(dvcSendList_XY);
    ScaLBL_FreeDeviceMemory(dvcSendList_xz);
    ScaLBL_FreeDeviceMemory(dvcSendList_xZ);
    ScaLBL_FreeDeviceMemory(dvcSendList_Xz);
    ScaLBL_FreeDeviceMemory(dvcSendList_XZ);
    ScaLBL_FreeDeviceMemory(dvcSendList_yz);
    ScaLBL_FreeDeviceMemory(dvcSendList_yZ);
    ScaLBL_FreeDeviceMemory(dvcSendList_Yz);
    ScaLBL_FreeDeviceMemory(dvcSendList_YZ);
    ScaLBL_FreeDeviceMemory(dvcRecvList_x);
    ScaLBL_FreeDeviceMemory(dvcRecvList_X);
    ScaLBL_FreeDeviceMemory(dvcRecvList_y);
    ScaLBL_FreeDeviceMemory(dvcRecvList_Y);
    ScaLBL_FreeDeviceMemory(dvcRecvList_z);
    ScaLBL_FreeDeviceMemory(dvcRecvList_Z);
    ScaLBL_FreeDeviceMemory(dvcRecvList_xy);
    ScaLBL_FreeDeviceMemory(dvcRecvList_xY);
    ScaLBL_FreeDeviceMemory(dvcRecvList_Xy);
    ScaLBL_FreeDeviceMemory(dvcRecvList_XY);
    ScaLBL_FreeDeviceMemory(dvcRecvList_xz);
    ScaLBL_FreeDeviceMemory(dvcRecvList_xZ);
    ScaLBL_FreeDeviceMemory(dvcRecvList_Xz);
    ScaLBL_FreeDeviceMemory(dvcRecvList_XZ);
    ScaLBL_FreeDeviceMemory(dvcRecvList_yz);
    ScaLBL_FreeDeviceMemory(dvcRecvList_yZ);
    ScaLBL_FreeDeviceMemory(dvcRecvList_Yz);
    ScaLBL_FreeDeviceMemory(dvcRecvList_YZ);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_x);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_X);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_y);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_Y);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_z);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_Z);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_xy);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_xY);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_Xy);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_XY);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_xz);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_xZ);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_Xz);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_XZ);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_yz);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_yZ);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_Yz);
    ScaLBL_FreeDeviceMemory(dvcRecvDist_YZ);
}

void ScaLBL_Communicator::start(
    std::vector<std::shared_ptr<MPI_Request>> &requests) {
    for (auto &req : requests)
        MPI_COMM_SCALBL.Start(*req);
}
void ScaLBL_Communicator::wait(
    std::vector<std::shared_ptr<MPI_Request>> &requests) {
    std::vector<MPI_Request> request2;
    for (auto &req : requests)
        request2.push_back(*req);
    MPI_COMM_SCALBL.waitAll(request2.size(), request2.data());
}

/********************************************************
 * Get send/recv lists                                   *
 ********************************************************/
int ScaLBL_Communicator::copyRecvList(const char *dir, int *buffer) {
    if (dir[0] == 'x') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[recvCount_x];
            ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_x,
                              recvCount_x * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  recvCount_x * sizeof(int));
            delete[] TempBuffer;
            return recvCount_x;
        } else if (dir[1] == 'y')
            return recvCount_xy;
        else if (dir[1] == 'Y')
            return recvCount_xY;
        else if (dir[1] == 'z')
            return recvCount_xz;
        else if (dir[1] == 'Z')
            return recvCount_xZ;
    } else if (dir[0] == 'y') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[recvCount_y];
            ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_y,
                              recvCount_y * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  recvCount_y * sizeof(int));
            delete[] TempBuffer;
            return recvCount_y;
        } else if (dir[1] == 'z')
            return recvCount_yz;
        else if (dir[1] == 'Z')
            return recvCount_yZ;
    } else if (dir[0] == 'z') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[recvCount_z];
            ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_z,
                              recvCount_z * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  recvCount_z * sizeof(int));
            delete[] TempBuffer;
            return recvCount_z;
        }
    } else if (dir[0] == 'X') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[recvCount_X];
            ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_X,
                              recvCount_X * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  recvCount_X * sizeof(int));
            delete[] TempBuffer;
            return recvCount_X;
        } else if (dir[1] == 'y')
            return recvCount_Xy;
        else if (dir[1] == 'Y')
            return recvCount_XY;
        else if (dir[1] == 'z')
            return recvCount_Xz;
        else if (dir[1] == 'Z')
            return recvCount_XZ;
    } else if (dir[0] == 'Y') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[recvCount_Y];
            ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_Y,
                              recvCount_Y * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  recvCount_Y * sizeof(int));
            delete[] TempBuffer;
            return recvCount_Y;
        } else if (dir[1] == 'z')
            return recvCount_Yz;
        else if (dir[1] == 'Z')
            return recvCount_YZ;
    } else if (dir[0] == 'Z') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[recvCount_Z];
            ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_Z,
                              recvCount_Z * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  recvCount_Z * sizeof(int));
            delete[] TempBuffer;
            return recvCount_Z;
        }
    }
    throw std::logic_error("Internal error");
}

int ScaLBL_Communicator::copySendList(const char *dir, int *buffer) {
    if (dir[0] == 'x') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[sendCount_x];
            ScaLBL_CopyToHost(TempBuffer, dvcSendList_x,
                              sendCount_x * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  sendCount_x * sizeof(int));
            delete[] TempBuffer;
            return sendCount_x;
        } else if (dir[1] == 'y')
            return sendCount_xy;
        else if (dir[1] == 'Y')
            return sendCount_xY;
        else if (dir[1] == 'z')
            return sendCount_xz;
        else if (dir[1] == 'Z')
            return sendCount_xZ;
    } else if (dir[0] == 'y') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[sendCount_y];
            ScaLBL_CopyToHost(TempBuffer, dvcSendList_y,
                              sendCount_y * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  sendCount_y * sizeof(int));
            delete[] TempBuffer;
            return sendCount_y;
        } else if (dir[1] == 'z')
            return sendCount_yz;
        else if (dir[1] == 'Z')
            return sendCount_yZ;
    } else if (dir[0] == 'z') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[sendCount_z];
            ScaLBL_CopyToHost(TempBuffer, dvcSendList_z,
                              sendCount_z * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  sendCount_z * sizeof(int));
            delete[] TempBuffer;
            return sendCount_z;
        }
    } else if (dir[0] == 'X') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[sendCount_X];
            ScaLBL_CopyToHost(TempBuffer, dvcSendList_X,
                              sendCount_X * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  sendCount_X * sizeof(int));
            delete[] TempBuffer;
            return sendCount_X;
        } else if (dir[1] == 'y')
            return sendCount_Xy;
        else if (dir[1] == 'Y')
            return sendCount_XY;
        else if (dir[1] == 'z')
            return sendCount_Xz;
        else if (dir[1] == 'Z')
            return sendCount_XZ;
    } else if (dir[0] == 'Y') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[sendCount_Y];
            ScaLBL_CopyToHost(TempBuffer, dvcSendList_Y,
                              sendCount_Y * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  sendCount_Y * sizeof(int));
            delete[] TempBuffer;
            return sendCount_Y;
        } else if (dir[1] == 'z')
            return sendCount_Yz;
        else if (dir[1] == 'Z')
            return sendCount_YZ;
    } else if (dir[0] == 'Z') {
        if (dir[1] == 0) {
            int *TempBuffer = new int[sendCount_Z];
            ScaLBL_CopyToHost(TempBuffer, dvcSendList_Z,
                              sendCount_Z * sizeof(int));
            ScaLBL_CopyToZeroCopy(buffer, TempBuffer,
                                  sendCount_Z * sizeof(int));
            delete[] TempBuffer;
            return sendCount_Z;
        }
    }
    throw std::logic_error("Internal error");
}

double ScaLBL_Communicator::GetPerformance(int *NeighborList, double *fq,
                                           int Np) {
    /* EACH MPI PROCESS GETS ITS OWN MEASUREMENT*/
    /* use MRT kernels to check performance without communication / synchronization */
    int TIMESTEPS = 500;
    double RLX_SETA = 1.0;
    double RLX_SETB = 8.f * (2.f - RLX_SETA) / (8.f - RLX_SETA);
    double FX = 0.0;
    double FY = 0.0;
    double FZ = 0.0;
    ScaLBL_D3Q19_Init(fq, Np);
    //.......create and start timer............
    Barrier();
    auto t1 = std::chrono::system_clock::now();
    for (int t = 0; t < TIMESTEPS; t++) {
        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, FirstInterior(),
                               LastInterior(), Np, RLX_SETA, RLX_SETB, FX, FY,
                               FZ);
        ScaLBL_D3Q19_AAodd_MRT(NeighborList, fq, 0, LastExterior(), Np,
                               RLX_SETA, RLX_SETB, FX, FY, FZ);
        ScaLBL_D3Q19_AAeven_MRT(fq, FirstInterior(), LastInterior(), Np,
                                RLX_SETA, RLX_SETB, FX, FY, FZ);
        ScaLBL_D3Q19_AAeven_MRT(fq, 0, LastExterior(), Np, RLX_SETA, RLX_SETB,
                                FX, FY, FZ);
    }
    auto t2 = std::chrono::system_clock::now();
    Barrier();
    // Compute the walltime per timestep
    double diff = std::chrono::duration<double>(t2 - t1).count();
    double cputime = 0.5 * diff / TIMESTEPS;
    // Performance obtained from each node
    double MLUPS = double(Np) / cputime / 1000000;
    return MLUPS;
}
int ScaLBL_Communicator::LastExterior() { return next; }
int ScaLBL_Communicator::FirstInterior() { return first_interior; }
int ScaLBL_Communicator::LastInterior() { return last_interior; }

void ScaLBL_Communicator::D3Q19_MapRecv(int Cqx, int Cqy, int Cqz,
                                        const int *list, int start, int count,
                                        int *d3q19_recvlist) {
    int i, j, k, n, nn, idx;
    int *ReturnDist;
    ReturnDist = new int[count];

    for (idx = 0; idx < count; idx++) {

        // Get the value from the list -- note that n is the index is from the send (non-local) process
        n = list[idx]; // if (rank == 0) printf("@ rank:%d n=%d\n",rank,n);
        // Get the 3-D indices from the send process
        k = n / (Nx * Ny);
        j = (n - Nx * Ny * k) / Nx;
        i = n - Nx * Ny * k - Nx * j;
        // if (rank ==0) printf("@ Get 3D indices from the send process: i=%d, j=%d, k=%d\n",i,j,k);

        // Streaming for the non-local distribution
        i += Cqx;
        j += Cqy;
        k += Cqz;
        // if (rank == 0) printf("@ Streaming for the non-local distribution: i=%d, j=%d, k=%d\n",i,j,k);

        // Compute 1D index for the neighbor and save
        nn = k * Nx * Ny + j * Nx + i;
        // if (rank == 0) printf("@ rank:%d: neighbor=%d\n",rank,nn);
        ReturnDist[idx] = nn;
    }

    // Return updated version to the device
    ScaLBL_CopyToDevice(&d3q19_recvlist[start], ReturnDist,
                        count * sizeof(int));

    // clean up the work arrays
    delete[] ReturnDist;
}

int ScaLBL_Communicator::MemoryOptimizedLayoutAA(IntArray &Map,
                                                 int *neighborList,
                                                 signed char *id, int Np,
                                                 int width) {
    /*
	 * Generate a memory optimized layout
	 *   id[n] == 0 implies that site n should be ignored (treat as a mask)
	 *   Map(i,j,k) = idx  <- this is the index for the memory optimized layout
	 *   neighborList(idx) <-stores the neighbors for the D3Q19 model
	 *   note that the number of communications remains the same
	 *   the index in the Send and Recv lists is also updated
	 *   this means that the commuincations are no longer valid for regular data structures
	 */
    int idx, i, j, k, n;

    // Check that Map has size matching sub-domain
    if ((int)Map.size(0) != Nx)
        ERROR("ScaLBL_Communicator::MemoryOptimizedLayout: Map array "
              "dimensions do not match! \n");

    // Initialize Map
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (id[n] > 0)
                    Map(i, j, k) =
                        -2; // this label is for parallel communication sites
                else
                    Map(i, j, k) =
                        -1; // this label is for solid bounce-back sites
            }
        }
    }

    //printf("Exterior... \n");

    // ********* Exterior **********
    // Step 1/2: Index the outer walls of the grid only
    idx = 0;
    next = 0;
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                // domain interior
                Map(i, j, k) = -1;
                // Local index
                n = k * Nx * Ny + j * Nx + i;
                if (id[n] > 0) {
                    // Counts for the six faces
                    if (i > 0 && i <= width)
                        Map(n) = idx++;
                    else if (j > 0 && j <= width)
                        Map(n) = idx++;
                    else if (k > 0 && k <= width)
                        Map(n) = idx++;
                    else if (i > Nx - width - 2 && i < Nx - 1)
                        Map(n) = idx++;
                    else if (j > Ny - width - 2 && j < Ny - 1)
                        Map(n) = idx++;
                    else if (k > Nz - width - 2 && k < Nz - 1)
                        Map(n) = idx++;
                }
            }
        }
    }
    next = idx;

    // ********* Interior **********
    // align the next read
    first_interior = (next / 16 + 1) * 16;
    idx = first_interior;
    // Step 2/2: Next loop over the domain interior in block-cyclic fashion
    for (k = width + 1; k < Nz - width - 1; k++) {
        for (j = width + 1; j < Ny - width - 1; j++) {
            for (i = width + 1; i < Nx - width - 1; i++) {
                // Local index (regular layout)
                n = k * Nx * Ny + j * Nx + i;
                if (id[n] > 0) {
                    Map(n) = idx++;
                    //neighborList[idx++] = n; // index of self in regular layout
                }
            }
        }
    }
    last_interior = idx;

    Np = (last_interior / 16 + 1) * 16;
    //printf("    Np=%i \n",Np);

    // Now use Map to determine the neighbors for each lattice direction
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                n = k * Nx * Ny + j * Nx + i;
                idx = Map(i, j, k);
                if (idx > Np)
                    printf("ScaLBL_Communicator::MemoryOptimizedLayout: "
                           "Map(%i,%i,%i) = %i > %i \n",
                           i, j, k, Map(i, j, k), Np);
                else if (!(idx < 0)) {
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
                    int neighbor; // cycle through the neighbors of lattice site idx
                    neighbor = Map(i - 1, j, k);
                    if (neighbor < 0)
                        neighborList[idx] = idx + 2 * Np;
                    else
                        neighborList[idx] = neighbor + 1 * Np;

                    neighbor = Map(i + 1, j, k);
                    if (neighbor < 0)
                        neighborList[Np + idx] = idx + 1 * Np;
                    else
                        neighborList[Np + idx] = neighbor + 2 * Np;

                    neighbor = Map(i, j - 1, k);
                    if (neighbor < 0)
                        neighborList[2 * Np + idx] = idx + 4 * Np;
                    else
                        neighborList[2 * Np + idx] = neighbor + 3 * Np;

                    neighbor = Map(i, j + 1, k);
                    if (neighbor < 0)
                        neighborList[3 * Np + idx] = idx + 3 * Np;
                    else
                        neighborList[3 * Np + idx] = neighbor + 4 * Np;

                    neighbor = Map(i, j, k - 1);
                    if (neighbor < 0)
                        neighborList[4 * Np + idx] = idx + 6 * Np;
                    else
                        neighborList[4 * Np + idx] = neighbor + 5 * Np;

                    neighbor = Map(i, j, k + 1);
                    if (neighbor < 0)
                        neighborList[5 * Np + idx] = idx + 5 * Np;
                    else
                        neighborList[5 * Np + idx] = neighbor + 6 * Np;

                    neighbor = Map(i - 1, j - 1, k);
                    if (neighbor < 0)
                        neighborList[6 * Np + idx] = idx + 8 * Np;
                    else
                        neighborList[6 * Np + idx] = neighbor + 7 * Np;

                    neighbor = Map(i + 1, j + 1, k);
                    if (neighbor < 0)
                        neighborList[7 * Np + idx] = idx + 7 * Np;
                    else
                        neighborList[7 * Np + idx] = neighbor + 8 * Np;

                    neighbor = Map(i - 1, j + 1, k);
                    if (neighbor < 0)
                        neighborList[8 * Np + idx] = idx + 10 * Np;
                    else
                        neighborList[8 * Np + idx] = neighbor + 9 * Np;

                    neighbor = Map(i + 1, j - 1, k);
                    if (neighbor < 0)
                        neighborList[9 * Np + idx] = idx + 9 * Np;
                    else
                        neighborList[9 * Np + idx] = neighbor + 10 * Np;

                    neighbor = Map(i - 1, j, k - 1);
                    if (neighbor < 0)
                        neighborList[10 * Np + idx] = idx + 12 * Np;
                    else
                        neighborList[10 * Np + idx] = neighbor + 11 * Np;

                    neighbor = Map(i + 1, j, k + 1);
                    if (neighbor < 0)
                        neighborList[11 * Np + idx] = idx + 11 * Np;
                    else
                        neighborList[11 * Np + idx] = neighbor + 12 * Np;

                    neighbor = Map(i - 1, j, k + 1);
                    if (neighbor < 0)
                        neighborList[12 * Np + idx] = idx + 14 * Np;
                    else
                        neighborList[12 * Np + idx] = neighbor + 13 * Np;

                    neighbor = Map(i + 1, j, k - 1);
                    if (neighbor < 0)
                        neighborList[13 * Np + idx] = idx + 13 * Np;
                    else
                        neighborList[13 * Np + idx] = neighbor + 14 * Np;

                    neighbor = Map(i, j - 1, k - 1);
                    if (neighbor < 0)
                        neighborList[14 * Np + idx] = idx + 16 * Np;
                    else
                        neighborList[14 * Np + idx] = neighbor + 15 * Np;

                    neighbor = Map(i, j + 1, k + 1);
                    if (neighbor < 0)
                        neighborList[15 * Np + idx] = idx + 15 * Np;
                    else
                        neighborList[15 * Np + idx] = neighbor + 16 * Np;

                    neighbor = Map(i, j - 1, k + 1);
                    if (neighbor < 0)
                        neighborList[16 * Np + idx] = idx + 18 * Np;
                    else
                        neighborList[16 * Np + idx] = neighbor + 17 * Np;

                    neighbor = Map(i, j + 1, k - 1);
                    if (neighbor < 0)
                        neighborList[17 * Np + idx] = idx + 17 * Np;
                    else
                        neighborList[17 * Np + idx] = neighbor + 18 * Np;
                }
            }
        }
    }

    //for (idx=0; idx<Np; idx++)	printf("%i: %i %i\n", idx, neighborList[Np],  neighborList[Np+idx]);
    //.......................................................................
    // Now map through  SendList and RecvList to update indices
    // First loop over the send lists

    int *TempBuffer;
    TempBuffer = new int[5 * RecvCount];

    //.......................................................................
    // Re-index the send lists
    ScaLBL_CopyToHost(TempBuffer, dvcSendList_x, sendCount_x * sizeof(int));
    for (i = 0; i < sendCount_x; i++) {
        n = TempBuffer[i];
        //if (rank==0) printf("s: n=%d ",n);
        idx = Map(n);
        //if (rank == 0) printf("s: mapped n=%d\n",idx);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_x, TempBuffer, sendCount_x * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_y, sendCount_y * sizeof(int));
    for (i = 0; i < sendCount_y; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_y, TempBuffer, sendCount_y * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_z, sendCount_z * sizeof(int));
    for (i = 0; i < sendCount_z; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_z, TempBuffer, sendCount_z * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_X, sendCount_X * sizeof(int));
    for (i = 0; i < sendCount_X; i++) {
        n = TempBuffer[i];
        //if (rank==0) printf("r: n=%d ",n);
        idx = Map(n);
        //if (rank == 0) printf("r: mapped n=%d\n",idx);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_X, TempBuffer, sendCount_X * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_Y, sendCount_Y * sizeof(int));
    for (i = 0; i < sendCount_Y; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_Y, TempBuffer, sendCount_Y * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_Z, sendCount_Z * sizeof(int));
    for (i = 0; i < sendCount_Z; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_Z, TempBuffer, sendCount_Z * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_xy, sendCount_xy * sizeof(int));
    for (i = 0; i < sendCount_xy; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_xy, TempBuffer, sendCount_xy * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_xY, sendCount_xY * sizeof(int));
    for (i = 0; i < sendCount_xY; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_xY, TempBuffer, sendCount_xY * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_Xy, sendCount_Xy * sizeof(int));
    for (i = 0; i < sendCount_Xy; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_Xy, TempBuffer, sendCount_Xy * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_XY, sendCount_XY * sizeof(int));
    for (i = 0; i < sendCount_XY; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_XY, TempBuffer, sendCount_XY * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_xz, sendCount_xz * sizeof(int));
    for (i = 0; i < sendCount_xz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_xz, TempBuffer, sendCount_xz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_xZ, sendCount_xZ * sizeof(int));
    for (i = 0; i < sendCount_xZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_xZ, TempBuffer, sendCount_xZ * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_Xz, sendCount_Xz * sizeof(int));
    for (i = 0; i < sendCount_Xz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_Xz, TempBuffer, sendCount_Xz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_XZ, sendCount_XZ * sizeof(int));
    for (i = 0; i < sendCount_XZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_XZ, TempBuffer, sendCount_XZ * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_yz, sendCount_yz * sizeof(int));
    for (i = 0; i < sendCount_yz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_yz, TempBuffer, sendCount_yz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_Yz, sendCount_Yz * sizeof(int));
    for (i = 0; i < sendCount_Yz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_Yz, TempBuffer, sendCount_Yz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_yZ, sendCount_yZ * sizeof(int));
    for (i = 0; i < sendCount_yZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_yZ, TempBuffer, sendCount_yZ * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcSendList_YZ, sendCount_YZ * sizeof(int));
    for (i = 0; i < sendCount_YZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcSendList_YZ, TempBuffer, sendCount_YZ * sizeof(int));
    //.......................................................................
    // Re-index the recieve lists for the D3Q19 distributions
    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_x, 5 * recvCount_x * sizeof(int));
    for (i = 0; i < 5 * recvCount_x; i++) {
        n = TempBuffer[i];
        //if (rank==0) printf("r: n=%d ",n);
        idx = Map(n);
        //if (rank == 0) printf("r: mapped n=%d\n",idx);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_x, TempBuffer,
                        5 * recvCount_x * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_y, 5 * recvCount_y * sizeof(int));
    for (i = 0; i < 5 * recvCount_y; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }

    ScaLBL_CopyToDevice(dvcRecvDist_y, TempBuffer,
                        5 * recvCount_y * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_z, 5 * recvCount_z * sizeof(int));
    for (i = 0; i < 5 * recvCount_z; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_z, TempBuffer,
                        5 * recvCount_z * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_X, 5 * recvCount_X * sizeof(int));
    for (i = 0; i < 5 * recvCount_X; i++) {
        n = TempBuffer[i];
        //if (rank==0) printf("r: n=%d ",n);
        idx = Map(n);
        //if (rank == 0) printf("r: mapped n=%d\n",idx);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_X, TempBuffer,
                        5 * recvCount_X * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_Y, 5 * recvCount_Y * sizeof(int));
    for (i = 0; i < 5 * recvCount_Y; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_Y, TempBuffer,
                        5 * recvCount_Y * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_Z, 5 * recvCount_Z * sizeof(int));
    for (i = 0; i < 5 * recvCount_Z; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_Z, TempBuffer,
                        5 * recvCount_Z * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_xy, recvCount_xy * sizeof(int));
    for (i = 0; i < recvCount_xy; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_xy, TempBuffer, recvCount_xy * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_xY, recvCount_xY * sizeof(int));
    for (i = 0; i < recvCount_xY; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_xY, TempBuffer, recvCount_xY * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_Xy, recvCount_Xy * sizeof(int));
    for (i = 0; i < recvCount_Xy; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_Xy, TempBuffer, recvCount_Xy * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_XY, recvCount_XY * sizeof(int));
    for (i = 0; i < recvCount_XY; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_XY, TempBuffer, recvCount_XY * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_xz, recvCount_xz * sizeof(int));
    for (i = 0; i < recvCount_xz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_xz, TempBuffer, recvCount_xz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_xZ, recvCount_xZ * sizeof(int));
    for (i = 0; i < recvCount_xZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_xZ, TempBuffer, recvCount_xZ * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_Xz, recvCount_Xz * sizeof(int));
    for (i = 0; i < recvCount_Xz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_Xz, TempBuffer, recvCount_Xz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_XZ, recvCount_XZ * sizeof(int));
    for (i = 0; i < recvCount_XZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_XZ, TempBuffer, recvCount_XZ * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_yz, recvCount_yz * sizeof(int));
    for (i = 0; i < recvCount_yz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_yz, TempBuffer, recvCount_yz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_Yz, recvCount_Yz * sizeof(int));
    for (i = 0; i < recvCount_Yz; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_Yz, TempBuffer, recvCount_Yz * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_yZ, recvCount_yZ * sizeof(int));
    for (i = 0; i < recvCount_yZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_yZ, TempBuffer, recvCount_yZ * sizeof(int));

    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_YZ, recvCount_YZ * sizeof(int));
    for (i = 0; i < recvCount_YZ; i++) {
        n = TempBuffer[i];
        idx = Map(n);
        TempBuffer[i] = idx;
    }
    ScaLBL_CopyToDevice(dvcRecvDist_YZ, TempBuffer, recvCount_YZ * sizeof(int));
    //.......................................................................

    // Reset the value of N to match the dense structure
    N = Np;

    // Clean up
    delete[] TempBuffer;
    return (Np);
}

void ScaLBL_Communicator::SetupBounceBackList(IntArray &Map, signed char *id,
                                              int Np, bool SlippingVelBC) {

    int idx, i, j, k;
    int neighbor;
    // save list of bounce-back distributions and interaction sites
    n_bb_d3q7 = 0;
    n_bb_d3q19 = 0;

    int local_count = 0;
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                n = k * Nx * Ny + j * Nx + i;
                idx = Map(i, j, k);
                if (!(idx < 0)) {

                    neighbor = Map(i - 1, j, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i + 1, j, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j - 1, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j + 1, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j, k - 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j, k + 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i - 1, j - 1, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i + 1, j + 1, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i - 1, j + 1, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i + 1, j - 1, k);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i - 1, j, k - 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i + 1, j, k + 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i - 1, j, k + 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i + 1, j, k - 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j - 1, k - 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j + 1, k + 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j - 1, k + 1);
                    if (neighbor == -1)
                        local_count++;

                    neighbor = Map(i, j + 1, k - 1);
                    if (neighbor == -1)
                        local_count++;
                }
            }
        }
    }
    if (local_count > 0) {

        int *bb_dist_tmp = new int[local_count];
        int *bb_interactions_tmp = new int[local_count];
        ScaLBL_AllocateDeviceMemory((void **)&bb_dist,
                                    sizeof(int) * local_count);
        ScaLBL_AllocateDeviceMemory((void **)&bb_interactions,
                                    sizeof(int) * local_count);
        int *fluid_boundary_tmp;
        double *lattice_weight_tmp;
        float *lattice_cx_tmp;
        float *lattice_cy_tmp;
        float *lattice_cz_tmp;
        /* allocate memory for bounce-back sites */
        fluid_boundary_tmp = new int[local_count];
        lattice_weight_tmp = new double[local_count];
        lattice_cx_tmp = new float[local_count];
        lattice_cy_tmp = new float[local_count];
        lattice_cz_tmp = new float[local_count];
        ScaLBL_AllocateDeviceMemory((void **)&fluid_boundary,
                                    sizeof(int) * local_count);
        ScaLBL_AllocateDeviceMemory((void **)&lattice_weight,
                                    sizeof(double) * local_count);
        ScaLBL_AllocateDeviceMemory((void **)&lattice_cx,
                                    sizeof(float) * local_count);
        ScaLBL_AllocateDeviceMemory((void **)&lattice_cy,
                                    sizeof(float) * local_count);
        ScaLBL_AllocateDeviceMemory((void **)&lattice_cz,
                                    sizeof(float) * local_count);

        local_count = 0;
        for (k = 1; k < Nz - 1; k++) {
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    n = k * Nx * Ny + j * Nx + i;
                    idx = Map(i, j, k);
                    if (!(idx < 0)) {

                        int neighbor; // cycle through the neighbors of lattice site idx
                        neighbor = Map(i - 1, j, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i - 1) + (j)*Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 18.0;
                            lattice_cx_tmp[local_count] = -1.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 2 * Np;
                        }

                        neighbor = Map(i + 1, j, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i + 1) + (j)*Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 18.0;
                            lattice_cx_tmp[local_count] = 1.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 1 * Np;
                        }

                        neighbor = Map(i, j - 1, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j - 1) * Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 18.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = -1.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 4 * Np;
                        }

                        neighbor = Map(i, j + 1, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j + 1) * Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 18.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = 1.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 3 * Np;
                        }

                        neighbor = Map(i, j, k - 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j)*Nx + (k - 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 18.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = -1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 6 * Np;
                        }

                        neighbor = Map(i, j, k + 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j)*Nx + (k + 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 18.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = 1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 5 * Np;
                        }
                    }
                }
            }
        }
        n_bb_d3q7 = local_count;
        for (k = 1; k < Nz - 1; k++) {
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    n = k * Nx * Ny + j * Nx + i;
                    idx = Map(i, j, k);
                    if (!(idx < 0)) {

                        neighbor = Map(i - 1, j - 1, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i - 1) + (j - 1) * Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = -1.0;
                            lattice_cy_tmp[local_count] = -1.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 8 * Np;
                        }

                        neighbor = Map(i + 1, j + 1, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i + 1) + (j + 1) * Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 1.0;
                            lattice_cy_tmp[local_count] = 1.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 7 * Np;
                        }

                        neighbor = Map(i - 1, j + 1, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i - 1) + (j + 1) * Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = -1.0;
                            lattice_cy_tmp[local_count] = 1.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 10 * Np;
                        }

                        neighbor = Map(i + 1, j - 1, k);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i + 1) + (j - 1) * Nx + (k)*Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 1.0;
                            lattice_cy_tmp[local_count] = -1.0;
                            lattice_cz_tmp[local_count] = 0.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 9 * Np;
                        }

                        neighbor = Map(i - 1, j, k - 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i - 1) + (j)*Nx + (k - 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = -1.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = -1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 12 * Np;
                        }

                        neighbor = Map(i + 1, j, k + 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i + 1) + (j)*Nx + (k + 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 1.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = 1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 11 * Np;
                        }

                        neighbor = Map(i - 1, j, k + 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i - 1) + (j)*Nx + (k + 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = -1.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = 1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 14 * Np;
                        }

                        neighbor = Map(i + 1, j, k - 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i + 1) + (j)*Nx + (k - 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 1.0;
                            lattice_cy_tmp[local_count] = 0.0;
                            lattice_cz_tmp[local_count] = -1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 13 * Np;
                        }

                        neighbor = Map(i, j - 1, k - 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j - 1) * Nx + (k - 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = -1.0;
                            lattice_cz_tmp[local_count] = -1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 16 * Np;
                        }

                        neighbor = Map(i, j + 1, k + 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j + 1) * Nx + (k + 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = 1.0;
                            lattice_cz_tmp[local_count] = 1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 15 * Np;
                        }

                        neighbor = Map(i, j - 1, k + 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j - 1) * Nx + (k + 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = -1.0;
                            lattice_cz_tmp[local_count] = 1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 18 * Np;
                        }

                        neighbor = Map(i, j + 1, k - 1);
                        if (neighbor == -1) {
                            bb_interactions_tmp[local_count] =
                                (i) + (j + 1) * Nx + (k - 1) * Nx * Ny;
                            //if(SlippingVelBC==true){
                            fluid_boundary_tmp[local_count] = idx;
                            lattice_weight_tmp[local_count] = 1.0 / 36.0;
                            lattice_cx_tmp[local_count] = 0.0;
                            lattice_cy_tmp[local_count] = 1.0;
                            lattice_cz_tmp[local_count] = -1.0;
                            //}
                            bb_dist_tmp[local_count++] = idx + 17 * Np;
                        }
                    }
                }
            }
        }
        n_bb_d3q19 =
            local_count; // this gives the d3q19 distributions not part of d3q7 model
        ScaLBL_CopyToDevice(bb_dist, bb_dist_tmp, local_count * sizeof(int));
        ScaLBL_CopyToDevice(bb_interactions, bb_interactions_tmp,
                            local_count * sizeof(int));
        ScaLBL_CopyToDevice(fluid_boundary, fluid_boundary_tmp,
                            local_count * sizeof(int));
        ScaLBL_CopyToDevice(lattice_weight, lattice_weight_tmp,
                            local_count * sizeof(double));
        ScaLBL_CopyToDevice(lattice_cx, lattice_cx_tmp,
                            local_count * sizeof(float));
        ScaLBL_CopyToDevice(lattice_cy, lattice_cy_tmp,
                            local_count * sizeof(float));
        ScaLBL_CopyToDevice(lattice_cz, lattice_cz_tmp,
                            local_count * sizeof(float));
        ScaLBL_DeviceBarrier();

        delete[] bb_dist_tmp;
        delete[] bb_interactions_tmp;
        delete[] fluid_boundary_tmp;
        delete[] lattice_weight_tmp;
        delete[] lattice_cx_tmp;
        delete[] lattice_cy_tmp;
        delete[] lattice_cz_tmp;
    } else {
        bb_dist = NULL;
        bb_interactions = NULL;
        fluid_boundary = NULL;
        lattice_weight = NULL;
        lattice_cx = NULL;
        lattice_cy = NULL;
        lattice_cz = NULL;
    }
}

void ScaLBL_Communicator::SolidDirichletD3Q7(double *fq,
                                             double *BoundaryValue) {
    // fq is a D3Q7 distribution
    // BoundaryValues is a list of values to assign at bounce-back solid sites
    ScaLBL_Solid_Dirichlet_D3Q7(fq, BoundaryValue, bb_dist, bb_interactions,
                                n_bb_d3q7);
}

void ScaLBL_Communicator::SolidNeumannD3Q7(double *fq, double *BoundaryValue) {
    // fq is a D3Q7 distribution
    // BoundaryValues is a list of values to assign at bounce-back solid sites
    ScaLBL_Solid_Neumann_D3Q7(fq, BoundaryValue, bb_dist, bb_interactions,
                              n_bb_d3q7);
}

void ScaLBL_Communicator::SolidDirichletAndNeumannD3Q7(double *fq,
                                                       double *BoundaryValue,
                                                       int *BoundaryLabel) {
    // fq is a D3Q7 distribution
    // BoundaryValues is a list of values to assign at bounce-back solid sites
    // BoundaryLabel: is a list of integer labels indicating the type of BCs
    //                1-> Dirichlet BC; 2-> Neumann BC.
    ScaLBL_Solid_DirichletAndNeumann_D3Q7(fq, BoundaryValue, BoundaryLabel,
                                          bb_dist, bb_interactions, n_bb_d3q7);
}

void ScaLBL_Communicator::SolidSlippingVelocityBCD3Q19(
    double *fq, double *zeta_potential, double *ElectricField,
    double *SolidGrad, double epsilon_LB, double tau, double rho0,
    double den_scale, double h, double time_conv) {
    // fq is a D3Q19 distribution
    // BoundaryValues is a list of values to assign at bounce-back solid sites
    ScaLBL_Solid_SlippingVelocityBC_D3Q19(
        fq, zeta_potential, ElectricField, SolidGrad, epsilon_LB, tau, rho0,
        den_scale, h, time_conv, bb_dist, bb_interactions, fluid_boundary,
        lattice_weight, lattice_cx, lattice_cy, lattice_cz, n_bb_d3q19, N);
}

void ScaLBL_Communicator::SendD3Q19AA(double *dist) {

    if (Lock == true) {
        ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did "
              "you forget to match Send/Recv calls?");
    } else {
        Lock = true;
    }
    ScaLBL_DeviceBarrier();
    // Pack the distributions
    //...Packing for x face(2,8,10,12,14)................................
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, 0, sendCount_x, sendbuf_x, dist, N);
    ScaLBL_D3Q19_Pack(8, dvcSendList_x, sendCount_x, sendCount_x, sendbuf_x,
                      dist, N);
    ScaLBL_D3Q19_Pack(10, dvcSendList_x, 2 * sendCount_x, sendCount_x,
                      sendbuf_x, dist, N);
    ScaLBL_D3Q19_Pack(12, dvcSendList_x, 3 * sendCount_x, sendCount_x,
                      sendbuf_x, dist, N);
    ScaLBL_D3Q19_Pack(14, dvcSendList_x, 4 * sendCount_x, sendCount_x,
                      sendbuf_x, dist, N);

    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, 0, sendCount_X, sendbuf_X, dist, N);
    ScaLBL_D3Q19_Pack(7, dvcSendList_X, sendCount_X, sendCount_X, sendbuf_X,
                      dist, N);
    ScaLBL_D3Q19_Pack(9, dvcSendList_X, 2 * sendCount_X, sendCount_X, sendbuf_X,
                      dist, N);
    ScaLBL_D3Q19_Pack(11, dvcSendList_X, 3 * sendCount_X, sendCount_X,
                      sendbuf_X, dist, N);
    ScaLBL_D3Q19_Pack(13, dvcSendList_X, 4 * sendCount_X, sendCount_X,
                      sendbuf_X, dist, N);

    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, 0, sendCount_y, sendbuf_y, dist, N);
    ScaLBL_D3Q19_Pack(8, dvcSendList_y, sendCount_y, sendCount_y, sendbuf_y,
                      dist, N);
    ScaLBL_D3Q19_Pack(9, dvcSendList_y, 2 * sendCount_y, sendCount_y, sendbuf_y,
                      dist, N);
    ScaLBL_D3Q19_Pack(16, dvcSendList_y, 3 * sendCount_y, sendCount_y,
                      sendbuf_y, dist, N);
    ScaLBL_D3Q19_Pack(18, dvcSendList_y, 4 * sendCount_y, sendCount_y,
                      sendbuf_y, dist, N);

    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, 0, sendCount_Y, sendbuf_Y, dist, N);
    ScaLBL_D3Q19_Pack(7, dvcSendList_Y, sendCount_Y, sendCount_Y, sendbuf_Y,
                      dist, N);
    ScaLBL_D3Q19_Pack(10, dvcSendList_Y, 2 * sendCount_Y, sendCount_Y,
                      sendbuf_Y, dist, N);
    ScaLBL_D3Q19_Pack(15, dvcSendList_Y, 3 * sendCount_Y, sendCount_Y,
                      sendbuf_Y, dist, N);
    ScaLBL_D3Q19_Pack(17, dvcSendList_Y, 4 * sendCount_Y, sendCount_Y,
                      sendbuf_Y, dist, N);

    //...Packing for z face(6,12,13,16,17)................................
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, 0, sendCount_z, sendbuf_z, dist, N);
    ScaLBL_D3Q19_Pack(12, dvcSendList_z, sendCount_z, sendCount_z, sendbuf_z,
                      dist, N);
    ScaLBL_D3Q19_Pack(13, dvcSendList_z, 2 * sendCount_z, sendCount_z,
                      sendbuf_z, dist, N);
    ScaLBL_D3Q19_Pack(16, dvcSendList_z, 3 * sendCount_z, sendCount_z,
                      sendbuf_z, dist, N);
    ScaLBL_D3Q19_Pack(17, dvcSendList_z, 4 * sendCount_z, sendCount_z,
                      sendbuf_z, dist, N);

    //...Packing for Z face(5,11,14,15,18)................................
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, 0, sendCount_Z, sendbuf_Z, dist, N);
    ScaLBL_D3Q19_Pack(11, dvcSendList_Z, sendCount_Z, sendCount_Z, sendbuf_Z,
                      dist, N);
    ScaLBL_D3Q19_Pack(14, dvcSendList_Z, 2 * sendCount_Z, sendCount_Z,
                      sendbuf_Z, dist, N);
    ScaLBL_D3Q19_Pack(15, dvcSendList_Z, 3 * sendCount_Z, sendCount_Z,
                      sendbuf_Z, dist, N);
    ScaLBL_D3Q19_Pack(18, dvcSendList_Z, 4 * sendCount_Z, sendCount_Z,
                      sendbuf_Z, dist, N);

    //...Pack the xy edge (8)................................
    ScaLBL_D3Q19_Pack(8, dvcSendList_xy, 0, sendCount_xy, sendbuf_xy, dist, N);
    //...Pack the Xy edge (9)................................
    ScaLBL_D3Q19_Pack(9, dvcSendList_Xy, 0, sendCount_Xy, sendbuf_Xy, dist, N);
    //...Pack the xY edge (10)................................
    ScaLBL_D3Q19_Pack(10, dvcSendList_xY, 0, sendCount_xY, sendbuf_xY, dist, N);
    //...Pack the XY edge (7)................................
    ScaLBL_D3Q19_Pack(7, dvcSendList_XY, 0, sendCount_XY, sendbuf_XY, dist, N);
    //...Pack the xz edge (12)................................
    ScaLBL_D3Q19_Pack(12, dvcSendList_xz, 0, sendCount_xz, sendbuf_xz, dist, N);

    //...Pack the xZ edge (14)................................
    ScaLBL_D3Q19_Pack(14, dvcSendList_xZ, 0, sendCount_xZ, sendbuf_xZ, dist, N);
    //...Pack the Xz edge (13)................................
    ScaLBL_D3Q19_Pack(13, dvcSendList_Xz, 0, sendCount_Xz, sendbuf_Xz, dist, N);

    //...Pack the XZ edge (11)................................
    ScaLBL_D3Q19_Pack(11, dvcSendList_XZ, 0, sendCount_XZ, sendbuf_XZ, dist, N);
    //...Pack the yz edge (16)................................
    ScaLBL_D3Q19_Pack(16, dvcSendList_yz, 0, sendCount_yz, sendbuf_yz, dist, N);
    //...Pack the yZ edge (18)................................
    ScaLBL_D3Q19_Pack(18, dvcSendList_yZ, 0, sendCount_yZ, sendbuf_yZ, dist, N);
    //...Pack the Yz edge (17)................................
    ScaLBL_D3Q19_Pack(17, dvcSendList_Yz, 0, sendCount_Yz, sendbuf_Yz, dist, N);
    //...Pack the YZ edge (15)................................
    ScaLBL_D3Q19_Pack(15, dvcSendList_YZ, 0, sendCount_YZ, sendbuf_YZ, dist, N);

    //...................................................................................

    ScaLBL_DeviceBarrier();
    start(req_D3Q19AA);
}

void ScaLBL_Communicator::RecvD3Q19AA(double *dist) {

    // NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
    //...................................................................................
    // Wait for completion of D3Q19 communication
    wait(req_D3Q19AA);
    ScaLBL_DeviceBarrier();

    //...................................................................................
    // NOTE: AA Routine writes to opposite
    // Unpack the distributions on the device
    //...................................................................................
    //...Unpacking for x face(2,8,10,12,14)................................
    ScaLBL_D3Q19_Unpack(2, dvcRecvDist_x, 0, recvCount_x, recvbuf_x, dist, N);
    ScaLBL_D3Q19_Unpack(8, dvcRecvDist_x, recvCount_x, recvCount_x, recvbuf_x,
                        dist, N);
    ScaLBL_D3Q19_Unpack(10, dvcRecvDist_x, 2 * recvCount_x, recvCount_x,
                        recvbuf_x, dist, N);
    ScaLBL_D3Q19_Unpack(12, dvcRecvDist_x, 3 * recvCount_x, recvCount_x,
                        recvbuf_x, dist, N);
    ScaLBL_D3Q19_Unpack(14, dvcRecvDist_x, 4 * recvCount_x, recvCount_x,
                        recvbuf_x, dist, N);
    //...................................................................................
    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q19_Unpack(1, dvcRecvDist_X, 0, recvCount_X, recvbuf_X, dist, N);
    ScaLBL_D3Q19_Unpack(7, dvcRecvDist_X, recvCount_X, recvCount_X, recvbuf_X,
                        dist, N);
    ScaLBL_D3Q19_Unpack(9, dvcRecvDist_X, 2 * recvCount_X, recvCount_X,
                        recvbuf_X, dist, N);
    ScaLBL_D3Q19_Unpack(11, dvcRecvDist_X, 3 * recvCount_X, recvCount_X,
                        recvbuf_X, dist, N);
    ScaLBL_D3Q19_Unpack(13, dvcRecvDist_X, 4 * recvCount_X, recvCount_X,
                        recvbuf_X, dist, N);
    //...................................................................................
    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q19_Unpack(4, dvcRecvDist_y, 0, recvCount_y, recvbuf_y, dist, N);
    ScaLBL_D3Q19_Unpack(8, dvcRecvDist_y, recvCount_y, recvCount_y, recvbuf_y,
                        dist, N);
    ScaLBL_D3Q19_Unpack(9, dvcRecvDist_y, 2 * recvCount_y, recvCount_y,
                        recvbuf_y, dist, N);
    ScaLBL_D3Q19_Unpack(16, dvcRecvDist_y, 3 * recvCount_y, recvCount_y,
                        recvbuf_y, dist, N);
    ScaLBL_D3Q19_Unpack(18, dvcRecvDist_y, 4 * recvCount_y, recvCount_y,
                        recvbuf_y, dist, N);
    //...................................................................................
    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q19_Unpack(3, dvcRecvDist_Y, 0, recvCount_Y, recvbuf_Y, dist, N);
    ScaLBL_D3Q19_Unpack(7, dvcRecvDist_Y, recvCount_Y, recvCount_Y, recvbuf_Y,
                        dist, N);
    ScaLBL_D3Q19_Unpack(10, dvcRecvDist_Y, 2 * recvCount_Y, recvCount_Y,
                        recvbuf_Y, dist, N);
    ScaLBL_D3Q19_Unpack(15, dvcRecvDist_Y, 3 * recvCount_Y, recvCount_Y,
                        recvbuf_Y, dist, N);
    ScaLBL_D3Q19_Unpack(17, dvcRecvDist_Y, 4 * recvCount_Y, recvCount_Y,
                        recvbuf_Y, dist, N);
    //...................................................................................

    //..................................................................................
    //...Pack the xy edge (8)................................
    ScaLBL_D3Q19_Unpack(8, dvcRecvDist_xy, 0, recvCount_xy, recvbuf_xy, dist,
                        N);
    //...Pack the Xy edge (9)................................
    ScaLBL_D3Q19_Unpack(9, dvcRecvDist_Xy, 0, recvCount_Xy, recvbuf_Xy, dist,
                        N);
    //...Pack the xY edge (10)................................
    ScaLBL_D3Q19_Unpack(10, dvcRecvDist_xY, 0, recvCount_xY, recvbuf_xY, dist,
                        N);
    //...Pack the XY edge (7)................................
    ScaLBL_D3Q19_Unpack(7, dvcRecvDist_XY, 0, recvCount_XY, recvbuf_XY, dist,
                        N);

    //if (BoundaryCondition == 0 || kproc != 0 ){
    ScaLBL_D3Q19_Unpack(6, dvcRecvDist_z, 0, recvCount_z, recvbuf_z, dist, N);
    ScaLBL_D3Q19_Unpack(12, dvcRecvDist_z, recvCount_z, recvCount_z, recvbuf_z,
                        dist, N);
    ScaLBL_D3Q19_Unpack(13, dvcRecvDist_z, 2 * recvCount_z, recvCount_z,
                        recvbuf_z, dist, N);
    ScaLBL_D3Q19_Unpack(16, dvcRecvDist_z, 3 * recvCount_z, recvCount_z,
                        recvbuf_z, dist, N);
    ScaLBL_D3Q19_Unpack(17, dvcRecvDist_z, 4 * recvCount_z, recvCount_z,
                        recvbuf_z, dist, N);

    //...Pack the xz edge (12)................................
    ScaLBL_D3Q19_Unpack(12, dvcRecvDist_xz, 0, recvCount_xz, recvbuf_xz, dist,
                        N);
    //...Pack the Xz edge (13)................................
    ScaLBL_D3Q19_Unpack(13, dvcRecvDist_Xz, 0, recvCount_Xz, recvbuf_Xz, dist,
                        N);
    //...Pack the yz edge (16)................................
    ScaLBL_D3Q19_Unpack(16, dvcRecvDist_yz, 0, recvCount_yz, recvbuf_yz, dist,
                        N);
    //...Pack the Yz edge (17)................................
    ScaLBL_D3Q19_Unpack(17, dvcRecvDist_Yz, 0, recvCount_Yz, recvbuf_Yz, dist,
                        N);

    //}
    //if (BoundaryCondition == 0 || kproc != nprocz-1){
    //...Packing for Z face(5,11,14,15,18)................................
    ScaLBL_D3Q19_Unpack(5, dvcRecvDist_Z, 0, recvCount_Z, recvbuf_Z, dist, N);
    ScaLBL_D3Q19_Unpack(11, dvcRecvDist_Z, recvCount_Z, recvCount_Z, recvbuf_Z,
                        dist, N);
    ScaLBL_D3Q19_Unpack(14, dvcRecvDist_Z, 2 * recvCount_Z, recvCount_Z,
                        recvbuf_Z, dist, N);
    ScaLBL_D3Q19_Unpack(15, dvcRecvDist_Z, 3 * recvCount_Z, recvCount_Z,
                        recvbuf_Z, dist, N);
    ScaLBL_D3Q19_Unpack(18, dvcRecvDist_Z, 4 * recvCount_Z, recvCount_Z,
                        recvbuf_Z, dist, N);

    //...Pack the xZ edge (14)................................
    ScaLBL_D3Q19_Unpack(14, dvcRecvDist_xZ, 0, recvCount_xZ, recvbuf_xZ, dist,
                        N);
    //...Pack the XZ edge (11)................................
    ScaLBL_D3Q19_Unpack(11, dvcRecvDist_XZ, 0, recvCount_XZ, recvbuf_XZ, dist,
                        N);
    //...Pack the yZ edge (18)................................
    ScaLBL_D3Q19_Unpack(18, dvcRecvDist_yZ, 0, recvCount_yZ, recvbuf_yZ, dist,
                        N);
    //...Pack the YZ edge (15)................................
    ScaLBL_D3Q19_Unpack(15, dvcRecvDist_YZ, 0, recvCount_YZ, recvbuf_YZ, dist,
                        N);
    //}

    //...................................................................................
    Lock = false; // unlock the communicator after communications complete
    //...................................................................................
}

void ScaLBL_Communicator::RecvGrad(double *phi, double *grad) {

    // Recieves halo and incorporates into D3Q19 based stencil gradient computation
    //...................................................................................
    // Wait for completion of D3Q19 communication
    MPI_COMM_SCALBL.waitAll(18, req1);
    MPI_COMM_SCALBL.waitAll(18, req2);
    ScaLBL_DeviceBarrier();

    //...................................................................................
    // Unpack the gradributions on the device
    //...................................................................................
    //...Unpacking for x face(2,8,10,12,14)................................
    ScaLBL_Gradient_Unpack(1.0, -1, 0, 0, dvcRecvDist_x, 0, recvCount_x,
                           recvbuf_x, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, -1, -1, 0, dvcRecvDist_x, recvCount_x,
                           recvCount_x, recvbuf_x, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, -1, 1, 0, dvcRecvDist_x, 2 * recvCount_x,
                           recvCount_x, recvbuf_x, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, -1, 0, 1, dvcRecvDist_x, 4 * recvCount_x,
                           recvCount_x, recvbuf_x, phi, grad, N);
    //...................................................................................
    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_Gradient_Unpack(1.0, 1, 0, 0, dvcRecvDist_X, 0, recvCount_X,
                           recvbuf_X, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, 1, 0, dvcRecvDist_X, recvCount_X,
                           recvCount_X, recvbuf_X, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, -1, 0, dvcRecvDist_X, 2 * recvCount_X,
                           recvCount_X, recvbuf_X, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, 0, 1, dvcRecvDist_X, 3 * recvCount_X,
                           recvCount_X, recvbuf_X, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, 0, -1, dvcRecvDist_X, 4 * recvCount_X,
                           recvCount_X, recvbuf_X, phi, grad, N);
    //...................................................................................
    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_Gradient_Unpack(1.0, 0, -1, 0, dvcRecvDist_y, 0, recvCount_y,
                           recvbuf_y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, -1, -1, 0, dvcRecvDist_y, recvCount_y,
                           recvCount_y, recvbuf_y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, -1, 0, dvcRecvDist_y, 2 * recvCount_y,
                           recvCount_y, recvbuf_y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, -1, -1, dvcRecvDist_y, 3 * recvCount_y,
                           recvCount_y, recvbuf_y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, -1, 1, dvcRecvDist_y, 4 * recvCount_y,
                           recvCount_y, recvbuf_y, phi, grad, N);
    //...................................................................................
    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_Gradient_Unpack(1.0, 0, 1, 0, dvcRecvDist_Y, 0, recvCount_Y,
                           recvbuf_Y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, 1, 0, dvcRecvDist_Y, recvCount_Y,
                           recvCount_Y, recvbuf_Y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, -1, 1, 0, dvcRecvDist_Y, 2 * recvCount_Y,
                           recvCount_Y, recvbuf_Y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, 1, 1, dvcRecvDist_Y, 3 * recvCount_Y,
                           recvCount_Y, recvbuf_Y, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, 1, -1, dvcRecvDist_Y, 4 * recvCount_Y,
                           recvCount_Y, recvbuf_Y, phi, grad, N);
    //...................................................................................
    //...Packing for z face(6,12,13,16,17)................................
    ScaLBL_Gradient_Unpack(1.0, 0, 0, -1, dvcRecvDist_z, 0, recvCount_z,
                           recvbuf_z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, -1, 0, -1, dvcRecvDist_z, recvCount_z,
                           recvCount_z, recvbuf_z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, 0, -1, dvcRecvDist_z, 2 * recvCount_z,
                           recvCount_z, recvbuf_z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, -1, -1, dvcRecvDist_z, 3 * recvCount_z,
                           recvCount_z, recvbuf_z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, 1, -1, dvcRecvDist_z, 4 * recvCount_z,
                           recvCount_z, recvbuf_z, phi, grad, N);
    //...Packing for Z face(5,11,14,15,18)................................
    ScaLBL_Gradient_Unpack(1.0, 0, 0, 1, dvcRecvDist_Z, 0, recvCount_Z,
                           recvbuf_Z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 1, 0, 1, dvcRecvDist_Z, recvCount_Z,
                           recvCount_Z, recvbuf_Z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, -1, 0, 1, dvcRecvDist_Z, 2 * recvCount_Z,
                           recvCount_Z, recvbuf_Z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, 1, 1, dvcRecvDist_Z, 3 * recvCount_Z,
                           recvCount_Z, recvbuf_Z, phi, grad, N);
    ScaLBL_Gradient_Unpack(0.5, 0, -1, 1, dvcRecvDist_Z, 4 * recvCount_Z,
                           recvCount_Z, recvbuf_Z, phi, grad, N);
    //..................................................................................
    //...Pack the xy edge (8)................................
    ScaLBL_Gradient_Unpack(0.5, -1, -1, 0, dvcRecvDist_xy, 0, recvCount_xy,
                           recvbuf_xy, phi, grad, N);
    //...Pack the Xy edge (9)................................
    ScaLBL_Gradient_Unpack(0.5, 1, -1, 0, dvcRecvDist_Xy, 0, recvCount_Xy,
                           recvbuf_Xy, phi, grad, N);
    //...Pack the xY edge (10)................................
    ScaLBL_Gradient_Unpack(0.5, -1, 1, 0, dvcRecvDist_xY, 0, recvCount_xY,
                           recvbuf_xY, phi, grad, N);
    //...Pack the XY edge (7)................................
    ScaLBL_Gradient_Unpack(0.5, 1, 1, 0, dvcRecvDist_XY, 0, recvCount_XY,
                           recvbuf_XY, phi, grad, N);
    //...Pack the xz edge (12)................................
    ScaLBL_Gradient_Unpack(0.5, -1, 0, -1, dvcRecvDist_xz, 0, recvCount_xz,
                           recvbuf_xz, phi, grad, N);
    //...Pack the xZ edge (14)................................
    ScaLBL_Gradient_Unpack(0.5, -1, 0, 1, dvcRecvDist_xZ, 0, recvCount_xZ,
                           recvbuf_xZ, phi, grad, N);
    //...Pack the Xz edge (13)................................
    ScaLBL_Gradient_Unpack(0.5, 1, 0, -1, dvcRecvDist_Xz, 0, recvCount_Xz,
                           recvbuf_Xz, phi, grad, N);
    //...Pack the XZ edge (11)................................
    ScaLBL_Gradient_Unpack(0.5, 1, 0, 1, dvcRecvDist_XZ, 0, recvCount_XZ,
                           recvbuf_XZ, phi, grad, N);
    //...Pack the yz edge (16)................................
    ScaLBL_Gradient_Unpack(0.5, 0, -1, -1, dvcRecvDist_yz, 0, recvCount_yz,
                           recvbuf_yz, phi, grad, N);
    //...Pack the yZ edge (18)................................
    ScaLBL_Gradient_Unpack(0.5, 0, -1, 1, dvcRecvDist_yZ, 0, recvCount_yZ,
                           recvbuf_yZ, phi, grad, N);
    //...Pack the Yz edge (17)................................
    ScaLBL_Gradient_Unpack(0.5, 0, 1, -1, dvcRecvDist_Yz, 0, recvCount_Yz,
                           recvbuf_Yz, phi, grad, N);
    //...Pack the YZ edge (15)................................
    ScaLBL_Gradient_Unpack(0.5, 0, 1, 1, dvcRecvDist_YZ, 0, recvCount_YZ,
                           recvbuf_YZ, phi, grad, N);
    //...................................................................................
    Lock = false; // unlock the communicator after communications complete
    //...................................................................................
}

void ScaLBL_Communicator::BiSendD3Q7AA(double *Aq, double *Bq) {

    // NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
    if (Lock == true) {
        ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did "
              "you forget to match Send/Recv calls?");
    } else {
        Lock = true;
    }
    // assign tag of 19 to D3Q19 communication
    sendtag = recvtag = 148;
    ScaLBL_DeviceBarrier();
    // Pack the distributions
    //...Packing for x face(2,8,10,12,14)................................
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, 0, sendCount_x, sendbuf_x, Aq, N);
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, sendCount_x, sendCount_x, sendbuf_x, Bq,
                      N);

    req1[0] =
        MPI_COMM_SCALBL.Isend(sendbuf_x, 2 * sendCount_x, rank_x, sendtag + 0);
    req2[0] =
        MPI_COMM_SCALBL.Irecv(recvbuf_X, 2 * recvCount_X, rank_X, recvtag + 0);

    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, 0, sendCount_X, sendbuf_X, Aq, N);
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, sendCount_X, sendCount_X, sendbuf_X, Bq,
                      N);

    req1[1] =
        MPI_COMM_SCALBL.Isend(sendbuf_X, 2 * sendCount_X, rank_X, sendtag + 1);
    req2[1] =
        MPI_COMM_SCALBL.Irecv(recvbuf_x, 2 * recvCount_x, rank_x, recvtag + 1);

    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, 0, sendCount_y, sendbuf_y, Aq, N);
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, sendCount_y, sendCount_y, sendbuf_y, Bq,
                      N);

    req1[2] =
        MPI_COMM_SCALBL.Isend(sendbuf_y, 2 * sendCount_y, rank_y, sendtag + 2);
    req2[2] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Y, 2 * recvCount_Y, rank_Y, recvtag + 2);

    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, 0, sendCount_Y, sendbuf_Y, Aq, N);
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, sendCount_Y, sendCount_Y, sendbuf_Y, Bq,
                      N);

    req1[3] =
        MPI_COMM_SCALBL.Isend(sendbuf_Y, 2 * sendCount_Y, rank_Y, sendtag + 3);
    req2[3] =
        MPI_COMM_SCALBL.Irecv(recvbuf_y, 2 * recvCount_y, rank_y, recvtag + 3);

    //...Packing for z face(6,12,13,16,17)................................
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, 0, sendCount_z, sendbuf_z, Aq, N);
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, sendCount_z, sendCount_z, sendbuf_z, Bq,
                      N);

    req1[4] =
        MPI_COMM_SCALBL.Isend(sendbuf_z, 2 * sendCount_z, rank_z, sendtag + 4);
    req2[4] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Z, 2 * recvCount_Z, rank_Z, recvtag + 4);

    //...Packing for Z face(5,11,14,15,18)................................
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, 0, sendCount_Z, sendbuf_Z, Aq, N);
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, sendCount_Z, sendCount_Z, sendbuf_Z, Bq,
                      N);

    //...................................................................................
    // Send all the distributions
    req1[5] =
        MPI_COMM_SCALBL.Isend(sendbuf_Z, 2 * sendCount_Z, rank_Z, sendtag + 5);
    req2[5] =
        MPI_COMM_SCALBL.Irecv(recvbuf_z, 2 * recvCount_z, rank_z, recvtag + 5);
}

void ScaLBL_Communicator::BiRecvD3Q7AA(double *Aq, double *Bq) {

    // NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
    //...................................................................................
    // Wait for completion of D3Q19 communication
    MPI_COMM_SCALBL.waitAll(6, req1);
    MPI_COMM_SCALBL.waitAll(6, req2);
    ScaLBL_DeviceBarrier();

    //...................................................................................
    // NOTE: AA Routine writes to opposite
    // Unpack the distributions on the device
    //...................................................................................
    //...Unpacking for x face(2,8,10,12,14)................................
    ScaLBL_D3Q7_Unpack(2, dvcRecvDist_x, 0, recvCount_x, recvbuf_x, Aq, N);
    ScaLBL_D3Q7_Unpack(2, dvcRecvDist_x, recvCount_x, recvCount_x, recvbuf_x,
                       Bq, N);
    //...................................................................................
    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q7_Unpack(1, dvcRecvDist_X, 0, recvCount_X, recvbuf_X, Aq, N);
    ScaLBL_D3Q7_Unpack(1, dvcRecvDist_X, recvCount_X, recvCount_X, recvbuf_X,
                       Bq, N);
    //...................................................................................
    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q7_Unpack(4, dvcRecvDist_y, 0, recvCount_y, recvbuf_y, Aq, N);
    ScaLBL_D3Q7_Unpack(4, dvcRecvDist_y, recvCount_y, recvCount_y, recvbuf_y,
                       Bq, N);
    //...................................................................................
    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q7_Unpack(3, dvcRecvDist_Y, 0, recvCount_Y, recvbuf_Y, Aq, N);
    ScaLBL_D3Q7_Unpack(3, dvcRecvDist_Y, recvCount_Y, recvCount_Y, recvbuf_Y,
                       Bq, N);
    //...................................................................................

    if (BoundaryCondition > 0 && kproc == 0) {
        // don't unpack little z
        //...Packing for Z face(5,11,14,15,18)................................
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 0, recvCount_Z, recvbuf_Z, Aq, N);
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, recvCount_Z, recvCount_Z,
                           recvbuf_Z, Bq, N);
    } else if (BoundaryCondition > 0 && kproc == nprocz - 1) {
        // don't unpack big z
        //...Packing for z face(6,12,13,16,17)................................
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 0, recvCount_z, recvbuf_z, Aq, N);
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, recvCount_z, recvCount_z,
                           recvbuf_z, Bq, N);
    } else {
        //...Packing for z face(6,12,13,16,17)................................
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 0, recvCount_z, recvbuf_z, Aq, N);
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, recvCount_z, recvCount_z,
                           recvbuf_z, Bq, N);
        //...Packing for Z face(5,11,14,15,18)................................
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 0, recvCount_Z, recvbuf_Z, Aq, N);
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, recvCount_Z, recvCount_Z,
                           recvbuf_Z, Bq, N);
    }

    //...................................................................................
    Lock = false; // unlock the communicator after communications complete
    //...................................................................................
}

void ScaLBL_Communicator::SendD3Q7AA(double *Aq, int Component) {

    // NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
    if (Lock == true) {
        ERROR("ScaLBL Error (SendD3Q7): ScaLBL_Communicator is locked -- did "
              "you forget to match Send/Recv calls?");
    } else {
        Lock = true;
    }
    // assign tag of 154 to D3Q19 communication
    sendtag = recvtag = 154;
    ScaLBL_DeviceBarrier();
    // Pack the distributions
    //...Packing for x face(2,8,10,12,14)................................
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, 0, sendCount_x, sendbuf_x,
                      &Aq[Component * 7 * N], N);
    req1[0] =
        MPI_COMM_SCALBL.Isend(sendbuf_x, sendCount_x, rank_x, sendtag + 0);
    req2[0] =
        MPI_COMM_SCALBL.Irecv(recvbuf_X, recvCount_X, rank_X, recvtag + 0);

    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, 0, sendCount_X, sendbuf_X,
                      &Aq[Component * 7 * N], N);
    req1[1] =
        MPI_COMM_SCALBL.Isend(sendbuf_X, sendCount_X, rank_X, sendtag + 1);
    req2[1] =
        MPI_COMM_SCALBL.Irecv(recvbuf_x, recvCount_x, rank_x, recvtag + 1);

    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, 0, sendCount_y, sendbuf_y,
                      &Aq[Component * 7 * N], N);
    req1[2] =
        MPI_COMM_SCALBL.Isend(sendbuf_y, sendCount_y, rank_y, sendtag + 2);
    req2[2] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Y, recvCount_Y, rank_Y, recvtag + 2);

    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, 0, sendCount_Y, sendbuf_Y,
                      &Aq[Component * 7 * N], N);
    req1[3] =
        MPI_COMM_SCALBL.Isend(sendbuf_Y, sendCount_Y, rank_Y, sendtag + 3);
    req2[3] =
        MPI_COMM_SCALBL.Irecv(recvbuf_y, recvCount_y, rank_y, recvtag + 3);

    //...Packing for z face(6,12,13,16,17)................................
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, 0, sendCount_z, sendbuf_z,
                      &Aq[Component * 7 * N], N);
    req1[4] =
        MPI_COMM_SCALBL.Isend(sendbuf_z, sendCount_z, rank_z, sendtag + 4);
    req2[4] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Z, recvCount_Z, rank_Z, recvtag + 4);

    //...Packing for Z face(5,11,14,15,18)................................
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, 0, sendCount_Z, sendbuf_Z,
                      &Aq[Component * 7 * N], N);
    req1[5] =
        MPI_COMM_SCALBL.Isend(sendbuf_Z, sendCount_Z, rank_Z, sendtag + 5);
    req2[5] =
        MPI_COMM_SCALBL.Irecv(recvbuf_z, recvCount_z, rank_z, recvtag + 5);
}

void ScaLBL_Communicator::RecvD3Q7AA(double *Aq, int Component) {

    // NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
    //...................................................................................
    // Wait for completion of D3Q19 communication
    MPI_COMM_SCALBL.waitAll(6, req1);
    MPI_COMM_SCALBL.waitAll(6, req2);
    ScaLBL_DeviceBarrier();

    //...................................................................................
    // NOTE: AA Routine writes to opposite
    // Unpack the distributions on the device
    //...................................................................................
    //...Unpacking for x face(2,8,10,12,14)................................
    ScaLBL_D3Q7_Unpack(2, dvcRecvDist_x, 0, recvCount_x, recvbuf_x,
                       &Aq[Component * 7 * N], N);
    //...................................................................................
    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q7_Unpack(1, dvcRecvDist_X, 0, recvCount_X, recvbuf_X,
                       &Aq[Component * 7 * N], N);
    //...................................................................................
    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q7_Unpack(4, dvcRecvDist_y, 0, recvCount_y, recvbuf_y,
                       &Aq[Component * 7 * N], N);
    //...................................................................................
    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q7_Unpack(3, dvcRecvDist_Y, 0, recvCount_Y, recvbuf_Y,
                       &Aq[Component * 7 * N], N);
    //...................................................................................

    if (BoundaryCondition > 0) {
        if (kproc != 0) {
            //...Packing for z face(6,12,13,16,17)................................
            ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 0, recvCount_z, recvbuf_z,
                               &Aq[Component * 7 * N], N);
        }
        if (kproc != nprocz - 1) {
            //...Packing for Z face(5,11,14,15,18)................................
            ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 0, recvCount_Z, recvbuf_Z,
                               &Aq[Component * 7 * N], N);
        }
    } else {
        //...Packing for z face(6,12,13,16,17)................................
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 0, recvCount_z, recvbuf_z,
                           &Aq[Component * 7 * N], N);
        //...Packing for Z face(5,11,14,15,18)................................
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 0, recvCount_Z, recvbuf_Z,
                           &Aq[Component * 7 * N], N);
    }

    //...................................................................................
    Lock = false; // unlock the communicator after communications complete
    //...................................................................................
}

void ScaLBL_Communicator::TriSendD3Q7AA(double *Aq, double *Bq, double *Cq) {

    // NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
    if (Lock == true) {
        ERROR("ScaLBL Error (SendD3Q19): ScaLBL_Communicator is locked -- did "
              "you forget to match Send/Recv calls?");
    } else {
        Lock = true;
    }
    // assign tag of 19 to D3Q19 communication
    sendtag = recvtag = 162;
    ScaLBL_DeviceBarrier();
    // Pack the distributions
    //...Packing for x face(2,8,10,12,14)................................
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, 0, sendCount_x, sendbuf_x, Aq, N);
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, sendCount_x, sendCount_x, sendbuf_x, Bq,
                      N);
    ScaLBL_D3Q19_Pack(2, dvcSendList_x, 2 * sendCount_x, sendCount_x, sendbuf_x,
                      Cq, N);
    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, 0, sendCount_X, sendbuf_X, Aq, N);
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, sendCount_X, sendCount_X, sendbuf_X, Bq,
                      N);
    ScaLBL_D3Q19_Pack(1, dvcSendList_X, 2 * sendCount_X, sendCount_X, sendbuf_X,
                      Cq, N);
    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, 0, sendCount_y, sendbuf_y, Aq, N);
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, sendCount_y, sendCount_y, sendbuf_y, Bq,
                      N);
    ScaLBL_D3Q19_Pack(4, dvcSendList_y, 2 * sendCount_y, sendCount_y, sendbuf_y,
                      Cq, N);
    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, 0, sendCount_Y, sendbuf_Y, Aq, N);
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, sendCount_Y, sendCount_Y, sendbuf_Y, Bq,
                      N);
    ScaLBL_D3Q19_Pack(3, dvcSendList_Y, 2 * sendCount_Y, sendCount_Y, sendbuf_Y,
                      Cq, N);
    //...Packing for z face(6,12,13,16,17)................................
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, 0, sendCount_z, sendbuf_z, Aq, N);
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, sendCount_z, sendCount_z, sendbuf_z, Bq,
                      N);
    ScaLBL_D3Q19_Pack(6, dvcSendList_z, 2 * sendCount_z, sendCount_z, sendbuf_z,
                      Cq, N);
    //...Packing for Z face(5,11,14,15,18)................................
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, 0, sendCount_Z, sendbuf_Z, Aq, N);
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, sendCount_Z, sendCount_Z, sendbuf_Z, Bq,
                      N);
    ScaLBL_D3Q19_Pack(5, dvcSendList_Z, 2 * sendCount_Z, sendCount_Z, sendbuf_Z,
                      Cq, N);

    //...................................................................................
    // Send all the distributions
    req1[0] =
        MPI_COMM_SCALBL.Isend(sendbuf_x, 3 * sendCount_x, rank_x, sendtag + 0);
    req2[0] =
        MPI_COMM_SCALBL.Irecv(recvbuf_X, 3 * recvCount_X, rank_X, recvtag + 0);
    req1[1] =
        MPI_COMM_SCALBL.Isend(sendbuf_X, 3 * sendCount_X, rank_X, sendtag + 1);
    req2[1] =
        MPI_COMM_SCALBL.Irecv(recvbuf_x, 3 * recvCount_x, rank_x, recvtag + 1);
    req1[2] =
        MPI_COMM_SCALBL.Isend(sendbuf_y, 3 * sendCount_y, rank_y, sendtag + 2);
    req2[2] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Y, 3 * recvCount_Y, rank_Y, recvtag + 2);
    req1[3] =
        MPI_COMM_SCALBL.Isend(sendbuf_Y, 3 * sendCount_Y, rank_Y, sendtag + 3);
    req2[3] =
        MPI_COMM_SCALBL.Irecv(recvbuf_y, 3 * recvCount_y, rank_y, recvtag + 3);
    req1[4] =
        MPI_COMM_SCALBL.Isend(sendbuf_z, 3 * sendCount_z, rank_z, sendtag + 4);
    req2[4] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Z, 3 * recvCount_Z, rank_Z, recvtag + 4);
    req1[5] =
        MPI_COMM_SCALBL.Isend(sendbuf_Z, 3 * sendCount_Z, rank_Z, sendtag + 5);
    req2[5] =
        MPI_COMM_SCALBL.Irecv(recvbuf_z, 3 * recvCount_z, rank_z, recvtag + 5);
}

void ScaLBL_Communicator::TriRecvD3Q7AA(double *Aq, double *Bq, double *Cq) {

    // NOTE: the center distribution f0 must NOT be at the start of feven, provide offset to start of f2
    //...................................................................................
    // Wait for completion of D3Q19 communication
    MPI_COMM_SCALBL.waitAll(6, req1);
    MPI_COMM_SCALBL.waitAll(6, req2);
    ScaLBL_DeviceBarrier();

    //...................................................................................
    // NOTE: AA Routine writes to opposite
    // Unpack the distributions on the device
    //...................................................................................
    //...Unpacking for x face(2,8,10,12,14)................................
    ScaLBL_D3Q7_Unpack(2, dvcRecvDist_x, 0, recvCount_x, recvbuf_x, Aq, N);
    ScaLBL_D3Q7_Unpack(2, dvcRecvDist_x, recvCount_x, recvCount_x, recvbuf_x,
                       Bq, N);
    ScaLBL_D3Q7_Unpack(2, dvcRecvDist_x, 2 * recvCount_x, recvCount_x,
                       recvbuf_x, Cq, N);
    //...................................................................................
    //...Packing for X face(1,7,9,11,13)................................
    ScaLBL_D3Q7_Unpack(1, dvcRecvDist_X, 0, recvCount_X, recvbuf_X, Aq, N);
    ScaLBL_D3Q7_Unpack(1, dvcRecvDist_X, recvCount_X, recvCount_X, recvbuf_X,
                       Bq, N);
    ScaLBL_D3Q7_Unpack(1, dvcRecvDist_X, 2 * recvCount_X, recvCount_X,
                       recvbuf_X, Cq, N);
    //...................................................................................
    //...Packing for y face(4,8,9,16,18).................................
    ScaLBL_D3Q7_Unpack(4, dvcRecvDist_y, 0, recvCount_y, recvbuf_y, Aq, N);
    ScaLBL_D3Q7_Unpack(4, dvcRecvDist_y, recvCount_y, recvCount_y, recvbuf_y,
                       Bq, N);
    ScaLBL_D3Q7_Unpack(4, dvcRecvDist_y, 2 * recvCount_y, recvCount_y,
                       recvbuf_y, Cq, N);
    //...................................................................................
    //...Packing for Y face(3,7,10,15,17).................................
    ScaLBL_D3Q7_Unpack(3, dvcRecvDist_Y, 0, recvCount_Y, recvbuf_Y, Aq, N);
    ScaLBL_D3Q7_Unpack(3, dvcRecvDist_Y, recvCount_Y, recvCount_Y, recvbuf_Y,
                       Bq, N);
    ScaLBL_D3Q7_Unpack(3, dvcRecvDist_Y, 2 * recvCount_Y, recvCount_Y,
                       recvbuf_Y, Cq, N);
    //...................................................................................

    if (BoundaryCondition > 0 && kproc == 0) {
        // don't unpack little z
        //...Packing for Z face(5,11,14,15,18)................................
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 0, recvCount_Z, recvbuf_Z, Aq, N);
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, recvCount_Z, recvCount_Z,
                           recvbuf_Z, Bq, N);
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 2 * recvCount_Z, recvCount_Z,
                           recvbuf_Z, Cq, N);
    } else if (BoundaryCondition > 0 && kproc == nprocz - 1) {
        // don't unpack big z
        //...Packing for z face(6,12,13,16,17)................................
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 0, recvCount_z, recvbuf_z, Aq, N);
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, recvCount_z, recvCount_z,
                           recvbuf_z, Bq, N);
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 2 * recvCount_z, recvCount_z,
                           recvbuf_z, Cq, N);
    } else {
        //...Packing for z face(6,12,13,16,17)................................
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 0, recvCount_z, recvbuf_z, Aq, N);
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, recvCount_z, recvCount_z,
                           recvbuf_z, Bq, N);
        ScaLBL_D3Q7_Unpack(6, dvcRecvDist_z, 2 * recvCount_z, recvCount_z,
                           recvbuf_z, Cq, N);
        //...Packing for Z face(5,11,14,15,18)................................
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 0, recvCount_Z, recvbuf_Z, Aq, N);
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, recvCount_Z, recvCount_Z,
                           recvbuf_Z, Bq, N);
        ScaLBL_D3Q7_Unpack(5, dvcRecvDist_Z, 2 * recvCount_Z, recvCount_Z,
                           recvbuf_Z, Cq, N);
    }

    //...................................................................................
    Lock = false; // unlock the communicator after communications complete
    //...................................................................................
}

void ScaLBL_Communicator::SendHalo(double *data) {
    //...................................................................................
    if (Lock == true) {
        ERROR("ScaLBL Error (SendHalo): ScaLBL_Communicator is locked -- did "
              "you forget to match Send/Recv calls?");
    } else {
        Lock = true;
    }
    ScaLBL_DeviceBarrier();
    //...................................................................................
    sendtag = recvtag = 168;
    //...................................................................................
    ScaLBL_Scalar_Pack(dvcSendList_x, sendCount_x, sendbuf_x, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_y, sendCount_y, sendbuf_y, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_z, sendCount_z, sendbuf_z, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_X, sendCount_X, sendbuf_X, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_Y, sendCount_Y, sendbuf_Y, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_Z, sendCount_Z, sendbuf_Z, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_xy, sendCount_xy, sendbuf_xy, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_xY, sendCount_xY, sendbuf_xY, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_Xy, sendCount_Xy, sendbuf_Xy, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_XY, sendCount_XY, sendbuf_XY, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_xz, sendCount_xz, sendbuf_xz, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_xZ, sendCount_xZ, sendbuf_xZ, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_Xz, sendCount_Xz, sendbuf_Xz, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_XZ, sendCount_XZ, sendbuf_XZ, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_yz, sendCount_yz, sendbuf_yz, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_yZ, sendCount_yZ, sendbuf_yZ, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_Yz, sendCount_Yz, sendbuf_Yz, data, N);
    ScaLBL_Scalar_Pack(dvcSendList_YZ, sendCount_YZ, sendbuf_YZ, data, N);
    //...................................................................................
    // Send / Recv all the phase indcator field values
    //...................................................................................

    req1[0] =
        MPI_COMM_SCALBL.Isend(sendbuf_x, sendCount_x, rank_x, sendtag + 0);
    req2[0] =
        MPI_COMM_SCALBL.Irecv(recvbuf_X, recvCount_X, rank_X, recvtag + 0);
    req1[1] =
        MPI_COMM_SCALBL.Isend(sendbuf_X, sendCount_X, rank_X, sendtag + 1);
    req2[1] =
        MPI_COMM_SCALBL.Irecv(recvbuf_x, recvCount_x, rank_x, recvtag + 1);
    req1[2] =
        MPI_COMM_SCALBL.Isend(sendbuf_y, sendCount_y, rank_y, sendtag + 2);
    req2[2] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Y, recvCount_Y, rank_Y, recvtag + 2);
    req1[3] =
        MPI_COMM_SCALBL.Isend(sendbuf_Y, sendCount_Y, rank_Y, sendtag + 3);
    req2[3] =
        MPI_COMM_SCALBL.Irecv(recvbuf_y, recvCount_y, rank_y, recvtag + 3);
    req1[4] =
        MPI_COMM_SCALBL.Isend(sendbuf_z, sendCount_z, rank_z, sendtag + 4);
    req2[4] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Z, recvCount_Z, rank_Z, recvtag + 4);
    req1[5] =
        MPI_COMM_SCALBL.Isend(sendbuf_Z, sendCount_Z, rank_Z, sendtag + 5);
    req2[5] =
        MPI_COMM_SCALBL.Irecv(recvbuf_z, recvCount_z, rank_z, recvtag + 5);
    req1[6] =
        MPI_COMM_SCALBL.Isend(sendbuf_xy, sendCount_xy, rank_xy, sendtag + 6);
    req2[6] =
        MPI_COMM_SCALBL.Irecv(recvbuf_XY, recvCount_XY, rank_XY, recvtag + 6);
    req1[7] =
        MPI_COMM_SCALBL.Isend(sendbuf_XY, sendCount_XY, rank_XY, sendtag + 7);
    req2[7] =
        MPI_COMM_SCALBL.Irecv(recvbuf_xy, recvCount_xy, rank_xy, recvtag + 7);
    req1[8] =
        MPI_COMM_SCALBL.Isend(sendbuf_Xy, sendCount_Xy, rank_Xy, sendtag + 8);
    req2[8] =
        MPI_COMM_SCALBL.Irecv(recvbuf_xY, recvCount_xY, rank_xY, recvtag + 8);
    req1[9] =
        MPI_COMM_SCALBL.Isend(sendbuf_xY, sendCount_xY, rank_xY, sendtag + 9);
    req2[9] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Xy, recvCount_Xy, rank_Xy, recvtag + 9);
    req1[10] =
        MPI_COMM_SCALBL.Isend(sendbuf_xz, sendCount_xz, rank_xz, sendtag + 10);
    req2[10] =
        MPI_COMM_SCALBL.Irecv(recvbuf_XZ, recvCount_XZ, rank_XZ, recvtag + 10);
    req1[11] =
        MPI_COMM_SCALBL.Isend(sendbuf_XZ, sendCount_XZ, rank_XZ, sendtag + 11);
    req2[11] =
        MPI_COMM_SCALBL.Irecv(recvbuf_xz, recvCount_xz, rank_xz, recvtag + 11);
    req1[12] =
        MPI_COMM_SCALBL.Isend(sendbuf_Xz, sendCount_Xz, rank_Xz, sendtag + 12);
    req2[12] =
        MPI_COMM_SCALBL.Irecv(recvbuf_xZ, recvCount_xZ, rank_xZ, recvtag + 12);
    req1[13] =
        MPI_COMM_SCALBL.Isend(sendbuf_xZ, sendCount_xZ, rank_xZ, sendtag + 13);
    req2[13] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Xz, recvCount_Xz, rank_Xz, recvtag + 13);
    req1[14] =
        MPI_COMM_SCALBL.Isend(sendbuf_yz, sendCount_yz, rank_yz, sendtag + 14);
    req2[14] =
        MPI_COMM_SCALBL.Irecv(recvbuf_YZ, recvCount_YZ, rank_YZ, recvtag + 14);
    req1[15] =
        MPI_COMM_SCALBL.Isend(sendbuf_YZ, sendCount_YZ, rank_YZ, sendtag + 15);
    req2[15] =
        MPI_COMM_SCALBL.Irecv(recvbuf_yz, recvCount_yz, rank_yz, recvtag + 15);
    req1[16] =
        MPI_COMM_SCALBL.Isend(sendbuf_Yz, sendCount_Yz, rank_Yz, sendtag + 16);
    req2[16] =
        MPI_COMM_SCALBL.Irecv(recvbuf_yZ, recvCount_yZ, rank_yZ, recvtag + 16);
    req1[17] =
        MPI_COMM_SCALBL.Isend(sendbuf_yZ, sendCount_yZ, rank_yZ, sendtag + 17);
    req2[17] =
        MPI_COMM_SCALBL.Irecv(recvbuf_Yz, recvCount_Yz, rank_Yz, recvtag + 17);
    //...................................................................................
}
void ScaLBL_Communicator::RecvHalo(double *data) {

    //...................................................................................
    MPI_COMM_SCALBL.waitAll(18, req1);
    MPI_COMM_SCALBL.waitAll(18, req2);
    ScaLBL_DeviceBarrier();
    //...................................................................................
    //...................................................................................
    ScaLBL_Scalar_Unpack(dvcRecvList_x, recvCount_x, recvbuf_x, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_y, recvCount_y, recvbuf_y, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_z, recvCount_z, recvbuf_z, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_X, recvCount_X, recvbuf_X, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_Y, recvCount_Y, recvbuf_Y, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_Z, recvCount_Z, recvbuf_Z, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_xy, recvCount_xy, recvbuf_xy, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_xY, recvCount_xY, recvbuf_xY, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_Xy, recvCount_Xy, recvbuf_Xy, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_XY, recvCount_XY, recvbuf_XY, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_xz, recvCount_xz, recvbuf_xz, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_xZ, recvCount_xZ, recvbuf_xZ, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_Xz, recvCount_Xz, recvbuf_Xz, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_XZ, recvCount_XZ, recvbuf_XZ, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_yz, recvCount_yz, recvbuf_yz, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_yZ, recvCount_yZ, recvbuf_yZ, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_Yz, recvCount_Yz, recvbuf_Yz, data, N);
    ScaLBL_Scalar_Unpack(dvcRecvList_YZ, recvCount_YZ, recvbuf_YZ, data, N);
    //...................................................................................
    Lock = false; // unlock the communicator after communications complete
    //...................................................................................
}

void ScaLBL_Communicator::RegularLayout(IntArray map, const double *data,
                                        DoubleArray &regdata) {
    // Gets data from the device and stores in regular layout
    int i, j, k, idx;
    int Nx = map.size(0);
    int Ny = map.size(1);
    int Nz = map.size(2);

    // initialize the array
    regdata.fill(0.f);

    double *TmpDat;
    double value;
    TmpDat = new double[N];
    ScaLBL_CopyToHost(&TmpDat[0], &data[0], N * sizeof(double));
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                idx = map(i, j, k);
                if (!(idx < 0)) {
                    value = TmpDat[idx];
                    regdata(i, j, k) = value;
                }
            }
        }
    }
    //printf("r=%i, value=%f   ",rank,value);

    delete[] TmpDat;
}

void ScaLBL_Communicator::Color_BC_z(int *Map, double *Phi, double *Den,
                                     double vA, double vB) {
    if (kproc == 0) {
        if (BoundaryCondition == 5) {
            //ScaLBL_CopySlice_z(Phi,Nx,Ny,Nz,1,0);
        } else {
            // Set the phase indicator field and density on the z inlet
            ScaLBL_Color_BC_z(dvcSendList_z, Map, Phi, Den, vA, vB, sendCount_z,
                              N);
        }
        //ScaLBL_SetSlice_z(Phi,Value,Nx,Ny,Nz,0);
    }
}

void ScaLBL_Communicator::Color_BC_Z(int *Map, double *Phi, double *Den,
                                     double vA, double vB) {
    if (kproc == nprocz - 1) {
        if (BoundaryCondition == 5) {
            //ScaLBL_CopySlice_z(Phi,Nx,Ny,Nz,Nz-2,Nz-1);
        } else {
            // Set the phase indicator field and density on the Z outlet
            ScaLBL_Color_BC_Z(dvcSendList_Z, Map, Phi, Den, vA, vB, sendCount_Z,
                              N);
        }
    }
}

void ScaLBL_Communicator::D3Q19_Pressure_BC_z(int *neighborList, double *fq,
                                              double din, int time) {
    //ScaLBL_D3Q19_Pressure_BC_z(int *LIST,fq,din,Nx,Ny,Nz);
    if (kproc == 0) {
        if (time % 2 == 0) {
            ScaLBL_D3Q19_AAeven_Pressure_BC_z(dvcSendList_z, fq, din,
                                              sendCount_z, N);
        } else {
            ScaLBL_D3Q19_AAodd_Pressure_BC_z(neighborList, dvcSendList_z, fq,
                                             din, sendCount_z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q19_Pressure_BC_Z(int *neighborList, double *fq,
                                              double dout, int time) {
    //ScaLBL_D3Q19_Pressure_BC_Z(int *LIST,fq,dout,Nx,Ny,Nz);
    if (kproc == nprocz - 1) {
        if (time % 2 == 0) {
            ScaLBL_D3Q19_AAeven_Pressure_BC_Z(dvcSendList_Z, fq, dout,
                                              sendCount_Z, N);
        } else {
            ScaLBL_D3Q19_AAodd_Pressure_BC_Z(neighborList, dvcSendList_Z, fq,
                                             dout, sendCount_Z, N);
        }
    }
}

double ScaLBL_Communicator::D3Q19_Flux_BC_z(int *neighborList, double *fq,
                                            double flux, int time) {
    double sum, locsum, din;
    double LocInletArea, InletArea;

    // Note that flux = rho_0 * Q

    // Compute the inlet area
    if (kproc == 0)
        LocInletArea = double(sendCount_z);
    else
        LocInletArea = 0.f;

    InletArea = MPI_COMM_SCALBL.sumReduce(LocInletArea);
    //printf("Inlet area = %f \n", InletArea);

    // Set the flux BC
    locsum = 0.f;
    if (time % 2 == 0) {
        if (kproc == 0)
            locsum = ScaLBL_D3Q19_AAeven_Flux_BC_z(dvcSendList_z, fq, flux,
                                                   InletArea, sendCount_z, N);

        sum = MPI_COMM_SCALBL.sumReduce(locsum);

        din = flux / InletArea + sum;
        //if (rank==0) printf("computed din (even) =%f \n",din);
        if (kproc == 0)
            ScaLBL_D3Q19_AAeven_Pressure_BC_z(dvcSendList_z, fq, din,
                                              sendCount_z, N);
    } else {
        if (kproc == 0)
            locsum =
                ScaLBL_D3Q19_AAodd_Flux_BC_z(neighborList, dvcSendList_z, fq,
                                             flux, InletArea, sendCount_z, N);

        sum = MPI_COMM_SCALBL.sumReduce(locsum);
        din = flux / InletArea + sum;

        //if (rank==0) printf("computed din (odd)=%f \n",din);
        if (kproc == 0)
            ScaLBL_D3Q19_AAodd_Pressure_BC_z(neighborList, dvcSendList_z, fq,
                                             din, sendCount_z, N);
    }
    //printf("Inlet pressure = %f \n", din);
    return din;
}

void ScaLBL_Communicator::D3Q19_Reflection_BC_z(double *fq) {
    if (kproc == 0)
        ScaLBL_D3Q19_Reflection_BC_z(dvcSendList_z, fq, sendCount_z, N);
}

void ScaLBL_Communicator::D3Q19_Reflection_BC_Z(double *fq) {
    if (kproc == nprocz - 1)
        ScaLBL_D3Q19_Reflection_BC_Z(dvcSendList_Z, fq, sendCount_Z, N);
}

void ScaLBL_Communicator::PrintD3Q19() {
    printf("Printing D3Q19 communication buffer contents \n");

    int i, n;
    double f;
    int *TempBuffer;
    TempBuffer = new int[5 * recvCount_x];

    //.......................................................................
    // Re-index the send lists
    ScaLBL_CopyToHost(TempBuffer, dvcRecvDist_x, 5 * recvCount_x * sizeof(int));
    for (i = 0; i < 5 * recvCount_x; i++) {
        n = TempBuffer[i];
        f = recvbuf_x[i];
        printf("Receive %f to %i from buffer index %i \n", f, n, i);
    }

    delete[] TempBuffer;
}

void ScaLBL_Communicator::D3Q7_Poisson_Potential_BC_z(int *neighborList,
                                                      double *fq, double Vin,
                                                      int time) {
    if (kproc == 0) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_z(dvcSendList_z, fq, Vin,
                                                      sendCount_z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_z(
                neighborList, dvcSendList_z, fq, Vin, sendCount_z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Poisson_Potential_BC_Z(int *neighborList,
                                                      double *fq, double Vout,
                                                      int time) {
    if (kproc == nprocz - 1) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Poisson_Potential_BC_Z(dvcSendList_Z, fq, Vout,
                                                      sendCount_Z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Poisson_Potential_BC_Z(
                neighborList, dvcSendList_Z, fq, Vout, sendCount_Z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q19_Poisson_Potential_BC_z(int *neighborList,
                                                       double *fq, double Vin,
                                                       int time) {
    if (kproc == 0) {
        if (time % 2 == 0) {
            ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_z(dvcSendList_z, fq, Vin,
                                                       sendCount_z, N);
        } else {
            ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_z(
                neighborList, dvcSendList_z, fq, Vin, sendCount_z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q19_Poisson_Potential_BC_Z(int *neighborList,
                                                       double *fq, double Vout,
                                                       int time) {
    if (kproc == nprocz - 1) {
        if (time % 2 == 0) {
            ScaLBL_D3Q19_AAeven_Poisson_Potential_BC_Z(dvcSendList_Z, fq, Vout,
                                                       sendCount_Z, N);
        } else {
            ScaLBL_D3Q19_AAodd_Poisson_Potential_BC_Z(
                neighborList, dvcSendList_Z, fq, Vout, sendCount_Z, N);
        }
    }
}

void ScaLBL_Communicator::Poisson_D3Q7_BC_z(int *Map, double *Psi, double Vin) {
    if (kproc == 0) {
        ScaLBL_Poisson_D3Q7_BC_z(dvcSendList_z, Map, Psi, Vin, sendCount_z);
    }
}

void ScaLBL_Communicator::Poisson_D3Q7_BC_Z(int *Map, double *Psi,
                                            double Vout) {
    if (kproc == nprocz - 1) {
        ScaLBL_Poisson_D3Q7_BC_Z(dvcSendList_Z, Map, Psi, Vout, sendCount_Z);
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Concentration_BC_z(int *neighborList,
                                                      double *fq, double Cin,
                                                      int time) {
    if (kproc == 0) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_z(dvcSendList_z, fq, Cin,
                                                      sendCount_z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_z(
                neighborList, dvcSendList_z, fq, Cin, sendCount_z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Concentration_BC_Z(int *neighborList,
                                                      double *fq, double Cout,
                                                      int time) {
    if (kproc == nprocz - 1) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Concentration_BC_Z(dvcSendList_Z, fq, Cout,
                                                      sendCount_Z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Concentration_BC_Z(
                neighborList, dvcSendList_Z, fq, Cout, sendCount_Z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Flux_Diff_BC_z(int *neighborList, double *fq,
                                                  double Cin, double tau,
                                                  double *VelocityZ, int time) {
    if (kproc == 0) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_z(dvcSendList_z, fq, Cin, tau,
                                                  VelocityZ, sendCount_z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_z(neighborList, dvcSendList_z,
                                                 fq, Cin, tau, VelocityZ,
                                                 sendCount_z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Flux_Diff_BC_Z(int *neighborList, double *fq,
                                                  double Cout, double tau,
                                                  double *VelocityZ, int time) {
    if (kproc == nprocz - 1) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Flux_Diff_BC_Z(dvcSendList_Z, fq, Cout, tau,
                                                  VelocityZ, sendCount_Z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Flux_Diff_BC_Z(neighborList, dvcSendList_Z,
                                                 fq, Cout, tau, VelocityZ,
                                                 sendCount_Z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Flux_DiffAdvc_BC_z(int *neighborList,
                                                      double *fq, double Cin,
                                                      double tau,
                                                      double *VelocityZ,
                                                      int time) {
    if (kproc == 0) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_z(
                dvcSendList_z, fq, Cin, tau, VelocityZ, sendCount_z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_z(
                neighborList, dvcSendList_z, fq, Cin, tau, VelocityZ,
                sendCount_z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Flux_DiffAdvc_BC_Z(int *neighborList,
                                                      double *fq, double Cout,
                                                      double tau,
                                                      double *VelocityZ,
                                                      int time) {
    if (kproc == nprocz - 1) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvc_BC_Z(
                dvcSendList_Z, fq, Cout, tau, VelocityZ, sendCount_Z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvc_BC_Z(
                neighborList, dvcSendList_Z, fq, Cout, tau, VelocityZ,
                sendCount_Z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Flux_DiffAdvcElec_BC_z(
    int *neighborList, double *fq, double Cin, double tau, double *VelocityZ,
    double *ElectricField_Z, double Di, double zi, double Vt, int time) {
    if (kproc == 0) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_z(
                dvcSendList_z, fq, Cin, tau, VelocityZ, ElectricField_Z, Di, zi,
                Vt, sendCount_z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_z(
                neighborList, dvcSendList_z, fq, Cin, tau, VelocityZ,
                ElectricField_Z, Di, zi, Vt, sendCount_z, N);
        }
    }
}

void ScaLBL_Communicator::D3Q7_Ion_Flux_DiffAdvcElec_BC_Z(
    int *neighborList, double *fq, double Cout, double tau, double *VelocityZ,
    double *ElectricField_Z, double Di, double zi, double Vt, int time) {
    if (kproc == nprocz - 1) {
        if (time % 2 == 0) {
            ScaLBL_D3Q7_AAeven_Ion_Flux_DiffAdvcElec_BC_Z(
                dvcSendList_Z, fq, Cout, tau, VelocityZ, ElectricField_Z, Di,
                zi, Vt, sendCount_Z, N);
        } else {
            ScaLBL_D3Q7_AAodd_Ion_Flux_DiffAdvcElec_BC_Z(
                neighborList, dvcSendList_Z, fq, Cout, tau, VelocityZ,
                ElectricField_Z, Di, zi, Vt, sendCount_Z, N);
        }
    }
}
