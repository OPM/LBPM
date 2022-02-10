#include <array>
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "common/Database.h"
#include "common/Utilities.h"
#include "common/MPI.h"


inline MPI_Request Isend( const Utilities::MPI &comm, const int *buf, int count, int rank, int tag )
{
    return comm.Isend( buf, count, rank, tag );
}
inline MPI_Request Irecv( const Utilities::MPI &comm, int *buf, int count, int rank, int tag )
{
    return comm.Irecv( buf, count, rank, tag );
}


struct RankInfoStruct2 {
    int nx;            //!<  The number of processors in the x direction
    int ny;            //!<  The number of processors in the y direction
    int nz;            //!<  The number of processors in the z direction
    int ix;            //!<  The index of the current process in the x direction
    int jy;            //!<  The index of the current process in the y direction
    int kz;            //!<  The index of the current process in the z direction
    int rank[3][3][3]; //!<  The rank for the neighbor [i][j][k]
    RankInfoStruct2() : RankInfoStruct2( 1, 0, 0, 0 ) {}
    RankInfoStruct2(int rank0, int nprocx, int nprocy, int nprocz) {
       memset(this, 0, sizeof(RankInfoStruct2));
       nx = nprocx;
       ny = nprocy;
       nz = nprocz;
       if (rank0 >= nprocx * nprocy * nprocz) {
           ix = -1;
           jy = -1;
           kz = -1;
           for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    for (int k = -1; k <= 1; k++) {
                       rank[i + 1][j + 1][k + 1] = -1;
                    }
                }
            }
        } else {
            ix = rank0 % nprocx;
            jy = (rank0 / nprocx) % nprocy;
            kz = rank0 / (nprocx * nprocy);
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    for (int k = -1; k <= 1; k++) {
                        rank[i + 1][j + 1][k + 1] =
                            getRankForBlock(ix + i, jy + j, kz + k);
                    }
                }
            }
            ASSERT(rank[1][1][1] == rank0);
        }
    }
    int getRankForBlock(int i, int j, int k) const {
        int i2 = (i + nx) % nx;
        int j2 = (j + ny) % ny;
        int k2 = (k + nz) % nz;
        return i2 + j2 * nx + k2 * nx * ny;
    }
};


class Domain2 {
public:
    Domain2(std::shared_ptr<Database> db, const Utilities::MPI &Communicator) {
        Comm = Communicator.dup();
        auto nproc = db->getVector<int>("nproc");
        auto n = db->getVector<int>("n");
        ASSERT(n.size() == 3u);
        ASSERT(nproc.size() == 3u);
        int nx = n[0];
        int ny = n[1];
        int nz = n[2];
        Nx = nx + 2;
        Ny = ny + 2;
        Nz = nz + 2;
        N = Nx * Ny * Nz;
        int myrank = Comm.getRank();
        rank_info = RankInfoStruct2(myrank, nproc[0], nproc[1], nproc[2]);
        int nprocs = Comm.getSize();
        INSIST(nprocs == nproc[0] * nproc[1] * nproc[2], "Fatal error in processor count!");
    }

    Domain2() = delete;
    Domain2(const Domain2 &) = delete;
    ~Domain2() = default;

public: // Public variables (need to create accessors instead)
    int Nx, Ny, Nz, N;
    RankInfoStruct2 rank_info;

    Utilities::MPI Comm; // MPI Communicator for this domain


    //**********************************
    // MPI ranks for all 18 neighbors
    //**********************************
    inline int iproc() const { return rank_info.ix; }
    inline int jproc() const { return rank_info.jy; }
    inline int kproc() const { return rank_info.kz; }
    inline int nprocx() const { return rank_info.nx; }
    inline int nprocy() const { return rank_info.ny; }
    inline int nprocz() const { return rank_info.nz; }
    inline int rank() const { return rank_info.rank[1][1][1]; }
    inline int rank_X() const { return rank_info.rank[2][1][1]; }
    inline int rank_x() const { return rank_info.rank[0][1][1]; }
    inline int rank_Y() const { return rank_info.rank[1][2][1]; }
    inline int rank_y() const { return rank_info.rank[1][0][1]; }
    inline int rank_Z() const { return rank_info.rank[1][1][2]; }
    inline int rank_z() const { return rank_info.rank[1][1][0]; }
    inline int rank_XY() const { return rank_info.rank[2][2][1]; }
    inline int rank_xy() const { return rank_info.rank[0][0][1]; }
    inline int rank_Xy() const { return rank_info.rank[2][0][1]; }
    inline int rank_xY() const { return rank_info.rank[0][2][1]; }
    inline int rank_XZ() const { return rank_info.rank[2][1][2]; }
    inline int rank_xz() const { return rank_info.rank[0][1][0]; }
    inline int rank_Xz() const { return rank_info.rank[2][1][0]; }
    inline int rank_xZ() const { return rank_info.rank[0][1][2]; }
    inline int rank_YZ() const { return rank_info.rank[1][2][2]; }
    inline int rank_yz() const { return rank_info.rank[1][0][0]; }
    inline int rank_Yz() const { return rank_info.rank[1][2][0]; }
    inline int rank_yZ() const { return rank_info.rank[1][0][2]; }

    // Initialize communication data structures within Domain object. 
    void CommInit() {
        int i, j, k, n;
        int sendtag = 21;
        int recvtag = 21;
        //......................................................................................
        int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y, sendCount_Z;
        int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz, sendCount_xZ;
        int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ, sendCount_XZ;
        sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
        sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
        sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
        //......................................................................................
        for (k = 1; k < Nz - 1; k++) {
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    // Counts for the six faces
                    if (i == 1)
                        sendCount_x++;
                    if (j == 1)
                        sendCount_y++;
                    if (k == 1)
                        sendCount_z++;
                    if (i == Nx - 2)
                        sendCount_X++;
                    if (j == Ny - 2)
                        sendCount_Y++;
                    if (k == Nz - 2)
                        sendCount_Z++;
                    // Counts for the twelve edges
                    if (i == 1 && j == 1)
                        sendCount_xy++;
                    if (i == 1 && j == Ny - 2)
                        sendCount_xY++;
                    if (i == Nx - 2 && j == 1)
                        sendCount_Xy++;
                    if (i == Nx - 2 && j == Ny - 2)
                        sendCount_XY++;

                    if (i == 1 && k == 1)
                        sendCount_xz++;
                    if (i == 1 && k == Nz - 2)
                        sendCount_xZ++;
                    if (i == Nx - 2 && k == 1)
                        sendCount_Xz++;
                    if (i == Nx - 2 && k == Nz - 2)
                        sendCount_XZ++;

                    if (j == 1 && k == 1)
                        sendCount_yz++;
                    if (j == 1 && k == Nz - 2)
                        sendCount_yZ++;
                    if (j == Ny - 2 && k == 1)
                        sendCount_Yz++;
                    if (j == Ny - 2 && k == Nz - 2)
                        sendCount_YZ++;
                }
            }
        }

        // allocate send lists
        sendList_x.resize(sendCount_x, 0);
        sendList_y.resize(sendCount_y, 0);
        sendList_z.resize(sendCount_z, 0);
        sendList_X.resize(sendCount_X, 0);
        sendList_Y.resize(sendCount_Y, 0);
        sendList_Z.resize(sendCount_Z, 0);
        sendList_xy.resize(sendCount_xy, 0);
        sendList_yz.resize(sendCount_yz, 0);
        sendList_xz.resize(sendCount_xz, 0);
        sendList_Xy.resize(sendCount_Xy, 0);
        sendList_Yz.resize(sendCount_Yz, 0);
        sendList_xZ.resize(sendCount_xZ, 0);
        sendList_xY.resize(sendCount_xY, 0);
        sendList_yZ.resize(sendCount_yZ, 0);
        sendList_Xz.resize(sendCount_Xz, 0);
        sendList_XY.resize(sendCount_XY, 0);
        sendList_YZ.resize(sendCount_YZ, 0);
        sendList_XZ.resize(sendCount_XZ, 0);
        // Populate the send list
        sendCount_x = sendCount_y = sendCount_z = sendCount_X = sendCount_Y = sendCount_Z = 0;
        sendCount_xy = sendCount_yz = sendCount_xz = sendCount_Xy = sendCount_Yz = sendCount_xZ = 0;
        sendCount_xY = sendCount_yZ = sendCount_Xz = sendCount_XY = sendCount_YZ = sendCount_XZ = 0;
        for (k = 1; k < Nz - 1; k++) {
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    // Local value to send
                    n = k * Nx * Ny + j * Nx + i;
                    // Counts for the six faces
                    if (i == 1)
                        sendList_x[sendCount_x++] = n;
                    if (j == 1)
                        sendList_y[sendCount_y++] = n;
                    if (k == 1)
                        sendList_z[sendCount_z++] = n;
                    if (i == Nx - 2)
                        sendList_X[sendCount_X++] = n;
                    if (j == Ny - 2)
                        sendList_Y[sendCount_Y++] = n;
                    if (k == Nz - 2)
                        sendList_Z[sendCount_Z++] = n;
                    // Counts for the twelve edges
                    if (i == 1 && j == 1)
                        sendList_xy[sendCount_xy++] = n;
                    if (i == 1 && j == Ny - 2)
                        sendList_xY[sendCount_xY++] = n;
                    if (i == Nx - 2 && j == 1)
                        sendList_Xy[sendCount_Xy++] = n;
                    if (i == Nx - 2 && j == Ny - 2)
                        sendList_XY[sendCount_XY++] = n;

                    if (i == 1 && k == 1)
                        sendList_xz[sendCount_xz++] = n;
                    if (i == 1 && k == Nz - 2)
                        sendList_xZ[sendCount_xZ++] = n;
                    if (i == Nx - 2 && k == 1)
                        sendList_Xz[sendCount_Xz++] = n;
                    if (i == Nx - 2 && k == Nz - 2)
                        sendList_XZ[sendCount_XZ++] = n;

                    if (j == 1 && k == 1)
                        sendList_yz[sendCount_yz++] = n;
                    if (j == 1 && k == Nz - 2)
                        sendList_yZ[sendCount_yZ++] = n;
                    if (j == Ny - 2 && k == 1)
                        sendList_Yz[sendCount_Yz++] = n;
                    if (j == Ny - 2 && k == Nz - 2)
                        sendList_YZ[sendCount_YZ++] = n;
                }
            }
        }

        //......................................................................................
        int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y, recvCount_Z;
        int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz, recvCount_xZ;
        int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ, recvCount_XZ;
        req1[0] = Isend( Comm, &sendCount_x, 1, rank_x(), 0 );
        req2[0] = Irecv( Comm, &recvCount_X, 1, rank_X(), 0 );
        req1[1] = Isend( Comm, &sendCount_X, 1, rank_X(), 1 );
        req2[1] = Irecv( Comm, &recvCount_x, 1, rank_x(), 1 );
        req1[2] = Isend( Comm, &sendCount_y, 1, rank_y(), 2 );
        req2[2] = Irecv( Comm, &recvCount_Y, 1, rank_Y(), 2 );
        req1[3] = Isend( Comm, &sendCount_Y, 1, rank_Y(), 3 );
        req2[3] = Irecv( Comm, &recvCount_y, 1, rank_y(), 3 );
        req1[4] = Isend( Comm, &sendCount_z, 1, rank_z(), 4 );
        req2[4] = Irecv( Comm, &recvCount_Z, 1, rank_Z(), 4 );
        req1[5] = Isend( Comm, &sendCount_Z, 1, rank_Z(), 5 );
        req2[5] = Irecv( Comm, &recvCount_z, 1, rank_z(), 5 );
        req1[6] = Isend( Comm, &sendCount_xy, 1, rank_xy(), 6 );
        req2[6] = Irecv( Comm, &recvCount_XY, 1, rank_XY(), 6 );
        req1[7] = Isend( Comm, &sendCount_XY, 1, rank_XY(), 7 );
        req2[7] = Irecv( Comm, &recvCount_xy, 1, rank_xy(), 7 );
        req1[8] = Isend( Comm, &sendCount_Xy, 1, rank_Xy(), 8 );
        req2[8] = Irecv( Comm, &recvCount_xY, 1, rank_xY(), 8 );
        req1[9] = Isend( Comm, &sendCount_xY, 1, rank_xY(), 9 );
        req2[9] = Irecv( Comm, &recvCount_Xy, 1, rank_Xy(), 9 );
        req1[10] = Isend( Comm, &sendCount_xz, 1, rank_xz(), 10 );
        req2[10] = Irecv( Comm, &recvCount_XZ, 1, rank_XZ(), 10 );
        req1[11] = Isend( Comm, &sendCount_XZ, 1, rank_XZ(), 11 );
        req2[11] = Irecv( Comm, &recvCount_xz, 1, rank_xz(), 11 );
        req1[12] = Isend( Comm, &sendCount_Xz, 1, rank_Xz(), 12 );
        req2[12] = Irecv( Comm, &recvCount_xZ, 1, rank_xZ(), 12 );
        req1[13] = Isend( Comm, &sendCount_xZ, 1, rank_xZ(), 13 );
        req2[13] = Irecv( Comm, &recvCount_Xz, 1, rank_Xz(), 13 );
        req1[14] = Isend( Comm, &sendCount_yz, 1, rank_yz(), 14 );
        req2[14] = Irecv( Comm, &recvCount_YZ, 1, rank_YZ(), 14 );
        req1[15] = Isend( Comm, &sendCount_YZ, 1, rank_YZ(), 15 );
        req2[15] = Irecv( Comm, &recvCount_yz, 1, rank_yz(), 15 );
        req1[16] = Isend( Comm, &sendCount_Yz, 1, rank_Yz(), 16 );
        req2[16] = Irecv( Comm, &recvCount_yZ, 1, rank_yZ(), 16 );
        req1[17] = Isend( Comm, &sendCount_yZ, 1, rank_yZ(), 17 );
        req2[17] = Irecv( Comm, &recvCount_Yz, 1, rank_Yz(), 17 );
        Comm.waitAll(18, req1);
        Comm.waitAll(18, req2);
        Comm.barrier();
        // allocate recv lists
        recvList_x.resize(recvCount_x, 0);
        recvList_y.resize(recvCount_y, 0);
        recvList_z.resize(recvCount_z, 0);
        recvList_X.resize(recvCount_X, 0);
        recvList_Y.resize(recvCount_Y, 0);
        recvList_Z.resize(recvCount_Z, 0);
        recvList_xy.resize(recvCount_xy, 0);
        recvList_yz.resize(recvCount_yz, 0);
        recvList_xz.resize(recvCount_xz, 0);
        recvList_Xy.resize(recvCount_Xy, 0);
        recvList_Yz.resize(recvCount_Yz, 0);
        recvList_xZ.resize(recvCount_xZ, 0);
        recvList_xY.resize(recvCount_xY, 0);
        recvList_yZ.resize(recvCount_yZ, 0);
        recvList_Xz.resize(recvCount_Xz, 0);
        recvList_XY.resize(recvCount_XY, 0);
        recvList_YZ.resize(recvCount_YZ, 0);
        recvList_XZ.resize(recvCount_XZ, 0);
        //......................................................................................
        req1[0] = Isend( Comm, sendList_x.data(), sendCount_x, rank_x(), 0 );
        req2[0] = Irecv( Comm, recvList_X.data(), recvCount_X, rank_X(), 0 );
        req1[1] = Isend( Comm, sendList_X.data(), sendCount_X, rank_X(), 1 );
        req2[1] = Irecv( Comm, recvList_x.data(), recvCount_x, rank_x(), 1 );
        req1[2] = Isend( Comm, sendList_y.data(), sendCount_y, rank_y(), 2 );
        req2[2] = Irecv( Comm, recvList_Y.data(), recvCount_Y, rank_Y(), 2 );
        req1[3] = Isend( Comm, sendList_Y.data(), sendCount_Y, rank_Y(), 3 );
        req2[3] = Irecv( Comm, recvList_y.data(), recvCount_y, rank_y(), 3 );
        req1[4] = Isend( Comm, sendList_z.data(), sendCount_z, rank_z(), 4 );
        req2[4] = Irecv( Comm, recvList_Z.data(), recvCount_Z, rank_Z(), 4 );
        req1[5] = Isend( Comm, sendList_Z.data(), sendCount_Z, rank_Z(), 5 );
        req2[5] = Irecv( Comm, recvList_z.data(), recvCount_z, rank_z(), 5 );
        req1[6] = Isend( Comm, sendList_xy.data(), sendCount_xy, rank_xy(), 6 );
        req2[6] = Irecv( Comm, recvList_XY.data(), recvCount_XY, rank_XY(), 6 );
        req1[7] = Isend( Comm, sendList_XY.data(), sendCount_XY, rank_XY(), 7 );
        req2[7] = Irecv( Comm, recvList_xy.data(), recvCount_xy, rank_xy(), 7 );
        req1[8] = Isend( Comm, sendList_Xy.data(), sendCount_Xy, rank_Xy(), 8 );
        req2[8] = Irecv( Comm, recvList_xY.data(), recvCount_xY, rank_xY(), 8 );
        req1[9] = Isend( Comm, sendList_xY.data(), sendCount_xY, rank_xY(), 9 );
        req2[9] = Irecv( Comm, recvList_Xy.data(), recvCount_Xy, rank_Xy(), 9 );
        req1[10] = Isend( Comm, sendList_xz.data(), sendCount_xz, rank_xz(), 10 );
        req2[10] = Irecv( Comm, recvList_XZ.data(), recvCount_XZ, rank_XZ(), 10 );
        req1[11] = Isend( Comm, sendList_XZ.data(), sendCount_XZ, rank_XZ(), 11 );
        req2[11] = Irecv( Comm, recvList_xz.data(), recvCount_xz, rank_xz(), 11 );
        req1[12] = Isend( Comm, sendList_Xz.data(), sendCount_Xz, rank_Xz(), 12 );
        req2[12] = Irecv( Comm, recvList_xZ.data(), recvCount_xZ, rank_xZ(), 12 );
        req1[13] = Isend( Comm, sendList_xZ.data(), sendCount_xZ, rank_xZ(), 13 );
        req2[13] = Irecv( Comm, recvList_Xz.data(), recvCount_Xz, rank_Xz(), 13 );
        req1[14] = Isend( Comm, sendList_yz.data(), sendCount_yz, rank_yz(), 14 );
        req2[14] = Irecv( Comm, recvList_YZ.data(), recvCount_YZ, rank_YZ(), 14 );
        req1[15] = Isend( Comm, sendList_YZ.data(), sendCount_YZ, rank_YZ(), 15 );
        req2[15] = Irecv( Comm, recvList_yz.data(), recvCount_yz, rank_yz(), 15 );
        req1[16] = Isend( Comm, sendList_Yz.data(), sendCount_Yz, rank_Yz(), 16 );
        req2[16] = Irecv( Comm, recvList_yZ.data(), recvCount_yZ, rank_yZ(), 16 );
        req1[17] = Isend( Comm, sendList_yZ.data(), sendCount_yZ, rank_yZ(), 17 );
        req2[17] = Irecv( Comm, recvList_Yz.data(), recvCount_Yz, rank_Yz(), 17 );
        Comm.waitAll(18, req1);
        Comm.waitAll(18, req2);
        //......................................................................................
        for (int idx = 0; idx < recvCount_x; idx++)
            recvList_x[idx] -= (Nx - 2);
        for (int idx = 0; idx < recvCount_X; idx++)
            recvList_X[idx] += (Nx - 2);
        for (int idx = 0; idx < recvCount_y; idx++)
            recvList_y[idx] -= (Ny - 2) * Nx;
        for (int idx = 0; idx < recvCount_Y; idx++)
            recvList_Y[idx] += (Ny - 2) * Nx;
        for (int idx = 0; idx < recvCount_z; idx++)
            recvList_z[idx] -= (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_Z; idx++)
            recvList_Z[idx] += (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_xy; idx++)
            recvList_xy[idx] -= (Nx - 2) + (Ny - 2) * Nx;
        for (int idx = 0; idx < recvCount_XY; idx++)
            recvList_XY[idx] += (Nx - 2) + (Ny - 2) * Nx;
        for (int idx = 0; idx < recvCount_xY; idx++)
            recvList_xY[idx] -= (Nx - 2) - (Ny - 2) * Nx;
        for (int idx = 0; idx < recvCount_Xy; idx++)
            recvList_Xy[idx] += (Nx - 2) - (Ny - 2) * Nx;
        for (int idx = 0; idx < recvCount_xz; idx++)
            recvList_xz[idx] -= (Nx - 2) + (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_XZ; idx++)
            recvList_XZ[idx] += (Nx - 2) + (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_xZ; idx++)
            recvList_xZ[idx] -= (Nx - 2) - (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_Xz; idx++)
            recvList_Xz[idx] += (Nx - 2) - (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_yz; idx++)
            recvList_yz[idx] -= (Ny - 2) * Nx + (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_YZ; idx++)
            recvList_YZ[idx] += (Ny - 2) * Nx + (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_yZ; idx++)
            recvList_yZ[idx] -= (Ny - 2) * Nx - (Nz - 2) * Nx * Ny;
        for (int idx = 0; idx < recvCount_Yz; idx++)
            recvList_Yz[idx] += (Ny - 2) * Nx - (Nz - 2) * Nx * Ny;
    }

private:

    MPI_Request req1[18], req2[18];
    std::vector<int> sendList_x, sendList_y, sendList_z, sendList_X, sendList_Y, sendList_Z;
    std::vector<int> sendList_xy, sendList_yz, sendList_xz, sendList_Xy, sendList_Yz, sendList_xZ;
    std::vector<int> sendList_xY, sendList_yZ, sendList_Xz, sendList_XY, sendList_YZ, sendList_XZ;
    std::vector<int> recvList_x, recvList_y, recvList_z, recvList_X, recvList_Y, recvList_Z;
    std::vector<int> recvList_xy, recvList_yz, recvList_xz, recvList_Xy, recvList_Yz, recvList_xZ;
    std::vector<int> recvList_xY, recvList_yZ, recvList_Xz, recvList_XY, recvList_YZ, recvList_XZ;
    const std::vector<int> &getRecvList(const char *dir) const;
    const std::vector<int> &getSendList(const char *dir) const;
};



int main(int argc, char **argv)
{
    Utilities::startup( argc, argv, true );
    Utilities::MPI comm( MPI_COMM_WORLD );
    {
        auto filename = argv[1];
        auto input_db = std::make_shared<Database>( filename );
        auto db = input_db->getDatabase( "Domain" );
        auto Dm  = std::make_shared<Domain2>(db,comm);
        Dm->CommInit();
        std::cout << "step 1" << std::endl << std::flush;
    }
    std::cout << "step 2" << std::endl << std::flush;
    comm.barrier();
    std::cout << "step 3" << std::endl << std::flush;
    Utilities::shutdown();
    std::cout << "step 4" << std::endl << std::flush;
    return 0;
}
