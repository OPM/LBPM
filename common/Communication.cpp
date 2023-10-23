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
#include "common/Communication.h"

/********************************************************
*  Structure to store the rank info                     *
********************************************************/
int RankInfoStruct::getRankForBlock(int i, int j, int k) const {
    int i2 = (i + nx) % nx;
    int j2 = (j + ny) % ny;
    int k2 = (k + nz) % nz;
    return i2 + j2 * nx + k2 * nx * ny;
}
RankInfoStruct::RankInfoStruct() {
    nx = 0;
    ny = 0;
    nz = 0;
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
}
RankInfoStruct::RankInfoStruct(int rank0, int nprocx, int nprocy, int nprocz) {
    memset(this, 0, sizeof(RankInfoStruct));
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

/********************************************************
*  Deprecated functions                                 *
********************************************************/
void InitializeRanks(const int rank, const int nprocx, const int nprocy,
                     const int nprocz, int &iproc, int &jproc, int &kproc,
                     int &rank_x, int &rank_y, int &rank_z, int &rank_X,
                     int &rank_Y, int &rank_Z, int &rank_xy, int &rank_XY,
                     int &rank_xY, int &rank_Xy, int &rank_xz, int &rank_XZ,
                     int &rank_xZ, int &rank_Xz, int &rank_yz, int &rank_YZ,
                     int &rank_yZ, int &rank_Yz) {
    const RankInfoStruct data(rank, nprocx, nprocy, nprocz);
    iproc = data.ix;
    jproc = data.jy;
    kproc = data.kz;
    rank_X = data.rank[2][1][1];
    rank_x = data.rank[0][1][1];
    rank_Y = data.rank[1][2][1];
    rank_y = data.rank[1][0][1];
    rank_Z = data.rank[1][1][2];
    rank_z = data.rank[1][1][0];
    rank_XY = data.rank[2][2][1];
    rank_xy = data.rank[0][0][1];
    rank_Xy = data.rank[2][0][1];
    rank_xY = data.rank[0][2][1];
    rank_XZ = data.rank[2][1][2];
    rank_xz = data.rank[0][1][0];
    rank_Xz = data.rank[2][1][0];
    rank_xZ = data.rank[0][1][2];
    rank_YZ = data.rank[1][2][2];
    rank_yz = data.rank[1][0][0];
    rank_Yz = data.rank[1][2][0];
    rank_yZ = data.rank[1][0][2];
}
