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
#ifndef COMMUNICATION_HPP_INC
#define COMMUNICATION_HPP_INC

#include "common/Communication.h"
#include "common/MPI.h"
#include "common/Utilities.h"

/********************************************************
*  Redistribute data between two grids                  *
********************************************************/
template <class TYPE>
Array<TYPE>
redistribute(const RankInfoStruct &src_rank, const Array<TYPE> &src_data,
             const RankInfoStruct &dst_rank, std::array<int, 3> dst_size,
             const Utilities::MPI &comm) {
    if (comm.getSize() == 1) {
        return src_data.subset({0, (size_t)dst_size[0] - 1, 0,
                                (size_t)dst_size[1] - 1, 0,
                                (size_t)dst_size[2] - 1});
    }
    // Get the src size
    std::array<int, 3> src_size;
    int size0[3] = {(int)src_data.size(0), (int)src_data.size(1),
                    (int)src_data.size(2)};
    comm.maxReduce(size0, src_size.data(), 3);
    if (!src_data.empty())
        ASSERT(src_size[0] == size0[0] && src_size[1] == size0[1] &&
               src_size[2] == size0[2]);
    // Check that dst_size matches on all ranks
    comm.maxReduce(dst_size.data(), size0, 3);
    ASSERT(dst_size[0] == size0[0] && dst_size[1] == size0[1] &&
           dst_size[2] == size0[2]);
    // Function to get overlap range
    auto calcOverlap = [](int i1[3], int i2[3], int j1[3], int j2[3]) {
        std::vector<size_t> index;
        if (i1[0] > j2[0] || i2[0] < j1[0] || i1[1] > j2[1] || i2[1] < j1[1] ||
            i1[2] > j2[2] || i2[2] < j1[2])
            return index;
        index.resize(6);
        index[0] = std::max(j1[0] - i1[0], 0);
        index[1] = std::min(j2[0] - i1[0], i2[0] - i1[0]);
        index[2] = std::max(j1[1] - i1[1], 0);
        index[3] = std::min(j2[1] - i1[1], i2[1] - i1[1]);
        index[4] = std::max(j1[2] - i1[2], 0);
        index[5] = std::min(j2[2] - i1[2], i2[2] - i1[2]);
        return index;
    };
    // Pack and send my data to the appropriate ranks (including myself)
    std::vector<int> send_rank;
    std::vector<Array<TYPE>> send_data;
    if (!src_data.empty()) {
        int i1[3] = {src_size[0] * src_rank.ix, src_size[1] * src_rank.jy,
                     src_size[2] * src_rank.kz};
        int i2[3] = {i1[0] + src_size[0] - 1, i1[1] + src_size[1] - 1,
                     i1[2] + src_size[2] - 1};
        for (int i = 0; i < dst_rank.nx; i++) {
            for (int j = 0; j < dst_rank.ny; j++) {
                for (int k = 0; k < dst_rank.nz; k++) {
                    int j1[3] = {i * dst_size[0], j * dst_size[1],
                                 k * dst_size[2]};
                    int j2[3] = {j1[0] + dst_size[0] - 1,
                                 j1[1] + dst_size[1] - 1,
                                 j1[2] + dst_size[2] - 1};
                    auto index = calcOverlap(i1, i2, j1, j2);
                    if (index.empty())
                        continue;
                    send_rank.push_back(dst_rank.getRankForBlock(i, j, k));
                    send_data.push_back(src_data.subset(index));
                }
            }
        }
    }
    std::vector<MPI_Request> send_request(send_rank.size());
    for (size_t i = 0; i < send_rank.size(); i++)
        send_request[i] = comm.Isend(send_data[i].data(), send_data[i].length(),
                                     send_rank[i], 5462);
    // Unpack data from the appropriate ranks (including myself)
    Array<TYPE> dst_data(dst_size[0], dst_size[1], dst_size[2]);
    int i1[3] = {dst_size[0] * dst_rank.ix, dst_size[1] * dst_rank.jy,
                 dst_size[2] * dst_rank.kz};
    int i2[3] = {i1[0] + dst_size[0] - 1, i1[1] + dst_size[1] - 1,
                 i1[2] + dst_size[2] - 1};
    for (int i = 0; i < src_rank.nx; i++) {
        for (int j = 0; j < src_rank.ny; j++) {
            for (int k = 0; k < src_rank.nz; k++) {
                int j1[3] = {i * src_size[0], j * src_size[1], k * src_size[2]};
                int j2[3] = {j1[0] + src_size[0] - 1, j1[1] + src_size[1] - 1,
                             j1[2] + src_size[2] - 1};
                auto index = calcOverlap(i1, i2, j1, j2);
                if (index.empty())
                    continue;
                int rank = src_rank.getRankForBlock(i, j, k);
                Array<TYPE> data(index[1] - index[0] + 1,
                                 index[3] - index[2] + 1,
                                 index[5] - index[4] + 1);
                comm.recv(data.data(), data.length(), rank, 5462);
                dst_data.copySubset(index, data);
            }
        }
    }
    // Free data
    comm.waitAll(send_request.size(), send_request.data());
    return dst_data;
}

/********************************************************
*  Structure to fill halo cells                         *
********************************************************/
template <class TYPE>
fillHalo<TYPE>::fillHalo(const Utilities::MPI &comm_,
                         const RankInfoStruct &info_, std::array<int, 3> n_,
                         std::array<int, 3> ng_, int tag0, int depth_,
                         std::array<bool, 3> fill, std::array<bool, 3> periodic)
    : comm(comm_), info(info_), n(n_), ng(ng_), depth(depth_) {
    // Set the fill pattern
    memset(fill_pattern, 0, sizeof(fill_pattern));
    if (fill[0]) {
        fill_pattern[0][1][1] = true;
        fill_pattern[2][1][1] = true;
        fill_pattern[1][0][1] = true;
        fill_pattern[1][2][1] = true;
        fill_pattern[1][1][0] = true;
        fill_pattern[1][1][2] = true;
    }
    if (fill[1]) {
        fill_pattern[0][0][1] = true;
        fill_pattern[0][2][1] = true;
        fill_pattern[2][0][1] = true;
        fill_pattern[2][2][1] = true;
        fill_pattern[0][1][0] = true;
        fill_pattern[0][1][2] = true;
        fill_pattern[2][1][0] = true;
        fill_pattern[2][1][2] = true;
        fill_pattern[1][0][0] = true;
        fill_pattern[1][0][2] = true;
        fill_pattern[1][2][0] = true;
        fill_pattern[1][2][2] = true;
    }
    if (fill[2]) {
        fill_pattern[0][0][0] = true;
        fill_pattern[0][0][2] = true;
        fill_pattern[0][2][0] = true;
        fill_pattern[0][2][2] = true;
        fill_pattern[2][0][0] = true;
        fill_pattern[2][0][2] = true;
        fill_pattern[2][2][0] = true;
        fill_pattern[2][2][2] = true;
    }
    // Remove communication for non-perioidic directions
    if (!periodic[0] && info.ix == 0) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++)
                fill_pattern[0][j][k] = false;
        }
    }
    if (!periodic[0] && info.ix == info.nx - 1) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++)
                fill_pattern[2][j][k] = false;
        }
    }
    if (!periodic[1] && info.jy == 0) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++)
                fill_pattern[i][0][k] = false;
        }
    }
    if (!periodic[1] && info.jy == info.ny - 1) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++)
                fill_pattern[i][2][k] = false;
        }
    }
    if (!periodic[2] && info.kz == 0) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++)
                fill_pattern[i][j][0] = false;
        }
    }
    if (!periodic[2] && info.kz == info.nz - 1) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++)
                fill_pattern[i][j][2] = false;
        }
    }
    // Determine the number of elements for each send/recv
    for (int i = 0; i < 3; i++) {
        int ni = (i - 1) == 0 ? n[0] : ng[0];
        for (int j = 0; j < 3; j++) {
            int nj = (j - 1) == 0 ? n[1] : ng[1];
            for (int k = 0; k < 3; k++) {
                int nk = (k - 1) == 0 ? n[2] : ng[2];
                if (fill_pattern[i][j][k])
                    N_send_recv[i][j][k] = ni * nj * nk;
                else
                    N_send_recv[i][j][k] = 0;
            }
        }
    }
    // Create send/recv buffers
    size_t N_mem = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++)
                N_mem += N_send_recv[i][j][k];
        }
    }
    mem = new TYPE[2 * depth * N_mem];
    size_t index = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                send[i][j][k] = &mem[index];
                index += depth * N_send_recv[i][j][k];
                recv[i][j][k] = &mem[index];
                index += depth * N_send_recv[i][j][k];
            }
        }
    }
    // Create the tags
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                tag[i][j][k] = tag0 + i + j * 3 + k * 9;
            }
        }
    }
}
template <class TYPE> fillHalo<TYPE>::~fillHalo() { delete[] mem; }
template <class TYPE> void fillHalo<TYPE>::fill(Array<TYPE> &data) {
    //PROFILE_START("fillHalo::fill",1);
    int depth2 = data.size(3);
    ASSERT((int)data.size(0) == n[0] + 2 * ng[0]);
    ASSERT((int)data.size(1) == n[1] + 2 * ng[1]);
    ASSERT((int)data.size(2) == n[2] + 2 * ng[2]);
    ASSERT(depth2 <= depth);
    ASSERT(data.ndim() == 3 || data.ndim() == 4);
    // Start the recieves
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (!fill_pattern[i][j][k])
                    continue;
                recv_req[i][j][k] =
                    comm.Irecv(recv[i][j][k], depth2 * N_send_recv[i][j][k],
                               info.rank[i][j][k], tag[2 - i][2 - j][2 - k]);
            }
        }
    }
    // Pack the src data and start the sends
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (!fill_pattern[i][j][k])
                    continue;
                pack(data, i - 1, j - 1, k - 1, send[i][j][k]);
                send_req[i][j][k] =
                    comm.Isend(send[i][j][k], depth2 * N_send_recv[i][j][k],
                               info.rank[i][j][k], tag[i][j][k]);
            }
        }
    }
    // Recv the dst data and unpack (we recive in reverse order to match the sends)
    for (int i = 2; i >= 0; i--) {
        for (int j = 2; j >= 0; j--) {
            for (int k = 2; k >= 0; k--) {
                if (!fill_pattern[i][j][k])
                    continue;
                comm.wait(recv_req[i][j][k]);
                unpack(data, i - 1, j - 1, k - 1, recv[i][j][k]);
            }
        }
    }
    // Wait until all sends have completed
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (!fill_pattern[i][j][k])
                    continue;
                comm.wait(send_req[i][j][k]);
            }
        }
    }
    //PROFILE_STOP("fillHalo::fill",1);
}
template <class TYPE>
void fillHalo<TYPE>::pack(const Array<TYPE> &data, int i0, int j0, int k0,
                          TYPE *buffer) {
    int depth2 = data.size(3);
    int ni = i0 == 0 ? n[0] : ng[0];
    int nj = j0 == 0 ? n[1] : ng[1];
    int nk = k0 == 0 ? n[2] : ng[2];
    int is = i0 == 0 ? ng[0] : ((i0 == -1) ? ng[0] : n[0]);
    int js = j0 == 0 ? ng[1] : ((j0 == -1) ? ng[1] : n[1]);
    int ks = k0 == 0 ? ng[2] : ((k0 == -1) ? ng[2] : n[2]);
    for (int d = 0; d < depth2; d++) {
        for (int k = 0; k < nk; k++) {
            for (int j = 0; j < nj; j++) {
                for (int i = 0; i < ni; i++) {
                    buffer[i + j * ni + k * ni * nj + d * ni * nj * nk] =
                        data(i + is, j + js, k + ks, d);
                }
            }
        }
    }
}
template <class TYPE>
void fillHalo<TYPE>::unpack(Array<TYPE> &data, int i0, int j0, int k0,
                            const TYPE *buffer) {
    int depth2 = data.size(3);
    int ni = i0 == 0 ? n[0] : ng[0];
    int nj = j0 == 0 ? n[1] : ng[1];
    int nk = k0 == 0 ? n[2] : ng[2];
    int is = i0 == 0 ? ng[0] : ((i0 == -1) ? 0 : n[0] + ng[0]);
    int js = j0 == 0 ? ng[1] : ((j0 == -1) ? 0 : n[1] + ng[1]);
    int ks = k0 == 0 ? ng[2] : ((k0 == -1) ? 0 : n[2] + ng[2]);
    for (int d = 0; d < depth2; d++) {
        for (int k = 0; k < nk; k++) {
            for (int j = 0; j < nj; j++) {
                for (int i = 0; i < ni; i++) {
                    data(i + is, j + js, k + ks, d) =
                        buffer[i + j * ni + k * ni * nj + d * ni * nj * nk];
                }
            }
        }
    }
}

/********************************************************
*  Function to remove the ghost halo                    *
********************************************************/
template <class TYPE>
template <class TYPE1, class TYPE2>
void fillHalo<TYPE>::copy(const Array<TYPE1> &src, Array<TYPE2> &dst) {
    //PROFILE_START("fillHalo::copy",1);
    ASSERT((int)src.size(0) == n[0] || (int)src.size(0) == n[0] + 2 * ng[0]);
    ASSERT((int)dst.size(0) == n[0] || (int)dst.size(0) == n[0] + 2 * ng[0]);
    bool src_halo = (int)src.size(0) == n[0] + 2 * ng[0];
    bool dst_halo = (int)dst.size(0) == n[0] + 2 * ng[0];
    if (src_halo) {
        ASSERT((int)src.size(0) == n[0] + 2 * ng[0]);
        ASSERT((int)src.size(1) == n[1] + 2 * ng[1]);
        ASSERT((int)src.size(2) == n[2] + 2 * ng[2]);
    } else {
        ASSERT((int)src.size(0) == n[0]);
        ASSERT((int)src.size(1) == n[1]);
        ASSERT((int)src.size(2) == n[2]);
    }
    if (dst_halo) {
        ASSERT((int)dst.size(0) == n[0] + 2 * ng[0]);
        ASSERT((int)dst.size(1) == n[1] + 2 * ng[1]);
        ASSERT((int)dst.size(2) == n[2] + 2 * ng[2]);
    } else {
        ASSERT((int)dst.size(0) == n[0]);
        ASSERT((int)dst.size(1) == n[1]);
        ASSERT((int)dst.size(2) == n[2]);
    }
    if (src_halo == dst_halo) {
        // Src and dst halos match
        for (size_t i = 0; i < src.length(); i++)
            dst(i) = src(i);
    } else if (src_halo && !dst_halo) {
        // Src has halos
        for (int k = 0; k < n[2]; k++) {
            for (int j = 0; j < n[1]; j++) {
                for (int i = 0; i < n[0]; i++) {
                    dst(i, j, k) = src(i + ng[0], j + ng[1], k + ng[2]);
                }
            }
        }
    } else if (!src_halo && dst_halo) {
        // Dst has halos
        for (int k = 0; k < n[2]; k++) {
            for (int j = 0; j < n[1]; j++) {
                for (int i = 0; i < n[0]; i++) {
                    dst(i + ng[0], j + ng[1], k + ng[2]) = src(i, j, k);
                }
            }
        }
        fill(dst);
    }
    //PROFILE_STOP("fillHalo::copy",1);
}

#endif
