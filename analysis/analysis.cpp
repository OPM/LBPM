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
#include "analysis/analysis.h"
#include "ProfilerApp.h"

#include <algorithm>
#include <iostream>

template <class TYPE> inline TYPE *getPtr(std::vector<TYPE> &x) {
    return x.empty() ? NULL : &x[0];
}
template <class TYPE> inline const TYPE *getPtr(const std::vector<TYPE> &x) {
    return x.empty() ? NULL : &x[0];
}

/******************************************************************
* Compute the blobs                                               *
******************************************************************/
int ComputeBlob(const Array<bool> &isPhase, BlobIDArray &LocalBlobID,
                bool periodic, int start_id) {
    PROFILE_START("ComputeBlob", 1);
    ASSERT(isPhase.size() == LocalBlobID.size());
    const int Nx = isPhase.size(0); // maxima for the meshes
    const int Ny = isPhase.size(1);
    const int Nz = isPhase.size(2);
    std::vector<int> map;
    map.reserve(128);
    // Get the list of neighbors we need to check
    int N_neighbors = 0;
    int d[26][3];
    bool include_corners =
        false; // Do we need to include cells that only touch at a corder/edge
    if (include_corners) {
        // Include corners/edges as neighbors, check all cells
        N_neighbors = 26;
        const int tmp[26][3] = {
            {1, 0, 0},   {-1, 0, 0},  {0, 1, 0},   {0, -1, 0},  {0, 0, 1},
            {0, 0, -1},  {1, 1, 0},   {1, -1, 0},  {-1, 1, 0},  {-1, -1, 0},
            {1, 0, 1},   {-1, 0, 1},  {1, 0, -1},  {-1, 0, -1}, {0, 1, 1},
            {0, -1, 1},  {0, 1, -1},  {0, -1, -1}, {1, 1, 1},   {1, 1, -1},
            {1, -1, 1},  {1, -1, -1}, {-1, 1, 1},  {-1, 1, -1}, {-1, -1, 1},
            {-1, -1, -1}}; // directions to neighbors
        memcpy(d, tmp, sizeof(tmp));
    } else {
        // Do not include corners/edges as neighbors
        if (periodic) {
            // Include all neighbors for periodic problems
            N_neighbors = 6;
            const int tmp[6][3] = {
                {1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                {0, -1, 0}, {0, 0, 1},  {0, 0, -1}}; // directions to neighbors
            memcpy(d, tmp, sizeof(tmp));
        } else {
            // We only need to include the lower neighbors for non-periodic problems
            N_neighbors = 3;
            const int tmp[3][3] = {
                {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}; // directions to neighbors
            memcpy(d, tmp, sizeof(tmp));
        }
    }
    // Loop through all the points
    int last = start_id - 1;
    std::vector<int> neighbor_ids;
    neighbor_ids.reserve(N_neighbors);
    const bool *isPhasePtr = isPhase.data();
    BlobIDType *LocalBlobIDPtr = LocalBlobID.data();
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int index = x + y * Nx + z * Nx * Ny;
                if (!isPhasePtr[index])
                    continue;
                // Get all neighbor indicies
                int N_list = 0, neighbor_ids[26];
                for (int p = 0; p < N_neighbors; p++) {
                    // Get the neighbor index
                    int x2 = x + d[p][0];
                    int y2 = y + d[p][1];
                    int z2 = z + d[p][2];
                    if (periodic) {
                        x2 = x2 < 0 ? Nx - 1 : x2; // Periodic BC for x
                        x2 = x2 > Nx - 1 ? 0 : x2;
                        y2 = y2 < 0 ? Ny - 1 : y2; // Periodic BC for x
                        y2 = y2 > Ny - 1 ? 0 : y2;
                        z2 = z2 < 0 ? Nz - 1 : z2; // Periodic BC for x
                        z2 = z2 > Nz - 1 ? 0 : z2;
                    } else {
                        if (x2 < 0 || x2 >= Nx || y2 < 0 || y2 >= Ny ||
                            z2 < 0 || z2 >= Nz)
                            continue;
                    }
                    // Check if a neighbor has a known blob id
                    size_t index2 = x2 + y2 * Nx + z2 * Nx * Ny;
                    int id = LocalBlobIDPtr[index2];
                    if (!isPhasePtr[index2] || id < 0)
                        continue;
                    neighbor_ids[N_list] = id;
                    N_list++;
                }
                if (N_list == 0) {
                    // No neighbors with a blob id, create a new one
                    LocalBlobIDPtr[index] = last + 1;
                    map.push_back(last + 1);
                    last++;
                } else if (N_list == 1) {
                    // We have one neighbor
                    LocalBlobIDPtr[index] = neighbor_ids[0];
                } else {
                    // We have multiple neighbors
                    int id = neighbor_ids[0];
                    for (int i = 1; i < N_list; i++)
                        id = std::min(id, neighbor_ids[i]);
                    LocalBlobIDPtr[index] = id;
                    for (int i = 0; i < N_list; i++)
                        map[neighbor_ids[i] - start_id] =
                            std::min(map[neighbor_ids[i] - start_id], id);
                }
            }
        }
    }
    // Collapse the ids that map to another id
    last = start_id - 1;
    for (int i = 0; i < (int)map.size(); i++) {
        if (map[i] == i + start_id) {
            map[i] = last + 1;
            last++;
        } else {
            ASSERT(map[i] < i + start_id);
            map[i] = map[map[i] - start_id];
        }
    }
    for (int i = 0; i < Nx * Ny * Nz; i++) {
        if (isPhasePtr[i]) {
            LocalBlobIDPtr[i] = map[LocalBlobIDPtr[i] - start_id];
        }
    }
    PROFILE_STOP("ComputeBlob", 1);
    return last - start_id + 1;
}

/******************************************************************
* Compute the local blob ids                                      *
******************************************************************/
int ComputeLocalBlobIDs(const DoubleArray &Phase, const DoubleArray &SignDist,
                        double vF, double vS, BlobIDArray &LocalBlobID,
                        bool periodic) {
    PROFILE_START("ComputeLocalBlobIDs");
    ASSERT(SignDist.size() == Phase.size());
    size_t Nx = Phase.size(0);
    size_t Ny = Phase.size(1);
    size_t Nz = Phase.size(2);
    // Initialize output
    LocalBlobID.resize(Nx, Ny, Nz);
    // Compute the local blob ids
    size_t N = Nx * Ny * Nz;
    Array<bool> isPhase(Nx, Ny, Nz);
    memset(isPhase.data(), 0, Nx * Ny * Nz * sizeof(bool));
    for (size_t i = 0; i < N; i++) {
        if (SignDist(i) <= vS) {
            // Solid phase
            LocalBlobID(i) = -2;
        } else {
            LocalBlobID(i) = -1;
            if (Phase(i) > vF && SignDist(i) > vS)
                isPhase(i) = true;
        }
    }
    int nblobs = ComputeBlob(isPhase, LocalBlobID, periodic, 0);
    PROFILE_STOP("ComputeLocalBlobIDs");
    return nblobs;
}
int ComputeLocalPhaseComponent(const IntArray &PhaseID, int &VALUE,
                               BlobIDArray &ComponentLabel, bool periodic) {
    PROFILE_START("ComputeLocalPhaseComponent");
    size_t Nx = PhaseID.size(0);
    size_t Ny = PhaseID.size(1);
    size_t Nz = PhaseID.size(2);
    size_t N = Nx * Ny * Nz;
    // Compute the local blob ids
    ComponentLabel.resize(Nx, Ny, Nz);
    Array<bool> isPhase(Nx, Ny, Nz);
    for (size_t i = 0; i < N; i++) {
        if (PhaseID(i) == VALUE) {
            ComponentLabel(i) = -1;
            isPhase(i) = true;
        } else {
            ComponentLabel(i) = -2;
            isPhase(i) = false;
        }
    }
    int ncomponents = ComputeBlob(isPhase, ComponentLabel, periodic, 0);
    PROFILE_STOP("ComputeLocalPhaseComponent");
    return ncomponents;
}

/******************************************************************
* Reorder the global blob ids                                     *
******************************************************************/
static int ReorderBlobIDs2(BlobIDArray &ID, int N_blobs, int ngx, int ngy,
                           int ngz, const Utilities::MPI &comm) {
    if (N_blobs == 0)
        return 0;
    PROFILE_START("ReorderBlobIDs2", 1);
    ASSERT(sizeof(long long int) == sizeof(int64_t));
    double *local_size = new double[N_blobs];
    double *global_size = new double[N_blobs];
    for (int i = 0; i < N_blobs; i++)
        local_size[i] = 0;
    for (int i = 0; i < N_blobs; i++)
        global_size[i] = 0;
    int max_id = -1;
    for (size_t k = ngz; k < ID.size(2) - ngz; k++) {
        for (size_t j = ngy; j < ID.size(1) - ngy; j++) {
            for (size_t i = ngx; i < ID.size(0) - ngx; i++) {
                int id = ID(i, j, k);
                if (id >= 0)
                    local_size[id] += 1;
                max_id = std::max(max_id, id);
            }
        }
    }
    ASSERT(max_id < N_blobs);
    comm.sumReduce(local_size, global_size, N_blobs);
    std::vector<std::pair<double, int>> map1(N_blobs);
    int N_blobs2 = 0;
    for (int i = 0; i < N_blobs; i++) {
        map1[i].first = global_size[i];
        map1[i].second = i;
        if (global_size[i] > 0)
            N_blobs2++;
    }
    std::sort(map1.begin(), map1.end());
    std::vector<int> map2(N_blobs, -1);
    for (int i = 0; i < N_blobs; i++) {
        map2[map1[N_blobs - i - 1].second] = i;
    }
    for (size_t i = 0; i < ID.length(); i++) {
        if (ID(i) >= 0)
            ID(i) = map2[ID(i)];
    }
    delete[] local_size;
    delete[] global_size;
    PROFILE_STOP("ReorderBlobIDs2", 1);
    return N_blobs2;
}
void ReorderBlobIDs(BlobIDArray &ID, const Utilities::MPI &comm) {
    PROFILE_START("ReorderBlobIDs");
    int tmp = ID.max() + 1;
    int N_blobs = 0;
    N_blobs = comm.maxReduce(tmp);
    ReorderBlobIDs2(ID, N_blobs, 1, 1, 1, comm);
    PROFILE_STOP("ReorderBlobIDs");
}

/******************************************************************
* Compute the global blob ids                                     *
******************************************************************/
struct global_id_info_struct {
    int64_t new_id;
    std::set<int64_t> remote_ids;
};
// Send the local ids and their new value to all neighbors
static void updateRemoteIds(const std::map<int64_t, global_id_info_struct> &map,
                            const std::vector<int> &neighbors, int N_send,
                            const std::vector<int> &N_recv, int64_t *send_buf,
                            std::vector<int64_t *> &recv_buf,
                            std::map<int64_t, int64_t> &remote_map,
                            const Utilities::MPI &comm) {
    std::vector<MPI_Request> send_req(neighbors.size());
    std::vector<MPI_Request> recv_req(neighbors.size());
    auto it = map.begin();
    ASSERT(N_send == (int)map.size());
    for (size_t i = 0; i < map.size(); i++, ++it) {
        send_buf[2 * i + 0] = it->first;
        send_buf[2 * i + 1] = it->second.new_id;
    }
    for (size_t i = 0; i < neighbors.size(); i++) {
        send_req[i] = comm.Isend(send_buf, 2 * N_send, neighbors[i], 0);
        recv_req[i] = comm.Irecv(recv_buf[i], 2 * N_recv[i], neighbors[i], 0);
    }
    for (it = map.begin(); it != map.end(); ++it) {
        remote_map[it->first] = it->second.new_id;
    }
    for (size_t i = 0; i < neighbors.size(); i++) {
        comm.wait(recv_req[i]);
        for (int j = 0; j < N_recv[i]; j++)
            remote_map[recv_buf[i][2 * j + 0]] = recv_buf[i][2 * j + 1];
    }
    comm.waitAll(neighbors.size(), getPtr(send_req));
}
// Compute a new local id for each local id
static bool updateLocalIds(const std::map<int64_t, int64_t> &remote_map,
                           std::map<int64_t, global_id_info_struct> &map) {
    bool changed = false;
    std::map<int64_t, global_id_info_struct>::iterator it;
    for (it = map.begin(); it != map.end(); ++it) {
        int64_t id = it->second.new_id;
        std::set<int64_t>::const_iterator it2;
        for (it2 = it->second.remote_ids.begin();
             it2 != it->second.remote_ids.end(); ++it2) {
            int64_t id2 = remote_map.find(*it2)->second;
            id = std::min(id, id2);
        }
        changed = changed || it->second.new_id != id;
        it->second.new_id = id;
    }
    return changed;
}
static int LocalToGlobalIDs(int nx, int ny, int nz,
                            const RankInfoStruct &rank_info, int nblobs,
                            BlobIDArray &IDs, const Utilities::MPI &comm) {
    PROFILE_START("LocalToGlobalIDs", 1);
    const int rank = rank_info.rank[1][1][1];
    int nprocs = comm.getSize();
    const int ngx = (IDs.size(0) - nx) / 2;
    const int ngy = (IDs.size(1) - ny) / 2;
    const int ngz = (IDs.size(2) - nz) / 2;
    // Get the number of blobs for each rank
    std::vector<int> N_blobs(nprocs, 0);
    PROFILE_START("LocalToGlobalIDs-Allgather", 1);
    comm.allGather(nblobs, getPtr(N_blobs));
    PROFILE_STOP("LocalToGlobalIDs-Allgather", 1);
    int64_t N_blobs_tot = 0;
    int offset = 0;
    for (int i = 0; i < rank; i++)
        offset += N_blobs[i];
    for (int i = 0; i < nprocs; i++)
        N_blobs_tot += N_blobs[i];
    INSIST(N_blobs_tot < 0x80000000, "Maximum number of blobs exceeded");
    // Compute temporary global ids
    for (size_t i = 0; i < IDs.length(); i++) {
        if (IDs(i) >= 0)
            IDs(i) += offset;
    }
    const BlobIDArray LocalIDs = IDs;
    // Copy the ids and get the neighbors through the halos
    fillHalo<BlobIDType> fillData(comm, rank_info, {nx, ny, nz}, {1, 1, 1}, 0,
                                  1, {true, true, true});
    fillData.fill(IDs);
    // Create a list of all neighbor ranks (excluding self)
    std::vector<int> neighbors;
    neighbors.push_back(rank_info.rank[0][1][1]);
    neighbors.push_back(rank_info.rank[2][1][1]);
    neighbors.push_back(rank_info.rank[1][0][1]);
    neighbors.push_back(rank_info.rank[1][2][1]);
    neighbors.push_back(rank_info.rank[1][1][0]);
    neighbors.push_back(rank_info.rank[1][1][2]);
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                    neighbors.end());
    // Create a map of all local ids to the neighbor ids
    std::map<int64_t, global_id_info_struct> map;
    std::set<int64_t> local;
    for (size_t i = 0; i < LocalIDs.length(); i++) {
        if (LocalIDs(i) >= 0) {
            local.insert(LocalIDs(i));
            if (LocalIDs(i) != IDs(i) && IDs(i) >= 0)
                map[LocalIDs(i)].remote_ids.insert(IDs(i));
        }
    }
    std::map<int64_t, global_id_info_struct>::iterator it;
    for (it = map.begin(); it != map.end(); ++it) {
        it->second.new_id = it->first;
        local.erase(it->first);
    }
    // Get the number of ids we will recieve from each rank
    int N_send = map.size();
    std::vector<int> N_recv(neighbors.size(), 0);
    std::vector<MPI_Request> send_req(neighbors.size());
    std::vector<MPI_Request> recv_req(neighbors.size());
    for (size_t i = 0; i < neighbors.size(); i++) {
        send_req[i] = comm.Isend(&N_send, 1, neighbors[i], 0);
        recv_req[i] = comm.Irecv(&N_recv[i], 1, neighbors[i], 0);
    }
    comm.waitAll(neighbors.size(), getPtr(send_req));
    comm.waitAll(neighbors.size(), getPtr(recv_req));
    // Allocate memory for communication
    int64_t *send_buf = new int64_t[2 * N_send];
    std::vector<int64_t *> recv_buf(neighbors.size());
    for (size_t i = 0; i < neighbors.size(); i++)
        recv_buf[i] = new int64_t[2 * N_recv[i]];
    // Compute a map for the remote ids, and new local id for each id
    std::map<int64_t, int64_t> remote_map;
    for (it = map.begin(); it != map.end(); ++it) {
        int64_t id = it->first;
        std::set<int64_t>::const_iterator it2;
        for (it2 = it->second.remote_ids.begin();
             it2 != it->second.remote_ids.end(); ++it2) {
            int64_t id2 = *it2;
            id = std::min(id, id2);
            remote_map.insert(std::pair<int64_t, int64_t>(id2, id2));
        }
        it->second.new_id = id;
    }
    // Iterate until we are done
    int iteration = 1;
    PROFILE_START("LocalToGlobalIDs-loop", 1);
    while (1) {
        iteration++;
        // Send the local ids and their new value to all neighbors
        updateRemoteIds(map, neighbors, N_send, N_recv, send_buf, recv_buf,
                        remote_map, comm);
        // Compute a new local id for each local id
        bool changed = updateLocalIds(remote_map, map);
        // Check if we are finished
        int test = changed ? 1 : 0;
        int result = comm.sumReduce(test);
        if (result == 0)
            break;
    }
    PROFILE_STOP("LocalToGlobalIDs-loop", 1);
    // Relabel the ids
    std::vector<int> final_map(nblobs, -1);
    for (it = map.begin(); it != map.end(); ++it)
        final_map[it->first - offset] = it->second.new_id;
    for (std::set<int64_t>::const_iterator it2 = local.begin();
         it2 != local.end(); ++it2)
        final_map[*it2 - offset] = *it2;
    for (size_t i = 0; i < final_map.size(); i++)
        ASSERT(final_map[i] >= 0);
    for (size_t k = ngz; k < IDs.size(2) - ngz; k++) {
        for (size_t j = ngy; j < IDs.size(1) - ngy; j++) {
            for (size_t i = ngx; i < IDs.size(0) - ngx; i++) {
                BlobIDType id = IDs(i, j, k);
                if (id >= 0)
                    IDs(i, j, k) = final_map[id - offset];
            }
        }
    }
    // Fill the ghosts
    fillHalo<BlobIDType> fillData2(comm, rank_info, {nx, ny, nz}, {1, 1, 1}, 0,
                                   1, {true, true, true});
    fillData2.fill(IDs);
    // Reorder based on size (and compress the id space
    int N_blobs_global = ReorderBlobIDs2(IDs, N_blobs_tot, ngx, ngy, ngz, comm);
    // Finished
    delete[] send_buf;
    for (size_t i = 0; i < neighbors.size(); i++)
        delete[] recv_buf[i];
    PROFILE_STOP("LocalToGlobalIDs", 1);
    return N_blobs_global;
}
int ComputeGlobalBlobIDs(int nx, int ny, int nz,
                         const RankInfoStruct &rank_info,
                         const DoubleArray &Phase, const DoubleArray &SignDist,
                         double vF, double vS, BlobIDArray &GlobalBlobID,
                         const Utilities::MPI &comm) {
    PROFILE_START("ComputeGlobalBlobIDs");
    // First compute the local ids
    int nblobs =
        ComputeLocalBlobIDs(Phase, SignDist, vF, vS, GlobalBlobID, false);
    // Compute the global ids
    int nglobal =
        LocalToGlobalIDs(nx, ny, nz, rank_info, nblobs, GlobalBlobID, comm);
    PROFILE_STOP("ComputeGlobalBlobIDs");
    return nglobal;
}
int ComputeGlobalPhaseComponent(int nx, int ny, int nz,
                                const RankInfoStruct &rank_info,
                                const IntArray &PhaseID, int &VALUE,
                                BlobIDArray &GlobalBlobID,
                                const Utilities::MPI &comm) {
    PROFILE_START("ComputeGlobalPhaseComponent");
    // First compute the local ids
    int nblobs =
        ComputeLocalPhaseComponent(PhaseID, VALUE, GlobalBlobID, false);
    // Compute the global ids
    int nglobal =
        LocalToGlobalIDs(nx, ny, nz, rank_info, nblobs, GlobalBlobID, comm);
    PROFILE_STOP("ComputeGlobalPhaseComponent");
    return nglobal;
}

/******************************************************************
* Compute the mapping of blob ids between timesteps               *
******************************************************************/
typedef std::map<BlobIDType, std::map<BlobIDType, int64_t>> map_type;
template <class TYPE>
void gatherSet(std::set<TYPE> &set, const Utilities::MPI &comm) {
    int nprocs = comm.getSize();
    std::vector<TYPE> send_data(set.begin(), set.end());
    int send_count = send_data.size();
    std::vector<int> recv_count(nprocs, 0), recv_disp(nprocs, 0);
    comm.allGather(send_count, getPtr(recv_count));
    for (int i = 1; i < nprocs; i++)
        recv_disp[i] = recv_disp[i - 1] + recv_count[i - 1];
    std::vector<TYPE> recv_data(recv_disp[nprocs - 1] + recv_count[nprocs - 1]);
    comm.allGather(getPtr(send_data), send_count, getPtr(recv_data),
                   getPtr(recv_count), getPtr(recv_disp), true);
    for (size_t i = 0; i < recv_data.size(); i++)
        set.insert(recv_data[i]);
}
void gatherSrcIDMap(map_type &src_map, const Utilities::MPI &comm) {
    int nprocs = comm.getSize();
    std::vector<int64_t> send_data;
    for (auto it = src_map.begin(); it != src_map.end(); ++it) {
        int id = it->first;
        const std::map<BlobIDType, int64_t> &src_ids = it->second;
        send_data.push_back(id);
        send_data.push_back(src_ids.size());
        std::map<BlobIDType, int64_t>::const_iterator it2;
        for (it2 = src_ids.begin(); it2 != src_ids.end(); ++it2) {
            send_data.push_back(it2->first);
            send_data.push_back(it2->second);
        }
    }
    int send_count = send_data.size();
    std::vector<int> recv_count(nprocs, 0), recv_disp(nprocs, 0);
    comm.allGather(send_count, getPtr(recv_count));
    for (int i = 1; i < nprocs; i++)
        recv_disp[i] = recv_disp[i - 1] + recv_count[i - 1];
    std::vector<int64_t> recv_data(recv_disp[nprocs - 1] +
                                   recv_count[nprocs - 1]);
    comm.allGather(getPtr(send_data), send_count, getPtr(recv_data),
                   getPtr(recv_count), getPtr(recv_disp), true);
    size_t i = 0;
    src_map.clear();
    while (i < recv_data.size()) {
        BlobIDType id = recv_data[i];
        size_t count = recv_data[i + 1];
        i += 2;
        auto &src_ids = src_map[id];
        for (size_t j = 0; j < count; j++, i += 2) {
            auto it = src_ids.find(recv_data[i]);
            if (it == src_ids.end())
                src_ids.insert(std::pair<BlobIDType, int64_t>(
                    recv_data[i], recv_data[i + 1]));
            else
                it->second += recv_data[i + 1];
        }
    }
}
void addSrcDstIDs(BlobIDType src_id, map_type &src_map, map_type &dst_map,
                  std::set<BlobIDType> &src, std::set<BlobIDType> &dst) {
    src.insert(src_id);
    const std::map<BlobIDType, int64_t> &dst_ids = dst_map[src_id];
    for (std::map<BlobIDType, int64_t>::const_iterator it = dst_ids.begin();
         it != dst_ids.end(); ++it) {
        if (dst.find(it->first) == dst.end())
            addSrcDstIDs(it->first, dst_map, src_map, dst, src);
    }
}
ID_map_struct computeIDMap(int nx, int ny, int nz, const BlobIDArray &ID1,
                           const BlobIDArray &ID2, const Utilities::MPI &comm) {
    ASSERT(ID1.size() == ID2.size());
    PROFILE_START("computeIDMap");
    const int ngx = (ID1.size(0) - nx) / 2;
    const int ngy = (ID1.size(1) - ny) / 2;
    const int ngz = (ID1.size(2) - nz) / 2;

    // Get a global list of all src/dst ids and the src map for each local blob
    std::set<BlobIDType> src_set, dst_set;
    map_type src_map; // Map of the src ids for each dst id
    for (int k = ngz; k < ngz + nz; k++) {
        for (int j = ngy; j < ngy + ny; j++) {
            for (int i = ngx; i < ngx + nx; i++) {
                int id1 = ID1(i, j, k);
                int id2 = ID2(i, j, k);
                if (id1 >= 0)
                    src_set.insert(id1);
                if (id2 >= 0)
                    dst_set.insert(id2);
                if (id1 >= 0 && id2 >= 0) {
                    std::map<BlobIDType, int64_t> &src_ids = src_map[id2];
                    std::map<BlobIDType, int64_t>::iterator it =
                        src_ids.find(id1);
                    if (it == src_ids.end()) {
                        src_ids.insert(std::pair<BlobIDType, int64_t>(id1, 0));
                        it = src_ids.find(id1);
                    }
                    it->second++;
                }
            }
        }
    }
    // Communicate the src/dst ids and src id map to all processors and reduce
    gatherSet(src_set, comm);
    gatherSet(dst_set, comm);
    gatherSrcIDMap(src_map, comm);
    // Compute the dst id map
    map_type dst_map; // Map of the dst ids for each src id
    for (map_type::const_iterator it = src_map.begin(); it != src_map.end();
         ++it) {
        BlobIDType id = it->first;
        const std::map<BlobIDType, int64_t> &src_ids = it->second;
        for (std::map<BlobIDType, int64_t>::const_iterator it2 =
                 src_ids.begin();
             it2 != src_ids.end(); ++it2) {
            std::map<BlobIDType, int64_t> &dst_ids = dst_map[it2->first];
            dst_ids.insert(std::pair<BlobIDType, int64_t>(id, it2->second));
        }
    }

    // Perform the mapping of ids
    ID_map_struct id_map;
    // Find new blobs
    for (std::set<BlobIDType>::const_iterator it = dst_set.begin();
         it != dst_set.end(); ++it) {
        if (src_map.find(*it) == src_map.end())
            id_map.created.push_back(*it);
    }
    // Fine blobs that disappeared
    for (std::set<BlobIDType>::const_iterator it = src_set.begin();
         it != src_set.end(); ++it) {
        if (dst_map.find(*it) == dst_map.end())
            id_map.destroyed.push_back(*it);
    }
    // Find blobs with a 1-to-1 mapping
    std::vector<BlobIDType> dst_list;
    dst_list.reserve(src_map.size());
    for (map_type::const_iterator it = src_map.begin(); it != src_map.end();
         ++it)
        dst_list.push_back(it->first);
    for (size_t i = 0; i < dst_list.size(); i++) {
        int dst_id = dst_list[i];
        const std::map<BlobIDType, int64_t> &src_ids = src_map[dst_id];
        if (src_ids.size() == 1) {
            int src_id = src_ids.begin()->first;
            const std::map<BlobIDType, int64_t> &dst_ids = dst_map[src_id];
            if (dst_ids.size() == 1) {
                ASSERT(dst_ids.begin()->first == dst_id);
                src_map.erase(dst_id);
                dst_map.erase(src_id);
                id_map.src_dst.push_back(
                    std::pair<BlobIDType, BlobIDType>(src_id, dst_id));
            }
        }
    }
    // Handle merge/splits
    while (!dst_map.empty()) {
        // Get a lit of the src-dst ids
        std::set<BlobIDType> src, dst;
        addSrcDstIDs(dst_map.begin()->first, src_map, dst_map, src, dst);
        if (src.size() == 1) {
            // Bubble split
            id_map.split.push_back(BlobIDSplitStruct(
                *src.begin(), std::vector<BlobIDType>(dst.begin(), dst.end())));
        } else if (dst.size() == 1) {
            // Bubble merge
            id_map.merge.push_back(BlobIDMergeStruct(
                std::vector<BlobIDType>(src.begin(), src.end()), *dst.begin()));
        } else {
            // Bubble split/merge
            id_map.merge_split.push_back(BlobIDMergeSplitStruct(
                std::vector<BlobIDType>(src.begin(), src.end()),
                std::vector<BlobIDType>(dst.begin(), dst.end())));
        }
        // Add the overlaps
        for (std::set<BlobIDType>::const_iterator it1 = src.begin();
             it1 != src.end(); ++it1) {
            const std::map<BlobIDType, int64_t> &dst_ids = dst_map[*it1];
            for (std::set<BlobIDType>::const_iterator it2 = dst.begin();
                 it2 != dst.end(); ++it2) {
                std::pair<BlobIDType, BlobIDType> id(*it1, *it2);
                int64_t overlap = 0;
                const std::map<BlobIDType, int64_t>::const_iterator it =
                    dst_ids.find(*it2);
                if (it != dst_ids.end()) {
                    overlap = it->second;
                }
                id_map.overlap.insert(
                    std::pair<OverlapID, int64_t>(id, overlap));
            }
        }
        // Clear the mapped entries
        for (std::set<BlobIDType>::const_iterator it = src.begin();
             it != src.end(); ++it)
            dst_map.erase(*it);
        for (std::set<BlobIDType>::const_iterator it = dst.begin();
             it != dst.end(); ++it)
            src_map.erase(*it);
    }
    ASSERT(src_map.empty());

    PROFILE_STOP("computeIDMap");
    return id_map;
}

/******************************************************************
* Renumber the ids                                                *
******************************************************************/
typedef std::vector<BlobIDType> IDvec;
inline void renumber(const std::vector<BlobIDType> &id1,
                     const std::vector<BlobIDType> &id2,
                     const std::map<OverlapID, int64_t> &overlap,
                     std::vector<BlobIDType> &new_ids, BlobIDType &id_max) {
    if (id2.empty()) {
        // No dst ids to set
    } else if (id1.empty()) {
        // No src ids
        for (size_t i = 0; i < id2.size(); i++) {
            id_max++;
            if ((BlobIDType)new_ids.size() < id2[i] + 1)
                new_ids.resize(id2[i] + 1, -1);
            new_ids[id2[i]] = id_max;
        }
    } else if (id1.size() == 1 && id2.size() == 1) {
        // Direct src-dst mapping
        if ((BlobIDType)new_ids.size() < id2[0] + 1)
            new_ids.resize(id2[0] + 1, -1);
        new_ids[id2[0]] = id1[0];
    } else {
        // General N to M mapping
        // Get the overlap weights
        Array<int64_t> cost(id1.size(), id2.size());
        for (size_t j = 0; j < id2.size(); j++) {
            for (size_t i = 0; i < id1.size(); i++) {
                cost(i, j) =
                    overlap
                        .find(std::pair<BlobIDType, BlobIDType>(id1[i], id2[j]))
                        ->second;
            }
        }
        // While we have not mapped all dst ids
        while (1) {
            size_t index = 1;
            int64_t cost2 = -1;
            for (size_t i = 0; i < cost.length(); i++) {
                if (cost(i) > cost2) {
                    cost2 = cost(i);
                    index = i;
                }
            }
            if (cost2 <= 0)
                break;
            // Use id1[i] for id2[j]
            int i = index % id1.size();
            int j = index / id1.size();
            if ((BlobIDType)new_ids.size() < id2[j] + 1)
                new_ids.resize(id2[j] + 1, -1);
            new_ids[id2[j]] = id1[i];
            for (size_t k = 0; k < id2.size(); k++)
                cost(i, k) = -1;
            for (size_t k = 0; k < id1.size(); k++)
                cost(k, j) = -1;
        }
        // No remaining src overlap with dst, create new ids for all remaining dst
        for (size_t i = 0; i < id2.size(); i++) {
            if ((BlobIDType)new_ids.size() < id2[i] + 1)
                new_ids.resize(id2[i] + 1, -1);
            if (new_ids[id2[i]] == -1) {
                id_max++;
                new_ids[id2[i]] = id_max;
            }
        }
    }
}
inline void renumberIDs(const std::vector<BlobIDType> &new_ids,
                        BlobIDType &id) {
    id = new_ids[id];
}
inline void renumberIDs(const std::vector<BlobIDType> &new_ids,
                        std::vector<BlobIDType> &ids) {
    for (size_t i = 0; i < ids.size(); i++)
        ids[i] = new_ids[ids[i]];
}
void getNewIDs(ID_map_struct &map, BlobIDType &id_max,
               std::vector<BlobIDType> &new_ids) {
    new_ids.resize(0);
    // Get the new id numbers for each map type
    for (size_t i = 0; i < map.src_dst.size(); i++)
        renumber(IDvec(1, map.src_dst[i].first),
                 IDvec(1, map.src_dst[i].second), map.overlap, new_ids, id_max);
    for (size_t i = 0; i < map.created.size(); i++)
        renumber(std::vector<BlobIDType>(), IDvec(1, map.created[i]),
                 map.overlap, new_ids, id_max);
    for (size_t i = 0; i < map.destroyed.size(); i++)
        renumber(IDvec(1, map.destroyed[i]), std::vector<BlobIDType>(),
                 map.overlap, new_ids, id_max);
    for (size_t i = 0; i < map.split.size(); i++)
        renumber(IDvec(1, map.split[i].first), map.split[i].second, map.overlap,
                 new_ids, id_max);
    for (size_t i = 0; i < map.merge.size(); i++)
        renumber(map.merge[i].first, IDvec(1, map.merge[i].second), map.overlap,
                 new_ids, id_max);
    for (size_t i = 0; i < map.merge_split.size(); i++)
        renumber(map.merge_split[i].first, map.merge_split[i].second,
                 map.overlap, new_ids, id_max);
    // Renumber the ids in the map
    for (size_t i = 0; i < map.src_dst.size(); i++)
        renumberIDs(new_ids, map.src_dst[i].second);
    renumberIDs(new_ids, map.created);
    for (size_t i = 0; i < map.split.size(); i++)
        renumberIDs(new_ids, map.split[i].second);
    for (size_t i = 0; i < map.merge.size(); i++)
        renumberIDs(new_ids, map.merge[i].second);
    for (size_t i = 0; i < map.merge_split.size(); i++)
        renumberIDs(new_ids, map.merge_split[i].second);
    std::map<OverlapID, int64_t> overlap2;
    for (std::map<OverlapID, int64_t>::const_iterator it = map.overlap.begin();
         it != map.overlap.begin(); ++it) {
        OverlapID id = it->first;
        renumberIDs(new_ids, id.second);
        overlap2.insert(std::pair<OverlapID, int64_t>(id, it->second));
    }
}
void renumberIDs(const std::vector<BlobIDType> &new_ids, BlobIDArray &IDs) {
    size_t N = IDs.length();
    BlobIDType *ids = IDs.data();
    for (size_t i = 0; i < N; i++) {
        BlobIDType id = ids[i];
        if (id >= 0)
            ids[i] = new_ids[id];
    }
}

/******************************************************************
* Write the id map for the given timestep                         *
******************************************************************/
void writeIDMap(const ID_map_struct &map, long long int timestep,
                const std::string &filename) {
    int rank = Utilities::MPI(MPI_COMM_WORLD).getRank();
    if (rank != 0)
        return;
    bool empty = map.created.empty() && map.destroyed.empty() &&
                 map.split.empty() && map.merge.empty() &&
                 map.merge_split.empty();
    for (size_t i = 0; i < map.src_dst.size(); i++)
        empty = empty && map.src_dst[i].first == map.src_dst[i].second;
    if (timestep != 0 && empty)
        return;
    FILE *fid = NULL;
    if (timestep == 0)
        fid = fopen(filename.c_str(), "wb");
    else
        fid = fopen(filename.c_str(), "ab");
    INSIST(fid != NULL, std::string("Error opening file: ") + filename);
    if (empty) {
        fclose(fid);
        return;
    }
    fprintf(fid, "%lli:", timestep);
    for (size_t i = 0; i < map.created.size(); i++)
        fprintf(fid, " -%lli", static_cast<long long int>(map.created[i]));
    for (size_t i = 0; i < map.destroyed.size(); i++)
        fprintf(fid, " %lli-", static_cast<long long int>(map.destroyed[i]));
    for (size_t i = 0; i < map.src_dst.size(); i++) {
        if (map.src_dst[i].first != map.src_dst[i].second)
            fprintf(fid, " %lli-%lli",
                    static_cast<long long int>(map.src_dst[i].first),
                    static_cast<long long int>(map.src_dst[i].second));
    }
    for (size_t i = 0; i < map.split.size(); i++) {
        fprintf(fid, " %lli-%lli",
                static_cast<long long int>(map.split[i].first),
                static_cast<long long int>(map.split[i].second[0]));
        for (size_t j = 1; j < map.split[i].second.size(); j++)
            fprintf(fid, "/%lli",
                    static_cast<long long int>(map.split[i].second[j]));
    }
    for (size_t i = 0; i < map.merge.size(); i++) {
        fprintf(fid, " %lli",
                static_cast<long long int>(map.merge[i].first[0]));
        for (size_t j = 1; j < map.merge[i].first.size(); j++)
            fprintf(fid, "/%lli",
                    static_cast<long long int>(map.merge[i].first[j]));
        fprintf(fid, "-%lli", static_cast<long long int>(map.merge[i].second));
    }
    for (size_t i = 0; i < map.merge_split.size(); i++) {
        fprintf(fid, " %lli",
                static_cast<long long int>(map.merge_split[i].first[0]));
        for (size_t j = 1; j < map.merge_split[i].first.size(); j++)
            fprintf(fid, "/%lli",
                    static_cast<long long int>(map.merge_split[i].first[j]));
        fprintf(fid, "-%lli",
                static_cast<long long int>(map.merge_split[i].second[0]));
        for (size_t j = 1; j < map.merge_split[i].second.size(); j++)
            fprintf(fid, "/%lli",
                    static_cast<long long int>(map.merge_split[i].second[j]));
    }
    fprintf(fid, "\n");
    fclose(fid);
}
