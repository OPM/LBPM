#include "analysis/analysis.h"
#include "ProfilerApp.h"
#include <iostream>


template<class TYPE>
inline TYPE* getPtr( std::vector<TYPE>& x ) { return x.empty() ? NULL:&x[0]; }
template<class TYPE>
inline const TYPE* getPtr( const std::vector<TYPE>& x ) { return x.empty() ? NULL:&x[0]; }


/******************************************************************
* Compute the blobs                                               *
******************************************************************/
int ComputeBlob( const Array<bool>& isPhase, IntArray& LocalBlobID, bool periodic, int start_id )
{
    PROFILE_START("ComputeBlob",1);
    ASSERT(isPhase.size(0)==LocalBlobID.size(0));
    ASSERT(isPhase.size(1)==LocalBlobID.size(1));
    ASSERT(isPhase.size(2)==LocalBlobID.size(2));
    const int Nx = isPhase.size(0);  // maxima for the meshes
    const int Ny = isPhase.size(1);
    const int Nz = isPhase.size(2);
    std::vector<int> map;
    map.reserve(128);
    // Get the list of neighbors we need to check
    int N_neighbors = 0;
    int d[26][3];
    bool include_corners = false;   // Do we need to include cells that only touch at a corder/edge
    if ( include_corners ) {
        // Include corners/edges as neighbors, check all cells
        N_neighbors = 26;
        const int tmp[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
                                {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
                                {1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
                                {1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
                                {-1,-1,1},{-1,-1,-1}};   // directions to neighbors
        memcpy(d,tmp,sizeof(tmp));
    } else {
        // Do not include corners/edges as neighbors
        if ( periodic ) {
            // Include all neighbors for periodic problems
            N_neighbors = 6;
            const int tmp[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};   // directions to neighbors
            memcpy(d,tmp,sizeof(tmp));
        } else {
            // We only need to include the lower neighbors for non-periodic problems
            N_neighbors = 3;
            const int tmp[3][3] = {{-1,0,0},{0,-1,0},{0,0,-1}};   // directions to neighbors
            memcpy(d,tmp,sizeof(tmp));
        }
    }
    // Loop through all the points
    int last = start_id-1;
    std::vector<int> neighbor_ids;
    neighbor_ids.reserve(N_neighbors);
    const bool *isPhasePtr = isPhase.get();
    int *LocalBlobIDPtr = LocalBlobID.get();
    for (int z=0; z<Nz; z++) {
        for (int y=0; y<Ny; y++) {
            for (int x=0; x<Nx; x++) {
                int index = x + y*Nx + z*Nx*Ny;
                if ( !isPhasePtr[index] )
                    continue;
                // Get all neighbor indicies
                int N_list=0, neighbor_ids[26];
                for (int p=0; p<N_neighbors; p++) {
                    // Get the neighbor index
                    int x2 = x+d[p][0];
                    int y2 = y+d[p][1];
                    int z2 = z+d[p][2];
                    if ( periodic ) {
                        x2 = x2<0 ? Nx-1:x2;         // Periodic BC for x
                        x2 = x2>Nx-1 ? 0:x2;
                        y2 = y2<0 ? Ny-1:y2;         // Periodic BC for x
                        y2 = y2>Ny-1 ? 0:y2;
                        z2 = z2<0 ? Nz-1:z2;         // Periodic BC for x
                        z2 = z2>Nz-1 ? 0:z2;
                    } else {
                        if ( x2<0 || x2>=Nx || y2<0 || y2>=Ny || z2<0 || z2>=Nz )
                            continue;
                    }
                    // Check if a neighbor has a known blob id
                    size_t index2 = x2 + y2*Nx + z2*Nx*Ny;
                    int id = LocalBlobIDPtr[index2];
                    if ( !isPhasePtr[index2] || id<0 )
                        continue;
                    neighbor_ids[N_list] = id;
                    N_list++;
                }
                if ( N_list==0 ) {
                    // No neighbors with a blob id, create a new one
                    LocalBlobIDPtr[index] = last+1;
                    map.push_back(last+1);
                    last++;
                } else if ( N_list==1 ) {
                    // We have one neighbor
                    LocalBlobIDPtr[index] = neighbor_ids[0];
                } else {
                    // We have multiple neighbors
                    int id = neighbor_ids[0];
                    for (int i=1; i<N_list; i++)
                        id = std::min(id,neighbor_ids[i]);
                    LocalBlobIDPtr[index] = id;
                    for (int i=0; i<N_list; i++)
                        map[neighbor_ids[i]-start_id] = std::min(map[neighbor_ids[i]-start_id],id);
                }
            }
        }
    }
    // Collapse the ids that map to another id
    last = start_id-1;
    for (int i=0; i<(int)map.size(); i++) {
        if ( map[i] == i+start_id ) {
            map[i] = last+1;
            last++;
        } else {
            ASSERT(map[i]<i+start_id);
            map[i] = map[map[i]-start_id];
        }
    }
    for (int i=0; i<Nx*Ny*Nz; i++) {
        if ( isPhasePtr[i] ) {
            LocalBlobIDPtr[i] = map[LocalBlobIDPtr[i]-start_id];
        }
    }
    PROFILE_STOP("ComputeBlob",1);
    return last-start_id+1;
}


/******************************************************************
* Compute the local blob ids                                      *
******************************************************************/
int ComputeLocalBlobIDs( const DoubleArray& Phase, const DoubleArray& SignDist, 
    double vF, double vS, IntArray& LocalBlobID, bool periodic )
{
    PROFILE_START("ComputeLocalBlobIDs");
    ASSERT(SignDist.size(0)==Phase.size(0));
    ASSERT(SignDist.size(1)==Phase.size(1));
    ASSERT(SignDist.size(2)==Phase.size(2));
    size_t Nx = Phase.size(0);
    size_t Ny = Phase.size(1);
    size_t Nz = Phase.size(2);
    // Initialize output
    LocalBlobID.resize(Nx,Ny,Nz);
    // Compute the local blob ids
    size_t N = Nx*Ny*Nz;
    Array<bool> isPhase(Nx,Ny,Nz);
    memset(isPhase.get(),0,Nx*Ny*Nz*sizeof(bool));
    for (size_t i=0; i<N; i++) {
        if ( SignDist(i) < 0.0) {
            // Solid phase 
            LocalBlobID(i) = -2;
        } else {
            LocalBlobID(i) = -1;
            if ( Phase(i)>vF && SignDist(i)>vS )
                isPhase(i) = true;
        }
    }
    int nblobs = ComputeBlob( isPhase, LocalBlobID, periodic, 0 );
    PROFILE_STOP("ComputeLocalBlobIDs");
    return nblobs;
}
int ComputeLocalPhaseComponent(const IntArray &PhaseID, int VALUE, IntArray &ComponentLabel, bool periodic )
{
    PROFILE_START("ComputeLocalPhaseComponent");
    size_t Nx = PhaseID.size(0);
    size_t Ny = PhaseID.size(1);
    size_t Nz = PhaseID.size(2);
    size_t N = Nx*Ny*Nz;
    // Compute the local blob ids
	ComponentLabel.resize(Nx,Ny,Nz);
    Array<bool> isPhase(Nx,Ny,Nz);
    for (size_t i=0; i<N; i++) {
        if ( PhaseID(i) == VALUE) {
            ComponentLabel(i) = -1;
            isPhase(i) = true;
        } else{
            ComponentLabel(i) = -2;
            isPhase(i) = false;
        }
    }
    int ncomponents = ComputeBlob( isPhase, ComponentLabel, periodic, 0 );
    PROFILE_STOP("ComputeLocalPhaseComponent");
    return ncomponents;
}


/******************************************************************
* Reorder the global blob ids                                     *
******************************************************************/
static int ReorderBlobIDs2( IntArray& ID, int N_blobs, int ngx, int ngy, int ngz )
{
    if ( N_blobs==0 )
        return 0;
    PROFILE_START("ReorderBlobIDs2",1);
    ASSERT(sizeof(long long int)==sizeof(int64_t));
    double *local_size = new double[N_blobs];
    double *global_size = new double[N_blobs];
    for (int i=0; i<N_blobs; i++)
        local_size[i] = 0;
    for (int i=0; i<N_blobs; i++)
        global_size[i] = 0;
    int max_id = -1;
    for (size_t k=ngz; k<ID.size(2)-ngz; k++) {
        for (size_t j=ngy; j<ID.size(1)-ngy; j++) {
            for (size_t i=ngx; i<ID.size(0)-ngx; i++) {
                int id = ID(i,j,k);
                if ( id >= 0 ) 
                    local_size[id] += 1;
                max_id = std::max(max_id,id);
            }
        }
    }
    ASSERT(max_id<N_blobs);
    MPI_Allreduce(local_size,global_size,N_blobs,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    std::vector<std::pair<double,int> > map1(N_blobs);
    int N_blobs2 = 0;
    for (int i=0; i<N_blobs; i++) {
        map1[i].first = global_size[i];
        map1[i].second = i;
        if ( global_size[i] > 0 )
            N_blobs2++;
    }
    std::sort( map1.begin(), map1.end() );
    std::vector<int> map2(N_blobs,-1);
    for (int i=0; i<N_blobs; i++) {
        map2[map1[N_blobs-i-1].second] = i;
    }
    for (size_t i=0; i<ID.length(); i++) {
        if ( ID(i) >= 0 )
            ID(i) = map2[ID(i)];
    }
    delete [] local_size;
    delete [] global_size;
    PROFILE_STOP("ReorderBlobIDs2",1);
    return N_blobs2;
}
void ReorderBlobIDs( IntArray& ID )
{
    PROFILE_START("ReorderBlobIDs");
    int tmp = ID.max()+1;
    int N_blobs = 0;
    MPI_Allreduce(&tmp,&N_blobs,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    ReorderBlobIDs2(ID,N_blobs,1,1,1);
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
static void updateRemoteIds( 
    const std::map<int64_t,global_id_info_struct>& map,
    const std::vector<int>& neighbors,
    int N_send, const std::vector<int>& N_recv,
    int64_t *send_buf, std::vector<int64_t*>& recv_buf,
    std::map<int64_t,int64_t>& remote_map )
{
    std::vector<MPI_Request> send_req(neighbors.size());
    std::vector<MPI_Request> recv_req(neighbors.size());
    std::vector<MPI_Status> status(neighbors.size());
    std::map<int64_t,global_id_info_struct>::const_iterator it = map.begin();
    ASSERT(N_send==(int)map.size());
    for (size_t i=0; i<map.size(); i++, ++it) {
        send_buf[2*i+0] = it->first;
        send_buf[2*i+1] = it->second.new_id;
    }
    for (size_t i=0; i<neighbors.size(); i++) {
        MPI_Isend( send_buf,    2*N_send,    MPI_LONG_LONG, neighbors[i], 0, MPI_COMM_WORLD, &send_req[i] );
        MPI_Irecv( recv_buf[i], 2*N_recv[i], MPI_LONG_LONG, neighbors[i], 0, MPI_COMM_WORLD, &recv_req[i] );
    }
    for (it=map.begin(); it!=map.end(); ++it) {
        remote_map[it->first] = it->second.new_id;
    }
    for (size_t i=0; i<neighbors.size(); i++) {
        MPI_Wait(&recv_req[i],&status[i]);
        for (int j=0; j<N_recv[i]; j++)
            remote_map[recv_buf[i][2*j+0]] = recv_buf[i][2*j+1];
    }
    MPI_Waitall(neighbors.size(),getPtr(send_req),getPtr(status));
}
// Compute a new local id for each local id
static bool updateLocalIds( const std::map<int64_t,int64_t>& remote_map, 
    std::map<int64_t,global_id_info_struct>& map )
{
    bool changed = false;
    std::map<int64_t,global_id_info_struct>::iterator it;
    for (it=map.begin(); it!=map.end(); ++it) {
        int64_t id = it->second.new_id;
        std::set<int64_t>::const_iterator it2;
        for (it2=it->second.remote_ids.begin(); it2!=it->second.remote_ids.end(); ++it2) {
            int64_t id2 = remote_map.find(*it2)->second;
            id = std::min(id,id2);
        }
        changed = changed || it->second.new_id!=id;
        it->second.new_id = id;
    }
    return changed;
}
static int LocalToGlobalIDs( int nx, int ny, int nz, const RankInfoStruct& rank_info, 
    int nblobs, IntArray& IDs )
{
    PROFILE_START("LocalToGlobalIDs",1);
    const int rank = rank_info.rank[1][1][1];
    int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    const int ngx = (IDs.size(0)-nx)/2;
    const int ngy = (IDs.size(1)-ny)/2;
    const int ngz = (IDs.size(2)-nz)/2;
    // Get the number of blobs for each rank
    std::vector<int> N_blobs(nprocs,0);
    PROFILE_START("LocalToGlobalIDs-Allgather",1);
    MPI_Allgather(&nblobs,1,MPI_INT,getPtr(N_blobs),1,MPI_INT,MPI_COMM_WORLD);
    PROFILE_STOP("LocalToGlobalIDs-Allgather",1);
    int64_t N_blobs_tot = 0;
    int offset = 0;
    for (int i=0; i<rank; i++)
        offset += N_blobs[i];
    for (int i=0; i<nprocs; i++)
        N_blobs_tot += N_blobs[i];
    INSIST(N_blobs_tot<0x80000000,"Maximum number of blobs exceeded");
    // Compute temporary global ids
    for (size_t i=0; i<IDs.length(); i++) {
        if ( IDs(i) >= 0 )
            IDs(i) += offset;
    }
    const Array<int> LocalIDs = IDs;
    // Copy the ids and get the neighbors through the halos
    fillHalo<int> fillData(rank_info,nx,ny,nz,1,1,1,0,1,true,true,true);
    fillData.fill(IDs);
    // Create a list of all neighbor ranks (excluding self)
    std::vector<int> neighbors;
    neighbors.push_back( rank_info.rank[0][1][1] );
    neighbors.push_back( rank_info.rank[2][1][1] );
    neighbors.push_back( rank_info.rank[1][0][1] );
    neighbors.push_back( rank_info.rank[1][2][1] );
    neighbors.push_back( rank_info.rank[1][1][0] );
    neighbors.push_back( rank_info.rank[1][1][2] );
    std::unique(neighbors.begin(),neighbors.end());
    if ( std::find(neighbors.begin(),neighbors.end(),rank) != neighbors.end() )
        neighbors.erase(std::find(neighbors.begin(),neighbors.end(),rank));
    // Create a map of all local ids to the neighbor ids
    std::map<int64_t,global_id_info_struct> map;
    std::set<int64_t> local;
    for (size_t i=0; i<LocalIDs.length(); i++) {
        if ( LocalIDs(i)>=0 ) {
            local.insert(LocalIDs(i));
            if ( LocalIDs(i)!=IDs(i) && IDs(i)>=0 )
                map[LocalIDs(i)].remote_ids.insert(IDs(i));
        }
    }
    std::map<int64_t,global_id_info_struct>::iterator it;
    for (it=map.begin(); it!=map.end(); ++it) {
        it->second.new_id = it->first;
        local.erase(it->first);
    }
    // Get the number of ids we will recieve from each rank
    int N_send = map.size();
    std::vector<int> N_recv(neighbors.size(),0);
    std::vector<MPI_Request> send_req(neighbors.size());
    std::vector<MPI_Request> recv_req(neighbors.size());
    std::vector<MPI_Status> status(neighbors.size());
    for (size_t i=0; i<neighbors.size(); i++) {
        MPI_Isend( &N_send,    1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &send_req[i] );
        MPI_Irecv( &N_recv[i], 1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &recv_req[i] );
    }
    MPI_Waitall(neighbors.size(),getPtr(send_req),getPtr(status));
    MPI_Waitall(neighbors.size(),getPtr(recv_req),getPtr(status));
    // Allocate memory for communication
    int64_t *send_buf = new int64_t[2*N_send];
    std::vector<int64_t*> recv_buf(neighbors.size());
    for (size_t i=0; i<neighbors.size(); i++)
        recv_buf[i] = new int64_t[2*N_recv[i]];
    // Compute a map for the remote ids, and new local id for each id
    std::map<int64_t,int64_t> remote_map;
    for (it=map.begin(); it!=map.end(); ++it) {
        int64_t id = it->first;
        std::set<int64_t>::const_iterator it2;
        for (it2=it->second.remote_ids.begin(); it2!=it->second.remote_ids.end(); ++it2) {
            int64_t id2 = *it2;
            id = std::min(id,id2);
            remote_map.insert(std::pair<int64_t,int64_t>(id2,id2));
        }
        it->second.new_id = id;
    }
    // Iterate until we are done
    int iteration = 1;
    PROFILE_START("LocalToGlobalIDs-loop",1);
    while ( 1 ) {
        iteration++;
        // Send the local ids and their new value to all neighbors
        updateRemoteIds( map, neighbors, N_send, N_recv,send_buf, recv_buf, remote_map );
        // Compute a new local id for each local id
        bool changed = updateLocalIds( remote_map, map );
        // Check if we are finished
        int test = changed ? 1:0;
        int result = 0;
        MPI_Allreduce(&test,&result,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        if ( result==0 )
            break;
    }
    PROFILE_STOP("LocalToGlobalIDs-loop",1);
    // Relabel the ids
    std::vector<int> final_map(nblobs,-1);
    for (it=map.begin(); it!=map.end(); ++it)
        final_map[it->first-offset] = it->second.new_id;
    for (std::set<int64_t>::const_iterator it2=local.begin(); it2!=local.end(); ++it2)
        final_map[*it2-offset] = *it2;
    for (size_t i=0; i<final_map.size(); i++)
        ASSERT(final_map[i]>=0);
    for (size_t k=ngz; k<IDs.size(2)-ngz; k++) {
        for (size_t j=ngy; j<IDs.size(1)-ngy; j++) {
            for (size_t i=ngx; i<IDs.size(0)-ngx; i++) {
                int id = IDs(i,j,k);
                if ( id >= 0 )
                    IDs(i,j,k) = final_map[id-offset];
            }
        }
    }
    // Fill the ghosts
    fillHalo<int> fillData2(rank_info,nx,ny,nz,1,1,1,0,1,true,true,true);
    fillData2.fill(IDs);
    // Reorder based on size (and compress the id space
    int N_blobs_global = ReorderBlobIDs2(IDs,N_blobs_tot,ngx,ngy,ngz);
    // Finished
    delete [] send_buf;
    for (size_t i=0; i<neighbors.size(); i++)
        delete [] recv_buf[i];
    PROFILE_STOP("LocalToGlobalIDs",1);
    return N_blobs_global;
}
int ComputeGlobalBlobIDs( int nx, int ny, int nz, const RankInfoStruct& rank_info, 
    const DoubleArray& Phase, const DoubleArray& SignDist, double vF, double vS,
    IntArray& GlobalBlobID )
{
    PROFILE_START("ComputeGlobalBlobIDs");
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    // First compute the local ids
    int nblobs = ComputeLocalBlobIDs(Phase,SignDist,vF,vS,GlobalBlobID,false);
    // Compute the global ids
    int nglobal = LocalToGlobalIDs( nx, ny, nz, rank_info, nblobs, GlobalBlobID );
    PROFILE_STOP("ComputeGlobalBlobIDs");
    return nglobal;
}
int ComputeGlobalPhaseComponent( int nx, int ny, int nz, const RankInfoStruct& rank_info,
    const IntArray &PhaseID, int VALUE, IntArray &GlobalBlobID )
{
    PROFILE_START("ComputeGlobalPhaseComponent");
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    // First compute the local ids
    int nblobs = ComputeLocalPhaseComponent(PhaseID,VALUE,GlobalBlobID,false);
    // Compute the global ids
    int nglobal = LocalToGlobalIDs( nx, ny, nz, rank_info, nblobs, GlobalBlobID );
    PROFILE_STOP("ComputeGlobalPhaseComponent");
    return nglobal;
}


/******************************************************************
* Compute the mapping of blob ids between timesteps               *
******************************************************************/
void gatherSet( std::set<int>& set )
{
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    std::vector<int> send_data(set.begin(),set.end());
    int send_count = send_data.size();
    std::vector<int> recv_count(nprocs,0), recv_disp(nprocs,0);
    MPI_Allgather(&send_count,1,MPI_INT,getPtr(recv_count),1,MPI_INT,MPI_COMM_WORLD);
    for (int i=1; i<nprocs; i++)
        recv_disp[i] = recv_disp[i-1] + recv_count[i-1];
    std::vector<int> recv_data(recv_disp[nprocs-1]+recv_count[nprocs-1]);
    MPI_Allgatherv(getPtr(send_data),send_count,MPI_INT,
        getPtr(recv_data),getPtr(recv_count),getPtr(recv_disp),MPI_INT,
        MPI_COMM_WORLD);
    for (size_t i=0; i<recv_data.size(); i++)
        set.insert(recv_data[i]);
}
void gatherSrcIDMap( std::map<int,std::set<int> >& src_map )
{
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    std::vector<int> send_data;
    for (std::map<int,std::set<int> >::const_iterator it=src_map.begin(); it!=src_map.end(); ++it) {
        int id = it->first;
        const std::set<int>& src_ids = it->second;
        send_data.push_back(id);
        send_data.push_back(src_ids.size());
        for (std::set<int>::const_iterator it2=src_ids.begin(); it2!=src_ids.end(); ++it2)
            send_data.push_back(*it2);
    }
    int send_count = send_data.size();
    std::vector<int> recv_count(nprocs,0), recv_disp(nprocs,0);
    MPI_Allgather(&send_count,1,MPI_INT,getPtr(recv_count),1,MPI_INT,MPI_COMM_WORLD);
    for (int i=1; i<nprocs; i++)
        recv_disp[i] = recv_disp[i-1] + recv_count[i-1];
    std::vector<int> recv_data(recv_disp[nprocs-1]+recv_count[nprocs-1]);
    MPI_Allgatherv(getPtr(send_data),send_count,MPI_INT,
        getPtr(recv_data),getPtr(recv_count),getPtr(recv_disp),MPI_INT,
        MPI_COMM_WORLD);
    size_t i=0;
    while ( i < recv_data.size() ) {
        int id = recv_data[i];
        int count = recv_data[i+1];
        i += 2;
        std::set<int>& src_ids = src_map[id];
        for (int j=0; j<count; j++,i++)
            src_ids.insert(recv_data[i]);
    }
}
void addSrcDstIDs( int src_id, std::map<int,std::set<int> >& src_map, 
    std::map<int,std::set<int> >& dst_map, std::set<int>& src, std::set<int>& dst )
{
    src.insert(src_id);
    const std::set<int>& dst_ids = dst_map[src_id];
    for (std::set<int>::const_iterator it=dst_ids.begin(); it!=dst_ids.end(); ++it) {
        if ( dst.find(*it)==dst.end() )
            addSrcDstIDs(*it,dst_map,src_map,dst,src);
    }
}
ID_map_struct computeIDMap( const IntArray& ID1, const IntArray& ID2 )
{
    ASSERT(ID1.size(0)==ID2.size(0)&&ID1.size(1)==ID2.size(1)&&ID1.size(2)==ID2.size(2));
    PROFILE_START("computeIDMap");

    // Get a global list of all src/dst ids and the src map for each local blob
    std::set<int> src_set, dst_set;
    std::map<int,std::set<int> > src_map;   // Map of the src ids for each dst id
    for (size_t i=0; i<ID1.length(); i++) {
        if ( ID1(i)>=0 )
            src_set.insert(ID1(i));
        if ( ID2(i)>=0 )
            dst_set.insert(ID2(i));
        if ( ID2(i)>=0 && ID1(i)>=0 ) {
            std::set<int>& src_ids = src_map[ID2(i)];
            src_ids.insert(ID1(i));
        }
    }
    // Communicate the src/dst ids and src id map to all processors and reduce
    gatherSet( src_set );
    gatherSet( dst_set );
    gatherSrcIDMap( src_map );
    // Compute the dst id map
    std::map<int,std::set<int> > dst_map;   // Map of the dst ids for each src id
    for (std::map<int,std::set<int> >::const_iterator it=src_map.begin(); it!=src_map.end(); ++it) {
        int id = it->first;
        const std::set<int>& src_ids = it->second;
        for (std::set<int>::const_iterator it2=src_ids.begin(); it2!=src_ids.end(); ++it2) {
            std::set<int>& dst_ids = dst_map[*it2];
            dst_ids.insert(id);
        }
    }

    // Perform the mapping of ids
    ID_map_struct id_map;
    // Find new blobs
    for (std::set<int>::const_iterator it=dst_set.begin(); it!=dst_set.end(); ++it) {
        if ( src_map.find(*it)==src_map.end() )
            id_map.created.push_back(*it);
    }
    // Fine blobs that disappeared
    for (std::set<int>::const_iterator it=src_set.begin(); it!=src_set.end(); ++it) {
        if ( dst_map.find(*it)==dst_map.end() )
            id_map.destroyed.push_back(*it);
    }
    // Find blobs with a 1-to-1 mapping
    std::vector<int> dst_list;
    dst_list.reserve(src_map.size());
    for (std::map<int,std::set<int> >::const_iterator it=src_map.begin(); it!=src_map.end(); ++it)
        dst_list.push_back(it->first);
    for (size_t i=0; i<dst_list.size(); i++) {
        int dst_id = dst_list[i];
        const std::set<int>& src_ids = src_map[dst_id];
        if ( src_ids.size()==1 ) {
            int src_id = *src_ids.begin();
            const std::set<int>& dst_ids = dst_map[src_id];
            if ( dst_ids.size()==1 ) {
                ASSERT(*dst_ids.begin()==dst_id);
                src_map.erase(dst_id);
                dst_map.erase(src_id);
                id_map.src_dst.push_back(std::pair<int,int>(src_id,dst_id));
            }
        }
    }
    // Handle merge/splits
    while ( !dst_map.empty() ) {
        // Get a lit of the src-dst ids
        std::set<int> src, dst;
        addSrcDstIDs( dst_map.begin()->first, src_map, dst_map, src, dst );
        for (std::set<int>::const_iterator it=src.begin(); it!=src.end(); ++it)
            dst_map.erase(*it);
        for (std::set<int>::const_iterator it=dst.begin(); it!=dst.end(); ++it)
            src_map.erase(*it);
        if ( src.size()==1 ) {
            // Bubble split
            id_map.split.push_back( BlobIDSplitStruct(*src.begin(),std::vector<int>(dst.begin(),dst.end())) );
        } else if ( dst.size()==1 ) {
            // Bubble merge
            id_map.merge.push_back( BlobIDMergeStruct(std::vector<int>(src.begin(),src.end()),*dst.begin()) );
        } else {
            // Bubble split/merge
            id_map.merge_split.push_back( BlobIDMergeSplitStruct(
                std::vector<int>(src.begin(),src.end()), std::vector<int>(dst.begin(),dst.end() ) ) );
        }
    }
    ASSERT(src_map.empty());

    PROFILE_STOP("computeIDMap");
    return id_map;
}



