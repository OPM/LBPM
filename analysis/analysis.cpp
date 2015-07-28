#include "analysis/analysis.h"
#include "ProfilerApp.h"
#include <iostream>


template<class TYPE>
inline TYPE* getPtr( std::vector<TYPE>& x ) { return x.empty() ? NULL:&x[0]; }
template<class TYPE>
inline const TYPE* getPtr( const std::vector<TYPE>& x ) { return x.empty() ? NULL:&x[0]; }


/******************************************************************
 * Compute the blobs                                              *
 ******************************************************************/
static int ComputePhaseComponent( IntArray &ComponentLabel,
		const IntArray &PhaseID, int VALUE, int &ncomponents,
		int startx, int starty, int startz, IntArray &temp, bool periodic)
{

	// Get the dimensions
	int Nx = PhaseID.size(0);  // maxima for the meshes
	int Ny = PhaseID.size(1);
	int Nz = PhaseID.size(2);

	int cubes_in_blob=0;
	int nrecent = 1;                    // number of nodes added at most recent sweep
	temp(0,0) = startx;                 // Set the initial point as a "seed" for the sweeps
	temp(1,0) = starty;
	temp(2,0) = startz;
	int ntotal = 1;                     // total number of nodes in blob
	ComponentLabel(startx,starty,startz) = ncomponents;

	int p,s,x,y,z,site,start,finish,nodx,nody,nodz;
	int imin=startx,imax=startx,jmin=starty,jmax=starty;    // initialize maxima / minima
	int kmin=startz,kmax=startz;
	int d[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
			{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
			{1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
			{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
			{-1,-1,1},{-1,-1,-1}};   // directions to neighbors
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	bool status = 1;                    // status == true => continue to look for points
	while (status == 1){
		start = ntotal - nrecent;
		finish = ntotal;
		nrecent = 0;                    // set recent points back to zero for next sweep through
		for (s=start;s<finish;s++){
			// Loop over recent points; look for new points
			x = temp(0,s);
			y = temp(1,s);
			z = temp(2,s);
			// Looop over the directions
			for (p=0;p<26;p++){
				nodx=x+d[p][0];
				nody=y+d[p][1];
				nodz=z+d[p][2];
				if ( periodic ) {
					if (nodx < 0 ){ nodx = Nx-1; }        // Periodic BC for x
					if (nodx > Nx-1 ){ nodx = 0; }
					if (nody < 0 ){ nody = Ny-1; }        // Periodic BC for y
					if (nody > Ny-1 ){ nody = 0; }
					if (nodz < 0 ){ nodz = Nz-1; }        // Periodic BC for z
					if (nodz > Nz-1 ){ nodz = 0; }
				} else {
					if ( nodx<0 || nodx>=Nx || nody<0 || nody>=Ny || nodz<0 || nodz>=Nz )
						continue;
				}
				site = nodz*Nx*Ny+nody*Nx+nodx;
				if ( PhaseID(nodx,nody,nodz) == VALUE && ComponentLabel(nodx,nody,nodz) == -1 ){
					// Node is a part of the blob - add it to the list
					temp(0,ntotal) = nodx;
					temp(1,ntotal) = nody;
					temp(2,ntotal) = nodz;
					ntotal++;
					nrecent++;
					// Update the component map for phase
					ComponentLabel(nodx,nody,nodz) = ncomponents;
					// Update the min / max for the cube loop
					if ( nodx < imin ){ imin = nodx; }
					if ( nodx > imax ){ imax = nodx; }
					if ( nody < jmin ){ jmin = nody; }
					if ( nody > jmax ){ jmax = nody; }
					if ( nodz < kmin ){ kmin = nodz; }
					if ( nodz > kmax ){ kmax = nodz; }
				}
			}

		}
		if ( nrecent == 0){
			status = 0;
		}
	}
	return ntotal;
}

int ComputeBlob(IntArray &blobs, int &nblobs, int &ncubes, IntArray &indicator,
    const DoubleArray &F, const DoubleArray &S, double vf, double vs, 
    int startx, int starty, int startz, IntArray &temp, bool periodic)
{
    // Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
    // F>vf => oil phase S>vs => in porespace
    // update the list of blobs, indicator mesh
    int m = F.size(0);  // maxima for the meshes
    int n = F.size(1);
    int o = F.size(2);

    int cubes_in_blob=0;
    int nrecent = 1;                    // number of nodes added at most recent sweep
    temp(0,0) = startx;                 // Set the initial point as a "seed" for the sweeps
    temp(1,0) = starty;
    temp(2,0) = startz;
    int ntotal = 1;                     // total number of nodes in blob
    indicator(startx,starty,startz) = nblobs;

    int p,s,x,y,z,start,finish,nodx,nody,nodz;
    int imin=startx,imax=startx,jmin=starty,jmax=starty;    // initialize maxima / minima
    int kmin=startz,kmax=startz;
    int d[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
                    {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
                    {1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
                    {1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
                    {-1,-1,1},{-1,-1,-1}};   // directions to neighbors
    int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
    bool status = 1;                    // status == true => continue to look for points
    while (status == 1){
        start = ntotal - nrecent;
        finish = ntotal;
        nrecent = 0;                    // set recent points back to zero for next sweep through
        for (s=start;s<finish;s++){
            // Loop over recent points; look for new points
            x = temp(0,s);
            y = temp(1,s);
            z = temp(2,s);
            // Looop over the directions
            for (p=0;p<26;p++){
                nodx=x+d[p][0];
                nody=y+d[p][1];
                nodz=z+d[p][2];
                if ( periodic ) {
                    if (nodx < 0 ){ nodx = m-1; }        // Periodic BC for x
                    if (nodx > m-1 ){ nodx = 0; }
                    if (nody < 0 ){ nody = n-1; }        // Periodic BC for y
                    if (nody > n-1 ){ nody = 0; }
                    if (nodz < 0 ){ nodz = o-1; }        // Periodic BC for z
                    if (nodz > o-1 ){ nodz = 0; }
                } else {
                    if ( nodx<0 || nodx>=m || nody<0 || nody>=n || nodz<0 || nodz>=o )
                        continue;
                }
                if ( F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
                     && indicator(nodx,nody,nodz) == -1 ){
                    // Node is a part of the blob - add it to the list
                    temp(0,ntotal) = nodx;
                    temp(1,ntotal) = nody;
                    temp(2,ntotal) = nodz;
                    ntotal++;
                    nrecent++;
                    // Update the indicator map
                    indicator(nodx,nody,nodz) = nblobs;
                    // Update the min / max for the cube loop
                    if ( nodx < imin ){ imin = nodx; }
                    if ( nodx > imax ){ imax = nodx; }
                    if ( nody < jmin ){ jmin = nody; }
                    if ( nody > jmax ){ jmax = nody; }
                    if ( nodz < kmin ){ kmin = nodz; }
                    if ( nodz > kmax ){ kmax = nodz; }
                }
                else if (F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
                         && indicator(nodx,nody,nodz) > -1 &&  indicator(nodx,nody,nodz) != nblobs){
                    // Some kind of error in algorithm
                    printf("Error in blob search algorithm!");
                }
            }

        }
        if ( nrecent == 0){
            status = 0;
        }
    }
    // Use points in temporary storage array to add cubes to the list of blobs
    if ( imin > 0) { imin = imin-1; }
//    if ( imax < m-1) { imax = imax+1; }
    if ( jmin > 0) { jmin = jmin-1; }
//    if ( jmax < n-1) { jmax = jmax+1; }
    if ( kmin > 0) { kmin = kmin-1; }
//    if ( kmax < o-1) { kmax = kmax+1; }
    int i,j,k;
    bool add;
    for (k=kmin;k<kmax;k++){
        for (j=jmin;j<jmax;j++){
            for (i=imin;i<imax;i++){
                // If cube(i,j,k) has any nodes in blob, add it to the list
                // Loop over cube edges
                add = 0;
                for (p=0;p<8;p++){
                    nodx = i+cube[p][0];
                    nody = j+cube[p][1];
                    nodz = k+cube[p][2];
                    if ( indicator(nodx,nody,nodz) == nblobs ){
                        // Cube corner is in this blob
                        add = 1;
                    }
                }
                if (add == 1){
                    // Add cube to the list
                    blobs(0,ncubes) = i;
                    blobs(1,ncubes) = j;
                    blobs(2,ncubes) = k;
                    ncubes++;
                    cubes_in_blob++;
                    // Loop again to check for overlap
                    for (p=0;p<8;p++){
                        nodx = i+cube[p][0];
                        nody = j+cube[p][1];
                        nodz = k+cube[p][2];
                        if (indicator(nodx,nody,nodz) > -1 && indicator(nodx,nody,nodz) != nblobs){
                            printf("Overlapping cube!");
                            std::cout << i << ", " << j << ", " << k << std::endl;
                        }
                    }
                }
            }
        }
    }

    return cubes_in_blob;
}


/******************************************************************
* Compute the local blob ids                                      *
******************************************************************/
int ComputeLocalBlobIDs( const DoubleArray& Phase, const DoubleArray& SignDist, 
    double vF, double vS, IntArray& LocalBlobID, bool periodic )
{
    PROFILE_START("ComputeLocalBlobIDs");
    size_t Nx = Phase.size(0);
    size_t Ny = Phase.size(1);
    size_t Nz = Phase.size(2);
    ASSERT(SignDist.size(0)==Nx&&SignDist.size(1)==Ny&&SignDist.size(2)==Nz);
    // Initialize output
    LocalBlobID.resize(Nx,Ny,Nz);
    // Compute the local blob ids
    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
    size_t N = Nx*Ny*Nz;
    int nblobs = 0;
    int ncubes = 0;         // total number of nodes in any blob
    IntArray blobs(3,N);    // store indices for blobs (cubes)
    IntArray temp(3,N);     // temporary storage array
    IntArray b(N);          // number of nodes in each blob
    for (size_t k=0; k<Nz; k++ ) {
        for (size_t j=0; j<Ny; j++) {
            for (size_t i=0; i<Nx; i++) {
                if ( SignDist(i,j,k) < 0.0) {
                    // Solid phase 
                    LocalBlobID(i,j,k) = -2;
                } else{
                    LocalBlobID(i,j,k) = -1;
                }
            }
        }
    }
    for (size_t k=0; k<Nz; k++ ) {
        for (size_t j=0; j<Ny; j++) {
            for (size_t i=0; i<Nx; i++) {
                if ( LocalBlobID(i,j,k)==-1 && Phase(i,j,k)>vF && SignDist(i,j,k)>vS ) {
                    // node i,j,k is in the porespace
                    b(nblobs) = ComputeBlob(blobs,nblobs,ncubes,LocalBlobID,
                        Phase,SignDist,vF,vS,i,j,k,temp,periodic);
                    nblobs++;                              
                }
                if ( nblobs > (int)b.length()-1){
                    printf("Increasing size of blob list \n");
                    b.resize(2*b.length());
                }
            }
        }
    }
    // Go over all cubes again -> add any that do not contain nw phase
    size_t count_in=0,count_out=0;
    size_t nodx,nody,nodz;
    for (size_t k=0; k<Nz-1; k++ ) {
        for (size_t j=0; j<Ny-1; j++) {
            for (size_t i=0; i<Nx-1; i++) {
                // Loop over cube corners
                int add=1;      // initialize to true - add unless corner occupied by nw-phase
                for (int p=0; p<8; p++) {
                    nodx=i+cube[p][0];
                    nody=j+cube[p][1];
                    nodz=k+cube[p][2];
                    if ( LocalBlobID(nodx,nody,nodz) > -1 ){
                        // corner occupied by nw-phase  -> do not add
                        add = 0;
                    }
                }
                if ( add == 1 ){
                    blobs(0,ncubes) = i;
                    blobs(1,ncubes) = j;
                    blobs(2,ncubes) = k;
                    ncubes++;
                    count_in++;
                }
                else { count_out++; }
            }
        }
    }
    b(nblobs) = count_in;
    PROFILE_STOP("ComputeLocalBlobIDs");
    return nblobs;
}

int ComputeLocalPhaseComponent(const IntArray &PhaseID, int VALUE, IntArray &ComponentLabel, bool periodic )
{
    PROFILE_START("ComputeLocalPhaseComponent");
    size_t Nx = PhaseID.size(0);
    size_t Ny = PhaseID.size(1);
    size_t Nz = PhaseID.size(2);
    // Compute the local blob ids

	ComponentLabel.resize(Nx,Ny,Nz);

    const int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
    size_t N = Nx*Ny*Nz;
    int ncomponents = 0;
    IntArray temp(3,N);     // temporary storage array

    for (size_t k=0; k<Nz; k++ ) {
        for (size_t j=0; j<Ny; j++) {
            for (size_t i=0; i<Nx; i++) {
                if ( PhaseID(i,j,k) == VALUE) {
                    // Solid phase
                    ComponentLabel(i,j,k) = -1;
                } else{
                    ComponentLabel(i,j,k) = -2;
                }
            }
        }
    }
    for (size_t k=0; k<Nz; k++ ) {
        for (size_t j=0; j<Ny; j++) {
            for (size_t i=0; i<Nx; i++) {
                if ( ComponentLabel(i,j,k)==-1 && PhaseID(i,j,k) == VALUE) {
                	// component is part of the phase we want
                    ComputePhaseComponent(ComponentLabel,PhaseID, VALUE, ncomponents,
                    		i,j,k,temp,periodic);
                    ncomponents++;
                }
            }
        }
    }
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
            if ( LocalIDs(i)!=IDs(i) )
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
    int ngx = (IDs.size(0)-nx)/2;
    int ngy = (IDs.size(1)-ny)/2;
    int ngz = (IDs.size(2)-nz)/2;
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
    const int rank = rank_info.rank[1][1][1];
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



