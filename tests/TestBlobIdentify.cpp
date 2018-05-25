// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "analysis/pmmc.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"



// Get a random number in [0 1]
inline double rand2()
{
    return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}


// Test if all ranks agree on a value
bool allAgree( int x, MPI_Comm comm ) {
    int x2 = x;
    MPI_Bcast(&x2,1,MPI_INT,0,comm);
    int diff = x==x2 ? 0:1;
    int diff2 = 0;
    MPI_Allreduce(&diff,&diff2,1,MPI_INT,MPI_SUM,comm);
    return diff2==0;
}
template<class T>
bool allAgree( const std::vector<T>& x, MPI_Comm comm ) {
    std::vector<T> x2 = x;
    MPI_Bcast(&x2[0],x.size()*sizeof(T)/sizeof(int),MPI_INT,0,comm);
    int diff = x==x2 ? 0:1;
    int diff2 = 0;
    MPI_Allreduce(&diff,&diff2,1,MPI_INT,MPI_SUM,comm);
    return diff2==0;
}


// Structure to hold a bubble
struct bubble_struct {
    Point center;
    double radius;
    // Get the distance to the bubble
    inline double dist( const Point& p, double Lx, double Ly, double Lz ) const {
        double x = fabs(p.x-center.x);
        double y = fabs(p.y-center.y);
        double z = fabs(p.z-center.z);
        x = std::min(x,fabs(Lx-x));
        y = std::min(y,fabs(Ly-y));
        z = std::min(z,fabs(Lz-z));
        return sqrt(x*x+y*y+z*z)-radius;
    }
    // Check if this bubble overlaps with rhs
    bool overlap( const bubble_struct& rhs, double Lx, double Ly, double Lz ) const {
        return dist(rhs.center,Lx,Ly,Lz) <= radius+rhs.radius;
    }
    // Create a random bubble
    static bubble_struct rand( double Lx, double Ly, double Lz, double R0 ) {
        bubble_struct bubble;
        bubble.center.x = Lx*rand2();
        bubble.center.y = Ly*rand2();
        bubble.center.z = Lz*rand2();
        bubble.radius = R0;
        bubble.radius *= 1 + 0.4*(1.0-2.0*rand2());
        return bubble;
    }
};


// Create a random set of bubles
std::vector<bubble_struct> create_bubbles( int N_bubbles, double Lx, double Ly, double Lz, MPI_Comm comm )
{
    int rank = comm_rank(comm);
    std::vector<bubble_struct> bubbles(N_bubbles);
    if ( rank == 0 ) {
        double R0 = 0.2*Lx*Ly*Lz/pow((double)N_bubbles,0.333);
        for (int i=0; i<N_bubbles; i++) {
            bool overlap = true;
            while ( overlap ) {
                bubbles[i] = bubble_struct::rand(Lx,Ly,Lz,R0);
                overlap = false;
                for (int k=0; k<i; k++)
                    overlap = overlap || bubbles[i].overlap(bubbles[k],Lx,Ly,Lz);
            }
        }
    }
    size_t N_bytes = N_bubbles*sizeof(bubble_struct);
    MPI_Bcast((char*)&bubbles[0],N_bytes,MPI_CHAR,0,comm);
    return bubbles;
}


// Fill the Phase/SignDist info from the bubble list
void fillBubbleData( const std::vector<bubble_struct>& bubbles, DoubleArray& Phase,
    DoubleArray& SignDist, double Lx, double Ly, double Lz, const RankInfoStruct rank_info )
{
    PROFILE_START("fillBubbleData");
    int nx = Phase.size(0)-2;
    int ny = Phase.size(1)-2;
    int nz = Phase.size(2)-2;
    Phase.fill(1);
    SignDist.fill(0);
    for (int k=0; k<nz; k++) {
        double z = Lz*(nz*rank_info.kz+k+0.5)/(nz*rank_info.nz);
        for (int j=0; j<ny; j++) {
            double y = Ly*(ny*rank_info.jy+j+0.5)/(ny*rank_info.ny);
            for (int i=0; i<nx; i++) {
                double x = Lx*(nx*rank_info.ix+i+0.5)/(nx*rank_info.nx);
                double dist = 1e100;
                for (size_t b=0; b<bubbles.size(); b++)
                    dist = std::min(dist,bubbles[b].dist(Point(x,y,z),Lx,Ly,Lz));
                SignDist(i+1,j+1,k+1) = -dist;
            }
        }
    }
    PROFILE_STOP("fillBubbleData");
}


// Shift all of the data by the given number of cells
void shift_data( DoubleArray& data, int sx, int sy, int sz, const RankInfoStruct& rank_info, MPI_Comm comm )
{
    int nx = data.size(0)-2;
    int ny = data.size(1)-2;
    int nz = data.size(2)-2;
    int ngx = nx+2*abs(sx);
    int ngy = ny+2*abs(sy);
    int ngz = nz+2*abs(sz);
    Array<double> tmp1(nx,ny,nz), tmp2(ngx,ngy,ngz), tmp3(ngx,ngy,ngz);
    fillHalo<double> fillData1(comm,rank_info,{nx,ny,nz},{1,1,1},0,1);
    fillHalo<double> fillData2(comm,rank_info,{nx,ny,nz},{abs(sx),abs(sy),abs(sz)},0,1);
    fillData1.copy(data,tmp1);
    fillData2.copy(tmp1,tmp2);
    fillData2.fill(tmp2);
    for (int k=abs(sz); k<nz+abs(sz); k++) {
        for (int j=abs(sy); j<ny+abs(sy); j++) {
            for (int i=abs(sx); i<nx+abs(sx); i++)
                tmp3(i,j,k) = tmp2(i-sx,j-sy,k-sz);
        }
    }
    fillData2.copy(tmp3,tmp1);
    fillData1.copy(tmp1,data);
    fillData1.fill(data);
}


// Main
int main(int argc, char **argv)
{
    // Initialize MPI
    int rank, nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    PROFILE_ENABLE(1);
    PROFILE_DISABLE_TRACE();
    PROFILE_SYNCHRONIZE();
    PROFILE_START("main");
    if ( rank==0 ) {
        printf("-----------------------------------------------------------\n");
        printf("Testing Blob Identification \n");
        printf("-----------------------------------------------------------\n");
    }

    // Set the domain information
    std::vector<int> factors = Utilities::factor(nprocs);
    int nproc[3]={1,1,1};
    for (size_t i=0; i<factors.size(); i++)
        nproc[i%3] *= factors[i];
    int nx = 100/nproc[0];
    int ny = 120/nproc[1];
    int nz = 110/nproc[2];
    double Lx=1.0, Ly=1.0, Lz=1.0;

    // Get the rank info
    const RankInfoStruct rank_info(rank,nproc[0],nproc[1],nproc[2]);

    // Create the dummy info
    DoubleArray Phase(nx+2,ny+2,nz+2);
    DoubleArray SignDist(nx+2,ny+2,nz+2);
    std::vector<bubble_struct> bubbles = create_bubbles(20,Lx,Ly,Lz,comm);
    fillBubbleData( bubbles, Phase, SignDist, Lx, Ly, Lz, rank_info );

    // Communication the halos
    fillHalo<double> fillData(comm,rank_info,{nx,ny,nz},{1,1,1},0,1);
    fillData.fill(Phase);
    fillData.fill(SignDist);

    // Find blob domains
    if ( rank==0 ) { printf("Finding blob domains\n"); }
    double vF=0.0;
    double vS=0.0;
    IntArray GlobalBlobID;
    int nblobs = ComputeGlobalBlobIDs(nx,ny,nz,rank_info,
        Phase,SignDist,vF,vS,GlobalBlobID,comm);
    if ( rank==0 ) { printf("Identified %i blobs\n",nblobs); }

    // Create the MeshDataStruct
    std::vector<IO::MeshDataStruct> meshData(1);
    meshData[0].meshName = "domain";
    meshData[0].mesh = std::make_shared<IO::DomainMesh>(rank_info,nx,ny,nz,Lx,Ly,Lz);
    std::shared_ptr<IO::Variable> PhaseVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> SignDistVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> BlobIDVar( new IO::Variable() );
    PhaseVar->name = "phase";
    PhaseVar->type = IO::VariableType::VolumeVariable;
    PhaseVar->dim = 1;
    PhaseVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(PhaseVar);
    SignDistVar->name = "SignDist";
    SignDistVar->type = IO::VariableType::VolumeVariable;
    SignDistVar->dim = 1;
    SignDistVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(SignDistVar);
    BlobIDVar->name = "BlobID";
    BlobIDVar->type = IO::VariableType::VolumeVariable;
    BlobIDVar->dim = 1;
    BlobIDVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(BlobIDVar);

    // Save the results
    fillData.copy(Phase,PhaseVar->data);
    fillData.copy(SignDist,SignDistVar->data);
    fillData.copy(GlobalBlobID,BlobIDVar->data);
    IO::writeData( 0, meshData, comm );
    writeIDMap(ID_map_struct(nblobs),0,"lbpm_id_map.txt");
    int save_it = 1;

    // Check the results
    int N_errors = 0;
    if ( nblobs != (int) bubbles.size() ) {
        printf("Error, detected number of bubbles %i does not match expected %i\n",nblobs,(int)bubbles.size());
        N_errors++;
    }


    // Move the blobs and connect them to the previous ids
    PROFILE_START("constant velocity test");
    if ( rank==0 ) { printf("Running constant velocity blob test\n"); }
    int id_max = nblobs-1;
    for (int i=0; i<20; i++, save_it++) {
        // Shift all the data
        shift_data( Phase,    3, -2, 1, rank_info, comm );
        shift_data( SignDist, 3, -2, 1, rank_info, comm );
        // Find blob domains
        IntArray GlobalBlobID2;
        int nblobs2 = ComputeGlobalBlobIDs(nx,ny,nz,rank_info,
            Phase,SignDist,vF,vS,GlobalBlobID2,comm);
        if ( nblobs2 != nblobs ) {
            printf("Error, number of blobs changed under constant velocity (%i,%i)\n",nblobs,nblobs2);
            N_errors++;
        }
        // Identify the blob maps and renumber the ids
        ID_map_struct map = computeIDMap(nx,ny,nz,GlobalBlobID,GlobalBlobID2,comm);
        std::swap(GlobalBlobID,GlobalBlobID2);
        std::vector<BlobIDType> new_list;
        getNewIDs(map,id_max,new_list);
        renumberIDs(new_list,GlobalBlobID);
        writeIDMap(map,save_it,"lbpm_id_map.txt");
        bool pass = (int)map.src_dst.size()==nblobs;
        pass = pass && map.created.empty();
        pass = pass && map.destroyed.empty();
        for (size_t j=0; j<map.src_dst.size(); j++)
            pass = pass && map.src_dst[j].first==map.src_dst[j].second;
        pass = pass && map.split.empty();
        pass = pass && map.merge.empty();
        pass = pass && map.merge_split.empty();
        pass = pass && id_max==nblobs-1;
        if ( !pass ) {
            printf("Error, blob ids do not match in constant velocity test\n");
            N_errors++;
        }
        // Save the results
        fillData.copy(Phase,PhaseVar->data);
        fillData.copy(SignDist,SignDistVar->data);
        fillData.copy(GlobalBlobID,BlobIDVar->data);
        IO::writeData( save_it, meshData, comm );
    }
    PROFILE_STOP("constant velocity test");


    // Move the bubbles in independent directions to test merge/splitting of bubbles
    PROFILE_START("moving bubble test");
    if ( rank==0 ) { printf("Running moving bubble test\n"); }
    std::vector<Point> velocity(bubbles.size());
    if ( rank==0 ) {
        for (size_t i=0; i<bubbles.size(); i++) {
            velocity[i].x = bubbles[i].radius*(2*rand2()-1);
            velocity[i].y = bubbles[i].radius*(2*rand2()-1);
            velocity[i].z = bubbles[i].radius*(2*rand2()-1);
        }
    }
    MPI_Bcast((char*)&velocity[0],bubbles.size()*sizeof(Point),MPI_CHAR,0,comm);
    fillBubbleData( bubbles, Phase, SignDist, Lx, Ly, Lz, rank_info );
    fillData.fill(Phase);
    fillData.fill(SignDist);
    ComputeGlobalBlobIDs(nx,ny,nz,rank_info,Phase,SignDist,vF,vS,GlobalBlobID,comm);
    fillData.copy(Phase,PhaseVar->data);
    fillData.copy(SignDist,SignDistVar->data);
    fillData.copy(GlobalBlobID,BlobIDVar->data);
    IO::writeData( save_it, meshData, comm );
    save_it++;
    id_max = nblobs-1;
    for (int i=0; i<25; i++, save_it++) {
        // Move the bubbles
        for (size_t j=0; j<bubbles.size(); j++) {
            bubbles[j].center.x = fmod(bubbles[j].center.x+velocity[j].x+Lx,Lx);
            bubbles[j].center.y = fmod(bubbles[j].center.y+velocity[j].y+Ly,Ly);
            bubbles[j].center.z = fmod(bubbles[j].center.z+velocity[j].z+Lz,Lz);
        }
        fillBubbleData( bubbles, Phase, SignDist, Lx, Ly, Lz, rank_info );
        fillData.fill(Phase);
        fillData.fill(SignDist);
        // Compute the ids
        IntArray GlobalBlobID2;
        int nblobs2 = ComputeGlobalBlobIDs(nx,ny,nz,rank_info,Phase,SignDist,vF,vS,GlobalBlobID2,comm);
        // Identify the blob maps and renumber the ids
        ID_map_struct map = computeIDMap(nx,ny,nz,GlobalBlobID,GlobalBlobID2,comm);
        std::swap(GlobalBlobID,GlobalBlobID2);
        std::vector<BlobIDType> new_list;
        getNewIDs(map,id_max,new_list);
        // Check id_max
        if ( !allAgree(id_max,comm) ) {
            if ( rank==0 )
                printf("All ranks do not agree on id_max\n");
            N_errors++;
            break;
        }
        // Check that the new id list matches on all ranks
        if ( !allAgree(new_list,comm) ) {
            if ( rank==0 )
                printf("All ranks do not agree on new_list\n");
            N_errors++;
            break;
        }
        // Renumber the ids and write the map
        renumberIDs(new_list,GlobalBlobID);
        writeIDMap(map,save_it,"lbpm_id_map.txt");
        // Check the number of blobs in the map
        int N1 = 0;
        int N2 = 0;
        if ( rank==0 ) {
            if ( !map.destroyed.empty() ) {
                N1 += map.destroyed.size();
                printf("  %i: destroyed:",i+1);
                for (size_t j=0; j<map.destroyed.size(); j++)
                    printf(" %i",map.destroyed[j]);
                printf("\n");
            }
            if ( !map.created.empty() ) {
                N2 += map.created.size();
                printf("  %i: created:",i+1);
                for (size_t j=0; j<map.created.size(); j++)
                    printf(" %i",map.created[j]);
                printf("\n");
            }
            N1 += map.src_dst.size();
            N2 += map.src_dst.size();
            for (size_t j=0; j<map.split.size(); j++ ) {
                N1 += 1;
                N2 += map.split[j].second.size();
                printf("  %i: split: %i -",i+1,map.split[j].first);
                for (size_t k=0; k<map.split[j].second.size(); k++)
                    printf(" %i",map.split[j].second[k]);
                printf("\n");
            }
            for (size_t j=0; j<map.merge.size(); j++ ) {
                N1 += map.merge[j].first.size();
                N2 += 1;
                printf("  %i: merged:",i+1);
                for (size_t k=0; k<map.merge[j].first.size(); k++)
                    printf(" %i",map.merge[j].first[k]);
                printf(" - %i\n",map.merge[j].second);
            }
            for (size_t j=0; j<map.merge_split.size(); j++ ) {
                N1 += map.merge_split[j].first.size();
                N2 += map.merge_split[j].second.size();
                printf("  %i: merged/split:",i+1);
                for (size_t k=0; k<map.merge_split[j].first.size(); k++)
                    printf(" %i",map.merge_split[j].first[k]);
                printf(" -");
                for (size_t k=0; k<map.merge_split[j].second.size(); k++)
                    printf(" %i",map.merge_split[j].second[k]);
                printf("\n");
            }
        }
        MPI_Bcast(&N1,1,MPI_INT,0,comm);
        MPI_Bcast(&N2,1,MPI_INT,0,comm);
        if ( N1!=nblobs || N2!=nblobs2 ) {
            if ( rank==0 )
                printf("Error, blob ids do not map in moving bubble test (%i,%i,%i,%i)\n",
                    nblobs,nblobs2,N1,N2);
            N_errors++;
        }
        nblobs = nblobs2;
        // Save the results
        fillData.copy(Phase,PhaseVar->data);
        fillData.copy(SignDist,SignDistVar->data);
        fillData.copy(GlobalBlobID,BlobIDVar->data);
        IO::writeData( save_it, meshData, comm );
    }
    PROFILE_STOP("moving bubble test");


    // Finished
    PROFILE_STOP("main");
    PROFILE_SAVE("TestBlobIdentify",false);
    MPI_Barrier(comm);
    MPI_Finalize();
    return N_errors;  
}

