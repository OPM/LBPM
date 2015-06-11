// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include <iostream>
#include <math.h>
#include "common/pmmc.h"
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


struct bubble_struct {
    Point center;
    double radius;
    // Get the distance to the bubble
    double dist( const Point& p, double Lx, double Ly, double Lz ) {
        double x = std::min(fabs(p.x-center.x),std::min(fabs(p.x-center.x-Lx),fabs(p.x-center.x+Lx)));
        double y = std::min(fabs(p.y-center.y),std::min(fabs(p.y-center.y-Ly),fabs(p.y-center.y+Ly)));
        double z = std::min(fabs(p.z-center.z),std::min(fabs(p.z-center.z-Lz),fabs(p.z-center.z+Lz)));
        return sqrt(x*x+y*y+z*z)-radius;
    }
    // Check if this bubble overlaps with rhs
    bool overlap( const bubble_struct& rhs, double Lx, double Ly, double Lz ) {
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
std::vector<bubble_struct> create_bubbles( int N_bubbles, double Lx, double Ly, double Lz )
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
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
    MPI_Bcast(&bubbles[0],N_bytes,MPI_CHAR,0,MPI_COMM_WORLD);
    return bubbles;
}


// Main
int main(int argc, char **argv)
{
    // Initialize MPI
    int rank, nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    PROFILE_ENABLE(0);
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
    Phase.fill(1);
    SignDist(0);
    std::vector<bubble_struct> bubbles = create_bubbles(20,Lx,Ly,Lz);
    for (int k=0; k<nz; k++) {
        double z = Lz*(nz*rank_info.kz+k+0.5)/(nz*nproc[2]);
        for (int j=0; j<ny; j++) {
            double y = Ly*(ny*rank_info.jy+j+0.5)/(ny*nproc[1]);
            for (int i=0; i<nx; i++) {
                double x = Lx*(nx*rank_info.ix+i+0.5)/(nx*nproc[0]);
                double dist = 1e100;
                for (size_t b=0; b<bubbles.size(); b++)
                    dist = std::min(dist,bubbles[b].dist(Point(x,y,z),Lx,Ly,Lz));
                SignDist(i+1,j+1,k+1) = -dist;
            }
        }
    }

    // Communication the halos
    fillHalo<double> fillData(rank_info,nx,ny,nz,1,1,1,0,1);
    fillData.fill(Phase);
    fillData.fill(SignDist);

    // Find blob domains
    if ( rank==0 ) { printf("Finding blob domains\n"); }
    double vF=0.0;
    double vS=0.0;
    IntArray GlobalBlobID;
    int nblobs = ComputeGlobalBlobIDs(nx,ny,nz,rank_info,
        Phase,SignDist,vF,vS,GlobalBlobID);
    if ( rank==0 ) { printf("Identified %i blobs\n",nblobs); }

    // Create the MeshDataStruct
    std::vector<IO::MeshDataStruct> meshData(1);
    meshData[0].meshName = "domain";
    meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(rank_info,nx,ny,nz,Lx,Ly,Lz) );
    std::shared_ptr<IO::Variable> PhaseVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> SignDistVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> BlobIDVar( new IO::Variable() );
    PhaseVar->name = "phase";
    PhaseVar->type = IO::VolumeVariable;
    PhaseVar->dim = 1;
    PhaseVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(PhaseVar);
    SignDistVar->name = "SignDist";
    SignDistVar->type = IO::VolumeVariable;
    SignDistVar->dim = 1;
    SignDistVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(SignDistVar);
    BlobIDVar->name = "BlobID";
    BlobIDVar->type = IO::VolumeVariable;
    BlobIDVar->dim = 1;
    BlobIDVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(BlobIDVar);

    // Save the results
    fillData.copy(Phase,PhaseVar->data);
    fillData.copy(SignDist,SignDistVar->data);
    fillData.copy(GlobalBlobID,BlobIDVar->data);
    IO::writeData( 0, meshData, 2 );

    // Check the results
    int N_errors = 0;
    if ( nblobs != (int) bubbles.size() ) {
        printf("Error, detected number of bubbles %i does not match expected %i\n",nblobs,(int)bubbles.size());
        N_errors++;
    }

    // Finished
    PROFILE_STOP("main");
    PROFILE_SAVE("TestBlobIdentify",false);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;  
}

