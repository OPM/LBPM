// Test of ComputeGlobalBlobIDs for special corner/edge cases
// Note: this is a short test, but requires 27 processors to run

#include <iostream>
#include <math.h>
#include "common/pmmc.h"
#include "common/Communication.h"
#include "analysis/analysis.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"




// Main
int main(int argc, char **argv)
{
    // Initialize MPI
    int rank, nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    /*if ( nprocs != 8 ) {
        printf("This tests requires 8 processors\n");
        return -1;
    }*/
    if ( rank==0 ) {
        printf("-----------------------------------------------------------\n");
        printf("Testing Blob Identification Special Cases\n");
        printf("-----------------------------------------------------------\n");
    }

    // Set the domain information
    std::vector<int> factors = Utilities::factor(nprocs);
    int nproc[3]={1,1,1};
    for (size_t i=0; i<factors.size(); i++)
        nproc[i%3] *= factors[i];
    const int nx = 4;
    const int ny = 4;
    const int nz = 4;
    double Lx=1.0, Ly=1.0, Lz=1.0;

    // Get the rank info
    const RankInfoStruct rank_info(rank,nproc[0],nproc[1],nproc[2]);
    const int rankx = rank_info.ix;
    const int ranky = rank_info.jy;
    const int rankz = rank_info.kz;

    // Test that local blobs that only touch at edges/corners do not share global ids
    DoubleArray Phase(nx+2,ny+2,nz+2);
    DoubleArray SignDist(nx+2,ny+2,nz+2);
    Phase.fill(1);
    SignDist.fill(1);
    // First block
    SignDist(1,1,1)=-1;  SignDist(2,1,1)=-1;  SignDist(1,2,1)=-1;  SignDist(2,2,1)=-1;
    SignDist(1,1,2)=-1;  SignDist(2,1,2)=-1;  SignDist(1,2,2)=-1;  SignDist(2,2,2)=-1;
    // Second block
    SignDist(3,3,1)=-1;  SignDist(4,3,1)=-1;  SignDist(3,4,1)=-1;  SignDist(4,4,1)=-1;
    SignDist(3,3,2)=-1;  SignDist(4,3,2)=-1;  SignDist(3,4,2)=-1;  SignDist(4,4,2)=-1;
    // Third block
    SignDist(3,1,3)=-1;  SignDist(4,1,3)=-1;  SignDist(3,2,3)=-1;  SignDist(4,2,3)=-1;
    SignDist(3,1,4)=-1;  SignDist(4,1,4)=-1;  SignDist(3,2,4)=-1;  SignDist(4,2,4)=-1;
    // Fourth block
    SignDist(1,3,3)=-1;  SignDist(1,4,3)=-1;  SignDist(2,3,3)=-1;  SignDist(2,4,3)=-1;
    SignDist(1,3,4)=-1;  SignDist(1,4,4)=-1;  SignDist(2,3,4)=-1;  SignDist(2,4,4)=-1;

    // Communication the halos
    fillHalo<double> fillData(rank_info,nx,ny,nz,1,1,1,0,1);
    fillData.fill(Phase);
    fillData.fill(SignDist);

    // Find blob domains
    if ( rank==0 ) { printf("Finding blob domains\n"); }
    double vF=0.0;
    double vS=0.0;
    IntArray LocalBlobID, GlobalBlobID;
    int nblobs0 = ComputeLocalBlobIDs(Phase,SignDist,vF,vS,LocalBlobID,false);
    int nblobs = ComputeGlobalBlobIDs(nx,ny,nz,rank_info,
        Phase,SignDist,vF,vS,GlobalBlobID);
    if ( rank==0 ) { printf("Identified %i blobs\n",nblobs); }
    
    // Create the MeshDataStruct
    std::vector<IO::MeshDataStruct> meshData(1);
    meshData[0].meshName = "domain";
    meshData[0].mesh = std::shared_ptr<IO::DomainMesh>( new IO::DomainMesh(rank_info,nx,ny,nz,Lx,Ly,Lz) );
    std::shared_ptr<IO::Variable> PhaseVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> SignDistVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> LocalBlobIDVar( new IO::Variable() );
    std::shared_ptr<IO::Variable> GlobalBlobIDVar( new IO::Variable() );
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
    LocalBlobIDVar->name = "LocalBlobID";
    LocalBlobIDVar->type = IO::VolumeVariable;
    LocalBlobIDVar->dim = 1;
    LocalBlobIDVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(LocalBlobIDVar);
    GlobalBlobIDVar->name = "GlobalBlobID";
    GlobalBlobIDVar->type = IO::VolumeVariable;
    GlobalBlobIDVar->dim = 1;
    GlobalBlobIDVar->data.resize(nx,ny,nz);
    meshData[0].vars.push_back(GlobalBlobIDVar);

    // Save the results
    fillData.copy(Phase,PhaseVar->data);
    fillData.copy(SignDist,SignDistVar->data);
    fillData.copy(LocalBlobID,LocalBlobIDVar->data);
    fillData.copy(GlobalBlobID,GlobalBlobIDVar->data);
    IO::writeData( 0, meshData );


    int N_errors = 0;
    return N_errors;  
}

