// Compute the signed distance from a digitized image 
// Two phases are present
// Phase 1 has value -1
// Phase 2 has value 1
// this code uses the segmented image to generate the signed distance 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "common/Array.h"
#include "common/Domain.h"
#include "IO/Writer.h"
#include "analysis/distance.h"


std::shared_ptr<Database> loadInputs( int nprocs )
{
    std::vector<int> nproc;
    if ( nprocs == 1 ) {
        nproc = { 1, 1, 1 };
    } else if ( nprocs == 8 ) {
        nproc = { 2, 2, 2 };
    } else {
        ERROR("TestSegDist: Unsupported number of processors");
    }
    auto db = std::make_shared<Database>( );
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", nproc );
    db->putVector<int>( "n", { 200, 200, 200 } );
    db->putScalar<int>( "nspheres", 0 );
    db->putVector<double>( "L", { 1, 1, 1 } );
    return db;
}


//***************************************************************************************
int main(int argc, char **argv)
{
    // Initialize MPI
    int rank, nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&nprocs);
    {


    // Load inputs
    auto db = loadInputs( nprocs );
    int Nx = db->getVector<int>( "n" )[0];
    int Ny = db->getVector<int>( "n" )[1];
    int Nz = db->getVector<int>( "n" )[2];


    // Get the rank info
    Domain Dm(db,comm);
    for (int k=0;k<Nz;k++){
        for (int j=0;j<Ny;j++){
            for (int i=0;i<Nx;i++){
                int n = k*Nx*Ny+j*Nx+i;
                Dm.id[n] = 1;
            }
        }
    }
    Dm.CommInit();

    int nx = Nx+2;
    int ny = Ny+2;
    int nz = Nz+2;

    // Initialize the bubble
    double BubbleRadius = 0.15*Nx*Dm.nprocx();
    double Cx = 0.40*Nx*Dm.nprocx();
    double Cy = 0.45*Nx*Dm.nprocy();
    double Cz = 0.50*Nx*Dm.nprocy();

    DoubleArray TrueDist(nx,ny,nz);
    Array<char> id(nx,ny,nz);
    id.fill(0);

    for (int k=1; k<nz-1; k++) {
        double z = k - 0.5 + Dm.kproc()*Nz;
        for (int j=1; j<ny-1; j++) {
            double y = j - 0.5 + Dm.jproc()*Ny;
            for (int i=1; i<nx-1; i++) {
                double x = i - 0.5 + Dm.iproc()*Nx;
                // True signed distance
                TrueDist(i,j,k) = sqrt((x-Cx)*(x-Cx)+(y-Cy)*(y-Cy)+(z-Cz)*(z-Cz)) - BubbleRadius;
                // Initialize phase positions
                if (TrueDist(i,j,k) < 0.0){
                    id(i,j,k) = 0;
                } else{
                    id(i,j,k)=1;
                }
            }
        }
    }

    MPI_Barrier(comm);
    if (rank==0) printf("Initialized! Converting to Signed Distance function \n");

    double t1 = MPI_Wtime();
    DoubleArray Distance(nx,ny,nz);
    CalcDist(Distance,id,Dm,{false,false,false});
    double t2 = MPI_Wtime();
    if (rank==0)
        printf("Total time: %f seconds \n",t2-t1);

    double err = 0.0;
    for (int i=1; i<Nx-1; i++) {
        for (int j=1; j<Ny-1; j++) {
            for (int k=1; k<Nz-1; k++) {
                err += (Distance(i,j,k)-TrueDist(i,j,k))*(Distance(i,j,k)-TrueDist(i,j,k));
            }
        }
    }
    err = sumReduce( Dm.Comm, err );
    err = sqrt( err / (nx*ny*nz*nprocs) );
    if (rank==0)
        printf("Mean error %0.4f \n", err);

    // Write the results
    Array<int> ID0(id.size());
    ID0.copy( id );
    Array<double> ID(Nx,Ny,Nz);
    Array<double> dist1(Nx,Ny,Nz);
    Array<double> dist2(Nx,Ny,Nz);
    fillHalo<double> fillData(Dm.Comm, Dm.rank_info,{Nx,Ny,Nz},{1,1,1},0,1);
    fillData.copy( ID0, ID );
    fillData.copy( TrueDist, dist1 );
    fillData.copy( Distance, dist2 );
    std::vector<IO::MeshDataStruct> data(1);
    data[0].meshName = "mesh";
    data[0].mesh.reset( new IO::DomainMesh( Dm.rank_info, Nx, Ny, Nz, Dm.Lx, Dm.Ly, Dm.Lz ) );
    data[0].vars.emplace_back( new IO::Variable( 1, IO::VariableType::VolumeVariable, "ID", ID ) );
    data[0].vars.emplace_back( new IO::Variable( 1, IO::VariableType::VolumeVariable, "TrueDist", dist1 ) );
    data[0].vars.emplace_back( new IO::Variable( 1, IO::VariableType::VolumeVariable, "Distance", dist2 ) );
    data[0].vars.emplace_back( new IO::Variable( 1, IO::VariableType::VolumeVariable, "error", dist2-dist1 ) );
    IO::initialize( "", "silo", false );
    IO::writeData( "testSegDist", data, MPI_COMM_WORLD );

    }
    MPI_Barrier(comm);
    MPI_Finalize();
    return 0;

}
