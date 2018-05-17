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
    INSIST(nprocs==8, "TestSegDist: Number of MPI processes must be equal to 8");
    auto db = std::make_shared<Database>( );
    db->putScalar<int>( "BC", 0 );
    db->putVector<int>( "nproc", { 2, 2, 2 } );
    db->putVector<int>( "n", { 100, 100, 100 } );
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
    Domain Dm(db);
    for (int k=0;k<Nz;k++){
        for (int j=0;j<Ny;j++){
            for (int i=0;i<Nx;i++){
                int n = k*Nx*Ny+j*Nx+i;
                Dm.id[n] = 1;
            }
        }
    }
    Dm.CommInit(comm);

    int nx = Nx+2;
    int ny = Ny+2;
    int nz = Nz+2;

    double BubbleRadius = 0.3*nx;

    // Initialize the bubble
    double Cx = 1.0*nx;
    double Cy = 1.0*ny;
    double Cz = 1.0*nz;

    DoubleArray TrueDist(nx,ny,nz);
    Array<char> id(nx,ny,nz);
    id.fill(0);

    for (int k=1;k<nz-1;k++){
        for (int j=1;j<ny-1;j++){
            for (int i=1;i<nx-1;i++){

                // True signed distance
                double x = (nx-2)*Dm.iproc()+i-1;
                double y = (ny-2)*Dm.jproc()+j-1;
                double z = (nz-2)*Dm.kproc()+k-1;
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

    double localError=0.0;
    int localCount = 0;
    for (int k=0;k<nz;k++){
        for (int j=0;j<ny;j++){
            for (int i=0;i<nx;i++){
                if (fabs(TrueDist(i,j,k)) < 3.0){
                    localError += (Distance(i,j,k)-TrueDist(i,j,k))*(Distance(i,j,k)-TrueDist(i,j,k));
                    localCount++;
                }
            }
        }
    }
    double globalError;
    int globalCount;
    MPI_Allreduce(&localError,&globalError,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&localCount,&globalCount,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    double err2 = sqrt(globalError)/(double (globalCount));
    if (rank==0) printf("Mean error %f \n", err2);

    // Write the results to visit
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
