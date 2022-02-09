#include <stdio.h>
#include <iostream>
#include "common/ScaLBL.h"
#include "common/MPI.h"


int main(int argc, char **argv)
{
    Utilities::startup( argc, argv, true );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    {
        auto filename = argv[1];
        auto input_db = std::make_shared<Database>( filename );
        auto db = input_db->getDatabase( "Domain" );
        auto Dm  = std::shared_ptr<Domain>(new Domain(db,comm));
        int Nx = db->getVector<int>( "n" )[0] + 2;
        int Ny = db->getVector<int>( "n" )[1] + 2;
        int Nz = db->getVector<int>( "n" )[2] + 2;
        char LocalRankString[8];
        sprintf(LocalRankString,"%05d",rank);
        char LocalRankFilename[40];
        sprintf(LocalRankFilename,"ID.%05i",rank);
        auto id = new char[Nx*Ny*Nz];
        for (int k=0;k<Nz;k++){
            for (int j=0;j<Ny;j++){
                for (int i=0;i<Nx;i++){
                    int n = k*Nx*Ny+j*Nx+i;
                    id[n] = 1;
                    Dm->id[n] = id[n];
                }
            }
        }
        Dm->CommInit();
        std::cout << "step 1" << std::endl << std::flush;
    }
    std::cout << "step 2" << std::endl << std::flush;
    comm.barrier();
    std::cout << "step 3" << std::endl << std::flush;
    Utilities::shutdown();
    std::cout << "step 4" << std::endl << std::flush;
    return 0;
}
