// Test reading high-resolution files from the microct database

#include "common/MPI.h"
#include "common/UnitTest.h"
#include "common/Database.h"
#include "common/Domain.h"
#include "common/ReadMicroCT.h"
#include "IO/Writer.h"

#include <memory>



void testReadMicroCT( const std::string& filename, UnitTest& ut )
{
    Utilities::MPI comm( MPI_COMM_WORLD );

    // Get the domain info
    auto db = std::make_shared<Database>( filename );
    auto domain_db = db->getDatabase( "Domain" );

    // Test reading microCT files
    auto data = readMicroCT( *domain_db, comm );
    
    // Check if we loaded the data correctly
    if ( data.size() == domain_db->getVector<size_t>( "n" ) )
        ut.passes( "Read data" );
    else
        ut.passes( "Data size does not match" );

    // Write the results to silo
    auto n = domain_db->getVector<int>( "n" );
    auto nproc = domain_db->getVector<int>( "nproc" );
    int N[3] = { n[0]*nproc[0], n[1]*nproc[1], n[2]*nproc[2] };
    int rank = comm.getRank();
    RankInfoStruct rankInfo( rank, nproc[0], nproc[1], nproc[2] );
    std::vector<IO::MeshDataStruct> meshData( 1 );
    auto Var = std::make_shared<IO::Variable>();
    Var->name = "Variable";
    Var->type = IO::VariableType::VolumeVariable;
    Var->dim = 1;
    Var->data.copy( data );
    meshData[0].meshName = "grid";
    meshData[0].mesh = std::make_shared<IO::DomainMesh>(rankInfo,n[0],n[1],n[2],N[0],N[1],N[2]);
    meshData[0].vars.push_back(Var);
    IO::writeData( 0, meshData, comm );
}


int main(int argc, char **argv)
{
    // Initialize MPI
    Utilities::startup( argc, argv );
    UnitTest ut;

    // Run the tests
    ASSERT( argc == 2 );
    testReadMicroCT( argv[1], ut );

    // Print the results
    ut.report();
    int N_errors = ut.NumFailGlobal();

    // Close MPI
    Utilities::shutdown();
    return N_errors;
}

