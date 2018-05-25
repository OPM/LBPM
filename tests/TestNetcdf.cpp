// Test reading/writing netcdf files

#include "IO/netcdf.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"
#include "common/UnitTest.h"

#include "ProfilerApp.h"


void load( const std::string& );


void test_NETCDF( UnitTest& ut )
{
    const int rank = comm_rank( MPI_COMM_WORLD );
    const int size = comm_size( MPI_COMM_WORLD );
    int nprocx = 2;
    int nprocy = 2;
    int nprocz = 2;
    RankInfoStruct info( rank, nprocx, nprocy, nprocz );
    Array<float> data( 4, 5, 6 );
    for (size_t i=0; i<data.length(); i++)
        data(i) = 120*rank + i;
    size_t x = info.ix*data.size(0);
    size_t y = info.jy*data.size(1);
    size_t z = info.kz*data.size(2);
    const char* filename = "test.nc";
    std::vector<int> dim = { (int) data.size(0)*nprocx, (int) data.size(1)*nprocy, (int) data.size(2)*nprocz };
    int fid = netcdf::open( filename, netcdf::CREATE, MPI_COMM_WORLD );
    auto dims =  netcdf::defDim( fid, {"X", "Y", "Z"}, dim );
    netcdf::write( fid, "tmp", dims, data, info );
    netcdf::close( fid );
    MPI_Barrier( MPI_COMM_WORLD );
    // Read the contents of the file we created
    fid = netcdf::open( filename, netcdf::READ );
    Array<float> tmp = netcdf::getVar<float>( fid, "tmp" );
    if ( (int)tmp.size(0)!=dim[0] || (int)tmp.size(1)!=dim[1] || (int)tmp.size(2)!=dim[2] ) {
        ut.failure("Array sizes do not match");
        return;
    }
    bool pass = true;
    for (size_t i=0; i<data.size(0); i++) {
        for (size_t j=0; j<data.size(1); j++) {
            for (size_t k=0; k<data.size(2); k++) {
                pass = pass && data(i,j,k) == tmp(i+x,j+y,k+z);
            }
        }
    }
    if ( pass ) {
        ut.passes("write/read simple parallel file");
    } else {
        ut.failure("write/read simple parallel file");
    }
}


inline void print( const std::string& name, const std::vector<size_t> size )
{
    printf("   Reading variable %s (%i",name.c_str(),(int)size[0]);
    for (size_t i=1; i<size.size(); i++)
        printf(",%i",(int)size[i]);
    printf(")\n");
}
void load( const std::string& filename )
{
    printf("Reading %s\n",filename.c_str());
    int fid = netcdf::open( filename, netcdf::READ );

    std::vector<std::string> vars = netcdf::getVarNames( fid );
    for (size_t i=0; i<vars.size(); i++) {
        netcdf::VariableType type = netcdf::getVarType( fid, vars[i] );
        print( vars[i], netcdf::getVarDim(fid,vars[i]) );
        if ( type == netcdf::STRING )
            Array<std::string> tmp = netcdf::getVar<std::string>( fid, vars[i] );
        else if ( type == netcdf::SHORT )
            Array<short> tmp = netcdf::getVar<short>( fid, vars[i] );
        else
            Array<double> tmp = netcdf::getVar<double>( fid, vars[i] );
    }

    std::vector<std::string> attr = netcdf::getAttNames( fid );
    for (size_t i=0; i<attr.size(); i++) {
        printf("   Reading attribute %s\n",attr[i].c_str());
        netcdf::VariableType type = netcdf::getAttType( fid, attr[i] );
        if ( type == netcdf::STRING )
            Array<std::string> tmp = netcdf::getAtt<std::string>( fid, attr[i] );
        else
            Array<double> tmp = netcdf::getAtt<double>( fid, attr[i] );
    }
    netcdf::close( fid );
}


int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc,&argv);
    int rank = comm_rank(MPI_COMM_WORLD);
    UnitTest ut;
    PROFILE_START("Main");

    // Test reading existing netcdf files
    if ( rank==0 ) {
        for (int i=1; i<argc; i++)
            load( argv[i] );
    }

    // Test writing/reading netcdf file
    test_NETCDF( ut );

    // Print the results
    ut.report();
    int N_errors = ut.NumFailGlobal();
    PROFILE_SAVE("TestNetcdf");

    // Close MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return N_errors;
}

