// Sequential blob analysis 
// Reads parallel simulation data and performs connectivity analysis
// and averaging on a blob-by-blob basis
// James E. McClure 2014

#include "IO/netcdf.h"

#include "ProfilerApp.h"


void load( const std::string filename )
{
    int fid = netcdf::open( filename );

    std::vector<std::string> vars = netcdf::getVarNames( fid );
    for (size_t i=0; i<vars.size(); i++) {
        printf("Reading variable %s\n",vars[i].c_str());
        netcdf::VariableType type = netcdf::getVarType( fid, vars[i] );
        if ( type == netcdf::STRING )
            Array<std::string> tmp = netcdf::getVar<std::string>( fid, vars[i] );
        else if ( type == netcdf::SHORT )
            Array<short> tmp = netcdf::getVar<short>( fid, vars[i] );
        else
            Array<double> tmp = netcdf::getVar<double>( fid, vars[i] );
    }

    std::vector<std::string> attr = netcdf::getAttNames( fid );
    for (size_t i=0; i<attr.size(); i++) {
        printf("Reading attribute %s\n",attr[i].c_str());
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
    PROFILE_START("Main");

    std::vector<std::string> filenames;
    
    if ( argc==0 ) {
        printf("At least one filename must be specified\n");
        return 1;
    }

    for (int i=1; i<argc; i++)
        load( argv[i] );

    PROFILE_SAVE("TestNetcdf");
    return 0;  
}

