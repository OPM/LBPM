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

    //std::vector<std::string> filenames;
    std::string filename;
    if ( argc==0 ) {
        printf("At least one filename must be specified\n");
        return 1;
    }
    else {
      filename=std::string(argv[1]);
      printf("Reading file: %s\n",filename.c_str());
    }
						

   int fid = netcdf::open(filename);

    // Read all of the attributes
    std::vector<std::string> attr = netcdf::getAttNames( fid );
    for (size_t i=0; i<attr.size(); i++) {
        printf("Reading attribute %s\n",attr[i].c_str());
        netcdf::VariableType type = netcdf::getAttType( fid, attr[i] );
        if ( type == netcdf::STRING ){
	  Array<std::string> tmp = netcdf::getAtt<std::string>( fid, attr[i] );
	}
        else{
	  //Array<double> tmp = netcdf::getAtt<double>( fid, attr[i] );
	}
    }

    // Read the data array
    std::string varname("VOLUME");
    printf("Reading %s\n",varname.c_str());
    Array<float> VOLUME = netcdf::getVar<float>( fid, varname);
    //printf("VOLUME size = %zu \n",VOLUME.length());	
    printf("VOLUME dims =  %zu x %zu x %zu \n",VOLUME.size(0),VOLUME.size(1),VOLUME.size(2));	
    printf("Sucess!! \n");
    netcdf::close( fid );

    PROFILE_SAVE("TestNetcdf");
    return 0;  
}

