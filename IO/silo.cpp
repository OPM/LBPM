#include "IO/silo.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"

#include "ProfilerApp.h"


#ifdef USE_SILO

#include <silo.h>



namespace silo {


/****************************************************
* Open/close a file                                 *
****************************************************/
DBfile* open( const std::string& filename, FileMode mode )
{
    DBfile *fid = nullptr;
    if ( mode == CREATE ) {
        fid = DBCreate( filename.c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_HDF5 );
    } else if ( mode == WRITE ) {
        fid = DBOpen( filename.c_str(), DB_HDF5, DB_APPEND );
    } else if ( mode == READ ) {
        fid = DBOpen( filename.c_str(), DB_HDF5, DB_READ );
    }
    return fid;
}
void close( DBfile* fid )
{
    DBClose( fid );
}


/****************************************************
* Write a uniform mesh to silo                      *
****************************************************/
template<int NDIM>
void writeUniformMesh( DBfile* fid, const std::string& meshname,
    const std::array<double,2*NDIM>& range, const std::array<int,NDIM>& N )
{
    PROFILE_START("writeUniformMesh",2);
    int dims[NDIM];
    for (size_t d=0; d<N.size(); d++)
        dims[d] = N[d]+1;
    float *x = nullptr;
    if ( NDIM >= 1 ) {
        x = new float[dims[0]];
        for (int i=0; i<N[0]; i++)
            x[i] = range[0] + i*(range[1]-range[0])/N[0];
        x[N[0]] = range[1];
    }
    float *y = nullptr;
    if ( NDIM >= 2 ) {
        y = new float[dims[1]];
        for (int i=0; i<N[1]; i++)
            y[i] = range[2] + i*(range[3]-range[2])/N[1];
        y[N[1]] = range[3];
    }
    float *z = nullptr;
    if ( NDIM >= 3 ) {
        z = new float[dims[2]];
        for (int i=0; i<N[2]; i++)
            z[i] = range[4] + i*(range[5]-range[4])/N[2];
        z[N[2]] = range[5];
    }
    float *coords[] = { x, y, z };
    int err = DBPutQuadmesh( fid, meshname.c_str(), nullptr, coords, dims, NDIM, DB_FLOAT, DB_COLLINEAR, nullptr );
    ASSERT( err == 0 );
    PROFILE_STOP("writeUniformMesh",2);
}


/****************************************************
* Write a vector/tensor quad variable               *
****************************************************/
std::vector<std::string> getVarSuffix( int ndim, int nvars )
{
    std::vector<std::string> suffix(nvars);
    if ( nvars == 1 ) {
        suffix[0] = "";
    } else if ( nvars == ndim ) {
        if ( ndim==2 ) {
            suffix[0] = "_x";
            suffix[1] = "_y";
        } else if ( ndim==3 ) {
            suffix[0] = "_x";
            suffix[1] = "_y";
            suffix[2] = "_z";
        } else {
            ERROR("Not finished");
        }
    } else if ( nvars == ndim*ndim ) {
        if ( ndim==2 ) {
            suffix[0] = "_xx";
            suffix[1] = "_xy";
            suffix[2] = "_yx";
            suffix[3] = "_yy";
        } else if ( ndim==3 ) {
            suffix[0] = "_xx";
            suffix[1] = "_xy";
            suffix[2] = "_xz";
            suffix[3] = "_yx";
            suffix[4] = "_yy";
            suffix[5] = "_yz";
            suffix[6] = "_zx";
            suffix[7] = "_zy";
            suffix[8] = "_zz";
        } else {
            ERROR("Not finished");
        }
    } else {
        for (int i=0; i<nvars; i++)
            suffix[i] = "_" + std::to_string(i+1);
    }
    return suffix;
}
template<int NDIM>
void writeUniformMeshVariable( DBfile* fid, const std::string& meshname, const std::array<int,NDIM>& N,
    const std::string& varname, const Array<double>& data, VariableType type )
{
    for (int d=0; d<NDIM; d++)
        ASSERT(N[d]==(int)data.size(d));
    PROFILE_START("writeUniformMeshVariable",2);
    int nvars=1, dims[NDIM]={1};
    const double *vars[NDIM] = { nullptr };
    if ( type == NodeVariable ) {
        ERROR("Not finished");
    } else if ( type == EdgeVariable ) {
        ERROR("Not finished");
    } else if ( type == SurfaceVariable ) {
        ERROR("Not finished");
    } else if ( type == VolumeVariable ) {
        ASSERT( data.ndim()==NDIM || data.ndim()==NDIM+1 );
        nvars = data.size(NDIM+1);
        size_t N = data.length()/nvars;
        for (int d=0; d<NDIM; d++)
            dims[d] = data.size(d);
        for (int i=0; i<nvars; i++)
            vars[i] = &data(i*N);
    } else {
        ERROR("Invalid variable type");
    }
    auto suffix = getVarSuffix( NDIM, nvars );
    std::vector<std::string> var_names(nvars);
    for (int i=0; i<nvars; i++)
        var_names[i] = varname + suffix[i];
    std::vector<char*> varnames(nvars,nullptr);
    for (int i=0; i<nvars; i++)
        varnames[i] = const_cast<char*>(var_names[i].c_str());
    int err = DBPutQuadvar( fid, varname.c_str(), meshname.c_str(), nvars,
        varnames.data(), vars, dims, NDIM, nullptr, 0, DB_DOUBLE, DB_ZONECENT, nullptr );
    ASSERT( err == 0 );
    PROFILE_STOP("writeUniformMeshVariable",2);
}


/****************************************************
* Write a multimesh                                 *
****************************************************/
void writeMultiMesh( DBfile* fid, const std::string& meshname,
    const std::vector<std::string>& meshNames,
    const std::vector<int>& meshTypes )
{
    std::vector<char*> meshnames(meshNames.size());
    for ( size_t i = 0; i < meshNames.size(); ++i )
        meshnames[i] = (char *) meshNames[i].c_str();
    std::string tree_name = meshname + "_tree";
    DBoptlist *optList    = DBMakeOptlist( 1 );
    DBAddOption( optList, DBOPT_MRGTREE_NAME, (char *) tree_name.c_str() );
    DBPutMultimesh( fid, meshname.c_str(), meshNames.size(), meshnames.data(), (int*) meshTypes.data(), nullptr );
    DBFreeOptlist( optList );
}


/****************************************************
* Write a multivariable                             *
****************************************************/
void writeMultiVar( DBfile* fid, const std::string& varname,
    const std::vector<std::string>& varNames,
    const std::vector<int>& varTypes, int ndim, int nvar )
{
    ASSERT(varNames.size()==varTypes.size());
    auto suffix = getVarSuffix( ndim, nvar );
    for (int i=0; i<nvar; i++) {
        auto varname2 = varname + suffix[i];
        auto varNames2 = varNames;
        for ( auto& tmp : varNames2 )
            tmp += suffix[i];
        std::vector<char*> varnames(varNames.size(),nullptr);
        for (size_t j=0; j<varNames.size(); j++)
            varnames[j] = const_cast<char*>(varNames2[j].c_str());
        DBPutMultivar( fid, varname2.c_str(), varNames.size(), varnames.data(), (int*) varTypes.data(), nullptr );
    }
}


}; // silo namespace


// Explicit instantiations
template void silo::writeUniformMesh<1>( DBfile*, const std::string&, const std::array<double,2>&, const std::array<int,1>& );
template void silo::writeUniformMesh<2>( DBfile*, const std::string&, const std::array<double,4>&, const std::array<int,2>& );
template void silo::writeUniformMesh<3>( DBfile*, const std::string&, const std::array<double,6>&, const std::array<int,3>& );
template void silo::writeUniformMeshVariable<1>( DBfile*, const std::string&, const std::array<int,1>&, const std::string&, const Array<double>&, silo::VariableType );
template void silo::writeUniformMeshVariable<2>( DBfile*, const std::string&, const std::array<int,2>&, const std::string&, const Array<double>&, silo::VariableType );
template void silo::writeUniformMeshVariable<3>( DBfile*, const std::string&, const std::array<int,3>&, const std::string&, const Array<double>&, silo::VariableType );

#else

#endif
