#include "IO/netcdf.h"
#include "common/Utilities.h"

#include "ProfilerApp.h"


#ifdef USE_NETCDF


#include <netcdf.h>


#define CHECK_NC_ERR( ERR )             \
    do {                                \
        if ( ERR != NC_NOERR ) {        \
            std::string msg = "Error calling netcdf routine: "; \
            msg += nc_strerror( ERR );  \
            ERROR( msg );               \
        }                               \
    } while (0)


namespace netcdf {


// Function to reverse an array
template<class TYPE>
inline std::vector<TYPE> reverse( const std::vector<TYPE>& x )
{
    std::vector<TYPE> y(x.size());
    for (size_t i=0; i<x.size(); i++)
        y[i] = x[x.size()-i-1];
    return y;
}


/****************************************************
* Open/close a file                                 *
****************************************************/
int open( const std::string& filename )
{
    int fid = 0;
    int err = nc_open( filename.c_str(), NC_NOWRITE, &fid );
    CHECK_NC_ERR( err );
    return fid;
}
void close( int fid )
{
    int err = nc_close( fid );
    if ( err != NC_NOERR )
        ERROR("Error closing file");
}


/****************************************************
* Query basic properties                            *
****************************************************/
static std::vector<size_t> getDimVar( int fid, int varid )
{
    int ndim = 0;
    int err = nc_inq_varndims( fid, varid, &ndim );
    CHECK_NC_ERR( err );
    std::vector<size_t> dims(ndim,0);
    int dimid[64] = {-1};
    err = nc_inq_vardimid( fid, varid, dimid );
    CHECK_NC_ERR( err );
    for (int i=0; i<ndim; i++) {
        err = nc_inq_dimlen( fid, dimid[i], &dims[i] );
        CHECK_NC_ERR( err );
    }
    return dims;
}
static int getVarID( int fid, const std::string& var )
{
    int id = -1;
    int err = nc_inq_varid( fid, var.c_str(), &id );
    CHECK_NC_ERR( err );
    return id;
}
std::vector<size_t> getVarDim( int fid, const std::string& var )
{
    return getDimVar( fid, getVarID( fid, var ) );
}
std::vector<size_t> getAttDim( int fid, const std::string& att )
{
    std::vector<size_t> dim(1,0);
    int err = nc_inq_attlen( fid, NC_GLOBAL, att.c_str(), dim.data() );
    return dim;
}
std::vector<std::string> getVarNames( int fid )
{
    int nvar;
    int err = nc_inq( fid, NULL, &nvar, NULL, NULL );
    CHECK_NC_ERR( err );
    std::vector<std::string> vars(nvar);
    for (int i=0; i<nvar; i++) {
        char name[NC_MAX_NAME+1];
        err = nc_inq_varname( fid, i, name );
        CHECK_NC_ERR( err );
        vars[i] = name;
    }
    return vars;
}
std::vector<std::string> getAttNames( int fid )
{
    int natt;
    int err = nc_inq( fid, NULL, NULL, &natt, NULL );
    CHECK_NC_ERR( err );
    std::vector<std::string> att(natt);
    for (int i=0; i<natt; i++) {
        char name[NC_MAX_NAME+1];
        err = nc_inq_attname( fid,  NC_GLOBAL, i, name );
        CHECK_NC_ERR( err );
        att[i] = name;
    }
    return att;
}
static inline VariableType convertType( nc_type type )
{
    VariableType type2;
    if ( type == NC_BYTE )
        type2 = BYTE;
    else if ( type == NC_CHAR )
        type2 = STRING;
    else if ( type == NC_SHORT )
        type2 = SHORT;
    else if ( type == NC_USHORT )
        type2 = USHORT;
    else if ( type == NC_INT )
        type2 = INT;
    else if ( type == NC_UINT )
        type2 = UINT;
    else if ( type == NC_INT64 )
        type2 = INT64;
    else if ( type == NC_UINT64 )
        type2 = UINT64;
    else if ( type == NC_FLOAT )
        type2 = FLOAT;
    else if ( type == NC_DOUBLE )
        type2 = DOUBLE;
    else
        ERROR("Unknown type");
    return type2;
}
VariableType getVarType( int fid, const std::string& var )
{
    int varid = -1;
    int err = nc_inq_varid( fid, var.c_str(), &varid );
    CHECK_NC_ERR( err );
    nc_type type;
    err = nc_inq_vartype( fid, varid, &type );
    CHECK_NC_ERR( err );
    return convertType(type);
}
VariableType getAttType( int fid, const std::string& att )
{
    nc_type type;
    int err = nc_inq_atttype( fid,  NC_GLOBAL, att.c_str(), &type );
    CHECK_NC_ERR( err );
    return convertType(type);
}



/****************************************************
* Read a variable                                   *
****************************************************/
template<>
Array<unsigned short> getVar<unsigned short>( int fid, const std::string& var )
{
    PROFILE_START("getVar<unsigned short>");
    Array<unsigned short> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_ushort( fid, getVarID(fid,var), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<unsigned short>");
    return x;
}
template<>
Array<short> getVar<short>( int fid, const std::string& var )
{
    PROFILE_START("getVar<short>");
    Array<short> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_short( fid, getVarID(fid,var), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<short>");
    return x;
}
template<>
Array<unsigned int> getVar<unsigned int>( int fid, const std::string& var )
{
    PROFILE_START("getVar<unsigned int>");
    Array<unsigned int> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_uint( fid, getVarID(fid,var), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<unsigned int>");
    return x;
}
template<>
Array<int> getVar<int>( int fid, const std::string& var )
{
    PROFILE_START("getVar<int>");
    Array<int> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_int( fid, getVarID(fid,var), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<int>");
    return x;
}
template<>
Array<float> getVar<float>( int fid, const std::string& var )
{
    PROFILE_START("getVar<float>");
    Array<float> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_float( fid, getVarID(fid,var), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<float>");
    return x;
}
template<>
Array<double> getVar<double>( int fid, const std::string& var )
{
    PROFILE_START("getVar<double>");
    Array<double> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_double( fid, getVarID(fid,var), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<double>");
    return x;
}
template<>
Array<char> getVar<char>( int fid, const std::string& var )
{ 
    PROFILE_START("getVar<char>");
    Array<char> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_text( fid, getVarID(fid,var), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<char>");
    return x;
}
template<>
Array<std::string> getVar<std::string>( int fid, const std::string& var )
{
    PROFILE_START("getVar<std::string>");
    Array<char> tmp = getVar<char>( fid, var );
    std::vector<size_t> dim = tmp.size();
    if ( dim.size() == 1 )
        dim[0] = 1;
    else
        dim.erase( dim.begin() );
    Array<std::string> text(dim);
    for (size_t i=0; i<text.length(); i++)
        text(i) = &(tmp(0,i));
    PROFILE_STOP("getVar<std::string>");
    return text;
}


/****************************************************
* Read an attribute                                 *
****************************************************/
template<>
Array<double> getAtt<double>( int fid, const std::string& att )
{
    PROFILE_START("getAtt<double>");
    Array<double> x( getAttDim(fid,att) );
    int err = nc_get_att_double( fid, NC_GLOBAL, att.c_str(), x.get() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getAtt<double>");
    return x;
}
template<>
Array<std::string> getAtt<std::string>( int fid, const std::string& att )
{
    PROFILE_START("getAtt<std::string>");
    char *tmp = new char[getAttDim(fid,att)[0]];
    Array<std::string> x(1);
    x(0) = tmp;
    delete [] tmp;
    PROFILE_STOP("getAtt<std::string>");
    return x;
}


}; // netcdf namespace

#else

#endif

