/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "IO/netcdf.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"

#include "ProfilerApp.h"


#ifdef USE_NETCDF


#include <netcdf.h>
#include <netcdf_par.h>


#define CHECK_NC_ERR( ERR )             \
    do {                                \
        if ( ERR != NC_NOERR ) {        \
            std::string msg = "Error calling netcdf routine: "; \
            msg += nc_strerror( ERR );  \
            ERROR( msg );               \
        }                               \
    } while (0)


namespace netcdf {


// Convert nc_type to VariableType
static inline VariableType convertType( nc_type type )
{
    VariableType type2 = UNKNOWN;
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


// Get nc_type from the template
template<class T> inline nc_type getType();
template<> inline nc_type getType<char>()   { return NC_CHAR; }
template<> inline nc_type getType<short>()  { return NC_SHORT; }
template<> inline nc_type getType<int>()    { return NC_INT; }
template<> inline nc_type getType<float>()  { return NC_FLOAT; }
template<> inline nc_type getType<double>() { return NC_DOUBLE; }


// Function to reverse an array
template<class TYPE>
inline std::vector<TYPE> reverse( const std::vector<TYPE>& x )
{
    std::vector<TYPE> y(x.size());
    for (size_t i=0; i<x.size(); i++)
        y[i] = x[x.size()-i-1];
    return y;
}
// Function to reverse an array
template<class TYPE1, class TYPE2>
inline std::vector<TYPE2> convert( const std::vector<TYPE1>& x )
{
    std::vector<TYPE2> y(x.size());
    for (size_t i=0; i<x.size(); i++)
        y[i] = static_cast<TYPE2>(x[i]);
    return y;
}


/****************************************************
* Convert the VariableType to a string              *
****************************************************/
std::string VariableTypeName( VariableType type )
{
    if ( type == BYTE )
        return "BYTE";
    else if ( type == SHORT )
        return "SHORT";
    else if ( type == USHORT )
        return "USHORT";
    else if ( type == INT )
        return "INT";
    else if ( type == UINT )
        return "UINT";
    else if ( type == INT64 )
        return "INT64";
    else if ( type == UINT64 )
        return "UINT64";
    else if ( type == FLOAT )
        return "FLOAT";
    else if ( type == DOUBLE )
        return "DOUBLE";
    else if ( type == STRING )
        return "STRING";
    return "Unknown";
}


/****************************************************
* Open/close a file                                 *
****************************************************/
int open( const std::string& filename, FileMode mode, MPI_Comm comm )
{
    int fid = 0;
    if ( comm == MPI_COMM_NULL ) {
        if ( mode == READ ) {
            int err = nc_open( filename.c_str(), NC_NOWRITE, &fid );
            CHECK_NC_ERR( err );
        } else if ( mode == WRITE ) {
            int err = nc_open( filename.c_str(), NC_WRITE, &fid );
            CHECK_NC_ERR( err );
        } else if ( mode == CREATE ) {
            int err = nc_create( filename.c_str(), NC_SHARE|NC_64BIT_OFFSET, &fid );
            CHECK_NC_ERR( err );
        } else {
            ERROR("Unknown file mode");
        }
    } else {
        if ( mode == READ ) {
            int err = nc_open_par( filename.c_str(), NC_MPIPOSIX, comm, MPI_INFO_NULL, &fid );
            CHECK_NC_ERR( err );
        } else if ( mode == WRITE ) {
            int err = nc_open_par( filename.c_str(), NC_WRITE|NC_MPIPOSIX, comm, MPI_INFO_NULL, &fid );
            CHECK_NC_ERR( err );
        } else if ( mode == CREATE ) {
            int err = nc_create_par( filename.c_str(), NC_NETCDF4|NC_MPIIO, comm, MPI_INFO_NULL, &fid );
            CHECK_NC_ERR( err );
        } else {
            ERROR("Unknown file mode");
        }
    }
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
VariableType getVarType( int fid, const std::string& var )
{
    int varid = -1;
    int err = nc_inq_varid( fid, var.c_str(), &varid );
    CHECK_NC_ERR( err );
    nc_type type=0;
    err = nc_inq_vartype( fid, varid, &type );
    CHECK_NC_ERR( err );
    return convertType(type);
}
VariableType getAttType( int fid, const std::string& att )
{
    nc_type type=0;
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
    int err = nc_get_var_ushort( fid, getVarID(fid,var), x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<unsigned short>");
    return x.reverseDim();
}
template<>
Array<short> getVar<short>( int fid, const std::string& var )
{
    PROFILE_START("getVar<short>");
    Array<short> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_short( fid, getVarID(fid,var), x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<short>");
    return x.reverseDim();
}
template<>
Array<unsigned int> getVar<unsigned int>( int fid, const std::string& var )
{
    PROFILE_START("getVar<unsigned int>");
    Array<unsigned int> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_uint( fid, getVarID(fid,var), x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<unsigned int>");
    return x.reverseDim();
}
template<>
Array<int> getVar<int>( int fid, const std::string& var )
{
    PROFILE_START("getVar<int>");
    Array<int> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_int( fid, getVarID(fid,var), x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<int>");
    return x.reverseDim();
}
template<>
Array<float> getVar<float>( int fid, const std::string& var )
{
    PROFILE_START("getVar<float>");
    Array<float> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_float( fid, getVarID(fid,var), x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<float>");
    return x.reverseDim();
}
template<>
Array<double> getVar<double>( int fid, const std::string& var )
{
    PROFILE_START("getVar<double>");
    Array<double> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_double( fid, getVarID(fid,var), x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<double>");
    return x.reverseDim();
}
template<>
Array<char> getVar<char>( int fid, const std::string& var )
{ 
    PROFILE_START("getVar<char>");
    Array<char> x( reverse(getVarDim(fid,var)) );
    int err = nc_get_var_text( fid, getVarID(fid,var), x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<char>");
    return x.reverseDim();
}
template<>
Array<std::string> getVar<std::string>( int fid, const std::string& var )
{
    PROFILE_START("getVar<std::string>");
    Array<char> tmp = getVar<char>( fid, var );
    std::vector<size_t> dim = {tmp.size(0), tmp.size(1), tmp.size(2) };
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
static inline void get_stride_args( const std::vector<int>& start,
    const std::vector<int>& count, const std::vector<int>& stride,
    size_t *startp, size_t *countp, ptrdiff_t *stridep )
{
    for (size_t i=0; i<start.size(); i++)
        startp[i] = start[i];
    for (size_t i=0; i<count.size(); i++)
        countp[i] = count[i];
    for (size_t i=0; i<stride.size(); i++)
        stridep[i] = stride[i];
}
template<class TYPE>
int nc_get_vars_TYPE( int fid, int varid, const size_t start[],
    const size_t count[], const ptrdiff_t stride[], TYPE *ptr );
template<>
int nc_get_vars_TYPE<short>( int fid, int varid, const size_t start[],
    const size_t count[], const ptrdiff_t stride[], short *ptr )
{
    return nc_get_vars_short( fid, varid, start, count, stride, ptr );
}
template<>
int nc_get_vars_TYPE<int>( int fid, int varid, const size_t start[],
    const size_t count[], const ptrdiff_t stride[], int *ptr )
{
    return nc_get_vars_int( fid, varid, start, count, stride, ptr );
}
template<>
int nc_get_vars_TYPE<float>( int fid, int varid, const size_t start[],
    const size_t count[], const ptrdiff_t stride[], float *ptr )
{
    return nc_get_vars_float( fid, varid, start, count, stride, ptr );
}
template<>
int nc_get_vars_TYPE<double>( int fid, int varid, const size_t start[],
    const size_t count[], const ptrdiff_t stride[], double *ptr )
{
    return nc_get_vars_double( fid, varid, start, count, stride, ptr );
}
template<class TYPE>
Array<TYPE> getVar( int fid, const std::string& var, const std::vector<int>& start,
    const std::vector<int>& count, const std::vector<int>& stride )
{
    PROFILE_START("getVar<> (strided)");
    std::vector<size_t> var_size = getVarDim( fid, var );
    for (int d=0; d<(int)var_size.size(); d++) {
        if ( start[d]<0 || start[d]+stride[d]*(count[d]-1)>(int)var_size[d] ) {
            int rank = comm_rank(MPI_COMM_WORLD);
            char tmp[1000];
            sprintf(tmp,"%i: Range exceeded array dimension:\n"
                "   start[%i]=%i, count[%i]=%i, stride[%i]=%i, var_size[%i]=%i",
                rank,d,start[d],d,count[d],d,stride[d],d,(int)var_size[d]);
            ERROR(tmp);
        }
    }
    Array<TYPE> x( reverse(convert<int,size_t>(count)) );
    size_t startp[10], countp[10];
    ptrdiff_t stridep[10];
    get_stride_args( start, count, stride, startp, countp, stridep );
    int err = nc_get_vars_TYPE<TYPE>( fid, getVarID(fid,var), startp, countp, stridep, x.data() );
    CHECK_NC_ERR( err );
    PROFILE_STOP("getVar<> (strided)");
    return x.reverseDim();
}
template Array<short>  getVar<short>(  int, const std::string&, const std::vector<int>&, const std::vector<int>&, const std::vector<int>& );
template Array<int>    getVar<int>(    int, const std::string&, const std::vector<int>&, const std::vector<int>&, const std::vector<int>& );
template Array<float>  getVar<float>(  int, const std::string&, const std::vector<int>&, const std::vector<int>&, const std::vector<int>& );
template Array<double> getVar<double>( int, const std::string&, const std::vector<int>&, const std::vector<int>&, const std::vector<int>& );


/****************************************************
* Read an attribute                                 *
****************************************************/
template<>
Array<double> getAtt<double>( int fid, const std::string& att )
{
    PROFILE_START("getAtt<double>");
    Array<double> x( getAttDim(fid,att) );
    int err = nc_get_att_double( fid, NC_GLOBAL, att.c_str(), x.data() );
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


/****************************************************
* Write an array to a file                          *
****************************************************/
std::vector<int> defDim( int fid, const std::vector<std::string>& names, const std::vector<int>& dims )
{
    std::vector<int> dimid(names.size(),0);
    for (size_t i=0; i<names.size(); i++) {
        int err = nc_def_dim( fid, names[i].c_str(), dims[i], &dimid[i]);
        CHECK_NC_ERR( err );
    }
    return dimid;
}
template<class TYPE>
void write( int fid, const std::string& var, const std::vector<int>& dimids,
    const Array<TYPE>& data, const RankInfoStruct& info )
{
    // Define the variable
    int varid = 0;
    int err = nc_def_var( fid, var.c_str(), getType<TYPE>(), data.ndim(), dimids.data(), &varid );
    CHECK_NC_ERR( err );
    // exit define mode 
    err = nc_enddef( fid );
    CHECK_NC_ERR( err );
    // set the access method to use MPI/PnetCDF collective I/O 
    err = nc_var_par_access( fid, varid, NC_INDEPENDENT );
    CHECK_NC_ERR( err );
    // parallel write: each process writes its subarray to the file
    auto x = data.reverseDim();
    std::vector<size_t> count = { data.size(0), data.size(1), data.size(2) };
    std::vector<size_t> start = { info.ix*data.size(0), info.jy*data.size(1), info.kz*data.size(2) };
    nc_put_vara( fid, varid, start.data(), count.data(), x.data() );
}
template void write<short>(  int fid, const std::string& var, const std::vector<int>& dimids, const Array<short>& data,  const RankInfoStruct& info );
template void write<int>(    int fid, const std::string& var, const std::vector<int>& dimids, const Array<int>& data,    const RankInfoStruct& info );
template void write<float>(  int fid, const std::string& var, const std::vector<int>& dimids, const Array<float>& data,  const RankInfoStruct& info );
template void write<double>( int fid, const std::string& var, const std::vector<int>& dimids, const Array<double>& data, const RankInfoStruct& info );



}; // netcdf namespace

#else

#endif


