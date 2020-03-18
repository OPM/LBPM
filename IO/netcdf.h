#ifndef NETCDF_READER
#define NETCDF_READER

#include <string>
#include <vector>

#include "common/Array.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"



namespace netcdf {


//! Enum to hold variable type
enum VariableType { BYTE, SHORT, USHORT, INT, UINT, INT64, UINT64, FLOAT, DOUBLE, STRING, UNKNOWN };

//! Enum to hold variable type
enum FileMode { READ, WRITE, CREATE };


//! Convert the VariableType to a string
std::string VariableTypeName( VariableType type );


/*!
 * @brief  Open netcdf file
 * @detailed  This function opens a netcdf file
 * @return This function returns a handle to the file
 * @param filename      File to open
 * @param mode          Open the file for reading or writing
 * @param comm          MPI communicator to use (MPI_COMM_WORLD: don't use parallel netcdf)
*/
int open( const std::string& filename, FileMode mode, MPI_Comm comm=MPI_COMM_NULL );


/*!
 * @brief  Close netcdf file
 * @detailed  This function closes a netcdf file
 * @param fid           Handle to the open file
*/
void close( int fid );


/*!
 * @brief  Read the variable names
 * @detailed  This function reads a list of the variable names in the file
 * @param fid           Handle to the open file
*/
std::vector<std::string> getVarNames( int fid );


/*!
 * @brief  Read the attribute names
 * @detailed  This function reads a list of the attribute names in the file
 * @param fid           Handle to the open file
*/
std::vector<std::string> getAttNames( int fid );


/*!
 * @brief  Return the variable type
 * @detailed  This function returns the type for a variable
 * @param fid           Handle to the open file
 * @param var           Variable to read
*/
VariableType getVarType( int fid, const std::string& var );


/*!
 * @brief  Return the attribute type
 * @detailed  This function returns the type for an attribute
 * @param fid           Handle to the open file
 * @param att           Attribute to read
*/
VariableType getAttType( int fid, const std::string& att );


/*!
 * @brief  Return the variable dimensions
 * @detailed  This function returns the die for a variable
 * @param fid           Handle to the open file
 * @param var           Variable to read
*/
std::vector<size_t> getVarDim( int fid, const std::string& var );


/*!
 * @brief  Read a variable
 * @detailed  This function reads a variable with the given name from the file
 * @param fid           Handle to the open file
 * @param var           Variable to read
*/
template<class TYPE>
Array<TYPE> getVar( int fid, const std::string& var );


/*!
 * @brief  Read a strided variable
 * @detailed  This function reads a strided variable with the given name from the file
 * @param fid           Handle to the open file
 * @param var           Variable to read
 * @param start         Starting corner for the read
 * @param count         Number of elements to read
 * @param stride        Stride size for the read
*/
template<class TYPE>
Array<TYPE> getVar( int fid, const std::string& var, const std::vector<int>& start,
    const std::vector<int>& count, const std::vector<int>& stride );


/*!
 * @brief  Read an attribute
 * @detailed  This function reads an attribute with the given name from the file
 * @param fid           Handle to the open file
 * @param att           Attribute to read
*/
template<class TYPE>
Array<TYPE> getAtt( int fid, const std::string& att );


/*!
 * @brief  Write the dimensions
 * @detailed  This function writes the grid dimensions to netcdf. 
 * @param fid           Handle to the open file
*/
std::vector<int> defDim( int fid, const std::vector<std::string>& names, const std::vector<int>& dims );


/*!
 * @brief  Write a variable
 * @detailed  This function writes a variable to netcdf. 
 * @param fid           Handle to the open file
*/
template<class TYPE>
void write( int fid, const std::string& var, const std::vector<int>& dimids, const Array<TYPE>& data, const RankInfoStruct& rank_info );


}; // netcdf namespace
#endif
