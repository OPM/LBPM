#ifndef NETCDF_READER
#define NETCDF_READER

#include <string>
#include <vector>
#include <array>

#include "common/Array.h"
#include "common/MPI_Helpers.h"
#include "common/Communication.h"


#ifdef USE_SILO
    #include <silo.h>
#else
    typedef int DBfile;
#endif



namespace silo {


enum FileMode { READ, WRITE, CREATE };

enum VariableType { NodeVariable=1, EdgeVariable=2, SurfaceVariable=2, VolumeVariable=3, NullVariable=0 };


/*!
 * @brief  Open silo file
 * @detailed  This function opens a silo file
 * @param[in] filename      File to open
 * @param[in] mode          Open the file for reading or writing
 * @return This function returns a handle to the file
*/
DBfile* open( const std::string& filename, FileMode mode );


/*!
 * @brief  Close silo file
 * @detailed  This function closes a silo file
 * @param[in] fid           Handle to the open file
*/
void close( DBfile* fid );


/*!
 * @brief  Write a uniform grid
 * @detailed  This function writes a uniform grid to silo as a Quadmesh
 * @return This function returns a handle to the file
 * @param[in] meshname      Mesh name
 * @param[in] range         Range of mesh { xmin, xmax, ymin, ymax, zmin, zmax }
 * @param[in] N             Number of cells in each direction
*/
template<int NDIM>
void writeUniformMesh( DBfile* fid, const std::string& meshname,
    const std::array<double,2*NDIM>& range, const std::array<int,NDIM>& N );


/*!
 * @brief  Write a multimesh
 * @detailed  This function writes a multimesh to silo
 * @return This function returns a handle to the file
 * @param[in] meshname      Mesh name
 * @param[in] subMeshNames  Names of the sub meshes in the form "filename:meshname"
 * @param[in] subMeshTypes  Type of each submesh
*/
void writeMultiMesh( DBfile* fid, const std::string& meshname,
    const std::vector<std::string>& subMeshNames,
    const std::vector<int>& subMeshTypes );


/*!
 * @brief  Write a uniform grid
 * @detailed  This function writes a uniform grid to silo as a Quadmesh
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] N             Number of cells in each direction
 * @param[in] varname       Variable name
 * @param[in] data          Variable data
 * @param[in] type          Variable type
*/
template<int NDIM>
void writeUniformMeshVariable( DBfile* fid, const std::string& meshname, const std::array<int,NDIM>& N,
    const std::string& varname, const Array<double>& data, VariableType type );


/*!
 * @brief  Write a multivariable
 * @detailed  This function writes a multivariable to silo
 * @return This function returns a handle to the file
 * @param[in] fid           Handle to the open file
 * @param[in] varname       Mesh name
 * @param[in] subVarNames   Names of the sub meshes in the form "filename:meshname"
 * @param[in] subVarTypes   Type of each submesh
 * @param[in] ndim          Dimension of variable (used to determine suffix)
 * @param[in] nvar          Number of subvariables (used to determine suffix)
*/
void writeMultiVar( DBfile* fid, const std::string& varname,
    const std::vector<std::string>& subVarNames,
    const std::vector<int>& subVarTypes, int ndim, int nvar );

}; // silo namespace
#endif
