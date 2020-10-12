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
#ifndef SILO_INTERFACE
#define SILO_INTERFACE

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

enum class VariableType : int { NodeVariable=1, EdgeVariable=2, SurfaceVariable=2, VolumeVariable=3, NullVariable=0 };

enum class VariableDataType { DOUBLE, FLOAT, INT, UNKNOWN };


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
 * @brief  Get the variable type
 * @detailed  This function returns the type of variable data
 * @param[in] fid           Handle to the open file
 * @param[in] name          Name of variable
*/
VariableDataType varDataType( DBfile *dbfile, const std::string& name );


/*!
 * @brief  Write data to silo
 * @detailed  This function writes an arbitrary array to silo
 * @param[in] fid           Handle to the open file
 * @param[in] varname       Variable name
 * @param[in] data          Data to write
*/
template<class TYPE>
void write( DBfile* fid, const std::string& varname, const std::vector<TYPE>& data );


/*!
 * @brief  Write data to silo
 * @detailed  This function writes an arbitrary array to silo
 * @param[in] fid           Handle to the open file
 * @param[in] varname       Variable name
 * @return                  Data read
*/
template<class TYPE>
std::vector<TYPE> read( DBfile* fid, const std::string& varname );


/*!
 * @brief  Write a uniform grid
 * @detailed  This function writes a uniform grid to silo as a Quadmesh
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] range         Range of mesh { xmin, xmax, ymin, ymax, zmin, zmax }
 * @param[in] N             Number of cells in each direction
*/
template<int NDIM>
void writeUniformMesh( DBfile* fid, const std::string& meshname,
    const std::array<double,2*NDIM>& range, const std::array<int,NDIM>& N );


/*!
 * @brief  Read a uniform grid
 * @detailed  This function reads a uniform grid from silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[out] range         Range of mesh { xmin, xmax, ymin, ymax, zmin, zmax }
 * @param[out] N             Number of cells in each direction
*/
void readUniformMesh( DBfile* fid, const std::string& meshname,
    std::vector<double>& range, std::vector<int>& N );


/*!
 * @brief  Write a uniform grid variable
 * @detailed  This function writes a uniform grid variable to silo as a Quadmesh
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] N             Number of cells in each direction
 * @param[in] varname       Variable name
 * @param[in] data          Variable data
 * @param[in] type          Variable type
*/
template< int NDIM, class TYPE >
void writeUniformMeshVariable( DBfile* fid, const std::string& meshname, const std::array<int,NDIM>& N,
    const std::string& varname, const Array<TYPE>& data, VariableType type );


/*!
 * @brief  Read a uniform mesh grid variable
 * @detailed  This function read a uniform mesh variable to silo
 * @param[in] fid           Handle to the open file
 * @param[in] varname       Variable name
 * @return                  Variable data
*/
template<class TYPE>
Array<TYPE> readUniformMeshVariable( DBfile* fid, const std::string& varname );


/*!
 * @brief  Write a pointmesh
 * @detailed  This function writes a pointmesh to silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] ndim          Number of dimensions
 * @param[in] N             Number of points
 * @param[in] coords        Coordinates of the points
*/
template<class TYPE>
void writePointMesh( DBfile* fid, const std::string& meshname,
    int ndim, int N, const TYPE *coords[] );


/*!
 * @brief  Read a pointmesh
 * @detailed  This function reads a pointmesh from silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @return                  Returns the coordinates as a N x ndim array 
*/
template<class TYPE>
Array<TYPE> readPointMesh( DBfile* fid, const std::string& meshname );


/*!
 * @brief  Write a pointmesh grid variable
 * @detailed  This function writes a pointmesh variable to silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] varname       Variable name
 * @param[in] data          Variable data
*/
template<class TYPE>
void writePointMeshVariable( DBfile* fid, const std::string& meshname,
    const std::string& varname, const Array<TYPE>& data );


/*!
 * @brief  Read a pointmesh grid variable
 * @detailed  This function reads a pointmesh variable from silo
 * @param[in] fid           Handle to the open file
 * @param[in] varname       Variable name
 * @return                  Variable data
*/
template<class TYPE>
Array<TYPE> readPointMeshVariable( DBfile* fid, const std::string& varname );


/*!
 * @brief  Write a triangle mesh
 * @detailed  This function writes a triangle (or simplex) based mesh to silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] ndim          Number of dimensions for the coordinates
 * @param[in] ndim_tri      Number of dimensions for the triangles (2: surface, 3: volume)
 * @param[in] N             Number of points
 * @param[in] coords        Coordinates of the points
 * @param[in] N_tri         Number of triangles
 * @param[in] tri           Coordinates of the points
*/
template<class TYPE>
void writeTriMesh( DBfile* fid, const std::string& meshname,
    int ndim, int ndim_tri, int N, const TYPE *coords[], int N_tri, const int *tri[] );


/*!
 * @brief  Read a triangle mesh
 * @detailed  This function reads a triangle (or simplex) based mesh to silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] coords        Coordinates of the points
 * @param[in] tri           Coordinates of the points
*/
template<class TYPE>
void readTriMesh( DBfile* fid, const std::string& meshname, Array<TYPE>& coords, Array<int>& tri );


/*!
 * @brief  Write a triangle mesh grid variable
 * @detailed  This function writes a triangle mesh variable to silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] varname       Variable name
 * @param[in] data          Variable data
 * @param[in] type          Variable type
*/
template<class TYPE>
void writeTriMeshVariable( DBfile* fid, int ndim, const std::string& meshname,
    const std::string& varname, const Array<TYPE>& data, VariableType type );


/*!
 * @brief  Read a triangle mesh grid variable
 * @detailed  This function read a triangle mesh variable to silo
 * @param[in] fid           Handle to the open file
 * @param[in] varname       Variable name
 * @return                  Variable data
*/
template<class TYPE>
Array<TYPE> readTriMeshVariable( DBfile* fid, const std::string& varname );


/*!
 * @brief  Write a multimesh
 * @detailed  This function writes a multimesh to silo
 * @param[in] fid           Handle to the open file
 * @param[in] meshname      Mesh name
 * @param[in] subMeshNames  Names of the sub meshes in the form "filename:meshname"
 * @param[in] subMeshTypes  Type of each submesh
*/
void writeMultiMesh( DBfile* fid, const std::string& meshname,
    const std::vector<std::string>& subMeshNames,
    const std::vector<int>& subMeshTypes );


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
    const std::vector<int>& subVarTypes );


}; // silo namespace
#endif

#include "IO/silo.hpp"

