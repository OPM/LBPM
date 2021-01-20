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
#ifndef WRITER_INC
#define WRITER_INC

#include <iostream>
#include <string.h>
#include <vector>

#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"


namespace IO {


/*!
 * @brief  Initialize the writer
 * @details  This function initializes the writer to the given path.  All subsequent
 *    writes will occur in this directory.  If this is not called, then it will default
 *    to the current path.
 * @param[in] path          The path to use for writes
 * @param[in] format        The data format to use:
 *                          old - Old mesh format (provided for backward compatibility, cannot write variables)
 *                          new - New format, 1 file/process
 *                          silo - Silo
 * @param[in] append        Append any existing data (default is false)
 */
void initialize( const std::string& path="", const std::string& format="silo", bool append=false );


/*!
 * @brief  Write the data for the timestep
 * @details  This function writes the mesh and variable data provided for the current timestep
 * @param[in] subdir        The subdirectory to use for the timestep
 * @param[in] meshData      The data to write
 * @param[in] comm          The comm to use for writing (usually MPI_COMM_WORLD or a dup thereof)
 */
void writeData( const std::string& subdir, const std::vector<IO::MeshDataStruct>& meshData, MPI_Comm comm );


/*!
 * @brief  Write the data for the timestep
 * @details  This function writes the mesh and variable data provided for the current timestep
 * @param[in] timestep      The timestep iteration
 * @param[in] meshData      The data to write
 * @param[in] comm          The comm to use for writing (usually MPI_COMM_WORLD or a dup thereof)
 */
inline void writeData( int timestep, const std::vector<IO::MeshDataStruct>& meshData, MPI_Comm comm )
{
    char subdir[100];
    sprintf(subdir,"vis%03i",timestep);
    writeData( subdir, meshData, comm );
}


} // IO namespace

#endif
