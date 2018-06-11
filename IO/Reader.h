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
#ifndef READER_INC
#define READER_INC

#include <iostream>
#include <string.h>
#include <vector>

#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"
#include "shared_ptr.h"


namespace IO {


//! Get the path to a file
std::string getPath( const std::string& filename );


//! List the timesteps in the given directors (dumps.LBPM)
std::vector<std::string> readTimesteps( const std::string& filename );


//! Read the list of mesh databases for the given timestep
std::vector<IO::MeshDatabase> getMeshList( const std::string& path, const std::string& timestep );


//! Read the given mesh domain
std::shared_ptr<IO::Mesh> getMesh( const std::string& path, const std::string& timestep, 
    const MeshDatabase& meshDatabase, int domain );


/*!
 * @brief Read the given variable
 * @details  This function reads the variable data provided for the current timestep
 * @param[in] path          The path to the file
 * @param[in] timestep      The timestep iteration
 * @param[in] meshDatabase  The mesh database (see getMeshList)
 * @param[in] domain        The index of the domain we want to read
 * @param[in] variable      The variable name to read
 * @return                  Returns the variable data as a linear array
 */
std::shared_ptr<IO::Variable> getVariable( const std::string& path, const std::string& timestep, 
    const MeshDatabase& meshDatabase, int domain, const std::string& variable );


/*!
 * @brief Reformat the variable to match the mesh
 * @details  This function modifies the dimensions of the array to match the mesh
 * @param[in] mesh          The underlying mesh
 * @param[in/out] variable  The variable name to read
 */
void reformatVariable( const IO::Mesh& mesh, IO::Variable& var );


} // IO namespace

#endif
