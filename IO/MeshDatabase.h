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
#ifndef MeshDatabase_INC
#define MeshDatabase_INC

#include "IO/Mesh.h" 
#include "common/MPI_Helpers.h"
#include "shared_ptr.h"

#include <iostream>
#include <string.h>
#include <vector>
#include <map>


namespace IO {

class Mesh;


//! Enum to identify mesh type
//enum class MeshType : char { PointMesh=1, SurfaceMesh=2, VolumeMesh=3, Unknown=-1 };
enum MeshType { PointMesh=1, SurfaceMesh=2, VolumeMesh=3, Unknown=-1 };


//! Helper struct for containing offsets for the mesh info
struct DatabaseEntry {
    std::string name;                   //!< Name of the entry
    std::string file;                   //!< Name of the file containing the entry
    size_t offset;                      //!< Offset in the file to start reading
    std::string write( ) const;         //!< Convert the data to a string
    void read( const char* line );      //!< Convert the string to data
    void read( const std::string& line ); //!< Convert the string to data
    DatabaseEntry( ) {}                 //!< Empty constructor
    DatabaseEntry( const char* line );  //!< Convert the string to data
    ~DatabaseEntry() {}                 //!< Destructor
};


//! Structure to hold the info about the variables
struct VariableDatabase {
    std::string name;                   //!< Name of the variable
    IO::VariableType type;              //!< Variable
    unsigned int dim;                   //!< Number of points per grid point (1: scalar, 3: vector, ...)
    // Overload key operators
    bool operator==(const VariableDatabase& rhs ) const;
    bool operator!=(const VariableDatabase& rhs ) const;
    bool operator>=(const VariableDatabase& rhs ) const;
    bool operator<=(const VariableDatabase& rhs ) const;
    bool operator> (const VariableDatabase& rhs ) const;
    bool operator< (const VariableDatabase& rhs ) const;
};


//! Structure to hold the info about the meshes
struct MeshDatabase {
    typedef  std::pair<std::string,std::string>  variable_id;
    std::string name;                   //!< Name of the mesh
    MeshType type;                      //!< Mesh type
    std::string meshClass;              //!< Mesh class
    unsigned char format;               //!< Data format (1: old, 2: new, 3: new (single), 4: silo)
    std::vector<DatabaseEntry> domains; //!< List of the domains
    std::vector<VariableDatabase> variables; //!< List of the variables
    std::map<variable_id,DatabaseEntry> variable_data; //!< Data for the variables
    VariableDatabase getVariableDatabase( const std::string& varname ) const;
public:
    MeshDatabase();
    ~MeshDatabase();
    MeshDatabase(const MeshDatabase&);
    MeshDatabase& operator=(const MeshDatabase&);
};


//! Gather the mesh databases from all processors
std::vector<MeshDatabase> gatherAll( const std::vector<MeshDatabase>& meshes, MPI_Comm comm );


//! Write the mesh databases to a file
void write( const std::vector<MeshDatabase>& meshes, const std::string& filename );


//! Read the mesh databases from a file
std::vector<MeshDatabase> read( const std::string& filename );


//! Return the mesh type
IO::MeshType meshType( const IO::Mesh& mesh );


} // IO namespace

#endif
