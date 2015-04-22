#ifndef READER_INC
#define READER_INC

#include <iostream>
#include <string.h>
#include <vector>

#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"
#include "shared_ptr.h"


namespace IO {


//! List the timesteps in the given directors (dumps.LBPM)
std::vector<std::string> readTimesteps( const std::string& filename );


//! Read the list of variables for the given timestep
std::vector<IO::MeshDatabase> getMeshList( const std::string& path, const std::string& timestep );


//! Read the given mesh domain
std::shared_ptr<IO::Mesh> getMesh( const std::string& path, const std::string& timestep, 
    const MeshDatabase& meshDatabase, int domain );


//! Read the given mesh domain
std::shared_ptr<IO::Variable> getVariable( const std::string& path, const std::string& timestep, 
    const MeshDatabase& meshDatabase, int domain, const std::string& variable );


} // IO namespace

#endif
