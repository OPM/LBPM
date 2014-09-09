#ifndef WRITER_INC
#define WRITER_INC

#include <iostream>
#include <string.h>
#include <memory>
#include <vector>

#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"


namespace IO {


void writeData( int timestep, const std::vector<IO::MeshDataStruct>& meshData );


} // IO namespace

#endif
