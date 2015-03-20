#ifndef WRITER_INC
#define WRITER_INC

#include <iostream>
#include <string.h>
#include <vector>

#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"


namespace IO {


/*!
 * @brief  Write the data for the timestep
 * @details  This function writes the mesh and variable data provided for the current timestep
 * @param[in] timestep      The timestep iteration
 * @param[in] meshData      The data to write
 * @param[in] format        The data format to use:
 *                          1 - Old mesh format (provided for backward compatibility, cannot write variables)
 *                          2 - New format, 1 file/process, double precision
 *                          3 - New format, 1 file/process, single precision (not finished)
 */
void writeData( int timestep, const std::vector<IO::MeshDataStruct>& meshData, int format );


} // IO namespace

#endif
