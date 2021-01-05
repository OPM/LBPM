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
