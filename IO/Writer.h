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
 * @details  This function initializes the writer to the given path.
 *    All subsequent writes will occur in this directory.
 *    If this is not called, then it will default to the current path.
 * @param[in] path          The path to use for writes
 * @param[in] format        The data format to use:
 *                              old - Old mesh format
 *                                    (provided for backward compatibility, cannot write variables)
 *                              new - New format, 1 file/process
 *                              silo - Silo
 *                              hdf5 - HDF5 + XMDF
 * @param[in] append        Append any existing data (default is false)
 */
void initialize(
    const std::string &path = "", const std::string &format = "hdf5", bool append = false );


/*!
 * @brief  Write the data for the timestep
 * @details  This function writes the mesh and variable data provided for the current timestep
 * @param[in] subdir        The subdirectory to use for the timestep
 * @param[in] meshData      The data to write
 * @param[in] comm          The comm to use for writing (usually MPI_COMM_WORLD or a dup thereof)
 */
void writeData( const std::string &subdir, const std::vector<IO::MeshDataStruct> &meshData,
    const Utilities::MPI &comm );


/*!
 * @brief  Write the data for the timestep
 * @details  This function writes the mesh and variable data provided for the current timestep
 * @param[in] timestep      The timestep iteration
 * @param[in] meshData      The data to write
 * @param[in] comm          The comm to use for writing (usually MPI_COMM_WORLD or a dup thereof)
 */
inline void writeData(
    int timestep, const std::vector<IO::MeshDataStruct> &meshData, const Utilities::MPI &comm )
{
    char subdir[100];
    sprintf( subdir, "vis%03i", timestep );
    writeData( subdir, meshData, comm );
}


// Create the database entry for the mesh data
IO::MeshDatabase getDatabase(
    const std::string &filename, const IO::MeshDataStruct &mesh, IO::FileFormat format, int rank );


} // namespace IO

#endif
