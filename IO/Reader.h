#ifndef READER_INC
#define READER_INC

#include <iostream>
#include <memory>
#include <string.h>
#include <vector>

#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"


namespace IO {


//! Get the path to a file
std::string getPath( const std::string &filename );


/*!
 * @brief  Get the maximum number of domains written
 * @details  This function reads the summary files to determine the maximum
 *    number of domains in the output.
 * @param[in] path          The path to use for reading
 * @param[in] format        The data format to use:
 *                              old - Old mesh format (provided for backward compatibility)
 *                              new - New format, 1 file/process
 *                              silo - Silo
 *                              auto - Auto-determin the format
 * @param[in] comm          Optional comm to use to reduce IO load by
 *                          reading on rank 0 and then communicating the result
 */
int maxDomains( const std::string &path, const std::string &format = "auto",
    const Utilities::MPI &comm = MPI_COMM_SELF );


/*!
 * @brief  Read the timestep list
 * @details  This function reads the timestep list from the summary file.
 * @param[in] path          The path to use for reading
 * @param[in] format        The data format to use:
 *                              old - Old mesh format (provided for backward compatibility)
 *                              new - New format, 1 file/process
 *                              silo - Silo
 *                              auto - Auto-determin the format
 * @return append           Append any existing data (default is false)
 */
std::vector<std::string> readTimesteps(
    const std::string &path, const std::string &format = "auto" );


/*!
 * @brief  Read the data for the timestep
 * @details  This function reads the mesh and variable data provided for the given timestep.
 * @param[in] path          The path to use for reading
 * @param[in] timestep      The timestep iteration
 * @param[in] domain        The desired domain to read
 */
std::vector<IO::MeshDataStruct> readData(
    const std::string &path, const std::string &timestep, int domain );


//! Read the list of mesh databases for the given timestep
std::vector<IO::MeshDatabase> getMeshList( const std::string &path, const std::string &timestep );


//! Read the given mesh domain
std::shared_ptr<IO::Mesh> getMesh( const std::string &path, const std::string &timestep,
    const MeshDatabase &meshDatabase, int domain );


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
std::shared_ptr<IO::Variable> getVariable( const std::string &path, const std::string &timestep,
    const MeshDatabase &meshDatabase, int domain, const std::string &variable );


/*!
 * @brief Reformat the variable to match the mesh
 * @details  This function modifies the dimensions of the array to match the mesh
 * @param[in] mesh          The underlying mesh
 * @param[in,out] var       The variable name to read
 */
void reformatVariable( const IO::Mesh &mesh, IO::Variable &var );


} // namespace IO

#endif
