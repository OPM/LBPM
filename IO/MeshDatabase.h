#ifndef MeshDatabase_INC
#define MeshDatabase_INC

#include <iostream>
#include <string.h>
#include <memory>
#include <vector>

#ifdef USE_MPI
    #include "mpi.h"
#else
    typedef int MPI_Comm;
#endif

namespace IO {


//! Enum to identify mesh type
enum class MeshType : char { PointMesh=1, SurfaceMesh=2, VolumeMesh=3 };


//! Structure to hold the info about the meshes
struct MeshDatabase {
    std::string name;                   //!< Name of the mesh
    MeshType type;                      //!< Mesh type
    unsigned char format;               //!< Data format
    std::vector<std::string> domains;   //!< List of the domains
    std::vector<std::string> file;      //!< File containing the mesh data for the given domain
    std::vector<std::string> variables; //!< List of the variables
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


} // IO namespace

#endif
