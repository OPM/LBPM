#ifndef MeshDatabase_INC
#define MeshDatabase_INC

#include <iostream>
#include <string.h>
#include <memory>
#include <vector>
#include <map>

#ifdef USE_MPI
    #include "mpi.h"
#else
    typedef int MPI_Comm;
#endif


namespace IO {

class Mesh;


//! Enum to identify mesh type
enum class MeshType : char { PointMesh=1, SurfaceMesh=2, VolumeMesh=3, Unknown=-1 };


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


//! Structure to hold the info about the meshes
struct MeshDatabase {
    std::string name;                   //!< Name of the mesh
    MeshType type;                      //!< Mesh type
    std::string meshClass;              //!< Mesh class
    unsigned char format;               //!< Data format
    std::vector<DatabaseEntry> domains; //!< List of the domains
    std::vector<std::string> variables; //!< List of the variables
    std::map<std::pair<std::string,std::string>,DatabaseEntry> variable_data; //!< Data for the variables
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
IO::MeshType meshType( std::shared_ptr<IO::Mesh> mesh );


} // IO namespace

#endif
