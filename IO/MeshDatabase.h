#ifndef MeshDatabase_INC
#define MeshDatabase_INC

#include "IO/Mesh.h"
#include "common/MPI.h"

#include <iostream>
#include <map>
#include <memory>
#include <string.h>
#include <vector>


namespace IO {


//! Helper struct for containing offsets for the mesh info
struct DatabaseEntry {
    std::string name;                     //!< Name of the entry
    std::string file;                     //!< Name of the file containing the entry
    size_t offset;                        //!< Offset in the file to start reading
    std::string write() const;            //!< Convert the data to a string
    void read( const char *line );        //!< Convert the string to data
    void read( const std::string &line ); //!< Convert the string to data
    DatabaseEntry() {}                    //!< Empty constructor
    DatabaseEntry( const char *line );    //!< Convert the string to data
    ~DatabaseEntry() {}                   //!< Destructor
};


//! Structure to hold the info about the variables
struct VariableDatabase {
    std::string name;      //!< Name of the variable
    IO::VariableType type; //!< Variable
    unsigned int dim;      //!< Number of points per grid point (1: scalar, 3: vector, ...)
    // Overload key operators
    bool operator==( const VariableDatabase &rhs ) const;
    bool operator!=( const VariableDatabase &rhs ) const;
    bool operator>=( const VariableDatabase &rhs ) const;
    bool operator<=( const VariableDatabase &rhs ) const;
    bool operator>( const VariableDatabase &rhs ) const;
    bool operator<( const VariableDatabase &rhs ) const;
};


//! Structure to hold the info about the meshes
struct MeshDatabase {
    typedef std::pair<std::string, std::string> variable_id;
    std::string name;                   //!< Name of the mesh
    MeshType type;                      //!< Mesh type
    std::string meshClass;              //!< Mesh class
    FileFormat format;                  //!< Data format (1: old, 2: new, 3: new (single), 4: silo)
    std::vector<DatabaseEntry> domains; //!< List of the domains
    std::vector<VariableDatabase> variables;            //!< List of the variables
    std::map<variable_id, DatabaseEntry> variable_data; //!< Data for the variables
    VariableDatabase getVariableDatabase( const std::string &varname ) const;

public:
    MeshDatabase();
    ~MeshDatabase();
    MeshDatabase( const MeshDatabase & );
    MeshDatabase &operator=( const MeshDatabase & );
};


//! Gather the mesh databases from all processors
std::vector<MeshDatabase> gatherAll(
    const std::vector<MeshDatabase> &meshes, const Utilities::MPI &comm );


//! Write the mesh databases to a file
void write( const std::vector<MeshDatabase> &meshes, const std::string &filename );


//! Read the mesh databases from a file
std::vector<MeshDatabase> read( const std::string &filename );


//! Return the mesh type
IO::MeshType meshType( const IO::Mesh &mesh );


} // namespace IO

#endif
