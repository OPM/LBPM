#include "IO/MeshDatabase.h"
#include "IO/IOHelpers.h"
#include "IO/Mesh.h"
#include "IO/PackData.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include <cstdio>
#include <map>
#include <set>
#include <vector>

#include <ProfilerApp.h>


// Default pack/unpack
// clang-format off
#define INSTANTIATE_PACK( TYPE )                            \
    template<>                                              \
    size_t packsize<TYPE>( const TYPE &rhs )                \
    {                                                       \
        return sizeof( TYPE );                              \
    }                                                       \
    template<>                                              \
    void pack<TYPE>( const TYPE &rhs, char *buffer )        \
    {                                                       \
        memcpy( buffer, &rhs, sizeof( IO::MeshType ) );     \
    }                                                       \
    template<>                                              \
    void unpack<TYPE>( TYPE &data, const char *buffer )     \
    {                                                       \
        memcpy( &data, buffer, sizeof( IO::MeshType ) );    \
    }
INSTANTIATE_PACK( IO::VariableType )
INSTANTIATE_PACK( IO::DataType )
INSTANTIATE_PACK( IO::MeshType )
INSTANTIATE_PACK( IO::FileFormat )
// clang-format on


// DatabaseEntry
template<>
size_t packsize<IO::DatabaseEntry>( const IO::DatabaseEntry &rhs )
{
    return packsize( rhs.name ) + packsize( rhs.file ) + packsize( rhs.offset );
}
template<>
void pack<IO::DatabaseEntry>( const IO::DatabaseEntry &rhs, char *buffer )
{
    size_t i = 0;
    pack( rhs.name, &buffer[i] );
    i += packsize( rhs.name );
    pack( rhs.file, &buffer[i] );
    i += packsize( rhs.file );
    pack( rhs.offset, &buffer[i] );
    i += packsize( rhs.offset );
}
template<>
void unpack<IO::DatabaseEntry>( IO::DatabaseEntry &data, const char *buffer )
{
    size_t i = 0;
    unpack( data.name, &buffer[i] );
    i += packsize( data.name );
    unpack( data.file, &buffer[i] );
    i += packsize( data.file );
    unpack( data.offset, &buffer[i] );
    i += packsize( data.offset );
}
// VariableDatabase
template<>
size_t packsize<IO::VariableDatabase>( const IO::VariableDatabase &rhs )
{
    return packsize( rhs.name ) + packsize( rhs.type ) + packsize( rhs.dim );
}
template<>
void pack<IO::VariableDatabase>( const IO::VariableDatabase &rhs, char *buffer )
{
    size_t i = 0;
    pack( rhs.name, &buffer[i] );
    i += packsize( rhs.name );
    pack( rhs.type, &buffer[i] );
    i += packsize( rhs.type );
    pack( rhs.dim, &buffer[i] );
    i += packsize( rhs.dim );
}
template<>
void unpack<IO::VariableDatabase>( IO::VariableDatabase &data, const char *buffer )
{
    size_t i = 0;
    unpack( data.name, &buffer[i] );
    i += packsize( data.name );
    unpack( data.type, &buffer[i] );
    i += packsize( data.type );
    unpack( data.dim, &buffer[i] );
    i += packsize( data.dim );
}
// MeshDatabase
template<>
size_t packsize<IO::MeshDatabase>( const IO::MeshDatabase &data )
{
    return packsize( data.name ) + packsize( data.type ) + packsize( data.meshClass ) +
           packsize( data.format ) + packsize( data.domains ) + packsize( data.variables ) +
           packsize( data.variable_data );
}
template<>
void pack<IO::MeshDatabase>( const IO::MeshDatabase &rhs, char *buffer )
{
    size_t i = 0;
    pack( rhs.name, &buffer[i] );
    i += packsize( rhs.name );
    pack( rhs.type, &buffer[i] );
    i += packsize( rhs.type );
    pack( rhs.meshClass, &buffer[i] );
    i += packsize( rhs.meshClass );
    pack( rhs.format, &buffer[i] );
    i += packsize( rhs.format );
    pack( rhs.domains, &buffer[i] );
    i += packsize( rhs.domains );
    pack( rhs.variables, &buffer[i] );
    i += packsize( rhs.variables );
    pack( rhs.variable_data, &buffer[i] );
    i += packsize( rhs.variable_data );
}
template<>
void unpack<IO::MeshDatabase>( IO::MeshDatabase &data, const char *buffer )
{
    size_t i = 0;
    unpack( data.name, &buffer[i] );
    i += packsize( data.name );
    unpack( data.type, &buffer[i] );
    i += packsize( data.type );
    unpack( data.meshClass, &buffer[i] );
    i += packsize( data.meshClass );
    unpack( data.format, &buffer[i] );
    i += packsize( data.format );
    unpack( data.domains, &buffer[i] );
    i += packsize( data.domains );
    unpack( data.variables, &buffer[i] );
    i += packsize( data.variables );
    unpack( data.variable_data, &buffer[i] );
    i += packsize( data.variable_data );
}


namespace IO {


/****************************************************
 * VariableDatabase                                  *
 ****************************************************/
bool VariableDatabase::operator==( const VariableDatabase &rhs ) const
{
    return type == rhs.type && dim == rhs.dim && name == rhs.name;
}
bool VariableDatabase::operator!=( const VariableDatabase &rhs ) const
{
    return type != rhs.type || dim != rhs.dim || name != rhs.name;
}
bool VariableDatabase::operator>=( const VariableDatabase &rhs ) const
{
    return operator>( rhs ) || operator==( rhs );
}
bool VariableDatabase::operator<=( const VariableDatabase &rhs ) const { return !operator>( rhs ); }
bool VariableDatabase::operator>( const VariableDatabase &rhs ) const
{
    if ( name > rhs.name )
        return true;
    else if ( name < rhs.name )
        return false;
    if ( type > rhs.type )
        return true;
    else if ( type < rhs.type )
        return false;
    if ( dim > rhs.dim )
        return true;
    else if ( dim < rhs.dim )
        return false;
    return false;
}
bool VariableDatabase::operator<( const VariableDatabase &rhs ) const
{
    return !operator>( rhs ) && operator!=( rhs );
}


/****************************************************
 * MeshDatabase                                      *
 ****************************************************/
MeshDatabase::MeshDatabase() {}
MeshDatabase::~MeshDatabase() {}
MeshDatabase::MeshDatabase( const MeshDatabase &rhs )
{
    name          = rhs.name;
    type          = rhs.type;
    meshClass     = rhs.meshClass;
    format        = rhs.format;
    domains       = rhs.domains;
    variables     = rhs.variables;
    variable_data = rhs.variable_data;
}
MeshDatabase &MeshDatabase::operator=( const MeshDatabase &rhs )
{
    this->name          = rhs.name;
    this->type          = rhs.type;
    this->meshClass     = rhs.meshClass;
    this->format        = rhs.format;
    this->domains       = rhs.domains;
    this->variables     = rhs.variables;
    this->variable_data = rhs.variable_data;
    return *this;
}
VariableDatabase MeshDatabase::getVariableDatabase( const std::string &varname ) const
{
    for ( size_t i = 0; i < variables.size(); i++ ) {
        if ( variables[i].name == varname )
            return variables[i];
    }
    return VariableDatabase();
}


/****************************************************
 * DatabaseEntry                                     *
 ****************************************************/
std::string DatabaseEntry::write() const
{
    char tmp[1000];
    sprintf( tmp, "%s; %s; %lu", name.c_str(), file.c_str(), offset );
    return std::string( tmp );
}
DatabaseEntry::DatabaseEntry( const char *line )
{
    auto list = splitList( line, ';' );
    name      = list[0];
    file      = list[1];
    offset    = atol( list[2].c_str() );
}
void DatabaseEntry::read( const char *line )
{
    auto list = splitList( line, ';' );
    name      = list[0];
    file      = list[1];
    offset    = atol( list[2].c_str() );
}
void DatabaseEntry::read( const std::string &line )
{
    auto list = splitList( line.c_str(), ';' );
    name      = list[0];
    file      = list[1];
    offset    = atol( list[2].c_str() );
}


// Gather the mesh databases from all processors
inline int tod( int N ) { return ( N + 7 ) / sizeof( double ); }
std::vector<MeshDatabase> gatherAll(
    const std::vector<MeshDatabase> &meshes, const Utilities::MPI &comm )
{
    if ( comm.getSize() == 1 )
        return meshes;
    PROFILE_START( "gatherAll" );
    PROFILE_START( "gatherAll-pack", 2 );
    int size = comm.getSize();
    // First pack the mesh data to local buffers
    int localsize = 0;
    for ( size_t i = 0; i < meshes.size(); i++ )
        localsize += tod( packsize( meshes[i] ) );
    auto localbuf = new double[localsize];
    int pos       = 0;
    for ( size_t i = 0; i < meshes.size(); i++ ) {
        pack( meshes[i], (char *) &localbuf[pos] );
        pos += tod( packsize( meshes[i] ) );
    }
    PROFILE_STOP( "gatherAll-pack", 2 );
    // Get the number of bytes each processor will be sending/recieving
    PROFILE_START( "gatherAll-send1", 2 );
    auto recvsize  = comm.allGather( localsize );
    int globalsize = recvsize[0];
    auto disp      = new int[size];
    disp[0]        = 0;
    for ( int i = 1; i < size; i++ ) {
        disp[i] = disp[i - 1] + recvsize[i];
        globalsize += recvsize[i];
    }
    PROFILE_STOP( "gatherAll-send1", 2 );
    // Send/recv the global data
    PROFILE_START( "gatherAll-send2", 2 );
    auto globalbuf = new double[globalsize];
    comm.allGather( localbuf, localsize, globalbuf, recvsize.data(), disp, true );
    PROFILE_STOP( "gatherAll-send2", 2 );
    // Unpack the data
    PROFILE_START( "gatherAll-unpack", 2 );
    std::map<std::string, MeshDatabase> data;
    pos = 0;
    while ( pos < globalsize ) {
        MeshDatabase tmp;
        unpack( tmp, (char *) &globalbuf[pos] );
        pos += tod( packsize( tmp ) );
        std::map<std::string, MeshDatabase>::iterator it = data.find( tmp.name );
        if ( it == data.end() ) {
            data[tmp.name] = tmp;
        } else {
            for ( size_t i = 0; i < tmp.domains.size(); i++ )
                it->second.domains.push_back( tmp.domains[i] );
            for ( size_t i = 0; i < tmp.variables.size(); i++ )
                it->second.variables.push_back( tmp.variables[i] );
            it->second.variable_data.insert( tmp.variable_data.begin(), tmp.variable_data.end() );
        }
    }
    for ( auto it = data.begin(); it != data.end(); ++it ) {
        // Get the unique variables
        std::set<VariableDatabase> data2(
            it->second.variables.begin(), it->second.variables.end() );
        it->second.variables = std::vector<VariableDatabase>( data2.begin(), data2.end() );
    }
    // Free temporary memory
    delete[] localbuf;
    delete[] disp;
    delete[] globalbuf;
    // Return the results
    std::vector<MeshDatabase> data2( data.size() );
    size_t i = 0;
    for ( auto it = data.begin(); it != data.end(); ++it, ++i )
        data2[i] = it->second;
    PROFILE_STOP( "gatherAll-unpack", 2 );
    PROFILE_STOP( "gatherAll" );
    return data2;
}


//! Write the mesh databases to a file
void write( const std::vector<MeshDatabase> &meshes, const std::string &filename )
{
    PROFILE_START( "write" );
    FILE *fid = fopen( filename.c_str(), "wb" );
    for ( size_t i = 0; i < meshes.size(); i++ ) {
        fprintf( fid, "%s\n", meshes[i].name.c_str() );
        fprintf( fid, "   type: %s\n", getString( meshes[i].type ).data() );
        fprintf( fid, "   meshClass: %s\n", meshes[i].meshClass.c_str() );
        fprintf( fid, "   format: %s\n", getString( meshes[i].format ).data() );
        for ( size_t j = 0; j < meshes[i].domains.size(); j++ )
            fprintf( fid, "   domain: %s\n", meshes[i].domains[j].write().c_str() );
        fprintf( fid, "   variables: " );
        for ( size_t j = 0; j < meshes[i].variables.size(); j++ ) {
            const VariableDatabase &var = meshes[i].variables[j];
            fprintf( fid, "%s|%s|%i; ", var.name.data(), getString( var.type ).data(), var.dim );
        }
        fprintf( fid, "\n" );
        for ( auto it = meshes[i].variable_data.begin(); it != meshes[i].variable_data.end();
              ++it ) {
            const char *domain   = it->first.first.c_str();
            const char *variable = it->first.second.c_str();
            fprintf(
                fid, "   variable(%s,%s): %s\n", domain, variable, it->second.write().c_str() );
        }
    }
    fclose( fid );
    PROFILE_STOP( "write" );
}


//! Read the mesh databases from a file
std::vector<MeshDatabase> read( const std::string &filename )
{
    std::vector<MeshDatabase> meshes;
    PROFILE_START( "read" );
    FILE *fid = fopen( filename.c_str(), "rb" );
    if ( fid == NULL )
        ERROR( "Error opening file: " + filename );
    char *line = new char[10000];
    while ( std::fgets( line, 1000, fid ) != NULL ) {
        if ( line[0] < 32 ) {
            // Empty line
            continue;
        } else if ( line[0] != ' ' ) {
            meshes.resize( meshes.size() + 1 );
            std::string name( line );
            name.resize( name.size() - 1 );
            meshes.back().name = name;
        } else if ( strncmp( line, "   format:", 10 ) == 0 ) {
            meshes.back().format = getFileFormat( &line[10] );
        } else if ( strncmp( line, "   type:", 8 ) == 0 ) {
            meshes.back().type = getMeshType( &line[8] );
        } else if ( strncmp( line, "   meshClass:", 13 ) == 0 ) {
            meshes.back().meshClass = deblank( std::string( &line[13] ) );
        } else if ( strncmp( line, "   domain:", 10 ) == 0 ) {
            DatabaseEntry data( &line[10] );
            meshes.back().domains.push_back( data );
        } else if ( strncmp( line, "   variables:", 13 ) == 0 ) {
            MeshDatabase &mesh                 = meshes.back();
            std::vector<std::string> variables = splitList( &line[13], ';' );
            mesh.variables.resize( variables.size() );
            for ( size_t i = 0; i < variables.size(); i++ ) {
                std::vector<std::string> tmp = splitList( variables[i].c_str(), '|' );
                ASSERT( tmp.size() == 3 );
                mesh.variables[i].name = tmp[0];
                mesh.variables[i].type = getVariableType( tmp[1] );
                mesh.variables[i].dim  = atoi( tmp[2].c_str() );
            }
        } else if ( strncmp( line, "   variable(", 12 ) == 0 ) {
            size_t i1            = find( line, ',' );
            size_t i2            = find( line, ':' );
            std::string domain   = deblank( std::string( line, 12, i1 - 12 ) );
            std::string variable = deblank( std::string( line, i1 + 1, i2 - i1 - 2 ) );
            std::pair<std::string, std::string> key( domain, variable );
            DatabaseEntry data( &line[i2 + 1] );
            meshes.back().variable_data.insert(
                std::pair<std::pair<std::string, std::string>, DatabaseEntry>( key, data ) );
        } else {
            ERROR( "Error reading line" );
        }
    }
    fclose( fid );
    delete[] line;
    PROFILE_STOP( "read" );
    return meshes;
}


// Return the mesh type
IO::MeshType meshType( const IO::Mesh &mesh )
{
    IO::MeshType type           = IO::MeshType::Unknown;
    const std::string meshClass = mesh.className();
    if ( meshClass == "PointList" ) {
        type = IO::MeshType::PointMesh;
    } else if ( meshClass == "TriList" || meshClass == "TriMesh" ) {
        type = IO::MeshType::SurfaceMesh;
    } else if ( meshClass == "DomainMesh" ) {
        type = IO::MeshType::VolumeMesh;
    } else {
        ERROR( "Unknown mesh" );
    }
    return type;
}


} // namespace IO
