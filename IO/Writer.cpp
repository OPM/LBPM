#include "IO/Writer.h"
#include "IO/HDF5_IO.h"
#include "IO/IOHelpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Xdmf.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <vector>


enum class Format { OLD, NEW, SILO, HDF5, UNKNOWN };


/****************************************************
 * External declerations                             *
 ****************************************************/
std::vector<IO::MeshDatabase> writeMeshesSilo(
    const std::vector<IO::MeshDataStruct> &, const std::string &, IO::FileFormat, int );
void writeSiloSummary( const std::vector<IO::MeshDatabase> &, const std::string & );
std::vector<IO::MeshDatabase> writeMeshesHDF5(
    const std::vector<IO::MeshDataStruct> &, const std::string &, IO::FileFormat, int, Xdmf & );


/****************************************************
 * Recursively create the subdirectory               *
 ****************************************************/
static void recursiveMkdir( const std::string &path, mode_t mode )
{
    // Iterate through the root directories until we create the desired path
    for ( size_t pos = 0; pos < path.size(); ) {
        // slide backwards in string until next slash found
        pos++;
        for ( ; pos < path.size(); pos++ ) {
            if ( path[pos] == '/' || path[pos] == 92 )
                break;
        }
        // Create the temporary path
        auto path2 = path.substr( 0, pos );
        // Check if the temporary path exists
        struct stat status;
        int result = stat( path2.data(), &status );
        if ( result == 0 ) {
            // if there is a part of the path that already exists make sure it is really a directory
            if ( !S_ISDIR( status.st_mode ) ) {
                ERROR(
                    "Error in recursiveMkdir...\n"
                    "    Cannot create directories in path = " +
                    path +
                    "\n    because some intermediate item in path exists and is NOT a directory" );
            }
            continue;
        }
        // Create the directory and test the result
        result = mkdir( path2.data(), mode );
        if ( result != 0 ) {
            // Maybe another rank created the directory, check
            int result = stat( path2.data(), &status );
            if ( result != 0 && !S_ISDIR( status.st_mode ) )
                ERROR( "Error in Utilities::recursiveMkdir...\n"
                       "    Cannot create directory  = " +
                       path2 );
        }
    }
}


/****************************************************
 * Initialize the writer                             *
 ****************************************************/
static std::string global_IO_path;
static Format global_IO_format = Format::UNKNOWN;
void IO::initialize( const std::string &path, const std::string &format, bool append )
{
    if ( path.empty() )
        global_IO_path = ".";
    else
        global_IO_path = path;
    if ( format == "old" )
        global_IO_format = Format::OLD;
    else if ( format == "new" )
        global_IO_format = Format::NEW;
    else if ( format == "silo" )
        global_IO_format = Format::SILO;
    else if ( format == "hdf5" )
        global_IO_format = Format::HDF5;
    else
        ERROR( "Unknown format" );
    int rank = Utilities::MPI( MPI_COMM_WORLD ).getRank();
    if ( !append && rank == 0 ) {
        recursiveMkdir( path, S_IRWXU | S_IRGRP );
        std::string filename;
        if ( global_IO_format == Format::OLD || global_IO_format == Format::NEW )
            filename = global_IO_path + "/summary.LBM";
        else if ( global_IO_format == Format::SILO || global_IO_format == Format::HDF5 )
            filename = global_IO_path + "/LBM.visit";
        else
            ERROR( "Unknown format" );
        auto fid = fopen( filename.c_str(), "wb" );
        fclose( fid );
    }
}


// Write the mesh data in the original format
static std::vector<IO::MeshDatabase> writeMeshesOrigFormat(
    const std::vector<IO::MeshDataStruct> &meshData, const std::string &path, int rank )
{
    std::vector<IO::MeshDatabase> meshes_written;
    for ( size_t i = 0; i < meshData.size(); i++ ) {
        char domainname[100], filename[100], fullpath[200];
        sprintf( domainname, "%05i", rank );
        sprintf( filename, "%s.%05i", meshData[i].meshName.c_str(), rank );
        sprintf( fullpath, "%s/%s", path.c_str(), filename );
        FILE *fid = fopen( fullpath, "wb" );
        INSIST( fid != NULL, std::string( "Error opening file: " ) + fullpath );
        std::shared_ptr<IO::Mesh> mesh = meshData[i].mesh;
        IO::MeshDatabase mesh_entry;
        mesh_entry.name      = meshData[i].meshName;
        mesh_entry.type      = meshType( *mesh );
        mesh_entry.meshClass = meshData[i].mesh->className();
        mesh_entry.format    = IO::FileFormat::OLD;
        IO::DatabaseEntry domain;
        domain.name   = domainname;
        domain.file   = filename;
        domain.offset = 0;
        mesh_entry.domains.push_back( domain );
        static bool printVariableSupportMsg = true;
        if ( !meshData[i].vars.empty() && printVariableSupportMsg ) {
            printVariableSupportMsg = false;
            printf( "Warning: variables are not supported with this format (original)\n" );
        }
        const std::string meshClass = mesh->className();
        if ( meshClass == "PointList" ) {
            // List of points
            std::shared_ptr<IO::PointList> pointlist =
                std::dynamic_pointer_cast<IO::PointList>( mesh );
            const std::vector<Point> &P = pointlist->points;
            for ( size_t i = 0; i < P.size(); i++ ) {
                double x[3];
                x[0] = P[i].x;
                x[1] = P[i].y;
                x[2] = P[i].z;
                fwrite( x, sizeof( double ), 3, fid );
            }
        } else if ( meshClass == "TriList" || meshClass == "TriMesh" ) {
            // Triangle mesh
            std::shared_ptr<IO::TriList> trilist = IO::getTriList( mesh );
            const std::vector<Point> &A          = trilist->A;
            const std::vector<Point> &B          = trilist->B;
            const std::vector<Point> &C          = trilist->C;
            for ( size_t i = 0; i < A.size(); i++ ) {
                double tri[9];
                tri[0] = A[i].x;
                tri[1] = A[i].y;
                tri[2] = A[i].z;
                tri[3] = B[i].x;
                tri[4] = B[i].y;
                tri[5] = B[i].z;
                tri[6] = C[i].x;
                tri[7] = C[i].y;
                tri[8] = C[i].z;
                fwrite( tri, sizeof( double ), 9, fid );
            }
        } else if ( meshClass == "DomainMesh" ) {
            // This format was never supported with the old format
        } else {
            ERROR( "Unknown mesh" );
        }
        fclose( fid );
        std::sort( mesh_entry.variables.begin(), mesh_entry.variables.end() );
        mesh_entry.variables.erase(
            std::unique( mesh_entry.variables.begin(), mesh_entry.variables.end() ),
            mesh_entry.variables.end() );
        meshes_written.push_back( mesh_entry );
    }
    return meshes_written;
}


// Create the database entry for the mesh data
IO::MeshDatabase IO::getDatabase(
    const std::string &filename, const IO::MeshDataStruct &mesh, IO::FileFormat format, int rank )
{
    char domainname[100];
    sprintf( domainname, "%s_%05i", mesh.meshName.c_str(), rank );
    // Create the MeshDatabase
    IO::MeshDatabase database;
    database.name      = mesh.meshName;
    database.type      = meshType( *( mesh.mesh ) );
    database.meshClass = mesh.mesh->className();
    database.format    = format;
    // Write the mesh
    IO::DatabaseEntry domain;
    domain.name   = domainname;
    domain.file   = filename;
    domain.offset = -1;
    database.domains.push_back( domain );
    // Write the variables
    for ( size_t i = 0; i < mesh.vars.size(); i++ ) {
        // Add basic variable info
        IO::VariableDatabase info;
        info.name = mesh.vars[i]->name;
        info.type = mesh.vars[i]->type;
        info.dim  = mesh.vars[i]->dim;
        database.variables.push_back( info );
        // Add domain variable info
        IO::DatabaseEntry variable;
        variable.name   = mesh.vars[i]->name;
        variable.file   = filename;
        variable.offset = -1;
        std::pair<std::string, std::string> key( domain.name, mesh.vars[i]->name );
        database.variable_data.insert(
            std::pair<std::pair<std::string, std::string>, IO::DatabaseEntry>( key, variable ) );
    }
    return database;
}


// Write a mesh (and variables) to a file
static IO::MeshDatabase write_domain( FILE *fid, const std::string &filename,
    const IO::MeshDataStruct &mesh, IO::FileFormat format, int rank )
{
    ASSERT( !mesh.meshName.empty() );
    const int level = 0;
    // Create the MeshDatabase
    IO::MeshDatabase database = getDatabase( filename, mesh, format, rank );
    // Write the mesh
    IO::DatabaseEntry &domain      = database.domains[0];
    domain.offset                  = ftell( fid );
    std::pair<size_t, void *> data = mesh.mesh->pack( level );
    fprintf( fid, "Mesh: %s-%05i: %lu\n", mesh.meshName.c_str(), rank, data.first );
    fwrite( data.second, 1, data.first, fid );
    fprintf( fid, "\n" );
    delete[]( char * ) data.second;
    // Write the variables
    for ( size_t i = 0; i < mesh.vars.size(); i++ ) {
        ASSERT( mesh.vars[i]->type != IO::VariableType::NullVariable );
        std::pair<std::string, std::string> key( domain.name, mesh.vars[i]->name );
        auto &variable  = database.variable_data[key];
        variable.offset = ftell( fid );
        int dim         = mesh.vars[i]->dim;
        auto type       = getString( mesh.vars[i]->type );
        size_t N        = mesh.vars[i]->data.length();
        size_t N_mesh   = mesh.mesh->numberPointsVar( mesh.vars[i]->type );
        ASSERT( N == dim * N_mesh );
        ASSERT( !type.empty() );
        ASSERT( !variable.name.empty() );
        fprintf( fid, "Var: %s-%05i-%s: %i, %s, %lu, %lu, double\n", database.name.c_str(), rank,
            variable.name.c_str(), dim, type.data(), N_mesh, N * sizeof( double ) );
        fwrite( mesh.vars[i]->data.data(), sizeof( double ), N, fid );
        fprintf( fid, "\n" );
    }
    return database;
}


// Write the mesh data in the new format
static std::vector<IO::MeshDatabase> writeMeshesNewFormat(
    const std::vector<IO::MeshDataStruct> &meshData, const std::string &path, IO::FileFormat format,
    int rank )
{
    std::vector<IO::MeshDatabase> meshes_written;
    char filename[100], fullpath[200];
    sprintf( filename, "%05i", rank );
    sprintf( fullpath, "%s/%s", path.c_str(), filename );
    FILE *fid = fopen( fullpath, "wb" );
    ASSERT( fid != nullptr );
    for ( size_t i = 0; i < meshData.size(); i++ )
        meshes_written.push_back( write_domain( fid, filename, meshData[i], format, rank ) );
    fclose( fid );
    return meshes_written;
}


/****************************************************
 * Write the mesh data                               *
 ****************************************************/
void IO::writeData( const std::string &subdir, const std::vector<IO::MeshDataStruct> &meshData,
    const Utilities::MPI &comm )
{
    if ( global_IO_path.empty() )
        IO::initialize();
    PROFILE_START( "writeData" );
    int rank = Utilities::MPI( MPI_COMM_WORLD ).getRank();
    // Check the meshData before writing
    for ( const auto &data : meshData )
        ASSERT( data.check() );
    // Create the output directory
    std::string path = global_IO_path + "/" + subdir;
    recursiveMkdir( path, S_IRWXU | S_IRGRP );
    // Write the mesh files
    Xdmf xmf;
    std::vector<IO::MeshDatabase> meshes_written;
    if ( global_IO_format == Format::OLD ) {
        // Write the original triangle format
        meshes_written = writeMeshesOrigFormat( meshData, path, rank );
    } else if ( global_IO_format == Format::NEW ) {
        // Write the new format (double precision)
        meshes_written = writeMeshesNewFormat( meshData, path, IO::FileFormat::NEW, rank );
    } else if ( global_IO_format == Format::SILO ) {
        // Write silo
        meshes_written = writeMeshesSilo( meshData, path, IO::FileFormat::SILO, rank );
    } else if ( global_IO_format == Format::HDF5 ) {
        // Write hdf5
        meshes_written = writeMeshesHDF5( meshData, path, IO::FileFormat::HDF5, rank, xmf );
    } else {
        ERROR( "Unknown format" );
    }
    // Gather a complete list of files on rank 0
    meshes_written = gatherAll( meshes_written, comm );
    // Gather xmf file (if applicable)
    if ( global_IO_format == Format::HDF5 ) {
        xmf.gather( comm );
    }
    // Write the summary files
    if ( rank == 0 ) {
        // Write the summary file for the current timestep
        write( meshes_written, path + "/LBM.summary" );
        // Write summary file if needed
        if ( global_IO_format == Format::SILO ) {
            writeSiloSummary( meshes_written, path + "/summary.silo" );
        } else if ( global_IO_format == Format::HDF5 ) {
            xmf.write( path + "/summary.xmf" );
        }
        // Add the timestep to the global summary file
        if ( global_IO_format == Format::OLD || global_IO_format == Format::NEW ) {
            auto filename = global_IO_path + "/summary.LBM";
            FILE *fid     = fopen( filename.c_str(), "ab" );
            fprintf( fid, "%s/\n", subdir.c_str() );
            fclose( fid );
        } else if ( global_IO_format == Format::SILO ) {
            auto filename = global_IO_path + "/LBM.visit";
            FILE *fid     = fopen( filename.c_str(), "ab" );
            fprintf( fid, "%s/summary.silo\n", subdir.c_str() );
            fclose( fid );
        } else if ( global_IO_format == Format::HDF5 ) {
            auto filename = global_IO_path + "/LBM.visit";
            FILE *fid     = fopen( filename.c_str(), "ab" );
            fprintf( fid, "%s/summary.xmf\n", subdir.c_str() );
            fclose( fid );
        } else {
            ERROR( "Unknown format" );
        }
    }
    PROFILE_STOP( "writeData" );
}
