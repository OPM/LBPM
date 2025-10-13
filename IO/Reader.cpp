#include "IO/Reader.h"
#include "IO/HDF5_IO.h"
#include "IO/IOHelpers.h"
#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"
#include "IO/silo.h"
#include "common/Utilities.h"

#include <ProfilerApp.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string.h>
#include <vector>


// Inline function to read line without a return argument
static inline void fgetl( char *str, int num, FILE *stream )
{
    char *ptr = fgets( str, num, stream );
    if ( 0 ) {
        char *temp = (char *) &ptr;
        temp++;
    }
}


// Check if the file exists
bool fileExists( const std::string &filename )
{
    std::ifstream ifile( filename.c_str() );
    return ifile.good();
}


// Get the path to a file
std::string IO::getPath( const std::string &filename )
{
    std::string file( filename );
    size_t k1 = file.rfind( 47 );
    size_t k2 = file.rfind( 92 );
    if ( k1 == std::string::npos ) {
        k1 = 0;
    }
    if ( k2 == std::string::npos ) {
        k2 = 0;
    }
    return file.substr( 0, std::max( k1, k2 ) );
}


// List the timesteps in the given directory (dumps.LBPM)
std::vector<std::string> IO::readTimesteps( const std::string &path, const std::string &format )
{
    // Get the name of the summary filename
    std::string filename = path + "/";
    if ( format == "old" || format == "new" ) {
        filename += "summary.LBM";
    } else if ( format == "silo" ) {
        filename += "LBM.visit";
    } else if ( format == "hdf5" ) {
        filename += "LBM.visit";
    } else if ( format == "auto" ) {
        bool test_old = fileExists( path + "/summary.LBM" );
        bool test_new = fileExists( path + "/LBM.visit" );
        if ( test_old && test_new ) {
            ERROR( "Unable to determine format (both summary.LBM and LBM.visit exist)" );
        } else if ( test_old ) {
            filename += "summary.LBM";
        } else if ( test_new ) {
            filename += "LBM.visit";
        } else {
            ERROR( "Unable to determine format (neither summary.LBM or LBM.visit exist)" );
        }
    } else {
        ERROR( "Unknown format: " + format );
    }
    PROFILE_START( "readTimesteps" );
    // Read the data
    FILE *fid = fopen( filename.c_str(), "rb" );
    if ( fid == NULL )
        ERROR( "Error opening file" );
    std::vector<std::string> timesteps;
    char buf[1000];
    while ( fgets( buf, sizeof( buf ), fid ) != NULL ) {
        std::string line( buf );
        line.resize( line.size() - 1 );
        auto pos = line.find( "summary.silo" );
        if ( pos != std::string::npos )
            line.resize( pos );
        pos = line.find( "summary.xmf" );
        if ( pos != std::string::npos )
            line.resize( pos );
        if ( line.empty() )
            continue;
        timesteps.push_back( line );
    }
    fclose( fid );
    PROFILE_STOP( "readTimesteps" );
    return timesteps;
    return timesteps;
}


// Get the maximum number of domains
int IO::maxDomains( const std::string &path, const std::string &format, const Utilities::MPI &comm )
{
    int rank      = comm.getRank();
    int n_domains = 0;
    if ( rank == 0 ) {
        // Get the timesteps
        auto timesteps = IO::readTimesteps( path, format );
        ASSERT( !timesteps.empty() );
        // Get the database for the first domain
        auto db = IO::getMeshList( path, timesteps[0] );
        for ( size_t i = 0; i < db.size(); i++ )
            n_domains = std::max<int>( n_domains, db[i].domains.size() );
    }
    return comm.bcast( n_domains, 0 );
}


// Read the data for the given timestep
std::vector<IO::MeshDataStruct> IO::readData(
    const std::string &path, const std::string &timestep, int rank )
{
    // Get the mesh databases
    auto db = IO::getMeshList( path, timestep );
    // Create the data
    std::vector<IO::MeshDataStruct> data( db.size() );
    for ( size_t i = 0; i < data.size(); i++ ) {
        data[i].precision = IO::DataType::Double;
        data[i].meshName  = db[i].name;
        data[i].mesh      = getMesh( path, timestep, db[i], rank );
        data[i].vars.resize( db[i].variables.size() );
        for ( size_t j = 0; j < db[i].variables.size(); j++ )
            data[i].vars[j] = getVariable( path, timestep, db[i], rank, db[i].variables[j].name );
        INSIST( data[i].check(), "Failed check of " + data[i].meshName );
    }
    return data;
}


// Read the list of variables for the given timestep
std::vector<IO::MeshDatabase> IO::getMeshList(
    const std::string &path, const std::string &timestep )
{
    std::string filename = path + "/" + timestep + "/LBM.summary";
    return IO::read( filename );
}


// Read the given mesh domain
std::shared_ptr<IO::Mesh> IO::getMesh( const std::string &path, const std::string &timestep,
    const IO::MeshDatabase &meshDatabase, int domain )
{
    PROFILE_START( "getMesh" );
    std::shared_ptr<IO::Mesh> mesh;
    if ( meshDatabase.format == FileFormat::OLD ) {
        // Old format (binary doubles)
        std::string filename = path + "/" + timestep + "/" + meshDatabase.domains[domain].file;
        FILE *fid            = fopen( filename.c_str(), "rb" );
        INSIST( fid != NULL, "Error opening file" );
        fseek( fid, 0, SEEK_END );
        size_t bytes = ftell( fid );
        size_t N_max = bytes / sizeof( double ) + 1000;
        double *data = new double[N_max];
        fseek( fid, 0, SEEK_SET );
        size_t count = fread( data, sizeof( double ), N_max, fid );
        fclose( fid );
        if ( count % 3 != 0 )
            ERROR( "Error reading file" );
        if ( meshDatabase.type == IO::MeshType::PointMesh ) {
            size_t N              = count / 3;
            auto pointlist        = std::make_shared<PointList>( N );
            std::vector<Point> &P = pointlist->points;
            for ( size_t i = 0; i < N; i++ ) {
                P[i].x = data[3 * i + 0];
                P[i].y = data[3 * i + 1];
                P[i].z = data[3 * i + 2];
            }
            mesh = pointlist;
        } else if ( meshDatabase.type == IO::MeshType::SurfaceMesh ) {
            if ( count % 9 != 0 )
                ERROR( "Error reading file (2)" );
            size_t N_tri          = count / 9;
            auto trilist          = std::make_shared<TriList>( N_tri );
            std::vector<Point> &A = trilist->A;
            std::vector<Point> &B = trilist->B;
            std::vector<Point> &C = trilist->C;
            for ( size_t i = 0; i < N_tri; i++ ) {
                A[i].x = data[9 * i + 0];
                A[i].y = data[9 * i + 1];
                A[i].z = data[9 * i + 2];
                B[i].x = data[9 * i + 3];
                B[i].y = data[9 * i + 4];
                B[i].z = data[9 * i + 5];
                C[i].x = data[9 * i + 6];
                C[i].y = data[9 * i + 7];
                C[i].z = data[9 * i + 8];
            }
            mesh = trilist;
        } else if ( meshDatabase.type == IO::MeshType::VolumeMesh ) {
            // this was never supported in the old format
            mesh = std::make_shared<DomainMesh>();
        } else {
            ERROR( "Unknown mesh type" );
        }
        delete[] data;
    } else if ( meshDatabase.format == FileFormat::NEW ||
                meshDatabase.format == FileFormat::NEW_SINGLE ) {
        const DatabaseEntry &database = meshDatabase.domains[domain];
        std::string filename          = path + "/" + timestep + "/" + database.file;
        FILE *fid                     = fopen( filename.c_str(), "rb" );
        fseek( fid, database.offset, SEEK_SET );
        char line[1000];
        fgetl( line, 1000, fid );
        size_t i1    = find( line, ':' );
        size_t i2    = find( &line[i1 + 1], ':' ) + i1 + 1;
        size_t bytes = atol( &line[i2 + 1] );
        char *data   = new char[bytes];
        size_t count = fread( data, 1, bytes, fid );
        fclose( fid );
        ASSERT( count == bytes );
        if ( meshDatabase.meshClass == "PointList" ) {
            mesh = std::make_shared<IO::PointList>();
        } else if ( meshDatabase.meshClass == "TriMesh" ) {
            mesh = std::make_shared<IO::TriMesh>();
        } else if ( meshDatabase.meshClass == "TriList" ) {
            mesh = std::make_shared<IO::TriList>();
        } else if ( meshDatabase.meshClass == "DomainMesh" ) {
            mesh = std::make_shared<IO::DomainMesh>();
        } else {
            ERROR( "Unknown mesh class" );
        }
        mesh->unpack( std::pair<size_t, void *>( bytes, data ) );
        delete[] data;
    } else if ( meshDatabase.format == FileFormat::SILO ) {
        // Reading a silo file
#ifdef USE_SILO
        const DatabaseEntry &database = meshDatabase.domains[domain];
        std::string filename          = path + "/" + timestep + "/" + database.file;
        auto fid                      = silo::open( filename, silo::READ );
        if ( meshDatabase.meshClass == "PointList" ) {
            Array<double> coords = silo::readPointMesh<double>( fid, database.name );
            ASSERT( coords.size( 1 ) == 3 );
            auto mesh2 = std::make_shared<IO::PointList>( coords.size( 0 ) );
            for ( size_t i = 0; i < coords.size( 1 ); i++ ) {
                mesh2->points[i].x = coords( i, 0 );
                mesh2->points[i].y = coords( i, 1 );
                mesh2->points[i].z = coords( i, 2 );
            }
            mesh = mesh2;
        } else if ( meshDatabase.meshClass == "TriMesh" || meshDatabase.meshClass == "TriList" ) {
            Array<double> coords;
            Array<int> tri;
            silo::readTriMesh( fid, database.name, coords, tri );
            ASSERT( tri.size( 1 ) == 3 && coords.size( 1 ) == 3 );
            int N_tri   = tri.size( 0 );
            int N_point = coords.size( 0 );
            auto mesh2  = std::make_shared<IO::TriMesh>( N_tri, N_point );
            for ( int i = 0; i < N_point; i++ ) {
                mesh2->vertices->points[i].x = coords( i, 0 );
                mesh2->vertices->points[i].y = coords( i, 1 );
                mesh2->vertices->points[i].z = coords( i, 2 );
            }
            for ( int i = 0; i < N_tri; i++ ) {
                mesh2->A[i] = tri( i, 0 );
                mesh2->B[i] = tri( i, 1 );
                mesh2->C[i] = tri( i, 2 );
            }
            if ( meshDatabase.meshClass == "TriMesh" ) {
                mesh = mesh2;
            } else if ( meshDatabase.meshClass == "TriList" ) {
                auto trilist = IO::getTriList( std::dynamic_pointer_cast<IO::Mesh>( mesh2 ) );
                mesh         = trilist;
            }
        } else if ( meshDatabase.meshClass == "DomainMesh" ) {
            std::vector<double> range;
            std::vector<int> N;
            silo::readUniformMesh( fid, database.name, range, N );
            auto rankinfo = silo::read<int>( fid, database.name + "_rankinfo" );
            RankInfoStruct rank_data( rankinfo[0], rankinfo[1], rankinfo[2], rankinfo[3] );
            mesh = std::make_shared<IO::DomainMesh>( rank_data, N[0], N[1], N[2],
                range[1] - range[0], range[3] - range[2], range[5] - range[4] );
        } else {
            ERROR( "Unknown mesh class" );
        }
        silo::close( fid );
#else
        ERROR( "Build without silo support" );
#endif
    } else if ( meshDatabase.format == FileFormat::HDF5 ) {
        // Reading an hdf5 file
#ifdef USE_HDF5
        auto &database = meshDatabase.domains[domain];
        auto filename  = path + "/" + timestep + "/" + database.file;
        auto fid       = IO::HDF5::openHDF5( filename, "r" );
        auto gid       = IO::HDF5::openGroup( fid, database.name );
        if ( meshDatabase.meshClass == "PointList" ) {
            std::vector<double> x, y, z;
            IO::HDF5::readHDF5( gid, "x", x );
            IO::HDF5::readHDF5( gid, "y", y );
            IO::HDF5::readHDF5( gid, "z", z );
            ASSERT( y.size() == x.size() && z.size() == x.size() );
            auto mesh2 = std::make_shared<IO::PointList>( x.size() );
            for ( size_t i = 0; i < x.size(); i++ ) {
                mesh2->points[i].x = x[i];
                mesh2->points[i].y = y[i];
                mesh2->points[i].z = z[i];
            }
            mesh = mesh2;
        } else if ( meshDatabase.meshClass == "TriMesh" || meshDatabase.meshClass == "TriList" ) {
            // Read the points
            std::vector<double> x, y, z;
            IO::HDF5::readHDF5( gid, "x", x );
            IO::HDF5::readHDF5( gid, "y", y );
            IO::HDF5::readHDF5( gid, "z", z );
            // Read the triangles
            Array<int> tri;
            IO::HDF5::readHDF5( gid, "tri", tri );
            ASSERT( tri.size( 0 ) == 3 );
            size_t N_tri   = tri.size( 1 );
            size_t N_point = x.size();
            auto mesh2     = std::make_shared<IO::TriMesh>( N_tri, N_point );
            for ( size_t i = 0; i < N_point; i++ ) {
                mesh2->vertices->points[i].x = x[i];
                mesh2->vertices->points[i].y = y[i];
                mesh2->vertices->points[i].z = z[i];
            }
            for ( size_t i = 0; i < N_tri; i++ ) {
                mesh2->A[i] = tri( 0, i );
                mesh2->B[i] = tri( 1, i );
                mesh2->C[i] = tri( 2, i );
            }
            if ( meshDatabase.meshClass == "TriMesh" ) {
                mesh = mesh2;
            } else if ( meshDatabase.meshClass == "TriList" ) {
                auto trilist = IO::getTriList( std::dynamic_pointer_cast<IO::Mesh>( mesh2 ) );
                mesh         = trilist;
            }
        } else if ( meshDatabase.meshClass == "DomainMesh" ) {
            std::vector<double> range;
            std::vector<int> N;
            std::vector<int> rankinfo;
            IO::HDF5::readHDF5( gid, "range", range );
            IO::HDF5::readHDF5( gid, "N", N );
            IO::HDF5::readHDF5( gid, "rankinfo", rankinfo );
            RankInfoStruct rank_data( rankinfo[0], rankinfo[1], rankinfo[2], rankinfo[3] );
            mesh = std::make_shared<IO::DomainMesh>( rank_data, N[0], N[1], N[2],
                range[1] - range[0], range[3] - range[2], range[5] - range[4] );
        } else {
            ERROR( "Unknown mesh class" );
        }
        IO::HDF5::closeGroup( gid );
        IO::HDF5::closeHDF5( fid );
#else
        ERROR( "Build without hdf5 support" );
#endif
    } else {
        ERROR( "Unknown format" );
    }
    PROFILE_STOP( "getMesh" );
    return mesh;
}


// Read the given variable for the given mesh domain
std::shared_ptr<IO::Variable> IO::getVariable( const std::string &path, const std::string &timestep,
    const MeshDatabase &meshDatabase, int domain, const std::string &variable )
{
    std::pair<std::string, std::string> key( meshDatabase.domains[domain].name, variable );
    auto it = meshDatabase.variable_data.find( key );
    if ( it == meshDatabase.variable_data.end() )
        return std::shared_ptr<IO::Variable>();
    std::shared_ptr<IO::Variable> var;
    if ( meshDatabase.format == FileFormat::NEW || meshDatabase.format == FileFormat::NEW_SINGLE ) {
        const DatabaseEntry &database = it->second;
        std::string filename          = path + "/" + timestep + "/" + database.file;
        FILE *fid                     = fopen( filename.c_str(), "rb" );
        fseek( fid, database.offset, SEEK_SET );
        char line[1000];
        fgetl( line, 1000, fid );
        size_t i1                       = find( line, ':' );
        size_t i2                       = find( &line[i1 + 1], ':' ) + i1 + 1;
        std::vector<std::string> values = splitList( &line[i2 + 1], ',' );
        ASSERT( values.size() == 5 );
        int dim               = atoi( values[0].c_str() );
        auto type             = values[1];
        size_t N              = atol( values[2].c_str() );
        size_t bytes          = atol( values[3].c_str() );
        std::string precision = values[4];
        var                   = std::make_shared<IO::Variable>();
        var->dim              = dim;
        var->type             = getVariableType( type );
        var->name             = variable;
        var->data.resize( N, dim );
        if ( precision == "double" ) {
            size_t count = fread( var->data.data(), sizeof( double ), N * dim, fid );
            ASSERT( count * sizeof( double ) == bytes );
        } else {
            ERROR( "Format not implimented" );
        }
        fclose( fid );
    } else if ( meshDatabase.format == FileFormat::SILO ) {
        // Reading a silo file
#ifdef USE_SILO
        const auto &database  = meshDatabase.domains[domain];
        auto variableDatabase = meshDatabase.getVariableDatabase( variable );
        std::string filename  = path + "/" + timestep + "/" + database.file;
        auto fid              = silo::open( filename, silo::READ );
        var = std::make_shared<Variable>( variableDatabase.dim, variableDatabase.type, variable );
        if ( meshDatabase.meshClass == "PointList" ) {
            var->data = silo::readPointMeshVariable<double>( fid, variable );
        } else if ( meshDatabase.meshClass == "TriMesh" || meshDatabase.meshClass == "TriList" ) {
            var->data = silo::readTriMeshVariable<double>( fid, variable );
        } else if ( meshDatabase.meshClass == "DomainMesh" ) {
            var->data = silo::readUniformMeshVariable<double>( fid, variable );
        } else {
            ERROR( "Unknown mesh class" );
        }
        silo::close( fid );
#else
        ERROR( "Build without silo support" );
#endif
    } else if ( meshDatabase.format == FileFormat::HDF5 ) {
        // Reading an hdf5 file
#ifdef USE_HDF5
        auto &database   = meshDatabase.domains[domain];
        auto varDatabase = meshDatabase.getVariableDatabase( variable );
        auto filename    = path + "/" + timestep + "/" + database.file;
        var      = std::make_shared<Variable>( varDatabase.dim, varDatabase.type, variable );
        auto fid = IO::HDF5::openHDF5( filename, "r" );
        auto gid = IO::HDF5::openGroup( fid, database.name );
        IO::HDF5::readHDF5( gid, var->name, var->data );
        IO::HDF5::closeHDF5( fid );
        if ( meshDatabase.meshClass == "PointList" || meshDatabase.meshClass == "TriMesh" ||
             meshDatabase.meshClass == "TriList" ) {
            if ( var->data.ndim() == 2 && var->data.size( 0 ) == 3 )
                var->data = var->data.permute( { 1, 0 } );
        } else if ( meshDatabase.meshClass == "DomainMesh" ) {
            if ( var->data.ndim() == 4 && var->data.size( 0 ) == 3 )
                var->data = var->data.permute( { 1, 2, 3, 0 } );
        } else {
            ERROR( "Unknown mesh class" );
        }
#else
        ERROR( "Build without silo support" );
#endif
    } else {
        ERROR( "Unknown format" );
    }
    return var;
}


/****************************************************
 * Reformat the variable to match the mesh           *
 ****************************************************/
void IO::reformatVariable( const IO::Mesh &mesh, IO::Variable &var )
{
    if ( mesh.className() == "DomainMesh" ) {
        const IO::DomainMesh &mesh2 = dynamic_cast<const IO::DomainMesh &>( mesh );
        if ( var.type == VariableType::NodeVariable ) {
            size_t N2 =
                var.data.length() / ( ( mesh2.nx + 1 ) * ( mesh2.ny + 1 ) * ( mesh2.nz + 1 ) );
            ASSERT(
                ( mesh2.nx + 1 ) * ( mesh2.ny + 1 ) * ( mesh2.nz + 1 ) * N2 == var.data.length() );
            var.data.reshape(
                { (size_t) mesh2.nx + 1, (size_t) mesh2.ny + 1, (size_t) mesh2.nz + 1, N2 } );
        } else if ( var.type == VariableType::EdgeVariable ) {
            ERROR( "Not finished" );
        } else if ( var.type == VariableType::SurfaceVariable ) {
            ERROR( "Not finished" );
        } else if ( var.type == VariableType::VolumeVariable ) {
            size_t N2 = var.data.length() / ( mesh2.nx * mesh2.ny * mesh2.nz );
            ASSERT( mesh2.nx * mesh2.ny * mesh2.nz * N2 == var.data.length() );
            var.data.reshape( { (size_t) mesh2.nx, (size_t) mesh2.ny, (size_t) mesh2.nz, N2 } );
        } else {
            ERROR( "Invalid variable type" );
        }
    } else if ( mesh.className() == "PointList" ) {
        const IO::PointList &mesh2 = dynamic_cast<const IO::PointList &>( mesh );
        size_t N                   = mesh2.points.size();
        size_t N_var               = var.data.length() / N;
        ASSERT( N * N_var == var.data.length() );
        var.data.reshape( { N, N_var } );
    } else if ( mesh.className() == "TriMesh" || mesh.className() == "TriList" ) {
        std::shared_ptr<Mesh> mesh_ptr( const_cast<Mesh *>( &mesh ), []( void * ) {} );
        std::shared_ptr<TriMesh> mesh2 = getTriMesh( mesh_ptr );
        if ( var.type == VariableType::NodeVariable ) {
            size_t N     = mesh2->vertices->points.size();
            size_t N_var = var.data.length() / N;
            ASSERT( N * N_var == var.data.length() );
            var.data.reshape( { N, N_var } );
        } else if ( var.type == VariableType::EdgeVariable ) {
            ERROR( "Not finished" );
        } else if ( var.type == VariableType::SurfaceVariable ) {
            ERROR( "Not finished" );
        } else if ( var.type == VariableType::VolumeVariable ) {
            size_t N     = mesh2->A.size();
            size_t N_var = var.data.length() / N;
            ASSERT( N * N_var == var.data.length() );
            var.data.reshape( { N, N_var } );
        } else {
            ERROR( "Invalid variable type" );
        }
    } else {
        ERROR( "Unknown mesh type" );
    }
}
