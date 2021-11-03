#include "IO/HDF5_IO.h"
#include "IO/IOHelpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Writer.h"
#include "IO/silo.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include <algorithm>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <vector>


#ifdef USE_SILO


// Write a PointList mesh (and variables) to a file
template<class TYPE>
static void writeSiloPointMesh(
    DBfile *fid, const IO::PointList &mesh, const std::string &meshname )
{
    const auto &points = mesh.getPoints();
    std::vector<TYPE> x( points.size() ), y( points.size() ), z( points.size() );
    for ( size_t i = 0; i < x.size(); i++ ) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        z[i] = points[i].z;
    }
    const TYPE *coords[] = { x.data(), y.data(), z.data() };
    IO::silo::writePointMesh<TYPE>( fid, meshname, 3, points.size(), coords );
}
static void writeSiloPointList(
    DBfile *fid, const IO::MeshDataStruct &meshData, IO::MeshDatabase database )
{
    const IO::PointList &mesh  = dynamic_cast<IO::PointList &>( *meshData.mesh );
    const std::string meshname = database.domains[0].name;
    if ( meshData.precision == IO::DataType::Double ) {
        writeSiloPointMesh<double>( fid, mesh, meshname );
    } else if ( meshData.precision == IO::DataType::Float ) {
        writeSiloPointMesh<float>( fid, mesh, meshname );
    } else {
        ERROR( "Unsupported format" );
    }
    const auto &points = mesh.getPoints();
    std::vector<double> x( points.size() ), y( points.size() ), z( points.size() );
    for ( size_t i = 0; i < x.size(); i++ ) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        z[i] = points[i].z;
    }
    const double *coords[] = { x.data(), y.data(), z.data() };
    IO::silo::writePointMesh( fid, meshname, 3, points.size(), coords );
    for ( size_t i = 0; i < meshData.vars.size(); i++ ) {
        const IO::Variable &var = *meshData.vars[i];
        if ( var.precision == IO::DataType::Double ) {
            IO::silo::writePointMeshVariable( fid, meshname, var.name, var.data );
        } else if ( var.precision == IO::DataType::Float ) {
            Array<float> data2( var.data.size() );
            data2.copy( var.data );
            IO::silo::writePointMeshVariable( fid, meshname, var.name, data2 );
        } else if ( var.precision == IO::DataType::Int ) {
            Array<int> data2( var.data.size() );
            data2.copy( var.data );
            IO::silo::writePointMeshVariable( fid, meshname, var.name, data2 );
        } else {
            ERROR( "Unsupported format" );
        }
    }
}
// Write a TriMesh mesh (and variables) to a file
template<class TYPE>
static void writeSiloTriMesh( DBfile *fid, const IO::TriMesh &mesh, const std::string &meshname )
{
    const auto &points = mesh.vertices->getPoints();
    std::vector<TYPE> x( points.size() ), y( points.size() ), z( points.size() );
    for ( size_t i = 0; i < x.size(); i++ ) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        z[i] = points[i].z;
    }
    const TYPE *coords[] = { x.data(), y.data(), z.data() };
    const int *tri[]     = { mesh.A.data(), mesh.B.data(), mesh.C.data() };
    IO::silo::writeTriMesh<TYPE>( fid, meshname, 3, 2, points.size(), coords, mesh.A.size(), tri );
}
static void writeSiloTriMesh2( DBfile *fid, const IO::MeshDataStruct &meshData,
    const IO::TriMesh &mesh, IO::MeshDatabase database )
{
    const std::string meshname = database.domains[0].name;
    if ( meshData.precision == IO::DataType::Double ) {
        writeSiloTriMesh<double>( fid, mesh, meshname );
    } else if ( meshData.precision == IO::DataType::Float ) {
        writeSiloTriMesh<float>( fid, mesh, meshname );
    } else {
        ERROR( "Unsupported format" );
    }
    for ( size_t i = 0; i < meshData.vars.size(); i++ ) {
        const IO::Variable &var = *meshData.vars[i];
        if ( var.precision == IO::DataType::Double ) {
            IO::silo::writeTriMeshVariable( fid, 3, meshname, var.name, var.data, var.type );
        } else if ( var.precision == IO::DataType::Float ) {
            Array<float> data2( var.data.size() );
            data2.copy( var.data );
            IO::silo::writeTriMeshVariable( fid, 3, meshname, var.name, data2, var.type );
        } else if ( var.precision == IO::DataType::Int ) {
            Array<int> data2( var.data.size() );
            data2.copy( var.data );
            IO::silo::writeTriMeshVariable( fid, 3, meshname, var.name, data2, var.type );
        } else {
            ERROR( "Unsupported format" );
        }
    }
}
static void writeSiloTriMesh(
    DBfile *fid, const IO::MeshDataStruct &meshData, IO::MeshDatabase database )
{
    const IO::TriMesh &mesh = dynamic_cast<IO::TriMesh &>( *meshData.mesh );
    writeSiloTriMesh2( fid, meshData, mesh, database );
}
static void writeSiloTriList(
    DBfile *fid, const IO::MeshDataStruct &meshData, IO::MeshDatabase database )
{
    auto mesh = getTriMesh( meshData.mesh );
    writeSiloTriMesh2( fid, meshData, *mesh, database );
}
// Write a DomainMesh mesh (and variables) to a file
static void writeSiloDomainMesh(
    DBfile *fid, const IO::MeshDataStruct &meshData, IO::MeshDatabase database )
{
    const IO::DomainMesh &mesh = dynamic_cast<IO::DomainMesh &>( *meshData.mesh );
    RankInfoStruct info( mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz );
    std::array<double, 6> range = { info.ix * mesh.Lx / info.nx,
        ( info.ix + 1 ) * mesh.Lx / info.nx, info.jy * mesh.Ly / info.ny,
        ( info.jy + 1 ) * mesh.Ly / info.ny, info.kz * mesh.Lz / info.nz,
        ( info.kz + 1 ) * mesh.Lz / info.nz };
    std::array<int, 3> N        = { mesh.nx, mesh.ny, mesh.nz };
    auto meshname               = database.domains[0].name;
    IO::silo::writeUniformMesh<3>( fid, meshname, range, N );
    IO::silo::write<int>(
        fid, meshname + "_rankinfo", { mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz } );
    for ( size_t i = 0; i < meshData.vars.size(); i++ ) {
        const auto &var = *meshData.vars[i];
        if ( var.precision == IO::DataType::Double ) {
            IO::silo::writeUniformMeshVariable<3>( fid, meshname, N, var.name, var.data, var.type );
        } else if ( var.precision == IO::DataType::Float ) {
            Array<float> data2( var.data.size() );
            data2.copy( var.data );
            IO::silo::writeUniformMeshVariable<3>( fid, meshname, N, var.name, data2, var.type );
        } else if ( var.precision == IO::DataType::Int ) {
            Array<int> data2( var.data.size() );
            data2.copy( var.data );
            IO::silo::writeUniformMeshVariable<3>( fid, meshname, N, var.name, data2, var.type );
        } else {
            ERROR( "Unsupported format" );
        }
    }
}
// Write a mesh (and variables) to a file
static IO::MeshDatabase write_domain_silo( DBfile *fid, const std::string &filename,
    const IO::MeshDataStruct &mesh, IO::FileFormat format, int rank )
{
    // Create the MeshDatabase
    auto database = getDatabase( filename, mesh, format, rank );
    if ( database.meshClass == "PointList" ) {
        writeSiloPointList( fid, mesh, database );
    } else if ( database.meshClass == "TriMesh" ) {
        writeSiloTriMesh( fid, mesh, database );
    } else if ( database.meshClass == "TriList" ) {
        writeSiloTriList( fid, mesh, database );
    } else if ( database.meshClass == "DomainMesh" ) {
        writeSiloDomainMesh( fid, mesh, database );
    } else {
        ERROR( "Unknown mesh class" );
    }
    return database;
}
// Write the summary file for silo
std::pair<int, int> getSiloMeshType( const std::string &meshClass )
{
    int meshType = 0;
    int varType  = 0;
    if ( meshClass == "PointList" ) {
        meshType = DB_POINTMESH;
        varType  = DB_POINTVAR;
    } else if ( meshClass == "TriMesh" ) {
        meshType = DB_UCDMESH;
        varType  = DB_UCDVAR;
    } else if ( meshClass == "TriList" ) {
        meshType = DB_UCDMESH;
        varType  = DB_UCDVAR;
    } else if ( meshClass == "DomainMesh" ) {
        meshType = DB_QUAD_RECT;
        varType  = DB_QUADVAR;
    } else {
        ERROR( "Unknown mesh class" );
    }
    return std::make_pair( meshType, varType );
}
void writeSiloSummary(
    const std::vector<IO::MeshDatabase> &meshes_written, const std::string &filename )
{
    auto fid = IO::silo::open( filename, IO::silo::CREATE );
    for ( const auto &data : meshes_written ) {
        auto type = getSiloMeshType( data.meshClass );
        std::vector<int> meshTypes( data.domains.size(), type.first );
        std::vector<int> varTypes( data.domains.size(), type.second );
        std::vector<std::string> meshNames;
        for ( const auto &tmp : data.domains )
            meshNames.push_back( tmp.file + ":" + tmp.name );
        IO::silo::writeMultiMesh( fid, data.name, meshNames, meshTypes );
        for ( const auto &variable : data.variables ) {
            std::vector<std::string> varnames;
            for ( const auto &tmp : data.domains )
                varnames.push_back( tmp.file + ":" + variable.name );
            IO::silo::writeMultiVar( fid, variable.name, varnames, varTypes );
        }
    }
    IO::silo::close( fid );
}
// Write the mesh data to silo
std::vector<IO::MeshDatabase> writeMeshesSilo( const std::vector<IO::MeshDataStruct> &meshData,
    const std::string &path, IO::FileFormat format, int rank )
{
    std::vector<IO::MeshDatabase> meshes_written;
    char filename[100], fullpath[200];
    sprintf( filename, "%05i.silo", rank );
    sprintf( fullpath, "%s/%s", path.c_str(), filename );
    auto fid = IO::silo::open( fullpath, IO::silo::CREATE );
    for ( size_t i = 0; i < meshData.size(); i++ ) {
        auto mesh = meshData[i].mesh;
        meshes_written.push_back( write_domain_silo( fid, filename, meshData[i], format, rank ) );
    }
    IO::silo::close( fid );
    return meshes_written;
}


#else


// Write the mesh data to silo
std::vector<IO::MeshDatabase> writeMeshesSilo(
    const std::vector<IO::MeshDataStruct> &, const std::string &, IO::FileFormat, int )
{
    return std::vector<IO::MeshDatabase>();
}
void writeSiloSummary( const std::vector<IO::MeshDatabase> &, const std::string & )
{
}


#endif
