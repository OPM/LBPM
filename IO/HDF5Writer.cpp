#include "IO/HDF5_IO.h"
#include "IO/IOHelpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Writer.h"
#include "IO/Xdmf.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include <algorithm>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <vector>


#ifdef USE_HDF5


std::string to_string( const ArraySize &s )
{
    std::string out = "[" + std::to_string( s[0] );
    for ( size_t i = 1; i < s.ndim(); i++ )
        out += "," + to_string( s[i] );
    out += "]";
    return out;
}


Xdmf::Center getXdmfType( IO::VariableType type )
{
    if ( type == IO::VariableType::NodeVariable ) {
        return Xdmf::Center::Node;
    } else if ( type == IO::VariableType::VolumeVariable ) {
        return Xdmf::Center::Cell;
    } else {
        ERROR( "Variable type not supported" );
    }
    return Xdmf::Center::Null;
}


// Write a PointList mesh (and variables) to a file
template<class TYPE>
static void writeCoordinates( hid_t fid, const std::vector<Point> &points )
{
    std::vector<TYPE> x( points.size() ), y( points.size() ), z( points.size() );
    for ( size_t i = 0; i < x.size(); i++ ) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        z[i] = points[i].z;
    }
    IO::HDF5::writeHDF5( fid, "x", x );
    IO::HDF5::writeHDF5( fid, "y", y );
    IO::HDF5::writeHDF5( fid, "z", z );
}
static void writeHDF5PointList( hid_t fid, const std::string &filename,
    const IO::MeshDataStruct &meshData, IO::MeshDatabase database, Xdmf &xmf )
{
    auto meshname    = database.domains[0].name;
    const auto &mesh = dynamic_cast<IO::PointList &>( *meshData.mesh );
    auto gid         = IO::HDF5::createGroup( fid, meshname );
    if ( meshData.precision == IO::DataType::Double ) {
        writeCoordinates<double>( gid, mesh.getPoints() );
    } else if ( meshData.precision == IO::DataType::Float ) {
        writeCoordinates<float>( gid, mesh.getPoints() );
    } else {
        ERROR( "Unsupported format" );
    }
    auto path   = filename + ":/" + meshname + "/";
    auto domain = Xdmf::createPointMesh(
        meshname, 3, mesh.getPoints().size(), path + "x", path + "y", path + "z" );
    for ( size_t i = 0; i < meshData.vars.size(); i++ ) {
        auto &var     = *meshData.vars[i];
        auto data     = var.data;
        auto rankType = Xdmf::RankType::Null;
        if ( data.ndim() == 1 ) {
            rankType = Xdmf::RankType::Scalar;
        } else if ( data.ndim() == 2 && data.size( 1 ) == 3 ) {
            // Vector data, need to permute for visit
            rankType = Xdmf::RankType::Vector;
            data     = data.permute( { 1, 0 } );
        } else {
            ERROR( "Unable to determine variable rank: " + to_string( var.data.size() ) );
        }
        if ( var.precision == IO::DataType::Double ) {
            IO::HDF5::writeHDF5( gid, var.name, data );
        } else if ( var.precision == IO::DataType::Float ) {
            IO::HDF5::writeHDF5( gid, var.name, data.cloneTo<float>() );
        } else if ( var.precision == IO::DataType::Int ) {
            IO::HDF5::writeHDF5( gid, var.name, data.cloneTo<int>() );
        } else {
            ERROR( "Unsupported format" );
        }
        domain.addVariable(
            meshname, var.name, data.size(), rankType, Xdmf::Center::Node, path + var.name );
    }
    xmf.addMesh( meshData.meshName, domain );
}
// Write a TriMesh mesh (and variables) to a file
static void writeHDF5TriMesh2( hid_t fid, const std::string &filename,
    const IO::MeshDataStruct &meshData, const IO::TriMesh &mesh, IO::MeshDatabase database,
    Xdmf &xmf )
{
    auto meshname = database.domains[0].name;
    auto gid      = IO::HDF5::createGroup( fid, meshname );
    auto path     = filename + ":/" + meshname + "/";
    // Write the verticies
    if ( meshData.precision == IO::DataType::Double ) {
        writeCoordinates<double>( gid, mesh.vertices->getPoints() );
    } else if ( meshData.precision == IO::DataType::Float ) {
        writeCoordinates<float>( gid, mesh.vertices->getPoints() );
    } else {
        ERROR( "Unsupported format" );
    }
    // Write the connectivity
    Array<int> tri( 3, mesh.A.size() );
    for ( size_t i = 0; i < mesh.A.size(); i++ ) {
        tri( 0, i ) = mesh.A[i];
        tri( 1, i ) = mesh.B[i];
        tri( 2, i ) = mesh.C[i];
    }
    IO::HDF5::writeHDF5( gid, "tri", tri );
    auto domain =
        Xdmf::createUnstructuredMesh( meshname, 3, Xdmf::TopologyType::Triangle, tri.size( 1 ),
            path + "tri", mesh.vertices->getPoints().size(), path + "x", path + "y", path + "z" );
    // Write the variables
    for ( size_t i = 0; i < meshData.vars.size(); i++ ) {
        auto &var     = *meshData.vars[i];
        auto data     = var.data;
        auto rankType = Xdmf::RankType::Null;
        if ( data.ndim() == 1 ) {
            rankType = Xdmf::RankType::Scalar;
        } else if ( data.ndim() == 2 && data.size( 1 ) == 3 ) {
            // Vector data, need to permute for visit
            rankType = Xdmf::RankType::Vector;
            data     = data.permute( { 1, 0 } );
        } else {
            ERROR( "Unable to determine variable rank: " + to_string( var.data.size() ) );
        }
        if ( var.precision == IO::DataType::Double ) {
            IO::HDF5::writeHDF5( gid, var.name, data );
        } else if ( var.precision == IO::DataType::Float ) {
            IO::HDF5::writeHDF5( gid, var.name, data.cloneTo<float>() );
        } else if ( var.precision == IO::DataType::Int ) {
            IO::HDF5::writeHDF5( gid, var.name, data.cloneTo<int>() );
        } else {
            ERROR( "Unsupported format" );
        }
        domain.addVariable(
            meshname, var.name, data.size(), rankType, getXdmfType( var.type ), path + var.name );
    }
    xmf.addMesh( meshData.meshName, domain );
}
static void writeHDF5TriMesh( hid_t fid, const std::string &filename,
    const IO::MeshDataStruct &meshData, IO::MeshDatabase database, Xdmf &xmf )
{
    const IO::TriMesh &mesh = dynamic_cast<IO::TriMesh &>( *meshData.mesh );
    writeHDF5TriMesh2( fid, filename, meshData, mesh, database, xmf );
}
static void writeHDF5TriList( hid_t fid, const std::string &filename,
    const IO::MeshDataStruct &meshData, IO::MeshDatabase database, Xdmf &xmf )
{
    auto mesh = getTriMesh( meshData.mesh );
    writeHDF5TriMesh2( fid, filename, meshData, *mesh, database, xmf );
}
// Write a DomainMesh mesh (and variables) to a file
static void writeHDF5DomainMesh( hid_t fid, const std::string &filename,
    const IO::MeshDataStruct &meshData, IO::MeshDatabase database, Xdmf &xmf )
{
    auto &mesh    = dynamic_cast<IO::DomainMesh &>( *meshData.mesh );
    auto meshname = database.domains[0].name;
    auto gid      = IO::HDF5::createGroup( fid, meshname );
    auto path     = filename + ":/" + meshname + "/";
    // Write the mesh
    RankInfoStruct info( mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz );
    std::vector<double> range = { info.ix * mesh.Lx / info.nx, ( info.ix + 1 ) * mesh.Lx / info.nx,
        info.jy * mesh.Ly / info.ny, ( info.jy + 1 ) * mesh.Ly / info.ny,
        info.kz * mesh.Lz / info.nz, ( info.kz + 1 ) * mesh.Lz / info.nz };
    std::vector<int> N        = { mesh.nx, mesh.ny, mesh.nz };
    std::vector<int> rankinfo = { mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz };
    IO::HDF5::writeHDF5( gid, "range", range );
    IO::HDF5::writeHDF5( gid, "N", N );
    IO::HDF5::writeHDF5( gid, "rankinfo", rankinfo );
    // xmf.addUniformMesh( meshname, range, ArraySize( N[0], N[1], N[2] ) );
    // Write a curvilinear mesh due to bug with vector data on nodes loading into visit
    Array<float> x( N[0] + 1, N[1] + 1, N[2] + 1 );
    Array<float> y( N[0] + 1, N[1] + 1, N[2] + 1 );
    Array<float> z( N[0] + 1, N[1] + 1, N[2] + 1 );
    double dx = ( range[1] - range[0] ) / N[0];
    double dy = ( range[3] - range[2] ) / N[1];
    double dz = ( range[5] - range[4] ) / N[2];
    for ( int k = 0; k <= N[2]; k++ ) {
        for ( int j = 0; j <= N[1]; j++ ) {
            for ( int i = 0; i <= N[0]; i++ ) {
                x( i, j, k ) = range[0] + dx * i;
                y( i, j, k ) = range[2] + dy * j;
                z( i, j, k ) = range[4] + dz * k;
            }
        }
    }
    IO::HDF5::writeHDF5( gid, "x", x );
    IO::HDF5::writeHDF5( gid, "y", y );
    IO::HDF5::writeHDF5( gid, "z", z );
    auto domain = Xdmf::createCurvilinearMesh(
        meshname, ArraySize( N[0], N[1], N[2] ), path + "x", path + "y", path + "z" );
    // Write the variables
    for ( size_t i = 0; i < meshData.vars.size(); i++ ) {
        auto &var     = *meshData.vars[i];
        auto data     = var.data;
        auto rankType = Xdmf::RankType::Null;
        if ( data.ndim() == 3 ) {
            rankType = Xdmf::RankType::Scalar;
        } else if ( data.ndim() == 4 && data.size( 3 ) == 3 ) {
            // Vector data, need to permute for visit
            rankType = Xdmf::RankType::Vector;
            data     = data.permute( { 3, 0, 1, 2 } );
        } else {
            ERROR( "Unable to determine variable rank: " + to_string( var.data.size() ) );
        }
        if ( var.precision == IO::DataType::Double ) {
            IO::HDF5::writeHDF5( gid, var.name, data );
        } else if ( var.precision == IO::DataType::Float ) {
            IO::HDF5::writeHDF5( gid, var.name, data.cloneTo<float>() );
        } else if ( var.precision == IO::DataType::Int ) {
            IO::HDF5::writeHDF5( gid, var.name, data.cloneTo<int>() );
        } else {
            ERROR( "Unsupported format" );
        }
        domain.addVariable(
            meshname, var.name, data.size(), rankType, getXdmfType( var.type ), path + var.name );
    }
    IO::HDF5::closeGroup( gid );
    xmf.addMesh( meshData.meshName, domain );
}
// Write a mesh (and variables) to a file
static IO::MeshDatabase write_domain_hdf5( hid_t fid, const std::string &filename,
    const IO::MeshDataStruct &mesh, IO::FileFormat format, int rank, Xdmf &xmf )
{
    // Create the MeshDatabase
    auto database = getDatabase( filename, mesh, format, rank );
    if ( database.meshClass == "PointList" ) {
        writeHDF5PointList( fid, filename, mesh, database, xmf );
    } else if ( database.meshClass == "TriMesh" ) {
        writeHDF5TriMesh( fid, filename, mesh, database, xmf );
    } else if ( database.meshClass == "TriList" ) {
        writeHDF5TriList( fid, filename, mesh, database, xmf );
    } else if ( database.meshClass == "DomainMesh" ) {
        writeHDF5DomainMesh( fid, filename, mesh, database, xmf );
    } else {
        ERROR( "Unknown mesh class" );
    }
    return database;
}
// Write the mesh data to hdf5
std::vector<IO::MeshDatabase> writeMeshesHDF5( const std::vector<IO::MeshDataStruct> &meshData,
    const std::string &path, IO::FileFormat format, int rank, Xdmf &xmf )
{

    std::vector<IO::MeshDatabase> meshes_written;
    char filename[100], fullpath[200];
    sprintf( filename, "%05i.h5", rank );
    sprintf( fullpath, "%s/%s", path.c_str(), filename );
    auto fid = IO::HDF5::openHDF5( fullpath, "w", IO::HDF5::Compression::GZIP );
    for ( size_t i = 0; i < meshData.size(); i++ ) {
        meshes_written.push_back(
            write_domain_hdf5( fid, filename, meshData[i], format, rank, xmf ) );
    }
    IO::HDF5::closeHDF5( fid );
    return meshes_written;
}


#else


std::vector<IO::MeshDatabase> writeMeshesHDF5(
    const std::vector<IO::MeshDataStruct> &, const std::string &, IO::FileFormat, int, Xdmf & )
{
    return std::vector<IO::MeshDatabase>();
}


#endif
