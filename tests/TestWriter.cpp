#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>

#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "ProfilerApp.h"
#include "common/MPI.h"
#include "common/UnitTest.h"
#include "common/Utilities.h"

inline bool approx_equal( const Point &A, const Point &B )
{
    double tol = 1e-7 * sqrt( A.x * A.x + A.y * A.y + A.z * A.z );
    return fabs( A.x - B.x ) <= tol && fabs( A.y - B.y ) <= tol && fabs( A.z - B.z ) <= tol;
}
inline bool approx_equal( const double &A, const double &B )
{
    return fabs( A - B ) <= std::max<double>( 1e-7 * fabs( A + B ), 1e-20 );
}


inline double distance( const Point &p ) { return sqrt( p.x * p.x + p.y * p.y + p.z * p.z ); }


bool checkMesh( const std::vector<IO::MeshDataStruct> &meshData, const std::string &format,
    std::shared_ptr<IO::Mesh> mesh )
{

    // Get direct access to the meshes used to test the reader
    const auto pointmesh = dynamic_cast<IO::PointList *>( meshData[0].mesh.get() );
    const auto trimesh   = dynamic_cast<IO::TriMesh *>( meshData[1].mesh.get() );
    const auto trilist   = dynamic_cast<IO::TriList *>( meshData[2].mesh.get() );
    const auto domain    = dynamic_cast<IO::DomainMesh *>( meshData[3].mesh.get() );
    const size_t N_tri   = trimesh->A.size();
    if ( mesh->className() == "pointmesh" ) {
        // Check the pointmesh
        auto pmesh = IO::getPointList( mesh );
        if ( pmesh.get() == NULL )
            return false;
        if ( pmesh->points.size() != pointmesh->points.size() )
            return false;
    }
    if ( mesh->className() == "trimesh" || mesh->className() == "trilist" ) {
        // Check the trimesh/trilist
        auto mesh1 = IO::getTriMesh( mesh );
        auto mesh2 = IO::getTriList( mesh );
        if ( mesh1.get() == NULL || mesh2.get() == NULL )
            return false;
        if ( mesh1->A.size() != N_tri || mesh1->B.size() != N_tri || mesh1->C.size() != N_tri ||
             mesh2->A.size() != N_tri || mesh2->B.size() != N_tri || mesh2->C.size() != N_tri )
            return false;
        const std::vector<Point> &P1 = mesh1->vertices->points;
        const std::vector<int> &A1   = mesh1->A;
        const std::vector<int> &B1   = mesh1->B;
        const std::vector<int> &C1   = mesh1->C;
        const std::vector<Point> &A2 = mesh2->A;
        const std::vector<Point> &B2 = mesh2->B;
        const std::vector<Point> &C2 = mesh2->C;
        const std::vector<Point> &A  = trilist->A;
        const std::vector<Point> &B  = trilist->B;
        const std::vector<Point> &C  = trilist->C;
        for ( size_t i = 0; i < N_tri; i++ ) {
            if ( !approx_equal( P1[A1[i]], A[i] ) || !approx_equal( P1[B1[i]], B[i] ) ||
                 !approx_equal( P1[C1[i]], C[i] ) )
                return false;
            if ( !approx_equal( A2[i], A[i] ) || !approx_equal( B2[i], B[i] ) ||
                 !approx_equal( C2[i], C[i] ) )
                return false;
        }
    }
    if ( mesh->className() == "domain" && format != "old" ) {
        // Check the domain mesh
        const IO::DomainMesh &mesh1 = *std::dynamic_pointer_cast<IO::DomainMesh>( mesh );
        if ( mesh1.nprocx != domain->nprocx || mesh1.nprocy != domain->nprocy ||
             mesh1.nprocz != domain->nprocz )
            return false;
        if ( mesh1.nx != domain->nx || mesh1.ny != domain->ny || mesh1.nz != domain->nz )
            return false;
        if ( mesh1.Lx != domain->Lx || mesh1.Ly != domain->Ly || mesh1.Lz != domain->Lz )
            return false;
    }
    return true;
}


bool checkVar( const std::string &format, std::shared_ptr<IO::Mesh> mesh,
    std::shared_ptr<IO::Variable> variable1, std::shared_ptr<IO::Variable> variable2 )
{
    if ( format == "new" )
        IO::reformatVariable( *mesh, *variable2 );
    bool pass        = true;
    const auto &var1 = *variable1;
    const auto &var2 = *variable2;
    pass             = var1.name == var2.name;
    pass             = pass && var1.dim == var2.dim;
    pass             = pass && var1.type == var2.type;
    pass             = pass && var1.data.length() == var2.data.length();
    if ( pass ) {
        for ( size_t m = 0; m < var1.data.length(); m++ )
            pass = pass && approx_equal( var1.data( m ), var2.data( m ) );
    }
    return pass;
}


// Test writing and reading the given format
void testWriter(
    const std::string &format, std::vector<IO::MeshDataStruct> &meshData, UnitTest &ut )
{
    PROFILE_SCOPED( path, 0, timer );

    Utilities::MPI comm( MPI_COMM_WORLD );
    int nprocs = comm.getSize();
    comm.barrier();


    // Set the path for the writer
    std::string path = "test_" + format;


    // Get the format
    std::string format2 = format;
    auto precision      = IO::DataType::Double;
    if ( format == "silo-double" ) {
        format2   = "silo";
        precision = IO::DataType::Double;
    } else if ( format == "silo-float" ) {
        format2   = "silo";
        precision = IO::DataType::Float;
    } else if ( format == "hdf5-double" ) {
        format2   = "hdf5";
        precision = IO::DataType::Double;
    } else if ( format == "hdf5-float" ) {
        format2   = "hdf5";
        precision = IO::DataType::Float;
    }


    // Set the precision for the variables
    for ( auto &data : meshData ) {
        data.precision = precision;
        for ( auto &var : data.vars )
            var->precision = precision;
    }

    // Write the data
    IO::initialize( path, format2, false );
    IO::writeData( 0, meshData, comm );
    IO::writeData( 3, meshData, comm );
    comm.barrier();


    // Get a list of the timesteps
    auto timesteps = IO::readTimesteps( path, format2 );
    if ( timesteps.size() == 2 )
        ut.passes( format + ": Corrent number of timesteps" );
    else
        ut.failure( format + ": Incorrent number of timesteps" );


    // Test the simple read interface
    bool pass = true;
    for ( const auto &timestep : timesteps ) {
        auto data = IO::readData( path, timestep, comm.getRank() );
        pass      = pass && data.size() == meshData.size();
        for ( size_t i = 0; i < data.size(); i++ ) {
            pass = pass && checkMesh( meshData, format, data[i].mesh );
        }
    }
    if ( pass )
        ut.passes( format + ": Simple read interface" );
    else
        ut.failure( format + ": Simple read interface" );


    // Test reading each mesh domain
    for ( const auto &timestep : timesteps ) {
        // Load the list of meshes and check its size
        auto databaseList = IO::getMeshList( path, timestep );
        if ( databaseList.size() == meshData.size() )
            ut.passes( format + ": Corrent number of meshes found" );
        else
            ut.failure( format + ": Incorrent number of meshes found" );
        // Check the number of domains for each mesh
        for ( const auto &database : databaseList ) {
            int N_domains = database.domains.size();
            if ( N_domains != nprocs ) {
                ut.failure( format + ": Incorrent number of domains for mesh" );
                continue;
            }
            // For each domain, load the mesh and check its data
            bool pass = true;
            for ( int k = 0; k < N_domains; k++ ) {
                auto mesh = IO::getMesh( path, timestep, database, k );
                if ( !mesh ) {
                    ut.failure( "Failed to load " + database.name );
                    pass = false;
                } else {
                    pass = pass && checkMesh( meshData, format, mesh );
                }
            }
            if ( pass ) {
                ut.passes( format + ": Mesh \"" + database.name + "\" loaded correctly" );
            } else {
                ut.failure( format + ": Mesh \"" + database.name + "\" did not load correctly" );
                continue;
            }
            // Load the variables and check their data
            if ( format == "old" )
                continue; // Old format does not support variables
            const IO::MeshDataStruct *mesh0 = nullptr;
            for ( size_t k = 0; k < meshData.size(); k++ ) {
                if ( meshData[k].meshName == database.name ) {
                    mesh0 = &meshData[k];
                    break;
                }
            }
            for ( int k = 0; k < N_domains; k++ ) {
                auto mesh = IO::getMesh( path, timestep, database, k );
                for ( size_t v = 0; v < mesh0->vars.size(); v++ ) {
                    PROFILE_START( format + "-read-getVariable" );
                    auto variable =
                        IO::getVariable( path, timestep, database, k, mesh0->vars[v]->name );
                    pass = checkVar( format, mesh, mesh0->vars[v], variable );
                    if ( pass ) {
                        ut.passes( format + ": Variable \"" + variable->name + "\" matched" );
                    } else {
                        ut.failure(
                            format + ": Variable \"" + variable->name + "\" did not match" );
                        break;
                    }
                }
            }
        }
    }
}


// Main
int main( int argc, char **argv )
{
    Utilities::startup( argc, argv );
    Utilities::MPI comm( MPI_COMM_WORLD );
    int rank   = comm.getRank();
    int nprocs = comm.getSize();
    Utilities::setAbortBehavior( true, 2 );
    Utilities::setErrorHandlers();
    UnitTest ut;

    // Create some points
    const int N_points = 8;
    const int N_tri    = 12;
    double x[8]        = { 0, 1, 0, 1, 0, 1, 0, 1 };
    double y[8]        = { 0, 0, 1, 1, 0, 0, 1, 1 };
    double z[8]        = { 0, 0, 0, 0, 1, 1, 1, 1 };
    int tri[N_tri][3]  = {
        { 0, 1, 3 }, { 0, 3, 2 }, { 4, 5, 7 }, { 4, 7, 6 }, // z faces
        { 0, 1, 4 }, { 1, 4, 5 }, { 2, 3, 6 }, { 3, 6, 7 }, // y faces
        { 0, 2, 4 }, { 2, 4, 6 }, { 1, 3, 5 }, { 3, 5, 7 }  // x faces
    };

    // Create the meshes
    auto set1 = std::make_shared<IO::PointList>( N_points );
    for ( int i = 0; i < N_points; i++ ) {
        set1->points[i].x = x[i];
        set1->points[i].y = y[i];
        set1->points[i].z = z[i];
    }
    auto trimesh = std::make_shared<IO::TriMesh>( N_tri, set1 );
    for ( int i = 0; i < N_tri; i++ ) {
        trimesh->A[i] = tri[i][0];
        trimesh->B[i] = tri[i][1];
        trimesh->C[i] = tri[i][2];
    }
    auto trilist = std::make_shared<IO::TriList>( *trimesh );
    for ( int i = 0; i < N_tri; i++ ) {
        Point A( x[tri[i][0]], y[tri[i][0]], z[tri[i][0]] );
        Point B( x[tri[i][1]], y[tri[i][1]], z[tri[i][1]] );
        Point C( x[tri[i][2]], y[tri[i][2]], z[tri[i][2]] );
        if ( !approx_equal( trilist->A[i], A ) || !approx_equal( trilist->B[i], B ) ||
             !approx_equal( trilist->C[i], C ) ) {
            printf( "Failed to create trilist\n" );
            return -1;
        }
    }
    RankInfoStruct rank_data( rank, nprocs, 1, 1 );
    auto domain = std::make_shared<IO::DomainMesh>( rank_data, 6, 7, 8, 1.0, 1.0, 1.0 );

    // Create the variables
    const auto NodeVar   = IO::VariableType::NodeVariable;
    const auto VolVar    = IO::VariableType::VolumeVariable;
    auto set_node_mag    = std::make_shared<IO::Variable>( 1, NodeVar, "Node_set_mag" );
    auto set_node_vec    = std::make_shared<IO::Variable>( 3, NodeVar, "Node_set_vec" );
    auto list_node_mag   = std::make_shared<IO::Variable>( 1, NodeVar, "Node_list_mag" );
    auto list_node_vec   = std::make_shared<IO::Variable>( 3, NodeVar, "Node_list_vec" );
    auto point_node_mag  = std::make_shared<IO::Variable>( 1, NodeVar, "Node_point_mag" );
    auto point_node_vec  = std::make_shared<IO::Variable>( 3, NodeVar, "Node_point_vec" );
    auto domain_node_mag = std::make_shared<IO::Variable>( 1, NodeVar, "Node_domain_mag" );
    auto domain_node_vec = std::make_shared<IO::Variable>( 3, NodeVar, "Node_domain_vec" );
    auto set_cell_mag    = std::make_shared<IO::Variable>( 1, VolVar, "Cell_set_mag" );
    auto set_cell_vec    = std::make_shared<IO::Variable>( 3, VolVar, "Cell_set_vec" );
    auto list_cell_mag   = std::make_shared<IO::Variable>( 1, VolVar, "Cell_list_mag" );
    auto list_cell_vec   = std::make_shared<IO::Variable>( 3, VolVar, "Cell_list_vec" );
    auto domain_cell_mag = std::make_shared<IO::Variable>( 1, VolVar, "Cell_domain_mag" );
    auto domain_cell_vec = std::make_shared<IO::Variable>( 3, VolVar, "Cell_domain_vec" );
    point_node_mag->data.resize( N_points );
    point_node_vec->data.resize( N_points, 3 );
    for ( int i = 0; i < N_points; i++ ) {
        point_node_mag->data( i )    = distance( set1->points[i] );
        point_node_vec->data( i, 0 ) = set1->points[i].x;
        point_node_vec->data( i, 1 ) = set1->points[i].y;
        point_node_vec->data( i, 2 ) = set1->points[i].z;
    }
    set_node_mag->data = point_node_mag->data;
    set_node_vec->data = point_node_vec->data;
    list_node_mag->data.resize( 3 * N_tri );
    list_node_vec->data.resize( 3 * N_tri, 3 );
    for ( int i = 0; i < N_tri; i++ ) {
        list_node_mag->data( 3 * i + 0 )    = distance( trilist->A[i] );
        list_node_mag->data( 3 * i + 1 )    = distance( trilist->B[i] );
        list_node_mag->data( 3 * i + 2 )    = distance( trilist->C[i] );
        list_node_vec->data( 3 * i + 0, 0 ) = trilist->A[i].x;
        list_node_vec->data( 3 * i + 0, 1 ) = trilist->A[i].y;
        list_node_vec->data( 3 * i + 0, 2 ) = trilist->A[i].z;
        list_node_vec->data( 3 * i + 1, 0 ) = trilist->B[i].x;
        list_node_vec->data( 3 * i + 1, 1 ) = trilist->B[i].y;
        list_node_vec->data( 3 * i + 1, 2 ) = trilist->B[i].z;
        list_node_vec->data( 3 * i + 2, 0 ) = trilist->C[i].x;
        list_node_vec->data( 3 * i + 2, 1 ) = trilist->C[i].y;
        list_node_vec->data( 3 * i + 2, 2 ) = trilist->C[i].z;
    }
    domain_node_mag->data.resize( domain->nx + 1, domain->ny + 1, domain->nz + 1 );
    domain_node_vec->data.resize(
        { (size_t) domain->nx + 1, (size_t) domain->ny + 1, (size_t) domain->nz + 1, 3 } );
    for ( int i = 0; i < domain->nx + 1; i++ ) {
        for ( int j = 0; j < domain->ny + 1; j++ ) {
            for ( int k = 0; k < domain->nz + 1; k++ ) {
                domain_node_mag->data( i, j, k )    = distance( Point( i, j, k ) );
                domain_node_vec->data( i, j, k, 0 ) = Point( i, j, k ).x;
                domain_node_vec->data( i, j, k, 1 ) = Point( i, j, k ).y;
                domain_node_vec->data( i, j, k, 2 ) = Point( i, j, k ).z;
            }
        }
    }
    set_cell_mag->data.resize( N_tri );
    set_cell_vec->data.resize( N_tri, 3 );
    for ( int i = 0; i < N_tri; i++ ) {
        set_cell_mag->data( i )    = i;
        set_cell_vec->data( i, 0 ) = 3 * i + 0;
        set_cell_vec->data( i, 1 ) = 3 * i + 1;
        set_cell_vec->data( i, 2 ) = 3 * i + 2;
    }
    list_cell_mag->data = set_cell_mag->data;
    list_cell_vec->data = set_cell_vec->data;
    domain_cell_mag->data.resize( domain->nx, domain->ny, domain->nz );
    domain_cell_vec->data.resize(
        { (size_t) domain->nx, (size_t) domain->ny, (size_t) domain->nz, 3 } );
    for ( int i = 0; i < domain->nx; i++ ) {
        for ( int j = 0; j < domain->ny; j++ ) {
            for ( int k = 0; k < domain->nz; k++ ) {
                domain_cell_mag->data( i, j, k )    = distance( Point( i, j, k ) );
                domain_cell_vec->data( i, j, k, 0 ) = Point( i, j, k ).x;
                domain_cell_vec->data( i, j, k, 1 ) = Point( i, j, k ).y;
                domain_cell_vec->data( i, j, k, 2 ) = Point( i, j, k ).z;
            }
        }
    }

    // Create the MeshDataStruct
    std::vector<IO::MeshDataStruct> meshData( 4 );
    meshData[0].meshName = "pointmesh";
    meshData[0].mesh     = set1;
    meshData[0].vars.push_back( point_node_mag );
    meshData[0].vars.push_back( point_node_vec );
    meshData[1].meshName = "trimesh";
    meshData[1].mesh     = trimesh;
    meshData[1].vars.push_back( set_node_mag );
    meshData[1].vars.push_back( set_node_vec );
    meshData[1].vars.push_back( set_cell_mag );
    meshData[1].vars.push_back( set_cell_vec );
    meshData[2].meshName = "trilist";
    meshData[2].mesh     = trilist;
    meshData[2].vars.push_back( list_node_mag );
    meshData[2].vars.push_back( list_node_vec );
    meshData[2].vars.push_back( list_cell_mag );
    meshData[2].vars.push_back( list_cell_vec );
    meshData[3].meshName = "domain";
    meshData[3].mesh     = domain;
    meshData[3].vars.push_back( domain_node_mag );
    meshData[3].vars.push_back( domain_node_vec );
    meshData[3].vars.push_back( domain_cell_mag );
    meshData[3].vars.push_back( domain_cell_vec );
    for ( const auto &data : meshData )
        ASSERT( data.check( true ) );

    // Run the tests
    testWriter( "old", meshData, ut );
    testWriter( "new", meshData, ut );
#ifdef USE_SILO
    testWriter( "silo-double", meshData, ut );
    testWriter( "silo-float", meshData, ut );
#endif
#ifdef USE_HDF5
    testWriter( "hdf5-double", meshData, ut );
    testWriter( "hdf5-float", meshData, ut );
#endif

    // Finished
    ut.report();
    PROFILE_SAVE( "TestWriter", true );
    int N_errors = ut.NumFailGlobal();
    comm.barrier();
    Utilities::shutdown();
    return N_errors;
}
