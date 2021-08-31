#include "IO/Xdmf.h"

#include "common/Array.h"
#include "common/UtilityMacros.h"


ArraySize squeeze( const ArraySize &x )
{
    int Nd      = 0;
    size_t N[5] = { 1 };
    for ( size_t i = 0; i < x.maxDim(); i++ ) {
        if ( x[i] != 1 )
            N[Nd++] = x[i];
    }
    return ArraySize( std::max( Nd, 1 ), N );
}


// Helper functions
static void addDataItem(
    FILE *xmf, const std::string &indent, ArraySize size, const std::string &location )
{
    size = squeeze( size );
    if ( size.ndim() == 1 ) {
        fprintf( xmf, "%s<DataItem Dimensions=\"%lu\"", indent.data(), size[0] );
    } else if ( size.ndim() == 2 ) {
        fprintf( xmf, "%s<DataItem Dimensions=\"%lu %lu\"", indent.data(), size[1], size[0] );
    } else if ( size.ndim() == 3 ) {
        fprintf( xmf, "%s<DataItem Dimensions=\"%lu %lu %lu\"", indent.data(), size[2], size[1],
            size[0] );
    } else if ( size.ndim() == 4 ) {
        fprintf( xmf, "%s<DataItem Dimensions=\"%lu %lu %lu %lu\"", indent.data(), size[3], size[2],
            size[1], size[0] );
    } else {
        ERROR( "Invalid number of dimensions" );
    }
    fprintf( xmf, " Format=\"HDF\">\n" );
    fprintf( xmf, "%s  %s\n", indent.data(), location.data() );
    fprintf( xmf, "%s</DataItem>\n", indent.data() );
}
template<class TYPE>
static void addVariable( FILE *xmf, const std::string &indent, const std::string &name,
    const std::string &type, const std::string &center, ArraySize size,
    const std::string &location )
{
    fprintf( xmf, "%s<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"%s\">\n", indent.data(),
        name.data(), type.data(), center.data() );
    addDataItem( xmf, indent + "  ", size, location );
    fprintf( xmf, "%s</Attribute>\n", indent.data() );
}


/****************************************************************
 * Enum functions                                                *
 ****************************************************************/
/*template<class TYPE>
static Xdmf::DataType getDataType()
{
    if ( std::is_same<TYPE, char>::value )
        return Xdmf::DataType::Char;
    else if ( std::is_same<TYPE, int32_t>::value )
        return Xdmf::DataType::Int32;
    else if ( std::is_same<TYPE, int64_t>::value )
        return Xdmf::DataType::Int64;
    else if ( std::is_same<TYPE, uint32_t>::value )
        return Xdmf::DataType::Uint32;
    else if ( std::is_same<TYPE, uint64_t>::value )
        return Xdmf::DataType::Uint64;
    else if ( std::is_same<TYPE, float>::value )
        return Xdmf::DataType::Float;
    else if ( std::is_same<TYPE, double>::value )
        return Xdmf::DataType::Double;
    else
        ERROR( "Invalid type" );
}*/
static const char *TopologyTypeNames[]  = { "", "Polyvertex", "Polyline", "Polygon", "Triangle",
    "Quadrilateral", "Tetrahedron", "Pyramid", "Wedge", "Hexahedron", "Edge_3", "Triangle_6",
    "Quadrilateral_8", "Tetrahedron_10", "Pyramid_13", "Wedge_15", "Hexahedron_20", "Mixed",
    "CurvilinearMesh2D", "CurvilinearMesh3D", "RectangularMesh2D", "RectangularMesh3D",
    "UniformMesh2D", "UniformMesh3D" };
static const uint8_t TopologyTypeDOFs[] = { 0, 1, 2, 0, 3, 4, 4, 5, 6, 8, 3, 6, 8, 10, 13, 15, 20,
    0, 0, 0, 0, 0, 0, 0 };


/****************************************************************
 * Create a mesh                                                 *
 ****************************************************************/
Xdmf::MeshData Xdmf::createPointMesh( const std::string &name, uint8_t NDIM, size_t N,
    const std::string &x, const std::string &y, const std::string &z )
{
    return createUnstructuredMesh( name, NDIM, TopologyType::Polyvertex, N, "", N, x, y, z );
}
Xdmf::MeshData Xdmf::createUniformMesh(
    const std::string &name, const std::vector<double> &range, ArraySize size )
{
    ASSERT( range.size() == 2 * size.ndim() );
    MeshData data;
    data.name = name;
    data.size = size;
    if ( size.ndim() == 2 )
        data.type = TopologyType::UniformMesh2D;
    else if ( size.ndim() == 3 )
        data.type = TopologyType::UniformMesh3D;
    else
        ERROR( "# of dimensions != 2 or 3" );
    for ( int i = 0; i < 2 * size.ndim(); i++ )
        data.range[i] = range[i];
    return data;
}
Xdmf::MeshData Xdmf::createCurvilinearMesh( const std::string &name, ArraySize size,
    const std::string &x, const std::string &y, const std::string &z )
{
    MeshData data;
    data.name = name;
    if ( size.ndim() == 2 )
        data.type = TopologyType::CurvilinearMesh2D;
    else if ( size.ndim() == 3 )
        data.type = TopologyType::CurvilinearMesh3D;
    else
        ERROR( "Invalid size for Curvilinear mesh" );
    data.size = size;
    data.x    = x;
    data.y    = y;
    data.z    = z;
    return data;
}
Xdmf::MeshData Xdmf::createUnstructuredMesh( const std::string &name, uint8_t NDIM,
    TopologyType type, size_t NumElements, const std::string &dofMap, size_t NumNodes,
    const std::string &x, const std::string &y, const std::string &z )
{
    ASSERT( type != TopologyType::Null );
    MeshData data;
    data.name   = name;
    data.type   = type;
    data.size   = { NDIM, NumElements, NumNodes };
    data.dofMap = dofMap;
    data.x      = x;
    data.y      = y;
    data.z      = z;
    return data;
}


/****************************************************************
 * Add a variable                                                *
 ****************************************************************/
void Xdmf::MeshData::addVariable( const std::string &meshName, const std::string &varName,
    ArraySize varSize, RankType rank, Center center, const std::string &varData )
{
    VarData var;
    var.name     = varName;
    var.size     = varSize;
    var.data     = varData;
    var.rankType = rank;
    var.center   = center;
    vars.push_back( std::move( var ) );
}


/****************************************************************
 * Add a mesh domain                                             *
 ****************************************************************/
void Xdmf::addMesh( const std::string &meshName, const MeshData &domain )
{
    auto &domains = d_meshData[meshName];
    for ( const auto &domain2 : domains )
        ASSERT( domain2.name != domain.name );
    domains.push_back( domain );
}


/****************************************************************
 * Write a variable                                              *
 ****************************************************************/
static void writeVariable( FILE *fid, const Xdmf::VarData &var, const std::string &indent )
{
    // Write the variable name
    fprintf( fid, "%s<Attribute Name=\"%s\"", indent.data(), var.name.data() );
    // Write the variable type
    if ( var.rankType == Xdmf::RankType::Scalar ) {
        fprintf( fid, " AttributeType=\"Scalar\"" );
    } else if ( var.rankType == Xdmf::RankType::Vector ) {
        fprintf( fid, " AttributeType=\"Vector\"" );
    } else if ( var.rankType == Xdmf::RankType::Tensor ) {
        fprintf( fid, " AttributeType=\"Tensor\"" );
    } else if ( var.rankType == Xdmf::RankType::Tensor6 ) {
        fprintf( fid, " AttributeType=\"Tensor6\"" );
    } else if ( var.rankType == Xdmf::RankType::Matrix ) {
        fprintf( fid, " AttributeType=\"Matrix\"" );
    } else if ( var.rankType == Xdmf::RankType::GlobalID ) {
        fprintf( fid, " AttributeType=\"GlobalID\"" );
    } else {
        ERROR( "Unknown center type" );
    }
    // Write the variable centering
    if ( var.center == Xdmf::Center::Node ) {
        fprintf( fid, " Center=\"Node\">\n" );
    } else if ( var.center == Xdmf::Center::Cell ) {
        fprintf( fid, " Center=\"Cell\">\n" );
    } else if ( var.center == Xdmf::Center::Grid ) {
        fprintf( fid, " Center=\"Grid\">\n" );
    } else if ( var.center == Xdmf::Center::Face ) {
        fprintf( fid, " Center=\"Face\">\n" );
    } else if ( var.center == Xdmf::Center::Edge ) {
        fprintf( fid, " Center=\"Edge\">\n" );
    } else if ( var.center == Xdmf::Center::Other ) {
        fprintf( fid, " Center=\"Other\">\n" );
    } else {
        ERROR( "Unknown center type" );
    }
    // Write the data item
    addDataItem( fid, indent + "  ", var.size, var.data );
    // Finished
    fprintf( fid, "%s</Attribute>\n", indent.data() );
}


/****************************************************************
 * Write the mesh grid                                           *
 ****************************************************************/
static void writeMeshGrid( FILE *fid, const Xdmf::MeshData &mesh, const std::string &indent )
{
    const char *s = indent.data();
    double x0[3]  = { mesh.range[0], mesh.range[2], mesh.range[4] };
    double dx[3]  = { ( mesh.range[1] - mesh.range[0] ) / mesh.size[0],
        ( mesh.range[3] - mesh.range[2] ) / mesh.size[1],
        ( mesh.range[5] - mesh.range[4] ) / mesh.size[2] };
    switch ( mesh.type ) {
    case Xdmf::TopologyType::UniformMesh2D:
        // Write a uniform 2d mesh
        fprintf( fid, "%s<Grid Name=\"%s\" GridType=\"Uniform\">\n", s, mesh.name.data() );
        fprintf( fid,
            "%s  <Topology TopologyType=\"2DCoRectMesh\" NumberOfElements=\"%lu %lu\"/>\n", s,
            mesh.size[1] + 1, mesh.size[0] + 1 );
        fprintf( fid, "%s  <Geometry GeometryType=\"ORIGIN_DXDY\">\n", s );
        fprintf(
            fid, "%s    <DataItem  Format=\"XML\" NumberType=\"float\" Dimensions=\"2\">\n", s );
        fprintf( fid, "%s      %0.12e  %0.12e\n", s, x0[0], x0[1] );
        fprintf( fid, "%s    </DataItem>\n", s );
        fprintf(
            fid, "%s    <DataItem  Format=\"XML\" NumberType=\"float\" Dimensions=\"2\">\n", s );
        fprintf( fid, "%s       %0.12e  %0.12e\n", s, dx[0], dx[1] );
        fprintf( fid, "%s    </DataItem>\n", s );
        fprintf( fid, "%s  </Geometry>\n", s );
        break;
    case Xdmf::TopologyType::UniformMesh3D:
        // Write a uniform 3d mesh
        fprintf( fid, "%s<Grid Name=\"%s\" GridType=\"Uniform\">\n", s, mesh.name.data() );
        fprintf( fid,
            "%s  <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%lu %lu\"/>\n", s,
            mesh.size[1] + 1, mesh.size[0] + 1 );
        fprintf( fid, "%s  <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n", s );
        fprintf(
            fid, "%s    <DataItem  Format=\"XML\" NumberType=\"float\" Dimensions=\"3\">\n", s );
        fprintf( fid, "%s      %0.12e  %0.12e  %0.12e\n", s, x0[0], x0[1], x0[2] );
        fprintf( fid, "%s    </DataItem>\n", s );
        fprintf(
            fid, "%s    <DataItem  Format=\"XML\" NumberType=\"float\" Dimensions=\"3\">\n", s );
        fprintf( fid, "%s       %0.12e  %0.12e  %0.12e\n", s, dx[0], dx[1], dx[2] );
        fprintf( fid, "%s    </DataItem>\n", s );
        fprintf( fid, "%s  </Geometry>\n", s );
        break;
    case Xdmf::TopologyType::CurvilinearMesh2D:
        // Write a 2D curvillinear mesh
        fprintf( fid, "%s<Grid Name=\"%s\" GridType=\"Uniform\">\n", s, mesh.name.data() );
        fprintf( fid, "%s  <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%lu %lu\"/>\n", s,
            mesh.size[1] + 1, mesh.size[0] + 1 );
        fprintf( fid, "%s  <Geometry GeometryType=\"X_Y\">\n", s );
        addDataItem( fid, indent + "    ", mesh.size + 1, mesh.x );
        addDataItem( fid, indent + "    ", mesh.size + 1, mesh.y );
        fprintf( fid, "%s  </Geometry>\n", s );
        break;
    case Xdmf::TopologyType::CurvilinearMesh3D:
        // Write a 3D curvillinear mesh
        fprintf( fid, "%s<Grid Name=\"%s\" GridType=\"Uniform\">\n", s, mesh.name.data() );
        fprintf( fid, "%s  <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%lu %lu %lu\"/>\n",
            s, mesh.size[2] + 1, mesh.size[1] + 1, mesh.size[0] + 1 );
        fprintf( fid, "%s  <Geometry GeometryType=\"X_Y_Z\">\n", s );
        addDataItem( fid, indent + "    ", mesh.size + 1, mesh.x );
        addDataItem( fid, indent + "    ", mesh.size + 1, mesh.y );
        addDataItem( fid, indent + "    ", mesh.size + 1, mesh.z );
        fprintf( fid, "%s  </Geometry>\n", s );
        break;
    case Xdmf::TopologyType::Polyvertex:
    case Xdmf::TopologyType::Polyline:
    case Xdmf::TopologyType::Polygon:
    case Xdmf::TopologyType::Triangle:
    case Xdmf::TopologyType::Quadrilateral:
    case Xdmf::TopologyType::Tetrahedron:
    case Xdmf::TopologyType::Pyramid:
    case Xdmf::TopologyType::Wedge:
    case Xdmf::TopologyType::Hexahedron:
    case Xdmf::TopologyType::Edge_3:
    case Xdmf::TopologyType::Triangle_6:
    case Xdmf::TopologyType::Quadrilateral_8:
    case Xdmf::TopologyType::Tetrahedron_10:
    case Xdmf::TopologyType::Pyramid_13:
    case Xdmf::TopologyType::Wedge_15:
    case Xdmf::TopologyType::Hexahedron_20:
        // Write an unstructured mesh
        {
            int NDIM      = mesh.size[0];
            size_t Nelem  = mesh.size[1];
            size_t Nnode  = mesh.size[2];
            uint8_t Ndofs = TopologyTypeDOFs[static_cast<int>( mesh.type )];
            auto type     = TopologyTypeNames[static_cast<int>( mesh.type )];
            fprintf( fid, "%s<Grid Name=\"%s\">\n", s, mesh.name.data() );
            fprintf( fid, "%s  <Topology TopologyType=\"%s\"", s, type );
            fprintf( fid, " NumberOfElements=\"%lu\">\n", Nelem );
            if ( !mesh.dofMap.empty() )
                addDataItem( fid, indent + "    ", { Ndofs, Nelem }, mesh.dofMap );
            fprintf( fid, "%s  </Topology>\n", s );
            if ( NDIM == 2 ) {
                if ( mesh.y.empty() ) {
                    fprintf( fid, "%s  <Geometry GeometryType=\"XY\">\n", s );
                    addDataItem( fid, indent + "    ", { 2, Nnode }, mesh.x );
                } else {
                    fprintf( fid, "%s  <Geometry GeometryType=\"X_Y\">\n", s );
                    addDataItem( fid, indent + "    ", Nnode, mesh.x );
                    addDataItem( fid, indent + "    ", Nnode, mesh.y );
                }
            } else if ( NDIM == 3 ) {
                if ( mesh.y.empty() ) {
                    fprintf( fid, "%s  <Geometry GeometryType=\"XYZ\">\n", s );
                    addDataItem( fid, indent + "    ", { 2, Nnode }, mesh.x );
                } else {
                    fprintf( fid, "%s  <Geometry GeometryType=\"X_Y_Z\">\n", s );
                    addDataItem( fid, indent + "    ", Nnode, mesh.x );
                    addDataItem( fid, indent + "    ", Nnode, mesh.y );
                    addDataItem( fid, indent + "    ", Nnode, mesh.z );
                }
            } else {
                ERROR( "Dimensions other than 2 or 3 are not supported" );
            }
            fprintf( fid, "%s  </Geometry>\n", s );
        }
        break;
    default: {
        auto msg = "Invalid mesh type: " + std::to_string( static_cast<int>( mesh.type ) ) + " - " +
                   mesh.name;
        ERROR( msg );
    }
    }
    // Write the variables
    for ( const auto &var : mesh.vars )
        writeVariable( fid, var, indent + "  " );
    fprintf( fid, "%s  </Grid>\n", s );
}


/****************************************************************
 * Write the XDMF xml file                                       *
 ****************************************************************/
void Xdmf::write( const std::string &filename ) const
{
    // Create XDMF file
    auto fid = fopen( filename.data(), "w" );
    fprintf( fid, "<?xml version=\"1.0\" ?>\n" );
    fprintf( fid, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" );
    fprintf( fid, "<Xdmf Version=\"2.0\">\n" );
    fprintf( fid, "<Domain>\n" );
    // Write an empty mesh to enable collections to work properly
    fprintf( fid, "  <Grid Name=\"\" GridType=\"Uniform\"></Grid>\n\n" );
    // Write each mesh
    for ( const auto &data : d_meshData ) {
        auto name    = data.first;
        auto domains = data.second;
        if ( domains.empty() )
            continue;
        if ( domains.size() == 1u && name == domains[0].name ) {
            writeMeshGrid( fid, domains[0], "  " );
        } else {
            fprintf( fid, "  <Grid Name=\"%s\" GridType=\"Collection\">\n", name.data() );
            for ( const auto &domain : domains )
                writeMeshGrid( fid, domain, "    " );
            fprintf( fid, "  </Grid>\n\n" );
        }
    }
    fprintf( fid, "</Domain>\n" );
    fprintf( fid, "</Xdmf>\n" );
    fclose( fid );
}


/****************************************************************
 * Pack/Unpack data                                              *
 ****************************************************************/
template<class T>
typename std::enable_if<std::is_trivially_copyable<T>::value, size_t>::type size( const T & )
{
    return sizeof( T );
}
template<class T>
typename std::enable_if<std::is_trivially_copyable<T>::value, char *>::type pack(
    char *ptr, const T &x )
{
    memcpy( ptr, &x, sizeof( T ) );
    return ptr + sizeof( T );
}
template<class T>
typename std::enable_if<std::is_trivially_copyable<T>::value, char *>::type unpack(
    char *ptr, T &x )
{
    memcpy( &x, ptr, sizeof( T ) );
    return ptr + sizeof( T );
}
static size_t size( const std::string &str ) { return sizeof( int ) + str.size(); }
static char *pack( char *ptr, const std::string &str )
{
    int N = str.size();
    memcpy( ptr, &N, sizeof( int ) );
    ptr += sizeof( int );
    memcpy( ptr, str.data(), str.size() );
    ptr += str.size();
    return ptr;
}
static char *unpack( char *ptr, std::string &str )
{
    int N = 0;
    memcpy( &N, ptr, sizeof( int ) );
    ASSERT( N >= 0 && N < 1000 );
    ptr += sizeof( int );
    str = std::string( ptr, N );
    ptr += N;
    return ptr;
}
static size_t size( const Xdmf::VarData &data )
{
    size_t bytes = 0;
    bytes += size( data.name );
    bytes += size( data.size );
    bytes += size( data.rankType );
    bytes += size( data.center );
    bytes += size( data.data );
    return bytes;
}
static char *pack( char *ptr, const Xdmf::VarData &data )
{
    ptr = pack( ptr, data.name );
    ptr = pack( ptr, data.size );
    ptr = pack( ptr, data.rankType );
    ptr = pack( ptr, data.center );
    ptr = pack( ptr, data.data );
    return ptr;
}
static char *unpack( char *ptr, Xdmf::VarData &data )
{
    int rankType = 0, center = 0;
    ptr           = unpack( ptr, data.name );
    ptr           = unpack( ptr, data.size );
    ptr           = unpack( ptr, rankType );
    ptr           = unpack( ptr, center );
    ptr           = unpack( ptr, data.data );
    data.rankType = static_cast<Xdmf::RankType>( rankType );
    data.center   = static_cast<Xdmf::Center>( center );
    return ptr;
}
static size_t size( const Xdmf::MeshData &data )
{
    int N_vars   = data.vars.size();
    size_t bytes = 0;
    bytes += size( data.name );
    bytes += size( data.type );
    bytes += size( data.size );
    bytes += size( data.range );
    bytes += size( data.x );
    bytes += size( data.y );
    bytes += size( data.z );
    bytes += size( N_vars );
    for ( int i = 0; i < N_vars; i++ )
        bytes += size( data.vars[i] );
    return bytes;
}
static char *pack( char *ptr, const Xdmf::MeshData &data )
{
    int N_vars = data.vars.size();
    ptr        = pack( ptr, data.name );
    ptr        = pack( ptr, data.type );
    ptr        = pack( ptr, data.size );
    ptr        = pack( ptr, data.range );
    ptr        = pack( ptr, data.x );
    ptr        = pack( ptr, data.y );
    ptr        = pack( ptr, data.z );
    ptr        = pack( ptr, N_vars );
    for ( int i = 0; i < N_vars; i++ )
        ptr = pack( ptr, data.vars[i] );
    return ptr;
}
static char *unpack( char *ptr, Xdmf::MeshData &data )
{
    int N_vars = 0;
    ptr        = unpack( ptr, data.name );
    ptr        = unpack( ptr, data.type );
    ptr        = unpack( ptr, data.size );
    ptr        = unpack( ptr, data.range );
    ptr        = unpack( ptr, data.x );
    ptr        = unpack( ptr, data.y );
    ptr        = unpack( ptr, data.z );
    ptr        = unpack( ptr, N_vars );
    data.vars.resize( N_vars );
    for ( int i = 0; i < N_vars; i++ )
        ptr = unpack( ptr, data.vars[i] );
    return ptr;
}
static size_t size( const std::vector<Xdmf::MeshData> &data )
{
    size_t bytes = 0;
    int N        = data.size();
    bytes += size( N );
    for ( int i = 0; i < N; i++ )
        bytes += size( data[i] );
    return bytes;
}
static char *pack( char *ptr, const std::vector<Xdmf::MeshData> &data )
{
    int N = data.size();
    ptr   = pack( ptr, N );
    for ( int i = 0; i < N; i++ )
        ptr = pack( ptr, data[i] );
    return ptr;
}
static char *unpack( char *ptr, std::vector<Xdmf::MeshData> &data )
{
    data.clear();
    int N = data.size();
    ptr   = unpack( ptr, N );
    data.resize( N );
    for ( int i = 0; i < N; i++ )
        ptr = unpack( ptr, data[i] );
    return ptr;
}
static size_t size( const std::map<std::string, std::vector<Xdmf::MeshData>> &data )
{
    size_t bytes = 0;
    int N_map    = data.size();
    bytes += size( N_map );
    for ( const auto &tmp : data ) {
        bytes += size( tmp.first );
        bytes += size( tmp.second );
    }
    return bytes;
}
static char *pack( char *ptr, const std::map<std::string, std::vector<Xdmf::MeshData>> &data )
{
    int N_map = data.size();
    ptr       = pack( ptr, N_map );
    for ( const auto &tmp : data ) {
        ptr = pack( ptr, tmp.first );
        ptr = pack( ptr, tmp.second );
    }
    return ptr;
}
static char *unpack( char *ptr, std::map<std::string, std::vector<Xdmf::MeshData>> &data )
{
    data.clear();
    int N_map = data.size();
    ptr       = unpack( ptr, N_map );
    for ( int i = 0; i < N_map; i++ ) {
        std::string name;
        std::vector<Xdmf::MeshData> data2;
        ptr        = unpack( ptr, name );
        ptr        = unpack( ptr, data2 );
        data[name] = std::move( data2 );
    }
    return ptr;
}


/****************************************************************
 * Gather all data to rank 0                                     *
 ****************************************************************/
void Xdmf::gather( const Utilities::MPI &comm )
{
    if ( comm.getRank() == 0 ) {
        for ( int i = 1; i < comm.getSize(); i++ ) {
            // Recieve the data
            size_t N_meshes = 0, N_bytes = 0;
            comm.recv( &N_meshes, 1, i, 717 );
            comm.recv( &N_bytes, 1, i, 718 );
            auto buf = new char[N_bytes];
            comm.recv( buf, N_bytes, i, 719 );
            // Unpack the data
            std::map<std::string, std::vector<MeshData>> data;
            unpack( buf, data );
            delete[] buf;
            // Add the meshes
            for ( auto tmp : data ) {
                const auto &name    = tmp.first;
                const auto &domains = tmp.second;
                if ( domains.size() == 1u && domains[0].name == name ) {
                    // We are dealing with a single mesh
                    ASSERT( d_meshData.find( name ) == d_meshData.end() );
                    d_meshData.insert( tmp );
                } else {
                    // Add the domains
                    auto &meshes = d_meshData[name];
                    for ( auto domain : domains ) {
                        for ( const auto &tmp : meshes )
                            ASSERT( tmp.name != domain.name );
                        meshes.push_back( domain );
                    }
                }
            }
        }
    } else {
        // Send the number of meshes
        size_t N_meshes = d_meshData.size();
        comm.send( &N_meshes, 1, 0, 717 );
        // Pack the send data
        size_t N_bytes = size( d_meshData );
        comm.send( &N_bytes, 1, 0, 718 );
        auto buf = new char[N_bytes];
        pack( buf, d_meshData );
        // Send the data to rank 0
        comm.send( buf, N_bytes, 0, 719 );
        delete[] buf;
        // Clear the internal data
        d_meshData.clear();
    }
}
