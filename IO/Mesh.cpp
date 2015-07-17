#include "Mesh.h"
#include "common/Utilities.h"
#include "shared_ptr.h"

#include <limits>
#include <stdint.h>

namespace IO {


inline Point nullPoint()
{
    Point tmp;
    tmp.x = std::numeric_limits<double>::quiet_NaN();
    tmp.y = std::numeric_limits<double>::quiet_NaN();
    tmp.z = std::numeric_limits<double>::quiet_NaN();
    return tmp;
}


/****************************************************
* Mesh                                              *
****************************************************/
Mesh::Mesh( )
{
}
Mesh::~Mesh( )
{
}


/****************************************************
* PointList                                         *
****************************************************/
PointList::PointList( )
{
}
PointList::PointList( size_t N )
{
    Point tmp = nullPoint();
    points.resize(N,tmp);
}
PointList::~PointList( )
{
}
size_t PointList::numberPointsVar( VariableType type ) const
{
    size_t N = 0;
    if ( type == NodeVariable )
        N = points.size();
    return N;
}
std::pair<size_t,void*> PointList::pack( int level ) const
{
    std::pair<size_t,void*> data_out(0,NULL);
    if ( level==0 ) {
        data_out.first = (2+3*points.size())*sizeof(double);
        double *data_ptr = new double[2+3*points.size()];
        data_out.second = data_ptr;
        uint64_t *data_int = reinterpret_cast<uint64_t*>(data_ptr);
        data_int[0] = level;
        data_int[1] = points.size();
        double *data = &data_ptr[2];
        for (size_t i=0; i<points.size(); i++) {
            data[3*i+0] = points[i].x;
            data[3*i+1] = points[i].y;
            data[3*i+2] = points[i].z;
        }
    }
    return data_out;
}
void PointList::unpack( const std::pair<size_t,void*>& data_in )
{
    uint64_t *data_int = reinterpret_cast<uint64_t*>(data_in.second);
    const double *data = reinterpret_cast<const double*>(data_in.second);
    int level = data_int[0];
    uint64_t N = data_int[1];
    data = &data[2];
    if ( level==0 ) {
        ASSERT((2+3*N)*sizeof(double)==data_in.first);
        points.resize(N);
        for (size_t i=0; i<points.size(); i++) {
            points[i].x = data[3*i+0];
            points[i].y = data[3*i+1];
            points[i].z = data[3*i+2];
        }
    }
}


/****************************************************
* TriList                                           *
****************************************************/
TriList::TriList( )
{
}
TriList::TriList( size_t N_tri )
{
    Point tmp = nullPoint();
    A.resize(N_tri,tmp);
    B.resize(N_tri,tmp);
    C.resize(N_tri,tmp);
}
TriList::TriList( const TriMesh& mesh )
{
    Point tmp = nullPoint();
    A.resize(mesh.A.size(),tmp);
    B.resize(mesh.B.size(),tmp);
    C.resize(mesh.C.size(),tmp);
    ASSERT(mesh.vertices.get()!=NULL);
    const std::vector<Point>& P = mesh.vertices->points;
    for (size_t i=0; i<A.size(); i++)
        A[i] = P[mesh.A[i]];
    for (size_t i=0; i<B.size(); i++)
        B[i] = P[mesh.B[i]];
    for (size_t i=0; i<C.size(); i++)
        C[i] = P[mesh.C[i]];
}
TriList::~TriList( )
{
}
size_t TriList::numberPointsVar( VariableType type ) const
{
    size_t N = 0;
    if ( type==NodeVariable )
        N = 3*A.size();
    else if ( type==SurfaceVariable ||  type==VolumeVariable )
        N = A.size();
    return N;
}
std::pair<size_t,void*> TriList::pack( int level ) const
{
    std::pair<size_t,void*> data_out(0,NULL);
    if ( level==0 ) {
        data_out.first = (2+9*A.size())*sizeof(double);
        double *data_ptr = new double[2+9*A.size()];
        data_out.second = data_ptr;
        uint64_t *data_int = reinterpret_cast<uint64_t*>(data_ptr);
        data_int[0] = level;
        data_int[1] = A.size();
        double *data = &data_ptr[2];
        for (size_t i=0; i<A.size(); i++) {
            data[9*i+0] = A[i].x;
            data[9*i+1] = A[i].y;
            data[9*i+2] = A[i].z;
            data[9*i+3] = B[i].x;
            data[9*i+4] = B[i].y;
            data[9*i+5] = B[i].z;
            data[9*i+6] = C[i].x;
            data[9*i+7] = C[i].y;
            data[9*i+8] = C[i].z;
        }
    }
    return data_out;
}
void TriList::unpack( const std::pair<size_t,void*>& data_in )
{
    uint64_t *data_int = reinterpret_cast<uint64_t*>(data_in.second);
    const double *data = reinterpret_cast<const double*>(data_in.second);
    int level = data_int[0];
    uint64_t N = data_int[1];
    data = &data[2];
    if ( level==0 ) {
        ASSERT((2+9*N)*sizeof(double)==data_in.first);
        A.resize(N);
        B.resize(N);
        C.resize(N);
        for (size_t i=0; i<A.size(); i++) {
            A[i].x = data[9*i+0];
            A[i].y = data[9*i+1];
            A[i].z = data[9*i+2];
            B[i].x = data[9*i+3];
            B[i].y = data[9*i+4];
            B[i].z = data[9*i+5];
            C[i].x = data[9*i+6];
            C[i].y = data[9*i+7];
            C[i].z = data[9*i+8];
        }
    }
}


/****************************************************
* TriMesh                                           *
****************************************************/
TriMesh::TriMesh(  )
{
}
TriMesh::TriMesh( size_t N_tri, size_t N_point )
{
    vertices.reset( new PointList(N_point) );
    A.resize(N_tri,-1);
    B.resize(N_tri,-1);
    C.resize(N_tri,-1);
}
TriMesh::TriMesh( size_t N_tri, std::shared_ptr<PointList> points )
{
    vertices = points;
    A.resize(N_tri,-1);
    B.resize(N_tri,-1);
    C.resize(N_tri,-1);
}
TriMesh::TriMesh( const TriList& mesh )
{
    // For simlicity we will just create a mesh with ~3x the verticies for now
    ASSERT(mesh.A.size()==mesh.B.size()&&mesh.A.size()==mesh.C.size());
    A.resize(mesh.A.size());
    B.resize(mesh.B.size());
    C.resize(mesh.C.size());
    vertices.reset( new PointList(3*mesh.A.size()) );
    for (size_t i=0; i<A.size(); i++) {
        A[i] = 3*i+0;
        B[i] = 3*i+1;
        C[i] = 3*i+2;
        vertices->points[A[i]] = mesh.A[i];
        vertices->points[B[i]] = mesh.B[i];
        vertices->points[C[i]] = mesh.C[i];
    }
}
TriMesh::~TriMesh( )
{
    vertices.reset();
    A.clear();
    B.clear();
    C.clear();
}
size_t TriMesh::numberPointsVar( VariableType type ) const
{
    size_t N = 0;
    if ( type==NodeVariable )
        N = vertices->points.size();
    else if ( type==SurfaceVariable ||  type==VolumeVariable )
        N = A.size();
    return N;
}
std::pair<size_t,void*> TriMesh::pack( int level ) const
{
    std::pair<size_t,void*> data_out(0,NULL);
    if ( level==0 ) {
        const std::vector<Point>& points = vertices->points;
        data_out.first = (3+3*points.size())*sizeof(double) + 3*A.size()*sizeof(int);
        double *data_ptr = new double[4+3*points.size()+(3*A.size()*sizeof(int))/sizeof(double)];
        data_out.second = data_ptr;
        uint64_t *data_int64 = reinterpret_cast<uint64_t*>(data_ptr);
        data_int64[0] = level;
        data_int64[1] = points.size();
        data_int64[2] = A.size();
        double *data = &data_ptr[3];
        for (size_t i=0; i<points.size(); i++) {
            data[3*i+0] = points[i].x;
            data[3*i+1] = points[i].y;
            data[3*i+2] = points[i].z;
        }
        int *data_int = reinterpret_cast<int*>(&data[3*points.size()]);
        for (size_t i=0; i<A.size(); i++) {
            data_int[3*i+0] = A[i];
            data_int[3*i+1] = B[i];
            data_int[3*i+2] = C[i];
        }
    }
    return data_out;
}
void TriMesh::unpack( const std::pair<size_t,void*>& data_in )
{
    uint64_t *data_int64 = reinterpret_cast<uint64_t*>(data_in.second);
    const double *data = reinterpret_cast<const double*>(data_in.second);
    int level = data_int64[0];
    uint64_t N_P = data_int64[1];
    uint64_t N_A = data_int64[2];
    data = &data[3];
    if ( level==0 ) {
        size_t size = (3+3*N_P)*sizeof(double)+3*N_A*sizeof(int);
        ASSERT(size==data_in.first);
        vertices.reset( new PointList(N_P) );
        std::vector<Point>& points = vertices->points;
        for (size_t i=0; i<points.size(); i++) {
            points[i].x = data[3*i+0];
            points[i].y = data[3*i+1];
            points[i].z = data[3*i+2];
        }
        const int *data_int = reinterpret_cast<const int*>(&data[3*N_P]);
        A.resize(N_A);
        B.resize(N_A);
        C.resize(N_A);
        for (size_t i=0; i<A.size(); i++) {
            A[i] = data_int[3*i+0];
            B[i] = data_int[3*i+1];
            C[i] = data_int[3*i+2];
        }
    }
}


/****************************************************
* Domain mesh                                       *
****************************************************/
DomainMesh::DomainMesh():
    nprocx(0), nprocy(0), nprocz(0), rank(0),
    nx(0), ny(0), nz(0),
    Lx(0), Ly(0), Lz(0)
{
}
DomainMesh::DomainMesh( RankInfoStruct data, 
    int nx2, int ny2, int nz2, double Lx2, double Ly2, double Lz2 ):
    nprocx(data.nx), nprocy(data.ny), nprocz(data.nz), rank(data.rank[1][1][1]),
    nx(nx2), ny(ny2), nz(nz2),
    Lx(Lx2), Ly(Ly2), Lz(Lz2)
{
}
DomainMesh::~DomainMesh()
{
}
size_t DomainMesh::numberPointsVar( VariableType type ) const
{
    size_t N = 0;
    if ( type==NodeVariable )
        N = (nx+1)*(ny+1)*(nz+1);
    else if ( type==SurfaceVariable )
        N = (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1);
    else if ( type==VolumeVariable )
        N = nx*ny*nz;
    return N;
}
std::pair<size_t,void*> DomainMesh::pack( int level ) const
{
    std::pair<size_t,void*> data(0,NULL);
    data.first = 7*sizeof(double);
    data.second = new double[7];
    memset(data.second,0,7*sizeof(double));
    int *data_int = reinterpret_cast<int*>(data.second);
    double *data_double = &reinterpret_cast<double*>(data.second)[4];
    data_int[0] = nprocx;
    data_int[1] = nprocy;
    data_int[2] = nprocz;
    data_int[3] = rank;
    data_int[4] = nx;
    data_int[5] = ny;
    data_int[6] = nz;
    data_double[0] = Lx;
    data_double[1] = Ly;
    data_double[2] = Lz;
    return data;
}
void DomainMesh::unpack( const std::pair<size_t,void*>& data )
{
    const int *data_int = reinterpret_cast<const int*>(data.second);
    const double *data_double = &reinterpret_cast<const double*>(data.second)[4];
    nprocx = data_int[0];
    nprocy = data_int[1];
    nprocz = data_int[2];
    rank = data_int[3];
    nx = data_int[4];
    ny = data_int[5];
    nz = data_int[6];
    Lx = data_double[0];
    Ly = data_double[1];
    Lz = data_double[2];
}


/****************************************************
* Converters                                        *
****************************************************/
std::shared_ptr<PointList> getPointList( std::shared_ptr<Mesh> mesh )
{
    return std::dynamic_pointer_cast<PointList>(mesh);
}
std::shared_ptr<TriMesh> getTriMesh( std::shared_ptr<Mesh> mesh )
{
    std::shared_ptr<TriMesh> mesh2;
    if ( std::dynamic_pointer_cast<TriMesh>(mesh).get() != NULL ) {
        mesh2 = std::dynamic_pointer_cast<TriMesh>(mesh);
    } else if ( std::dynamic_pointer_cast<TriList>(mesh).get() != NULL ) {
        std::shared_ptr<TriList> trilist = std::dynamic_pointer_cast<TriList>(mesh);
        ASSERT(trilist.get()!=NULL);
        mesh2.reset( new TriMesh(*trilist) );
    }
    return mesh2;
}
std::shared_ptr<TriList> getTriList( std::shared_ptr<Mesh> mesh )
{
    std::shared_ptr<TriList> mesh2;
    if ( std::dynamic_pointer_cast<TriList>(mesh).get() != NULL ) {
        mesh2 = std::dynamic_pointer_cast<TriList>(mesh);
    } else if ( std::dynamic_pointer_cast<TriMesh>(mesh).get() != NULL ) {
        std::shared_ptr<TriMesh> trimesh = std::dynamic_pointer_cast<TriMesh>(mesh);
        ASSERT(trimesh.get()!=NULL);
        mesh2.reset( new TriList(*trimesh) );
    }
    return mesh2;
}
std::shared_ptr<const PointList> getPointList( std::shared_ptr<const Mesh> mesh )
{
    return getPointList( std::const_pointer_cast<Mesh>(mesh) );
}
std::shared_ptr<const TriMesh> getTriMesh( std::shared_ptr<const Mesh> mesh )
{
    return getTriMesh( std::const_pointer_cast<Mesh>(mesh) );
}
std::shared_ptr<const TriList> getTriList( std::shared_ptr<const Mesh> mesh )
{
    return getTriList( std::const_pointer_cast<Mesh>(mesh) );
}


} // IO namespace

