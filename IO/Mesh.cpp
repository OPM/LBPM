#include "Mesh.h"
#include "common/Utilities.h"

#include <limits>


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
    if ( std::dynamic_pointer_cast<TriMesh>(mesh) != NULL ) {
        mesh2 = std::dynamic_pointer_cast<TriMesh>(mesh);
    } else if ( std::dynamic_pointer_cast<TriList>(mesh) != NULL ) {
        std::shared_ptr<TriList> trilist = std::dynamic_pointer_cast<TriList>(mesh);
        ASSERT(trilist!=NULL);
        mesh2.reset( new TriMesh(*trilist) );
    }
    return mesh2;
}
std::shared_ptr<TriList> getTriList( std::shared_ptr<Mesh> mesh )
{
    std::shared_ptr<TriList> mesh2;
    if ( std::dynamic_pointer_cast<TriList>(mesh) != NULL ) {
        mesh2 = std::dynamic_pointer_cast<TriList>(mesh);
    } else if ( std::dynamic_pointer_cast<TriMesh>(mesh) != NULL ) {
        std::shared_ptr<TriMesh> trimesh = std::dynamic_pointer_cast<TriMesh>(mesh);
        ASSERT(trimesh!=NULL);
        mesh2.reset( new TriList(*trimesh) );
    }
    return mesh2;
}


} // IO namespace

