#ifndef MESH_INC
#define MESH_INC

#include <iostream>
#include <string.h>
#include <memory>
#include <vector>

#include "common/PointList.h"


namespace IO {


/*! \class Mesh
    \brief A base class for meshes
*/
class Mesh
{
public:
    //! Destructor
    virtual ~Mesh();
protected:
    //! Empty constructor
    Mesh();
    Mesh(const Mesh&);
    Mesh& operator=(const Mesh&);
};


/*! \class PointList
    \brief A class used to hold a list of verticies
*/
class PointList: public Mesh
{
public:
    //! Empty constructor
    PointList();
    //! Constructor for N points
    PointList( size_t N );
    //! Destructor
    virtual ~PointList();
public:
    std::vector<Point>  points;  //!< List of points vertex
};


/*! \class TriMesh
    \brief A class used to hold a list of trianges specified by their vertex number and list of coordiantes
*/
class TriList;
class TriMesh: public Mesh
{
public:
    //! TriMesh constructor
    TriMesh();
    //! Constructor for Nt triangles and Np points
    TriMesh( size_t N_tri, size_t N_point );
    //! Constructor for Nt triangles and the given points
    TriMesh( size_t N_tri, std::shared_ptr<PointList> points );
    //! Constructor from TriList
    TriMesh( const TriList& );
    //! Destructor
    virtual ~TriMesh();
public:
    std::shared_ptr<PointList> vertices;    //!< List of verticies
    std::vector<int>    A;                  //!< First vertex
    std::vector<int>    B;                  //!< Second vertex
    std::vector<int>    C;                  //!< Third vertex
};


/*! \class TriList
    \brief A class used to hold a list of triangles specified by their vertex coordinates
*/
class TriList: public Mesh
{
public:
    //! Empty constructor
    TriList();
    //! Constructor for N triangles
    TriList( size_t N_tri );
    //! Constructor from TriMesh
    TriList( const TriMesh& );
    //! Destructor
    virtual ~TriList();
public:
    std::vector<Point>  A;      //!< First vertex
    std::vector<Point>  B;      //!< Second vertex
    std::vector<Point>  C;      //!< Third vertex
};



/*! \class MeshDataStruct
    \brief A class used to hold database info for saving a mesh
*/
struct MeshDataStruct {
    std::string             meshName;
    std::shared_ptr<Mesh>   mesh;
    //std::vector<std::string> dataName;
    //std::vector<int>        dataType;
    //std::vector<double*>    data;
};


//! Convert the mesh to a TriMesh (will return NULL if this is invalid)
std::shared_ptr<PointList> getPointList( std::shared_ptr<Mesh> mesh );
std::shared_ptr<TriMesh> getTriMesh( std::shared_ptr<Mesh> mesh );
std::shared_ptr<TriList> getTriList( std::shared_ptr<Mesh> mesh );


} // IO namespace

#endif

