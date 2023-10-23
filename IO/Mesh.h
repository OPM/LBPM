#ifndef MESH_INC
#define MESH_INC

#include <iostream>
#include <memory>
#include <string.h>
#include <vector>

#include "analysis/PointList.h"
#include "common/Array.h"
#include "common/Communication.h"


namespace IO {


//! Enums to define types
enum class VariableType {
    NodeVariable,
    EdgeVariable,
    SurfaceVariable,
    VolumeVariable,
    NullVariable
};
enum class DataType { Double, Float, Int, Null };
enum class MeshType { PointMesh, SurfaceMesh, VolumeMesh, Unknown };
enum class FileFormat { OLD, NEW, NEW_SINGLE, SILO, HDF5 };


//! Convert enums to/from strings (more future-proof than static_cast<int>)
std::string getString( VariableType );
std::string getString( DataType );
std::string getString( MeshType );
std::string getString( FileFormat );
VariableType getVariableType( const std::string & );
DataType getDataType( const std::string & );
MeshType getMeshType( const std::string & );
FileFormat getFileFormat( const std::string & );


/*! \class Mesh
    \brief A base class for meshes
*/
class Mesh
{
public:
    //! Destructor
    virtual ~Mesh();
    //! Mesh class name (eg. PointList)
    virtual std::string className() const = 0;
    //! Number of points for the given variable type
    virtual size_t numberPointsVar( VariableType type ) const = 0;
    //! Pack the data
    virtual std::pair<size_t, void *> pack( int level ) const = 0;
    //! Unpack the data
    virtual void unpack( const std::pair<size_t, void *> &data ) = 0;

protected:
    //! Empty constructor
    Mesh();
    Mesh( const Mesh & );
    Mesh &operator=( const Mesh & );
};


/*! \class PointList
    \brief A class used to hold a list of verticies
*/
class PointList : public Mesh
{
public:
    //! Empty constructor
    PointList();
    //! Constructor for N points
    PointList( size_t N );
    //! Destructor
    virtual ~PointList();
    //! Mesh class name
    virtual std::string className() const { return "PointList"; }
    //! Number of points for the given variable type
    virtual size_t numberPointsVar( VariableType type ) const;
    //! Pack the data
    virtual std::pair<size_t, void *> pack( int level ) const;
    //! Unpack the data
    virtual void unpack( const std::pair<size_t, void *> &data );
    //! Access the points
    const std::vector<Point> &getPoints() const { return points; }

public:
    std::vector<Point> points; //!< List of points vertex
};


/*! \class TriList
    \brief A class used to hold a list of triangles specified by their vertex coordinates
*/
class TriMesh;
class TriList : public Mesh
{
public:
    //! Empty constructor
    TriList();
    //! Constructor for N triangles
    TriList( size_t N_tri );
    //! Constructor from TriMesh
    TriList( const TriMesh & );
    //! Destructor
    virtual ~TriList();
    //! Mesh class name
    virtual std::string className() const { return "TriList"; }
    //! Number of points for the given variable type
    virtual size_t numberPointsVar( VariableType type ) const;
    //! Pack the data
    virtual std::pair<size_t, void *> pack( int level ) const;
    //! Unpack the data
    virtual void unpack( const std::pair<size_t, void *> &data );

public:
    std::vector<Point> A; //!< First vertex
    std::vector<Point> B; //!< Second vertex
    std::vector<Point> C; //!< Third vertex
};


/*! \class TriMesh
    \brief A class used to hold a list of trianges specified by their vertex number and list of
   coordiantes
*/
class TriMesh : public Mesh
{
public:
    //! TriMesh constructor
    TriMesh();
    //! Constructor for Nt triangles and Np points
    TriMesh( size_t N_tri, size_t N_point );
    //! Constructor for Nt triangles and the given points
    TriMesh( size_t N_tri, std::shared_ptr<PointList> points );
    //! Constructor from TriList
    TriMesh( const TriList & );
    //! Destructor
    virtual ~TriMesh();
    //! Mesh class name
    virtual std::string className() const { return "TriMesh"; }
    //! Number of points for the given variable type
    virtual size_t numberPointsVar( VariableType type ) const;
    //! Pack the data
    virtual std::pair<size_t, void *> pack( int level ) const;
    //! Unpack the data
    virtual void unpack( const std::pair<size_t, void *> &data );

public:
    std::shared_ptr<PointList> vertices; //!< List of verticies
    std::vector<int> A;                  //!< First vertex
    std::vector<int> B;                  //!< Second vertex
    std::vector<int> C;                  //!< Third vertex
};


/*! \class Domain
    \brief A class used to hold the domain
*/
class DomainMesh : public Mesh
{
public:
    //! Empty constructor
    DomainMesh();
    //! Default constructor
    DomainMesh( RankInfoStruct rank_data, int nx, int ny, int nz, double Lx, double Ly, double Lz );
    //! Destructor
    virtual ~DomainMesh();
    //! Mesh class name
    virtual std::string className() const { return "DomainMesh"; }
    //! Number of points for the given variable type
    virtual size_t numberPointsVar( VariableType type ) const;
    //! Pack the data
    virtual std::pair<size_t, void *> pack( int level ) const;
    //! Unpack the data
    virtual void unpack( const std::pair<size_t, void *> &data );

public:
    int nprocx, nprocy, nprocz, rank;
    int nx, ny, nz;
    double Lx, Ly, Lz;
};


/*! \class Variable
    \brief A base class for variables
*/
struct Variable {
public:
    // Internal variables
    unsigned char dim;  //!< Number of points per grid point (1: scalar, 3: vector, ...)
    VariableType type;  //!< Variable type
    DataType precision; //!< Variable precision to use for IO
    std::string name;   //!< Variable name
    Array<double> data; //!< Variable data
    //! Empty constructor
    Variable() : dim( 0 ), type( VariableType::NullVariable ), precision( DataType::Double ) {}
    //! Constructor
    Variable( int dim_, IO::VariableType type_, const std::string &name_ )
        : dim( dim_ ), type( type_ ), precision( DataType::Double ), name( name_ )
    {
    }
    //! Constructor
    Variable(
        int dim_, IO::VariableType type_, const std::string &name_, const Array<double> &data_ )
        : dim( dim_ ), type( type_ ), precision( DataType::Double ), name( name_ ), data( data_ )
    {
    }
    //! Destructor
    virtual ~Variable() {}

protected:
    //! Empty constructor
    Variable( const Variable & );
    Variable &operator=( const Variable & );
};


/*! \class MeshDataStruct
    \brief A class used to hold database info for saving a mesh
*/
struct MeshDataStruct {
    DataType precision;         //!< Precision to use for IO (mesh)
    std::string meshName;       //!< Mesh name
    std::shared_ptr<Mesh> mesh; //!< Mesh data
    std::vector<std::shared_ptr<Variable>> vars;
    //! Empty constructor
    MeshDataStruct() : precision( DataType::Double ) {}
    //! Check the data
    bool check( bool abort = true ) const;
};


//! Convert the mesh to a TriMesh (will return NULL if this is invalid)
std::shared_ptr<PointList> getPointList( std::shared_ptr<Mesh> mesh );
std::shared_ptr<TriMesh> getTriMesh( std::shared_ptr<Mesh> mesh );
std::shared_ptr<TriList> getTriList( std::shared_ptr<Mesh> mesh );
std::shared_ptr<const PointList> getPointList( std::shared_ptr<const Mesh> mesh );
std::shared_ptr<const TriMesh> getTriMesh( std::shared_ptr<const Mesh> mesh );
std::shared_ptr<const TriList> getTriList( std::shared_ptr<const Mesh> mesh );


} // namespace IO

#endif
