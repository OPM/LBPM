// vtk headers
#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkCellType.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkStructuredGrid.h>


// LBPM headers
#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"


extern std::ostream DebugStream1;
extern std::ostream DebugStream2;


// Convert a point array to vtkFloatArray
inline vtkFloatArray* pointToVTK( const std::vector<Point>& points )
{
    vtkFloatArray *coords = vtkFloatArray::New();
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(points.size());
    for (size_t i=0; i<points.size(); i++)
        coords->SetTuple3(i,points[i].x,points[i].y,points[i].z);
    return coords;
}


// Return the mesh type
inline avtMeshType vtkMeshType( IO::MeshType meshType )
{
    avtMeshType vtkType = AVT_UNKNOWN_MESH;
    if ( meshType==IO::PointMesh )
        vtkType = AVT_POINT_MESH;
    else if ( meshType==IO::SurfaceMesh )
        vtkType = AVT_SURFACE_MESH;
    else if ( meshType==IO::VolumeMesh )
        vtkType = AVT_RECTILINEAR_MESH;
    return vtkType;
}


// Return the topological dimension
inline int vtkTopDim( IO::MeshType meshType )
{
    int dim = -1;
    if ( meshType==IO::PointMesh )
        dim = 1;
    else if ( meshType==IO::SurfaceMesh )
        dim = 2;
    else if ( meshType==IO::VolumeMesh )
        dim = 3;
    return dim;
}


// Convert a PointList to vtkDataSet
inline vtkDataSet* PointListToVTK( std::shared_ptr<const IO::PointList> mesh )
{
    vtkFloatArray *coords = pointToVTK(mesh->points);
    vtkPoints *points = vtkPoints::New(VTK_FLOAT);
    points->SetData(coords);
    points->ComputeBounds();
    size_t N = mesh->points.size();
    vtkPolyData *vtkMesh = vtkPolyData::New();
    vtkMesh->SetPoints(points);
    vtkMesh->Allocate(N);
    for (int i=0; i<(int)N; i++) {
        vtkIdType ids[1] = {i};
        vtkMesh->InsertNextCell(VTK_VERTEX,1,ids);
    }
    // vtkMesh->BuildCells();
    vtkMesh->ComputeBounds();
    points->Delete();
    coords->Delete();
    return vtkMesh;
}


// Convert a TriMesh to vtkDataSet
inline vtkDataSet* TriMeshToVTK( std::shared_ptr<const IO::TriMesh> mesh )
{
    vtkFloatArray *coords = pointToVTK(mesh->vertices->points);
    ASSERT(coords!=NULL);
    vtkPoints *points = vtkPoints::New(VTK_FLOAT);
    points->SetData(coords);
    points->ComputeBounds();
    const std::vector<int>& A = mesh->A;
    const std::vector<int>& B = mesh->B;
    const std::vector<int>& C = mesh->C;
    size_t N_tri = A.size();
    vtkPolyData *vtkMesh = vtkPolyData::New();
    vtkMesh->SetPoints(points);
    vtkMesh->Allocate(N_tri);
    for (size_t i=0; i<N_tri; i++) {
        vtkIdType ids[3] = {A[i],B[i],C[i]};
        vtkMesh->InsertNextCell(VTK_TRIANGLE,3,ids);
    }
    vtkMesh->BuildCells();
    vtkMesh->ComputeBounds();
    points->Delete();
    coords->Delete();
    return vtkMesh;
}


// Convert a TriList to vtkDataSet
inline vtkDataSet* TriListToVTK( std::shared_ptr<const IO::TriList> mesh )
{
    std::vector<LBPM_Point> point_set(3*mesh->A.size());
    for (size_t i=0; i<mesh->A.size(); i++) {
        point_set[3*i+0] = mesh->A[i];
        point_set[3*i+1] = mesh->B[i];
        point_set[3*i+2] = mesh->C[i];
    }
    vtkFloatArray *coords = pointToVTK(point_set);
    ASSERT(coords!=NULL);
    vtkPoints *points = vtkPoints::New(VTK_FLOAT);
    points->SetData(coords);
    points->ComputeBounds();
    size_t N_tri = mesh->A.size();
    vtkPolyData *vtkMesh = vtkPolyData::New();
    vtkMesh->SetPoints(points);
    vtkMesh->Allocate(N_tri);
    for (int i=0; i<(int)N_tri; i++) {
        vtkIdType ids[3] = {3*i+0,3*i+1,3*i+2};
        vtkMesh->InsertNextCell(VTK_TRIANGLE,3,ids);
    }
    vtkMesh->BuildCells();
    vtkMesh->ComputeBounds();
    points->Delete();
    coords->Delete();
    return vtkMesh;
}


// Convert a DomainMesh to vtkDataSet
inline vtkDataSet* DomainToVTK( std::shared_ptr<const IO::DomainMesh> mesh )
{
    int nx = mesh->nx;
    int ny = mesh->ny;
    int nz = mesh->nz;
    if ( nx==0 || ny==0 || nz==0 ) {
        DebugStream::Stream1() << "   Domain mesh is empty" << std::endl;
        return NULL;
    }
    RankInfoStruct rank_data(mesh->rank,mesh->nprocx,mesh->nprocy,mesh->nprocz);
    vtkFloatArray *x = vtkFloatArray::New();
    vtkFloatArray *y = vtkFloatArray::New();
    vtkFloatArray *z = vtkFloatArray::New();
    for(int i=0; i<=nx; i++)
        x->InsertNextValue(static_cast<double>(nx*rank_data.ix+i)/static_cast<double>(nx*mesh->nprocx));
    for(int j=0; j<=ny; j++)
        y->InsertNextValue(static_cast<double>(ny*rank_data.jy+j)/static_cast<double>(ny*mesh->nprocy));
    for(int k=0; k<=nz; k++)
        z->InsertNextValue(static_cast<double>(nz*rank_data.kz+k)/static_cast<double>(nz*mesh->nprocz));
    vtkRectilinearGrid *vtkMesh = vtkRectilinearGrid::New();
    vtkMesh->SetDimensions(nx+1,ny+1,nz+1);
    vtkMesh->SetXCoordinates(x);
    vtkMesh->SetYCoordinates(y);
    vtkMesh->SetZCoordinates(z);
    vtkMesh->ComputeBounds();
    x->Delete();
    y->Delete();
    z->Delete();
    return vtkMesh;
}


// Convert a mesh to vtkDataSet
inline vtkDataSet* meshToVTK( std::shared_ptr<const IO::Mesh> mesh )
{
    if ( mesh.get() == NULL ){
        DebugStream::Stream1() << "   Mesh is NULL" << std::endl;
        return NULL;
    }
    vtkDataSet* mesh2 = NULL;
    const std::string className = mesh->className();
    if ( className == "PointList" ) {
        // We are dealing with a point mesh
        mesh2 = PointListToVTK( std::dynamic_pointer_cast<const IO::PointList>(mesh) );
        DebugStream1 << "   Point mesh created" << std::endl;
    } else if ( className == "TriMesh" ) {
        mesh2 = TriMeshToVTK( std::dynamic_pointer_cast<const IO::TriMesh>(mesh) );
        DebugStream1 << "   Surface mesh created" << std::endl;
    } else if ( className == "TriList" ) {
        mesh2 = TriListToVTK( std::dynamic_pointer_cast<const IO::TriList>(mesh) );
        DebugStream1 << "   Surface mesh created" << std::endl;
    } else if ( className == "DomainMesh" ) {
        mesh2 = DomainToVTK( std::dynamic_pointer_cast<const IO::DomainMesh>(mesh) );
        DebugStream1 << "   Volume mesh created" << std::endl;
    } else {
        //DebugStream1 << "   Error, unknown mesh type" << std::endl;
        return NULL;
    }
    return mesh2;
}


// Convert a variable to vtkDataSet
inline vtkDataArray* varToVTK( std::shared_ptr<const IO::Variable> var )
{
    vtkFloatArray* var2 = NULL;
    if ( var->dim==1 ) {
        // Scalar variable
        var2 = vtkFloatArray::New();
        var2->SetNumberOfTuples(var->data.length());
        for (size_t i=0; i<var->data.length(); i++)
            var2->SetTuple1(i,var->data(i));
    } else {
        DebugStream::Stream1() << "   Error, variable not yet supported" << std::endl;
        return NULL;
    }
    return var2;
}


