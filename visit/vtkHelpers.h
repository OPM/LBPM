// vtk headers
#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkCellType.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>


// LBPM headers
#include "IO/Mesh.h"


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
    if ( meshType==IO::MeshType::PointMesh )
        vtkType = AVT_POINT_MESH;
    else if ( meshType==IO::MeshType::SurfaceMesh )
        vtkType = AVT_SURFACE_MESH;
    else if ( meshType==IO::MeshType::VolumeMesh )
        vtkType = AVT_UNSTRUCTURED_MESH;
    return vtkType;
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


// Convert a mesh to vtkDataSet
inline vtkDataSet* meshToVTK( std::shared_ptr<const IO::Mesh> mesh )
{
    vtkDataSet* mesh2 = NULL;
    if ( std::dynamic_pointer_cast<const IO::PointList>(mesh) != NULL ) {
        // We are dealing with a point mesh
        mesh2 = PointListToVTK( std::dynamic_pointer_cast<const IO::PointList>(mesh) );
        DebugStream::Stream1() << "   Point mesh created" << std::endl;
    } else if ( std::dynamic_pointer_cast<const IO::TriMesh>(mesh) != NULL ) {
        mesh2 = TriMeshToVTK( std::dynamic_pointer_cast<const IO::TriMesh>(mesh) );
        DebugStream::Stream1() << "   Surface mesh created" << std::endl;
    } else if ( std::dynamic_pointer_cast<const IO::TriList>(mesh) != NULL ) {
        mesh2 = TriListToVTK( std::dynamic_pointer_cast<const IO::TriList>(mesh) );
        DebugStream::Stream1() << "   Surface mesh created" << std::endl;
    } else {
        DebugStream::Stream1() << "   Error, unknown mesh type" << std::endl;
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
        var2->SetNumberOfTuples(var->data.size());
        for (size_t i=0; i<var->data.size(); i++)
            var2->SetTuple1(i,var->data[i]);
    } else {
        DebugStream::Stream1() << "   Error, variable not yet supported" << std::endl;
        return NULL;
    }
    return var2;
}


