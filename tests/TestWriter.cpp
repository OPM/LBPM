#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <memory>

#include "common/UnitTest.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"


inline bool approx_equal( const Point& A, const Point& B )
{
    double tol = 1e-8*sqrt(A.x*A.x+A.y*A.y+A.z*A.z); 
    return fabs(A.x-B.x)<=tol && fabs(A.y-B.y)<=tol && fabs(A.z-B.z)<=tol;
}
inline bool approx_equal( const double& A, const double& B )
{
    return fabs(A-B)<=1e-8*fabs(A+B);
}


inline double distance( const Point& p )
{
    return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}



int main(int argc, char **argv)
{
    int rank,nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    Utilities::setAbortBehavior(true,true,false);
    Utilities::setErrorHandlers();
    UnitTest ut;

    // Create some points
    const int N_points = 8;
    const int N_tri = 12;
    double x[8] = { 0, 1, 0, 1, 0, 1, 0, 1 };
    double y[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };
    double z[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
    int tri[N_tri][3] = { 
        {0,1,3}, {0,3,2}, {4,5,7}, {4,7,6}, // z faces
        {0,1,4}, {1,4,5}, {2,3,6}, {3,6,7}, // y faces
        {0,2,4}, {2,4,6}, {1,3,5}, {3,5,7}  // x faces
    };

    // Create the meshes
    std::shared_ptr<IO::PointList> set1( new IO::PointList(N_points) );
    for (int i=0; i<N_points; i++) {
        set1->points[i].x = x[i];
        set1->points[i].y = y[i];
        set1->points[i].z = z[i];
    }
    std::shared_ptr<IO::TriMesh> trimesh( new IO::TriMesh(N_tri,set1) );
    for (int i=0; i<N_tri; i++) {
        trimesh->A[i] = tri[i][0];
        trimesh->B[i] = tri[i][1];
        trimesh->C[i] = tri[i][2];
    }
    std::shared_ptr<IO::TriList> trilist( new IO::TriList(*trimesh) );
    for (int i=0; i<N_tri; i++) {
        Point A(x[tri[i][0]],y[tri[i][0]],z[tri[i][0]]);
        Point B(x[tri[i][1]],y[tri[i][1]],z[tri[i][1]]);
        Point C(x[tri[i][2]],y[tri[i][2]],z[tri[i][2]]);
        if ( !approx_equal(trilist->A[i],A) || !approx_equal(trilist->B[i],B) || !approx_equal(trilist->C[i],C) )
        {
            printf("Failed to create trilist\n");
            return -1;
        }
    }

    // Create the variables
    std::shared_ptr<IO::Variable> dist_set1( new IO::Variable() );
    std::shared_ptr<IO::Variable> dist_list( new IO::Variable() );
    dist_set1->dim = 1;
    dist_list->dim = 1;
    dist_set1->name = "Distance";
    dist_list->name = "Distance";
    dist_set1->type = IO::VariableType::NodeVariable;
    dist_list->type = IO::VariableType::NodeVariable;
    dist_set1->data.resize( N_points );
    for (int i=0; i<N_points; i++)
        dist_set1->data[i] = distance(set1->points[i]);
    dist_list->data.resize( 3*N_tri );
    for (int i=0; i<N_tri; i++) {
        dist_list->data[3*i+0] = distance(trilist->A[i]);
        dist_list->data[3*i+1] = distance(trilist->B[i]);
        dist_list->data[3*i+2] = distance(trilist->C[i]);
    }

    // Create the MeshDataStruct
    std::vector<IO::MeshDataStruct> meshData(3);
    meshData[0].meshName = "pointmesh";
    meshData[0].mesh = set1;
    meshData[0].vars.push_back(dist_set1);
    meshData[1].meshName = "trimesh";
    meshData[1].mesh = trimesh;
    meshData[1].vars.push_back(dist_set1);
    meshData[2].meshName = "trilist";
    meshData[2].mesh = trilist;
    meshData[2].vars.push_back(dist_list);

    // Write the data
    IO::writeData( 0, meshData, 1 );
    IO::writeData( 3, meshData, 2 );
    MPI_Barrier(MPI_COMM_WORLD);

    // Get a list of the timesteps 
    std::vector<std::string> timesteps = IO::readTimesteps("summary.LBM");
    if ( timesteps.size()==2 )
        ut.passes("Corrent number of timesteps");
    else
        ut.failure("Incorrent number of timesteps");

    // Check the mesh lists
    for (size_t i=0; i<timesteps.size(); i++) {
        // Load the list of meshes and check its size
        std::vector<IO::MeshDatabase> list = IO::getMeshList(".",timesteps[i]);
        if ( list.size()==meshData.size() )
            ut.passes("Corrent number of meshes found");
        else
            ut.failure("Incorrent number of meshes found");
        // Check the number of domains for each mesh
        bool pass = true;
        for (size_t j=0; j<list.size(); j++)
            pass = pass && (int)list[j].domains.size()==nprocs;
        if ( pass ) {
            ut.passes("Corrent number of domains for mesh");
        } else {
            ut.failure("Incorrent number of domains for mesh");
            continue;
        }
        // For each domain, load the mesh and check its data
        for (size_t j=0; j<list.size(); j++) {
            for (size_t k=0; k<list[i].domains.size(); k++) {
                std::shared_ptr<IO::Mesh> mesh = IO::getMesh(".",timesteps[i],list[j],k);
                if ( mesh==NULL ) {
                    printf("Failed to load %s\n",meshData[i].meshName.c_str());
                    pass = false;
                    break;
                }
                if ( meshData[j].meshName=="pointmesh" ) {
                    // Check the pointmesh
                    std::shared_ptr<IO::PointList> pmesh = IO::getPointList(mesh);
                    if ( pmesh==NULL ) {
                        pass = false;
                        break;
                    }
                    if ( pmesh->points.size()!=N_points ) {
                        pass = false;
                        break;
                    }                    
                }
                if ( meshData[j].meshName=="trimesh" || meshData[j].meshName=="trilist" ) {
                    // Check the trimesh/trilist
                    std::shared_ptr<IO::TriMesh> mesh1 = IO::getTriMesh(mesh);
                    std::shared_ptr<IO::TriList> mesh2 = IO::getTriList(mesh);
                    if ( mesh1==NULL || mesh2==NULL ) {
                        pass = false;
                        break;
                    }
                    if ( mesh1->A.size()!=N_tri || mesh1->B.size()!=N_tri || mesh1->C.size()!=N_tri || 
                         mesh2->A.size()!=N_tri || mesh2->B.size()!=N_tri || mesh2->C.size()!=N_tri ) 
                    {
                        pass = false;
                        break;
                    }
                    const std::vector<Point>& P1 = mesh1->vertices->points;
                    const std::vector<int>& A1 = mesh1->A;
                    const std::vector<int>& B1 = mesh1->B;
                    const std::vector<int>& C1 = mesh1->C;
                    const std::vector<Point>& A2 = mesh2->A;
                    const std::vector<Point>& B2 = mesh2->B;
                    const std::vector<Point>& C2 = mesh2->C;
                    const std::vector<Point>& A = trilist->A;
                    const std::vector<Point>& B = trilist->B;
                    const std::vector<Point>& C = trilist->C;
                    for (size_t i=0; i<N_tri; i++) {
                        if ( !approx_equal(P1[A1[i]],A[i]) || !approx_equal(P1[B1[i]],B[i]) || !approx_equal(P1[C1[i]],C[i]) )
                            pass = false;
                        if ( !approx_equal(A2[i],A[i]) || !approx_equal(B2[i],B[i]) || !approx_equal(C2[i],C[i]) )
                            pass = false;
                    }
                }
            }
        }
        if ( pass ) {
            ut.passes("Meshes loaded correctly");
        } else {
            ut.failure("Meshes did not load correctly");
            continue;
        }
        // For each domain, load the variables and check their data
        if ( i==0 )
            continue;   // Format 1 does not support variables
        for (size_t j=0; j<list.size(); j++) {
            const IO::MeshDataStruct* mesh0 = NULL;
            for (size_t k=0; k<meshData.size(); k++) {
                if ( meshData[k].meshName == list[j].name ) {
                    mesh0 = &meshData[k];
                    break;
                }
            }
            for (size_t v=0; v<list[i].variables.size(); v++) {
                for (size_t k=0; k<list[i].domains.size(); k++) {
                    std::shared_ptr<const IO::Variable> variable = 
                        IO::getVariable(".",timesteps[i],list[j],k,list[j].variables[v].name);
                    const IO::Variable& var1 = *mesh0->vars[v];
                    const IO::Variable& var2 = *variable;
                    pass = var1.name == var2.name;
                    pass = pass && var1.dim == var2.dim;
                    pass = pass && var1.type == var2.type;
                    pass = pass && var1.data.size() == var2.data.size();
                    if ( pass ) {
                        for (size_t m=0; m<var1.data.size(); m++)
                            pass = pass && approx_equal(var1.data[m],var2.data[m]);
                    }
                    if ( pass ) {
                        ut.passes("Variables match");
                    } else {
                        ut.failure("Variables did not match");
                        break;
                    }
                }
            }
        }
    }

    // Finished
    ut.report();
    int N_errors = ut.NumFailGlobal();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return N_errors;
}


