#include "IO/Writer.h"
#include "IO/MeshDatabase.h"
#include "common/Utilities.h"

#include "mpi.h"
#include <sys/stat.h>
#include <algorithm>
#include <vector>
#include <set>


static bool global_summary_created = false;

// Write the mesh data in the original format
static std::vector<IO::MeshDatabase> writeMeshesOrigFormat( const std::vector<IO::MeshDataStruct>& meshData, const char* path )
{
    int rank = -1;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    std::vector<IO::MeshDatabase> meshes_written;
    for (size_t i=0; i<meshData.size(); i++) {
        char domain[100], filename[100], fullpath[200];
        sprintf(domain,"%05i",rank);
        sprintf(filename,"%s.%05i",meshData[i].meshName.c_str(),rank);
	    sprintf(fullpath,"%s/%s",path,filename);
        FILE *fid = fopen(fullpath,"wb");
        std::shared_ptr<IO::Mesh> mesh = meshData[i].mesh;
        IO::MeshDatabase mesh_entry;
        mesh_entry.name = meshData[i].meshName;
        mesh_entry.format = 1;
        mesh_entry.domains.push_back(domain);
        mesh_entry.file.push_back(filename);
        if ( std::dynamic_pointer_cast<IO::PointList>(mesh)!=NULL ) {
            // List of points
            mesh_entry.type = IO::MeshType::PointMesh;
            std::shared_ptr<IO::PointList> pointlist = std::dynamic_pointer_cast<IO::PointList>(mesh);
            const std::vector<Point>& P = pointlist->points;
            for (size_t i=0; i<P.size(); i++) {
                double x[3];
                x[0] = P[i].x;  x[1] = P[i].y;  x[2] = P[i].z;
                fwrite(x,sizeof(double),3,fid);
            }
        } else if ( std::dynamic_pointer_cast<IO::TriList>(mesh)!=NULL || std::dynamic_pointer_cast<IO::TriMesh>(mesh)!=NULL ) {
            // Triangle mesh
            mesh_entry.type = IO::MeshType::SurfaceMesh;
            std::shared_ptr<IO::TriList> trilist = IO::getTriList(mesh);
            const std::vector<Point>& A = trilist->A;
            const std::vector<Point>& B = trilist->B;
            const std::vector<Point>& C = trilist->C;
            for (size_t i=0; i<A.size(); i++) {
                double tri[9];
                tri[0] = A[i].x;  tri[1] = A[i].y;  tri[2] = A[i].z;
                tri[3] = B[i].x;  tri[4] = B[i].y;  tri[5] = B[i].z;
                tri[6] = C[i].x;  tri[7] = C[i].y;  tri[8] = C[i].z;
                fwrite(tri,sizeof(double),9,fid);
            }
        } else {
            ERROR("Unknown mesh");
        }
        fclose(fid);
        meshes_written.push_back(mesh_entry);
    }
    return meshes_written;
}


// Write the mesh data
void IO::writeData( int timestep, const std::vector<IO::MeshDataStruct>& meshData )
{
    const int format = 1;
    int rank = -1;
    int size = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    // Create the output directory
    char path[100];
    sprintf(path,"vis%03i",timestep);
    if ( rank == 0 )
        mkdir(path,S_IRWXU|S_IRGRP);
    // Write the mesh files
    std::vector<IO::MeshDatabase> meshes_written;
    if ( format == 1 ) {
        // Write the original triangle format
        meshes_written = writeMeshesOrigFormat( meshData, path );
    } else {
        ERROR("Unknown format");
    }
    // Gather a complete list of files on rank 0
    meshes_written = gatherAll(meshes_written,MPI_COMM_WORLD);
    // Write the summary files
    if ( rank == 0 ) {
        // Write the summary file for the current timestep
        char filename[200];
        sprintf(filename,"%s/LBM.summary",path);
        write(meshes_written,filename);
        // Add the timestep to the global summary file
        FILE *fid = NULL;
        if ( !global_summary_created ) {
            fid = fopen("summary.LBM","wb");
            global_summary_created = true;
        } else {
            fid = fopen("summary.LBM","ab");
        }
        fprintf(fid,"%s\n",path);
        fclose(fid);
    }
}


