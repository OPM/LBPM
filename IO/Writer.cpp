#include "IO/Writer.h"
#include "IO/MeshDatabase.h"
#include "IO/IOHelpers.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include <sys/stat.h>
#include <algorithm>
#include <vector>
#include <set>


static bool global_summary_created = false;


// Write the mesh data in the original format
static std::vector<IO::MeshDatabase> writeMeshesOrigFormat( const std::vector<IO::MeshDataStruct>& meshData, const char* path )
{
    int rank = MPI_WORLD_RANK();
    std::vector<IO::MeshDatabase> meshes_written;
    for (size_t i=0; i<meshData.size(); i++) {
        char domainname[100], filename[100], fullpath[200];
        sprintf(domainname,"%05i",rank);
        sprintf(filename,"%s.%05i",meshData[i].meshName.c_str(),rank);
	    sprintf(fullpath,"%s/%s",path,filename);
        FILE *fid = fopen(fullpath,"wb");
        INSIST(fid!=NULL,std::string("Error opening file: ")+fullpath);
        std::shared_ptr<IO::Mesh> mesh = meshData[i].mesh;
        IO::MeshDatabase mesh_entry;
        mesh_entry.name = meshData[i].meshName;
        mesh_entry.type = meshType(mesh);
        mesh_entry.meshClass = meshData[i].mesh->className();
        mesh_entry.format = 1;
        IO::DatabaseEntry domain;
        domain.name = domainname;
        domain.file = filename;
        domain.offset = 0;
        mesh_entry.domains.push_back(domain);
        if ( !meshData[i].vars.empty() ) {
            printf("Warning: variables are not supported with this format\n");
            //for (size_t j=0; j<meshData[i].vars.size(); j++)
            //    mesh_entry.variables.push_back( meshData[i].vars[j]->name );
        }
        if ( std::dynamic_pointer_cast<IO::PointList>(mesh)!=NULL ) {
            // List of points
            std::shared_ptr<IO::PointList> pointlist = std::dynamic_pointer_cast<IO::PointList>(mesh);
            const std::vector<Point>& P = pointlist->points;
            for (size_t i=0; i<P.size(); i++) {
                double x[3];
                x[0] = P[i].x;  x[1] = P[i].y;  x[2] = P[i].z;
                fwrite(x,sizeof(double),3,fid);
            }
        } else if ( std::dynamic_pointer_cast<IO::TriList>(mesh)!=NULL || std::dynamic_pointer_cast<IO::TriMesh>(mesh)!=NULL ) {
            // Triangle mesh
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
        std::unique(mesh_entry.variables.begin(),mesh_entry.variables.end());
        meshes_written.push_back(mesh_entry);
    }
    return meshes_written;
}


// Write a mesh (and variables) to a file
static IO::MeshDatabase write_domain( FILE *fid, const std::string& filename,
    const IO::MeshDataStruct& mesh, int format )
{
    int rank = MPI_WORLD_RANK();
    char domainname[10];
    sprintf(domainname,"%05i",rank);
    int level = 0;
    // Create the MeshDatabase
    IO::MeshDatabase database;
    database.name = mesh.meshName;
    database.type = meshType(mesh.mesh);
    database.meshClass = mesh.mesh->className();
    database.format = format;
    // Write the mesh
    IO::DatabaseEntry domain;
    domain.name = domainname;
    domain.file = filename;
    domain.offset = ftell(fid);
    database.domains.push_back(domain);
    std::pair<size_t,void*> data = mesh.mesh->pack(level);
    fprintf(fid,"Mesh: %s-%05i: %lu\n",mesh.meshName.c_str(),rank,data.first);
    fwrite(data.second,1,data.first,fid);
    fprintf(fid,"\n");
    delete [] (char*) data.second;
    // Write the variables
    for (size_t i=0; i<mesh.vars.size(); i++) {
        IO::DatabaseEntry variable;
        variable.name = mesh.vars[i]->name;
        variable.file = filename;
        variable.offset = ftell(fid);
        IO::VariableDatabase info;
        info.name = variable.name;
        info.type = mesh.vars[i]->type;
        info.dim = mesh.vars[i]->dim;
        database.variables.push_back(info);
        std::pair<std::string,std::string> key(domainname,variable.name);
        database.variable_data.insert( 
            std::pair<std::pair<std::string,std::string>,IO::DatabaseEntry>(key,variable) );
        int dim = mesh.vars[i]->dim;
        int type = static_cast<int>(mesh.vars[i]->type);
        size_t N = mesh.vars[i]->data.size();
        const void* data = N==0 ? 0:&mesh.vars[i]->data[0];
        if ( type == static_cast<int>(IO::VariableType::Null) ) {
            ERROR("Variable type not set");
        }
        size_t N_mesh = mesh.mesh->numberPointsVar(mesh.vars[i]->type);
        ASSERT(N==dim*N_mesh);
        fprintf(fid,"Var: %s-%05i-%s: %i, %i, %lu, %lu, double\n",
            database.name.c_str(), rank, variable.name.c_str(),
            dim, type, N, dim*N*sizeof(double) );
        fwrite(data,sizeof(double),dim*N,fid);
        fprintf(fid,"\n");
    }
    return database;
}


// Write the mesh data in the new format
static std::vector<IO::MeshDatabase> writeMeshesNewFormat( 
    const std::vector<IO::MeshDataStruct>& meshData, const char* path, int format )
{
    int rank = MPI_WORLD_RANK();
    std::vector<IO::MeshDatabase> meshes_written;
    char filename[100], fullpath[200];
    sprintf(filename,"%05i",rank);
    sprintf(fullpath,"%s/%s",path,filename);
    FILE *fid = fopen(fullpath,"wb");
    for (size_t i=0; i<meshData.size(); i++) {
        std::shared_ptr<IO::Mesh> mesh = meshData[i].mesh;
        meshes_written.push_back( write_domain(fid,filename,meshData[i],format) );
    }
    fclose(fid);
    return meshes_written;
}


// Write the mesh data
void IO::writeData( int timestep, const std::vector<IO::MeshDataStruct>& meshData, int format )
{
    int rank = MPI_WORLD_RANK();
    int size = MPI_WORLD_SIZE();
    // Create the output directory
    char path[100];
    sprintf(path,"vis%03i",timestep);
    if ( rank == 0 )
        mkdir(path,S_IRWXU|S_IRGRP);
    MPI_Barrier(MPI_COMM_WORLD);
    // Write the mesh files
    std::vector<IO::MeshDatabase> meshes_written;
    if ( format == 1 ) {
        // Write the original triangle format
        meshes_written = writeMeshesOrigFormat( meshData, path );
    } else if ( format == 2 ) {
        // Write the new format
        meshes_written = writeMeshesNewFormat( meshData, path, format );
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


