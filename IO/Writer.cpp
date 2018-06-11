/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "IO/Writer.h"
#include "IO/MeshDatabase.h"
#include "IO/IOHelpers.h"
#include "IO/silo.h"
#include "common/MPI_Helpers.h"
#include "common/Utilities.h"
#include "shared_ptr.h"

#include <sys/stat.h>
#include <algorithm>
#include <vector>
#include <set>



enum class Format { OLD, NEW, SILO, UNKNOWN }; 



/****************************************************
* Initialize the writer                             *
****************************************************/
static std::string global_IO_path;
static Format global_IO_format = Format::UNKNOWN;
void IO::initialize( const std::string& path, const std::string& format, bool append )
{
    if ( path.empty() )
        global_IO_path = ".";
    else
        global_IO_path = path;
    if ( format == "old" )
        global_IO_format = Format::OLD;
    else if ( format == "new" )
        global_IO_format = Format::NEW;
    else if ( format == "silo" )
        global_IO_format = Format::SILO;
    else
        ERROR("Unknown format");
    int rank = comm_rank(MPI_COMM_WORLD);
    if ( !append && rank==0 ) {
        mkdir(path.c_str(),S_IRWXU|S_IRGRP);
        std::string filename;
        if ( global_IO_format==Format::OLD || global_IO_format==Format::NEW )
            filename = global_IO_path + "/summary.LBM";
        else if ( global_IO_format==Format::SILO )
            filename = global_IO_path + "/LBM.visit";
        else
            ERROR("Unknown format");
        auto fid = fopen(filename.c_str(),"wb");
        fclose(fid);
    }
}


// Write the mesh data in the original format
static std::vector<IO::MeshDatabase> writeMeshesOrigFormat( const std::vector<IO::MeshDataStruct>& meshData, const std::string& path )
{
    int rank = MPI_WORLD_RANK();
    std::vector<IO::MeshDatabase> meshes_written;
    for (size_t i=0; i<meshData.size(); i++) {
        char domainname[100], filename[100], fullpath[200];
        sprintf(domainname,"%05i",rank);
        sprintf(filename,"%s.%05i",meshData[i].meshName.c_str(),rank);
	    sprintf(fullpath,"%s/%s",path.c_str(),filename);
        FILE *fid = fopen(fullpath,"wb");
        INSIST(fid!=NULL,std::string("Error opening file: ")+fullpath);
        std::shared_ptr<IO::Mesh> mesh = meshData[i].mesh;
        IO::MeshDatabase mesh_entry;
        mesh_entry.name = meshData[i].meshName;
        mesh_entry.type = meshType(*mesh);
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
        const std::string meshClass = mesh->className();
        if ( meshClass=="PointList" ) {
            // List of points
            std::shared_ptr<IO::PointList> pointlist = std::dynamic_pointer_cast<IO::PointList>(mesh);
            const std::vector<Point>& P = pointlist->points;
            for (size_t i=0; i<P.size(); i++) {
                double x[3];
                x[0] = P[i].x;  x[1] = P[i].y;  x[2] = P[i].z;
                fwrite(x,sizeof(double),3,fid);
            }
        } else if ( meshClass=="TriList" || meshClass=="TriMesh" ) {
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
        } else if ( meshClass=="DomainMesh" ) {
            // This format was never supported with the old format
        } else {
            ERROR("Unknown mesh");
        }
        fclose(fid);
        std::sort( mesh_entry.variables.begin(), mesh_entry.variables.end() );
        mesh_entry.variables.erase( std::unique( mesh_entry.variables.begin(), mesh_entry.variables.end() ), mesh_entry.variables.end() );
        meshes_written.push_back(mesh_entry);
    }
    return meshes_written;
}


// Create the database entry for the mesh data
static IO::MeshDatabase getDatabase( const std::string& filename, const IO::MeshDataStruct& mesh, int format )
{
    int rank = MPI_WORLD_RANK();
    char domainname[100];
    sprintf(domainname,"%s_%05i",mesh.meshName.c_str(),rank);
    // Create the MeshDatabase
    IO::MeshDatabase database;
    database.name = mesh.meshName;
    database.type = meshType(*(mesh.mesh));
    database.meshClass = mesh.mesh->className();
    database.format = format;
    // Write the mesh
    IO::DatabaseEntry domain;
    domain.name = domainname;
    domain.file = filename;
    domain.offset = -1;
    database.domains.push_back(domain);
    // Write the variables
    for (size_t i=0; i<mesh.vars.size(); i++) {
        // Add basic variable info
        IO::VariableDatabase info;
        info.name = mesh.vars[i]->name;
        info.type = mesh.vars[i]->type;
        info.dim = mesh.vars[i]->dim;
        database.variables.push_back(info);
        // Add domain variable info
        IO::DatabaseEntry variable;
        variable.name = mesh.vars[i]->name;
        variable.file = filename;
        variable.offset = -1;
        std::pair<std::string,std::string> key(domain.name,mesh.vars[i]->name);
        database.variable_data.insert( 
            std::pair<std::pair<std::string,std::string>,IO::DatabaseEntry>(key,variable) );
    }
    return database;
}


// Write a mesh (and variables) to a file
static IO::MeshDatabase write_domain( FILE *fid, const std::string& filename,
    const IO::MeshDataStruct& mesh, int format )
{
    const int level = 0;
    int rank = MPI_WORLD_RANK();
    // Create the MeshDatabase
    IO::MeshDatabase database = getDatabase( filename, mesh, format );
    // Write the mesh
    IO::DatabaseEntry& domain = database.domains[0];
    domain.offset = ftell(fid);
    std::pair<size_t,void*> data = mesh.mesh->pack(level);
    fprintf(fid,"Mesh: %s-%05i: %lu\n",mesh.meshName.c_str(),rank,data.first);
    fwrite(data.second,1,data.first,fid);
    fprintf(fid,"\n");
    delete [] (char*) data.second;
    // Write the variables
    for (size_t i=0; i<mesh.vars.size(); i++) {
        std::pair<std::string,std::string> key(domain.name,mesh.vars[i]->name);
        IO::DatabaseEntry& variable = database.variable_data[key];
        variable.offset = ftell(fid);
        int dim = mesh.vars[i]->dim;
        int type = static_cast<int>(mesh.vars[i]->type);
        size_t N = mesh.vars[i]->data.length();
        if ( type == static_cast<int>(IO::VariableType::NullVariable) ) {
            ERROR("Variable type not set");
        }
        size_t N_mesh = mesh.mesh->numberPointsVar(mesh.vars[i]->type);
        ASSERT(N==dim*N_mesh);
        fprintf(fid,"Var: %s-%05i-%s: %i, %i, %lu, %lu, double\n",
            database.name.c_str(), rank, variable.name.c_str(),
            dim, type, N_mesh, N*sizeof(double) );
        fwrite(mesh.vars[i]->data.data(),sizeof(double),N,fid);
        fprintf(fid,"\n");
    }
    return database;
}


#ifdef USE_SILO
// Write a PointList mesh (and variables) to a file
template<class TYPE>
static void writeSiloPointMesh( DBfile *fid, const IO::PointList& mesh, const std::string& meshname )
{
    const auto& points = mesh.getPoints();
    std::vector<TYPE> x(points.size()), y(points.size()), z(points.size());
    for (size_t i=0; i<x.size(); i++) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        z[i] = points[i].z;
    }
    const TYPE *coords[] = { x.data(), y.data(), z.data() };
    silo::writePointMesh<TYPE>( fid, meshname, 3, points.size(), coords );
}
static void writeSiloPointList( DBfile *fid, const IO::MeshDataStruct& meshData, IO::MeshDatabase database )
{
    const IO::PointList& mesh = dynamic_cast<IO::PointList&>( *meshData.mesh );
    const std::string meshname = database.domains[0].name;
    if ( meshData.precision == IO::DataType::Double ) {
        writeSiloPointMesh<double>( fid, mesh, meshname );
    } else if ( meshData.precision == IO::DataType::Float ) {
        writeSiloPointMesh<float>( fid, mesh, meshname );
    } else {
        ERROR("Unsupported format");
    }
    const auto& points = mesh.getPoints();
    std::vector<double> x(points.size()), y(points.size()), z(points.size());
    for (size_t i=0; i<x.size(); i++) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        z[i] = points[i].z;
    }
    const double *coords[] = { x.data(), y.data(), z.data() };
    silo::writePointMesh( fid, meshname, 3, points.size(), coords );
    for (size_t i=0; i<meshData.vars.size(); i++) {
        const IO::Variable& var = *meshData.vars[i];
        if ( var.precision == IO::DataType::Double ) {
            silo::writePointMeshVariable( fid, meshname, var.name, var.data );
        } else if ( var.precision == IO::DataType::Float ) {
            Array<float> data2( var.data.size() );
            data2.copy( var.data );
            silo::writePointMeshVariable( fid, meshname, var.name, data2 );
        } else if ( var.precision == IO::DataType::Int ) {
            Array<int> data2( var.data.size() );
            data2.copy( var.data );
            silo::writePointMeshVariable( fid, meshname, var.name, data2 );
        } else {
            ERROR("Unsupported format");
        }
    }
}
// Write a TriMesh mesh (and variables) to a file
template<class TYPE>
static void writeSiloTriMesh( DBfile *fid, const IO::TriMesh& mesh, const std::string& meshname )
{
    const auto& points = mesh.vertices->getPoints();
    std::vector<TYPE> x(points.size()), y(points.size()), z(points.size());
    for (size_t i=0; i<x.size(); i++) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        z[i] = points[i].z;
    }
    const TYPE *coords[] = { x.data(), y.data(), z.data() };
    const int *tri[] = { mesh.A.data(), mesh.B.data(), mesh.C.data() };
    silo::writeTriMesh<TYPE>( fid, meshname, 3, 2, points.size(), coords, mesh.A.size(), tri );
}
static void writeSiloTriMesh2( DBfile *fid, const IO::MeshDataStruct& meshData,
    const IO::TriMesh& mesh, IO::MeshDatabase database )
{
    const std::string meshname = database.domains[0].name;
    if ( meshData.precision == IO::DataType::Double ) {
        writeSiloTriMesh<double>( fid, mesh, meshname );
    } else if ( meshData.precision == IO::DataType::Float ) {
        writeSiloTriMesh<float>( fid, mesh, meshname );
    } else {
        ERROR("Unsupported format");
    }
    for (size_t i=0; i<meshData.vars.size(); i++) {
        const IO::Variable& var = *meshData.vars[i];
        auto type = static_cast<silo::VariableType>( var.type );
        if ( var.precision == IO::DataType::Double ) {
            silo::writeTriMeshVariable( fid, 3, meshname, var.name, var.data, type );
        } else if ( var.precision == IO::DataType::Float ) {
            Array<float> data2( var.data.size() );
            data2.copy( var.data );
            silo::writeTriMeshVariable( fid, 3, meshname, var.name, data2, type );
        } else if ( var.precision == IO::DataType::Int ) {
            Array<int> data2( var.data.size() );
            data2.copy( var.data );
            silo::writeTriMeshVariable( fid, 3, meshname, var.name, data2, type );
        } else {
            ERROR("Unsupported format");
        }
    }
}
static void writeSiloTriMesh( DBfile *fid, const IO::MeshDataStruct& meshData, IO::MeshDatabase database )
{
    const IO::TriMesh& mesh = dynamic_cast<IO::TriMesh&>( *meshData.mesh );
    writeSiloTriMesh2( fid, meshData, mesh, database );
}
static void writeSiloTriList( DBfile *fid, const IO::MeshDataStruct& meshData, IO::MeshDatabase database )
{
    auto mesh = getTriMesh( meshData.mesh );
    writeSiloTriMesh2( fid, meshData, *mesh, database );
}
// Write a DomainMesh mesh (and variables) to a file
static void writeSiloDomainMesh( DBfile *fid, const IO::MeshDataStruct& meshData, IO::MeshDatabase database )
{
    const IO::DomainMesh& mesh = dynamic_cast<IO::DomainMesh&>( *meshData.mesh );
    RankInfoStruct info( mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz );
    std::array<double,6> range = { info.ix*mesh.Lx/info.nx, (info.ix+1)*mesh.Lx/info.nx, 
                                   info.jy*mesh.Ly/info.ny, (info.jy+1)*mesh.Ly/info.ny,
                                   info.kz*mesh.Lz/info.nz, (info.kz+1)*mesh.Lz/info.nz };
    std::array<int,3> N = { mesh.nx, mesh.ny, mesh.nz };
    auto meshname = database.domains[0].name;
    silo::writeUniformMesh<3>( fid, meshname, range, N );
    silo::write<int>( fid, meshname+"_rankinfo", { mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz } );
    for (size_t i=0; i<meshData.vars.size(); i++) {
        const auto& var = *meshData.vars[i];
        auto type = static_cast<silo::VariableType>( var.type );
        if ( var.precision == IO::DataType::Double ) {
            silo::writeUniformMeshVariable<3>( fid, meshname, N, var.name, var.data, type );
        } else if ( var.precision == IO::DataType::Float ) {
            Array<float> data2( var.data.size() );
            data2.copy( var.data );
            silo::writeUniformMeshVariable<3>( fid, meshname, N, var.name, data2, type );
        } else if ( var.precision == IO::DataType::Int ) {
            Array<int> data2( var.data.size() );
            data2.copy( var.data );
            silo::writeUniformMeshVariable<3>( fid, meshname, N, var.name, data2, type );
        } else {
            ERROR("Unsupported format");
        }
    }
}
// Write a mesh (and variables) to a file
static IO::MeshDatabase write_domain_silo( DBfile *fid, const std::string& filename,
    const IO::MeshDataStruct& mesh, int format )
{
    // Create the MeshDatabase
    auto database = getDatabase( filename, mesh, format );
    if ( database.meshClass=="PointList" ) {
        writeSiloPointList( fid, mesh, database );
    } else if ( database.meshClass=="TriMesh" ) {
        writeSiloTriMesh( fid, mesh, database );
    } else if ( database.meshClass=="TriList" ) {
        writeSiloTriList( fid, mesh, database );
    } else if ( database.meshClass=="DomainMesh" ) {
        writeSiloDomainMesh( fid, mesh, database );
    } else {
        ERROR("Unknown mesh class");
    }
    return database;
}
// Write the summary file for silo
std::pair<int,int> getSiloMeshType( const std::string& meshClass )
{
    int meshType = 0;
    int varType = 0;
    if ( meshClass=="PointList" ) {
        meshType = DB_POINTMESH;
        varType  = DB_POINTVAR;
    } else if ( meshClass=="TriMesh" ) {
        meshType = DB_UCDMESH;
        varType  = DB_UCDVAR;
    } else if ( meshClass=="TriList" ) {
        meshType = DB_UCDMESH;
        varType  = DB_UCDVAR;
    } else if ( meshClass=="DomainMesh" ) {
        meshType = DB_QUAD_RECT;
        varType  = DB_QUADVAR;
    } else {
        ERROR("Unknown mesh class");
    }
    return std::make_pair( meshType, varType );
}
void writeSiloSummary( const std::vector<IO::MeshDatabase>& meshes_written, const std::string& filename )
{
    auto fid = silo::open( filename, silo::CREATE );
    for ( const auto& data : meshes_written ) {
        auto type = getSiloMeshType( data.meshClass );
        std::vector<int> meshTypes( data.domains.size(), type.first );
        std::vector<int> varTypes( data.domains.size(), type.second );
        std::vector<std::string> meshNames;
        for ( const auto& tmp : data.domains )
            meshNames.push_back( tmp.file + ":" + tmp.name );
        silo::writeMultiMesh( fid, data.name, meshNames, meshTypes );
        for (const auto& variable : data.variables ) {
            std::vector<std::string> varnames;
            for ( const auto& tmp : data.domains )
                varnames.push_back( tmp.file + ":" + variable.name );
            silo::writeMultiVar( fid, variable.name, varnames, varTypes );
        }
    }
    silo::close( fid );
}
#endif


// Write the mesh data in the new format
static std::vector<IO::MeshDatabase> writeMeshesNewFormat( 
    const std::vector<IO::MeshDataStruct>& meshData, const std::string& path, int format )
{
    int rank = MPI_WORLD_RANK();
    std::vector<IO::MeshDatabase> meshes_written;
    char filename[100], fullpath[200];
    sprintf(filename,"%05i",rank);
    sprintf(fullpath,"%s/%s",path.c_str(),filename);
    FILE *fid = fopen(fullpath,"wb");
    for (size_t i=0; i<meshData.size(); i++) {
        std::shared_ptr<IO::Mesh> mesh = meshData[i].mesh;
        meshes_written.push_back( write_domain(fid,filename,meshData[i],format) );
    }
    fclose(fid);
    return meshes_written;
}


// Write the mesh data to silo
static std::vector<IO::MeshDatabase> writeMeshesSilo( 
    const std::vector<IO::MeshDataStruct>& meshData, const std::string& path, int format )
{
#ifdef USE_SILO
    int rank = MPI_WORLD_RANK();
    std::vector<IO::MeshDatabase> meshes_written;
    char filename[100], fullpath[200];
    sprintf(filename,"%05i.silo",rank);
    sprintf(fullpath,"%s/%s",path.c_str(),filename);
    auto fid = silo::open( fullpath, silo::CREATE );
    for (size_t i=0; i<meshData.size(); i++) {
        auto mesh = meshData[i].mesh;
        meshes_written.push_back( write_domain_silo(fid,filename,meshData[i],format) );
    }
    silo::close( fid );
    return meshes_written;
#else
    ERROR("Application built without silo support");
    return std::vector<IO::MeshDatabase>();
#endif
}        


/****************************************************
* Write the mesh data                               *
****************************************************/
void IO::writeData( const std::string& subdir, const std::vector<IO::MeshDataStruct>& meshData, MPI_Comm comm )
{
    if ( global_IO_path.empty() )
        IO::initialize( );
    PROFILE_START("writeData");
    int rank = comm_rank(comm);
    // Check the meshData before writing
    for ( const auto& data : meshData ) {
        if ( !data.check() )
            ERROR("Error in meshData");
    }
    // Create the output directory
    std::string path = global_IO_path + "/" + subdir;
    if ( rank == 0 ) {
        mkdir(path.c_str(),S_IRWXU|S_IRGRP);
    }
    MPI_Barrier(comm);
    // Write the mesh files
    std::vector<IO::MeshDatabase> meshes_written;
    if ( global_IO_format == Format::OLD ) {
        // Write the original triangle format
        meshes_written = writeMeshesOrigFormat( meshData, path );
    } else if ( global_IO_format == Format::NEW ) {
        // Write the new format (double precision)
        meshes_written = writeMeshesNewFormat( meshData, path, 2 );
    } else if ( global_IO_format == Format::SILO ) {
        // Write silo
        meshes_written = writeMeshesSilo( meshData, path, 4 );
    } else {
        ERROR("Unknown format");
    }
    // Gather a complete list of files on rank 0
    meshes_written = gatherAll(meshes_written,comm);
    // Write the summary files
    if ( rank == 0 ) {
        // Write the summary file for the current timestep
        char filename[200];
        sprintf(filename,"%s/LBM.summary",path.c_str());
        write(meshes_written,filename);
        // Write summary silo file if needed
        #ifdef USE_SILO
        if ( global_IO_format == Format::SILO ) {
            sprintf(filename,"%s/summary.silo",path.c_str());
            writeSiloSummary(meshes_written,filename);
        }
        #endif
        // Add the timestep to the global summary file
        if ( global_IO_format == Format::OLD || global_IO_format == Format::NEW ) {
            auto filename = global_IO_path+"/summary.LBM";
            FILE *fid = fopen(filename.c_str(),"ab");
            fprintf(fid,"%s/\n",subdir.c_str());
            fclose(fid);
        } else if ( global_IO_format == Format::SILO ) {
            auto filename = global_IO_path+"/LBM.visit";
            FILE *fid = fopen(filename.c_str(),"ab");
            fprintf(fid,"%s/summary.silo\n",subdir.c_str());
            fclose(fid);
        } else {
            ERROR("Unknown format");
        }
    }
    PROFILE_STOP("writeData");
}


