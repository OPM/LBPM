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
#include "IO/Reader.h"
#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"
#include "IO/IOHelpers.h"
#include "common/Utilities.h"

#ifdef USE_SILO
#include "IO/silo.h"
#endif


#include <ProfilerApp.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <cstdio>


// Inline function to read line without a return argument
static inline void fgetl( char * str, int num, FILE * stream )
{
    char* ptr = fgets( str, num, stream );
    if ( 0 ) {char *temp = (char *)&ptr; temp++;}
}


// Get the path to a file
std::string IO::getPath( const std::string& filename )
{
    std::string file(filename);
    size_t k1 = file.rfind(47);
    size_t k2 = file.rfind(92);
    if ( k1==std::string::npos ) { k1=0; }
    if ( k2==std::string::npos ) { k2=0; }
    return file.substr(0,std::max(k1,k2));
}


// List the timesteps in the given directors (dumps.LBPM)
std::vector<std::string> IO::readTimesteps( const std::string& filename )
{
    PROFILE_START("readTimesteps");
    FILE *fid= fopen(filename.c_str(),"rb");
    if ( fid==NULL )
        ERROR("Error opening file");
    std::vector<std::string> timesteps;
    char buf[1000];
    while (fgets(buf,sizeof(buf),fid) != NULL) {
        std::string line(buf);
        line.resize(line.size()-1);
        auto pos = line.find( "summary.silo" );
        if ( pos != std::string::npos )
            line.resize(pos);
        if ( line.empty() )
            continue;
        timesteps.push_back(line);
    }
    fclose(fid);
    PROFILE_STOP("readTimesteps");
    return timesteps;
}


// Read the list of variables for the given timestep
std::vector<IO::MeshDatabase> IO::getMeshList( const std::string& path, const std::string& timestep )
{
    std::string filename = path + "/" + timestep + "/LBM.summary";
    return IO::read( filename );
}


// Read the given mesh domain
std::shared_ptr<IO::Mesh> IO::getMesh( const std::string& path, const std::string& timestep, 
    const IO::MeshDatabase& meshDatabase, int domain )
{
    PROFILE_START("getMesh");
    std::shared_ptr<IO::Mesh> mesh;
    if ( meshDatabase.format==1 ) {
        // Old format (binary doubles)
        std::string filename = path + "/" + timestep + "/" + meshDatabase.domains[domain].file;
        FILE *fid = fopen(filename.c_str(),"rb");
        INSIST(fid!=NULL,"Error opening file");
        fseek( fid, 0, SEEK_END );
        size_t bytes = ftell(fid);
        size_t N_max = bytes/sizeof(double)+1000;
        double *data = new double[N_max];
        fseek(fid,0,SEEK_SET);
        size_t count = fread(data,sizeof(double),N_max,fid);
        fclose(fid);
        if ( count%3 != 0 )
            ERROR("Error reading file");
        if ( meshDatabase.type==IO::PointMesh ) {
            size_t N = count/3;
            std::shared_ptr<PointList> pointlist( new PointList(N) );
            std::vector<Point>& P = pointlist->points;
            for (size_t i=0; i<N; i++) {
                P[i].x = data[3*i+0];
                P[i].y = data[3*i+1];
                P[i].z = data[3*i+2];
            }
            mesh = pointlist;
        } else if ( meshDatabase.type==IO::SurfaceMesh ) {
            if ( count%9 != 0 )
                ERROR("Error reading file (2)");
            size_t N_tri = count/9;
            std::shared_ptr<TriList> trilist( new TriList(N_tri) );
            std::vector<Point>& A = trilist->A;
            std::vector<Point>& B = trilist->B;
            std::vector<Point>& C = trilist->C;
            for (size_t i=0; i<N_tri; i++) {
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
            mesh = trilist;
        } else if ( meshDatabase.type==IO::VolumeMesh ) {
            // this was never supported in the old format
            mesh = std::shared_ptr<DomainMesh>( new DomainMesh() );
        } else {
            ERROR("Unknown mesh type");
        }
        delete [] data;
    } else if ( meshDatabase.format==2 ) {
        const DatabaseEntry& database = meshDatabase.domains[domain];
        std::string filename = path + "/" + timestep + "/" + database.file;
        FILE *fid = fopen(filename.c_str(),"rb");
        fseek(fid,database.offset,SEEK_SET);
        char line[1000];
        fgetl(line,1000,fid);
        size_t i1 = find(line,':');
        size_t i2 = find(&line[i1+1],':')+i1+1;
        size_t bytes = atol(&line[i2+1]);
        char *data = new char[bytes];
        size_t count = fread(data,1,bytes,fid);
        fclose(fid);
        ASSERT(count==bytes);
        if ( meshDatabase.meshClass=="PointList" ) {
            mesh.reset( new IO::PointList() );
        } else if ( meshDatabase.meshClass=="TriMesh" ) {
            mesh.reset( new IO::TriMesh() );
        } else if ( meshDatabase.meshClass=="TriList" ) {
            mesh.reset( new IO::TriList() );
        } else if ( meshDatabase.meshClass=="DomainMesh" ) {
            mesh.reset( new IO::DomainMesh() );
        } else {
            ERROR("Unknown mesh class");
        }
        mesh->unpack( std::pair<size_t,void*>(bytes,data) );
        delete [] data;
    } else if ( meshDatabase.format==4 ) {
        // Reading a silo file
#ifdef USE_SILO
        const DatabaseEntry& database = meshDatabase.domains[domain];
        std::string filename = path + "/" + timestep + "/" + database.file;
        auto fid = silo::open( filename, silo::READ );
        if ( meshDatabase.meshClass=="PointList" ) {
            Array<double> coords = silo::readPointMesh<double>( fid, database.name );
            ASSERT(coords.size(1)==3);
            std::shared_ptr<IO::PointList> mesh2( new IO::PointList( coords.size(0) ) );
            for (size_t i=0; i<coords.size(1); i++) {
                mesh2->points[i].x = coords(i,0);
                mesh2->points[i].y = coords(i,1);
                mesh2->points[i].z = coords(i,2);
            }
            mesh = mesh2;
        } else if ( meshDatabase.meshClass=="TriMesh" || meshDatabase.meshClass=="TriList" ) {
            Array<double> coords;
            Array<int> tri;
            silo::readTriMesh( fid, database.name, coords, tri );
            ASSERT( tri.size(1)==3 && coords.size(1)==3 );
            int N_tri = tri.size(0);
            int N_point = coords.size(0);
            std::shared_ptr<IO::TriMesh> mesh2( new IO::TriMesh( N_tri, N_point ) );
            for (int i=0; i<N_point; i++) {
                mesh2->vertices->points[i].x = coords(i,0);
                mesh2->vertices->points[i].y = coords(i,1);
                mesh2->vertices->points[i].z = coords(i,2);
            }
            for (int i=0; i<N_tri; i++) {
                mesh2->A[i] = tri(i,0);
                mesh2->B[i] = tri(i,1);
                mesh2->C[i] = tri(i,2);
            }
            if ( meshDatabase.meshClass=="TriMesh" ) {
                mesh = mesh2;
            } else if ( meshDatabase.meshClass=="TriList" ) {
                auto trilist = IO::getTriList( std::dynamic_pointer_cast<IO::Mesh>( mesh2 ) );
                mesh = trilist;
            }
        } else if ( meshDatabase.meshClass=="DomainMesh" ) {
            std::vector<double> range;
            std::vector<int> N;
            silo::readUniformMesh( fid, database.name, range, N );
            auto rankinfo = silo::read<int>( fid, database.name+"_rankinfo" );
            RankInfoStruct rank_data( rankinfo[0], rankinfo[1], rankinfo[2], rankinfo[3] );
            mesh.reset( new IO::DomainMesh( rank_data, N[0], N[1], N[2], range[1]-range[0], range[3]-range[2], range[5]-range[4] ) );
        } else {
            ERROR("Unknown mesh class");
        }
        silo::close( fid );
#else
        ERROR("Build without silo support");
#endif
    } else {
        ERROR("Unknown format");
    }
    PROFILE_STOP("getMesh");
    return mesh;
}


// Read the given variable for the given mesh domain
std::shared_ptr<IO::Variable> IO::getVariable( const std::string& path, const std::string& timestep, 
    const MeshDatabase& meshDatabase, int domain, const std::string& variable )
{
    std::pair<std::string,std::string> key(meshDatabase.domains[domain].name,variable);
    std::map<std::pair<std::string,std::string>,DatabaseEntry>::const_iterator it;
    it = meshDatabase.variable_data.find(key);
    if ( it==meshDatabase.variable_data.end() )
        return std::shared_ptr<IO::Variable>();
    std::shared_ptr<IO::Variable> var;
    if ( meshDatabase.format == 2 ) {
        const DatabaseEntry& database = it->second;
        std::string filename = path + "/" + timestep + "/" + database.file;
        FILE *fid = fopen(filename.c_str(),"rb");
        fseek(fid,database.offset,SEEK_SET);
        char line[1000];
        fgetl(line,1000,fid);
        size_t i1 = find(line,':');
        size_t i2 = find(&line[i1+1],':')+i1+1;
        std::vector<std::string> values = splitList(&line[i2+1],',');
        ASSERT(values.size()==5);
        int dim = atoi(values[0].c_str());
        int type = atoi(values[1].c_str());
        size_t N = atol(values[2].c_str());
        size_t bytes = atol(values[3].c_str());
        std::string precision  = values[4];
        var = std::shared_ptr<IO::Variable>( new IO::Variable() );
        var->dim = dim;
        var->type = static_cast<IO::VariableType>(type);
        var->name = variable;
        var->data.resize(N*dim);
        if ( precision=="double" ) {
            size_t count = fread(var->data.data(),sizeof(double),N*dim,fid);
            ASSERT(count*sizeof(double)==bytes);
        } else {
            ERROR("Format not implimented");
        }
        fclose(fid);
    } else if ( meshDatabase.format == 4 ) {
        // Reading a silo file
#ifdef USE_SILO
        const auto& database = meshDatabase.domains[domain];
        auto variableDatabase = meshDatabase.getVariableDatabase( variable );
        std::string filename = path + "/" + timestep + "/" + database.file;
        auto fid = silo::open( filename, silo::READ );
        var.reset( new Variable( variableDatabase.dim, variableDatabase.type, variable ) );
        if ( meshDatabase.meshClass=="PointList" ) {
            var->data = silo::readPointMeshVariable<double>( fid, variable );
        } else if ( meshDatabase.meshClass=="TriMesh" || meshDatabase.meshClass=="TriList" ) {
            var->data = silo::readTriMeshVariable<double>( fid, variable );
        } else if ( meshDatabase.meshClass=="DomainMesh" ) {
            var->data = silo::readUniformMeshVariable<double>( fid, variable );
        } else {
            ERROR("Unknown mesh class");
        }
        silo::close( fid );
#else
        ERROR("Build without silo support");
#endif

    } else {
        ERROR("Unknown format");
    }
    return var;
}


/****************************************************
* Reformat the variable to match the mesh           *
****************************************************/
void IO::reformatVariable( const IO::Mesh& mesh, IO::Variable& var )
{
    if ( mesh.className() == "DomainMesh" ) {
        const IO::DomainMesh& mesh2 = dynamic_cast<const IO::DomainMesh&>( mesh );
        if ( var.type == VariableType::NodeVariable ) {
            size_t N2 = var.data.length() / ((mesh2.nx+1)*(mesh2.ny+1)*(mesh2.nz+1));
            ASSERT( (mesh2.nx+1)*(mesh2.ny+1)*(mesh2.nz+1)*N2 == var.data.length() );
            var.data.reshape( { (size_t) mesh2.nx+1, (size_t) mesh2.ny+1, (size_t) mesh2.nz+1, N2 } );
        } else if ( var.type == VariableType::EdgeVariable ) {
            ERROR("Not finished");
        } else if ( var.type == VariableType::SurfaceVariable ) {
            ERROR("Not finished");
        } else if ( var.type == VariableType::VolumeVariable ) {
            size_t N2 = var.data.length() / (mesh2.nx*mesh2.ny*mesh2.nz);
            ASSERT( mesh2.nx*mesh2.ny*mesh2.nz*N2 == var.data.length() );
            var.data.reshape( { (size_t) mesh2.nx, (size_t) mesh2.ny, (size_t) mesh2.nz, N2 } );
        } else {
            ERROR("Invalid variable type");
        }
    } else if ( mesh.className() == "PointList" ) {
        const IO::PointList& mesh2 = dynamic_cast<const IO::PointList&>( mesh );
        size_t N = mesh2.points.size();
        size_t N_var = var.data.length()/N;
        ASSERT( N*N_var == var.data.length() );
        var.data.reshape( { N, N_var } );
    } else if ( mesh.className()=="TriMesh" || mesh.className() == "TriList" ) {
        std::shared_ptr<Mesh> mesh_ptr( const_cast<Mesh*>(&mesh), []( void* ) {} );
        std::shared_ptr<TriMesh> mesh2 = getTriMesh( mesh_ptr );
        if ( var.type == VariableType::NodeVariable ) {
            size_t N = mesh2->vertices->points.size();
            size_t N_var = var.data.length()/N;
            ASSERT( N*N_var == var.data.length() );
            var.data.reshape( { N, N_var } );
        } else if ( var.type == VariableType::EdgeVariable ) {
            ERROR("Not finished");
        } else if ( var.type == VariableType::SurfaceVariable ) {
            ERROR("Not finished");
        } else if ( var.type == VariableType::VolumeVariable ) {
            size_t N = mesh2->A.size();
            size_t N_var = var.data.length()/N;
            ASSERT( N*N_var == var.data.length() );
            var.data.reshape( { N, N_var } );
        } else {
            ERROR("Invalid variable type");
        }
    } else {
        ERROR("Unknown mesh type");
    }
}



