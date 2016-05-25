#include "IO/Reader.h"
#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"
#include "IO/IOHelpers.h"
#include "common/Utilities.h"

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
    char *data = new char[bytes];
    size_t count = fread(data,1,bytes,fid);
    fclose(fid);
    ASSERT(count==bytes);
    std::shared_ptr<IO::Variable> var( new IO::Variable() );
    var->dim = dim;
    var->type = static_cast<IO::VariableType>(type);
    var->name = variable;
    var->data.resize(N);
    if ( precision=="double" ) {
        memcpy(var->data.data(),data,bytes);
    } else {
        ERROR("Format not implimented");
    }
    delete [] data;
    return var;
}



