#include "IO/Reader.h"
#include "IO/Mesh.h"
#include "IO/MeshDatabase.h"
#include "common/Utilities.h"

#include <ProfilerApp.h>
#include <iostream>
#include <string.h>
#include <memory>
#include <vector>


// Helper function
static inline size_t find( const char *line, char key )
{
    size_t i=0;
    while ( 1 ) {
        if ( line[i]==key || line[i]<32 || line[i]==0 )
            break;
        ++i;
    }
    return i;
}


// Remove preceeding/trailing whitespace
inline std::string deblank( const std::string& str )
{
    size_t i1 = str.size();
    size_t i2 = 0;
    for (size_t i=0; i<str.size(); i++) {
        if ( str[i]!=' ' && str[i]>=32 ) {
            i1 = std::min(i1,i);
            i2 = std::max(i2,i);
        }
    }
    return str.substr(i1,i2-i1+1);
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
        if ( meshDatabase.type==IO::MeshType::PointMesh ) {
            size_t N = count/3;
            std::shared_ptr<PointList> pointlist( new PointList(N) );
            std::vector<Point>& P = pointlist->points;
            for (size_t i=0; i<N; i++) {
                P[i].x = data[3*i+0];
                P[i].y = data[3*i+1];
                P[i].z = data[3*i+2];
            }
            mesh = pointlist;
        } else if ( meshDatabase.type==IO::MeshType::SurfaceMesh ) {
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
        std::fgets(line,0x100000,fid);
        size_t i1 = find(line,':');
        size_t i2 = find(&line[i1+1],':')+i1+1;
        size_t bytes = atol(&line[i2+1]);
        char *data = new char[bytes];
        size_t count = fread(data,1,bytes,fid);
        fclose(fid);
        if ( meshDatabase.meshClass=="PointList" ) {
            mesh.reset( new IO::PointList() );
        } else if ( meshDatabase.meshClass=="TriMesh" ) {
            mesh.reset( new IO::TriMesh() );
        } else if ( meshDatabase.meshClass=="TriList" ) {
            mesh.reset( new IO::TriList() );
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
    return std::shared_ptr<IO::Variable>();
}



