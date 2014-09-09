#include "IO/MeshDatabase.h"
#include "common/Utilities.h"

#include <vector>
#include <map>
#include <set>
#include <ProfilerApp.h>


namespace IO {


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


// MeshDatabase::MeshDatabase
MeshDatabase::MeshDatabase()
{
}
MeshDatabase::~MeshDatabase()
{
}
MeshDatabase::MeshDatabase(const MeshDatabase& rhs)
{
    name = rhs.name;
    type = rhs.type;
    format = rhs.format;
    domains = rhs.domains;
    file = rhs.file;
    variables = rhs.variables;
}
MeshDatabase& MeshDatabase::operator=(const MeshDatabase& rhs)
{
    this->name = rhs.name;
    this->type = rhs.type;
    this->format = rhs.format;
    this->domains = rhs.domains;
    this->file = rhs.file;
    this->variables = rhs.variables;
}


// Pack/Unpack the MeshDatabase to/from buffer
static size_t MeshDatabaseSize( const MeshDatabase& data )
{
    size_t bytes = 2;
    bytes += data.name.size()+1;
    bytes += sizeof(int);
    for (size_t i=0; i<data.domains.size(); i++)
        bytes += data.domains[i].size()+1;
    bytes += sizeof(int);
    for (size_t i=0; i<data.file.size(); i++)
        bytes += data.file[i].size()+1;
    bytes += sizeof(int);
    for (size_t i=0; i<data.variables.size(); i++)
        bytes += data.variables[i].size()+1;
    return bytes;
}
static void pack( const MeshDatabase& data, char *buffer )
{
    buffer[0] = static_cast<char>(data.type);
    buffer[1] = data.format;
    size_t pos = 2;
    memcpy(&buffer[pos],data.name.c_str(),data.name.size()+1);
    pos += data.name.size()+1;
    *reinterpret_cast<int*>(&buffer[pos]) = data.domains.size();
    pos += sizeof(int);
    for (size_t i=0; i<data.domains.size(); i++) {
        memcpy(&buffer[pos],data.domains[i].c_str(),data.domains[i].size()+1);
        pos += data.domains[i].size()+1;
    }
    *reinterpret_cast<int*>(&buffer[pos]) = data.file.size();
    pos += sizeof(int);
    for (size_t i=0; i<data.file.size(); i++) {
        memcpy(&buffer[pos],data.file[i].c_str(),data.file[i].size()+1);
        pos += data.file[i].size()+1;
    }
    *reinterpret_cast<int*>(&buffer[pos]) = data.variables.size();
    pos += sizeof(int);
    for (size_t i=0; i<data.variables.size(); i++) {
        memcpy(&buffer[pos],data.variables[i].c_str(),data.variables[i].size()+1);
        pos += data.variables[i].size()+1;
    }
}
static MeshDatabase unpack( const char *buffer )
{
    MeshDatabase data;
    data.type = static_cast<MeshType>(buffer[0]);
    data.format = buffer[1];
    size_t pos = 2;
    data.name = std::string(&buffer[pos]);
    pos += data.name.size()+1;
    int N_domains = *reinterpret_cast<const int*>(&buffer[pos]);
    pos += sizeof(int);
    data.domains.resize(N_domains);
    for (size_t i=0; i<data.domains.size(); i++) {
        data.domains[i] = std::string(&buffer[pos]);
        pos += data.domains[i].size()+1;
    }
    int N_file = *reinterpret_cast<const int*>(&buffer[pos]);
    pos += sizeof(int);
    data.file.resize(N_file);
    for (size_t i=0; i<data.file.size(); i++) {
        data.file[i] = std::string(&buffer[pos]);
        pos += data.file[i].size()+1;
    }
    int N_variables = *reinterpret_cast<const int*>(&buffer[pos]);
    pos += sizeof(int);
    data.variables.resize(N_variables);
    for (size_t i=0; i<data.variables.size(); i++) {
        data.variables[i] = std::string(&buffer[pos]);
        pos += data.variables[i].size()+1;
    }
    return data;
}


// Gather the mesh databases from all processors
std::vector<MeshDatabase> gatherAll( const std::vector<MeshDatabase>& meshes, MPI_Comm comm )
{
    PROFILE_START("gatherAll");
    int rank = -1;
    int size = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    // First pack the mesh data to local buffers
    size_t localsize = 0;
    for (size_t i=0; i<meshes.size(); i++)
        localsize += MeshDatabaseSize(meshes[i]);
    char *localbuf = new char[localsize];
    size_t pos = 0;
    for (size_t i=0; i<meshes.size(); i++) {
        pack( meshes[i], &localbuf[pos] );
        pos += MeshDatabaseSize(meshes[i]);
    }
    // Get the number of bytes each processor will be sending/recieving
    int sendsize = static_cast<int>(localsize);
    int *recvsize = new int[size];
    MPI_Allgather(&sendsize,1,MPI_INT,recvsize,1,MPI_INT,comm);
    size_t globalsize = recvsize[0];
    int *disp = new int[size];
    disp[0] = 0;    
    for (size_t i=1; i<size; i++) {
        disp[i] = disp[i-1] + recvsize[i];
        globalsize += recvsize[i];
    }
    // Send/recv the global data
    char *globalbuf = new char[globalsize];
    MPI_Allgatherv(localbuf,sendsize,MPI_CHAR,globalbuf,recvsize,disp,MPI_CHAR,comm);
    // Unpack the data
    std::map<std::string,MeshDatabase> data;
    pos = 0;
    while ( pos < globalsize ) {
        MeshDatabase tmp = unpack(&globalbuf[pos]);
        pos += MeshDatabaseSize(tmp);
        std::map<std::string,MeshDatabase>::iterator it = data.find(tmp.name);
        if ( it==data.end() ) {
            data[tmp.name] = tmp;
        } else {
            for (size_t i=0; i<tmp.domains.size(); i++)
                it->second.domains.push_back(tmp.domains[i]);
            for (size_t i=0; i<tmp.domains.size(); i++)
                it->second.file.push_back(tmp.file[i]);
            for (size_t i=0; i<tmp.variables.size(); i++)
                it->second.variables.push_back(tmp.variables[i]);
        }
    }
    for (std::map<std::string,MeshDatabase>::iterator it=data.begin(); it!=data.end(); ++it) {
        // Get the unique variables
        std::set<std::string> data2(it->second.variables.begin(),it->second.variables.end());
        it->second.variables = std::vector<std::string>(data2.begin(),data2.end());
    }
    // Free temporary memory
    delete [] localbuf;
    delete [] recvsize;
    delete [] disp;
    delete [] globalbuf;
    // Return the results
    std::vector<MeshDatabase> data2(data.size());
    size_t i=0; 
    for (std::map<std::string,MeshDatabase>::iterator it=data.begin(); it!=data.end(); ++it, ++i)
        data2[i] = it->second;
    PROFILE_STOP("gatherAll");
    return data2;
}


//! Write the mesh databases to a file
void write( const std::vector<MeshDatabase>& meshes, const std::string& filename )
{
    PROFILE_START("write");
    FILE *fid = fopen(filename.c_str(),"wb");
    for (size_t i=0; i<meshes.size(); i++) {
        fprintf(fid,"%s\n",meshes[i].name.c_str());
        fprintf(fid,"   format: %i\n",static_cast<int>(meshes[i].format));
        fprintf(fid,"   type: %i\n",static_cast<int>(meshes[i].type));
        fprintf(fid,"   domains: ");
        for (size_t j=0; j<meshes[i].domains.size(); j++) {
            fprintf(fid,"%s; ",meshes[i].domains[j].c_str());
        }
        fprintf(fid,"\n");
        fprintf(fid,"   file: ");
        for (size_t j=0; j<meshes[i].file.size(); j++) {
            fprintf(fid,"%s; ",meshes[i].file[j].c_str());
        }
        fprintf(fid,"\n");
        fprintf(fid,"   variables: ");
        for (size_t j=0; j<meshes[i].variables.size(); j++) {
            fprintf(fid,"%s; ",meshes[i].variables[j].c_str());
        }
        fprintf(fid,"\n");
    }
    fclose(fid);
    PROFILE_STOP("write");
}


// Helper function to split a line into the sub variables
static inline std::vector<std::string> splitList( const char *line )
{
    std::vector<std::string> list;
    size_t i1 = 0;
    size_t i2 = 0;
    while ( 1 ) {
        if ( line[i2]==';' || line[i2]<32 ) {
            std::string tmp(&line[i1],i2-i1);
            tmp = deblank(tmp);
            if ( !tmp.empty() )
                list.push_back(tmp);
            i1 = i2+1;
        }
        if ( line[i2]==0 )
            break;
        i2++;
    }
    return list;
}


//! Read the mesh databases from a file
std::vector<MeshDatabase> read( const std::string& filename )
{
    std::vector<MeshDatabase> meshes;
    PROFILE_START("read");
    FILE *fid = fopen(filename.c_str(),"rb");
    if ( fid==NULL )
        ERROR("Error opening file");
    char *line = new char[0x100000];
    while ( std::fgets(line,0x100000,fid) != NULL ) {
        if ( line[0]<32 ) {
            // Empty line
            continue;
        } else if ( line[0] != ' ' ) {
            meshes.resize(meshes.size()+1);
            std::string name(line);
            name.resize(name.size()-1);
            meshes.back().name = name;
        } else if ( strncmp(line,"   format:",10)==0 ) {
            meshes.back().format = static_cast<unsigned char>(atoi(&line[10]));
        } else if ( strncmp(line,"   type:",8)==0 ) {
            meshes.back().type = static_cast<MeshType>(atoi(&line[8]));
        } else if ( strncmp(line,"   domains:",11)==0 ) {
            meshes.back().domains = splitList(&line[11]);
        } else if ( strncmp(line,"   file:",8)==0 ) {
            meshes.back().file = splitList(&line[8]);
        } else if ( strncmp(line,"   variables:",13)==0 ) {
            meshes.back().variables = splitList(&line[13]);
        } else {
            ERROR("Error reading line");
        }
    }
    fclose(fid);
    delete [] line;
    PROFILE_STOP("read");
    return meshes;
}



} // IO namespace

