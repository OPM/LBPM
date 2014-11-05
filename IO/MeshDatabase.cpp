#include "IO/MeshDatabase.h"
#include "IO/Mesh.h"
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



/****************************************************
* Pack/unpack data from a buffer                    *
****************************************************/
template<class TYPE>
size_t packsize( const TYPE& rhs );
template<class TYPE>
void pack( const TYPE& rhs, char *buffer );
template<class TYPE>
void unpack( TYPE& data, const char *buffer );
// std::vector
template<class TYPE>
size_t packsize( const std::vector<TYPE>& rhs )
{
    size_t bytes = sizeof(size_t);
    for (size_t i=0; i<rhs.size(); i++)
        bytes += packsize(rhs[i]);
    return bytes;
}
template<class TYPE>
void pack( const std::vector<TYPE>& rhs, char *buffer )
{
    size_t size = rhs.size();
    memcpy(buffer,&size,sizeof(size_t));
    size_t pos = sizeof(size_t);
    for (int i=0; i<rhs.size(); i++) {
        pack(rhs[i],&buffer[pos]);
        pos += packsize(rhs[i]);
    }
}
template<class TYPE>
void unpack( std::vector<TYPE>& data, const char *buffer )
{
    size_t size;
    memcpy(&size,buffer,sizeof(size_t));
    data.clear();
    data.resize(size);
    size_t pos = sizeof(size_t);
    for (int i=0; i<data.size(); i++) {
        unpack(data[i],&buffer[pos]);
        pos += packsize(data[i]);
    }
}
// std::pair
template<class TYPE1, class TYPE2>
size_t packsize( const std::pair<TYPE1,TYPE2>& rhs )
{
    return packsize(rhs.first)+packsize(rhs.second);
}
template<class TYPE1, class TYPE2>
void pack( const std::pair<TYPE1,TYPE2>& rhs, char *buffer )
{
    pack(rhs.first,buffer);
    pack(rhs.second,&buffer[packsize(rhs.first)]);
}
template<class TYPE1, class TYPE2>
void unpack( std::pair<TYPE1,TYPE2>& data, const char *buffer )
{
    unpack(data.first,buffer);
    unpack(data.second,&buffer[packsize(data.first)]);
}
// std::map
template<class TYPE1, class TYPE2>
size_t packsize( const std::map<TYPE1,TYPE2>& rhs )
{
    size_t bytes = sizeof(size_t);
    typename std::map<TYPE1,TYPE2>::const_iterator it;
    for (it=rhs.begin(); it!=rhs.end(); ++it) {
        bytes += packsize(it->first);
        bytes += packsize(it->second);
    }
    return bytes;
}
template<class TYPE1, class TYPE2>
void pack( const std::map<TYPE1,TYPE2>& rhs, char *buffer )
{
    size_t N = rhs.size();
    pack(N,buffer);
    size_t pos = sizeof(size_t);
    typename std::map<TYPE1,TYPE2>::const_iterator it;
    for (it=rhs.begin(); it!=rhs.end(); ++it) {
        pack(it->first,&buffer[pos]);   pos+=packsize(it->first);
        pack(it->second,&buffer[pos]);  pos+=packsize(it->second);
    }
}
template<class TYPE1, class TYPE2>
void unpack( std::map<TYPE1,TYPE2>& data, const char *buffer )
{
    size_t N = 0;
    unpack(N,buffer);
    size_t pos = sizeof(size_t);
    data.clear();
    for (size_t i=0; i<N; i++) {
        std::pair<TYPE1,TYPE2> tmp;
        unpack(tmp.first,&buffer[pos]);   pos+=packsize(tmp.first);
        unpack(tmp.second,&buffer[pos]);  pos+=packsize(tmp.second);
        data.insert(tmp);
    }
}
// size_t
template<>
size_t packsize<size_t>( const size_t& rhs )
{
    return sizeof(size_t);
}
template<>
void pack<size_t>( const size_t& rhs, char *buffer )
{
    memcpy(buffer,&rhs,sizeof(size_t));
}
template<>
void unpack<size_t>( size_t& data, const char *buffer )
{
    memcpy(&data,buffer,sizeof(size_t));
}
// unsigned char
template<>
size_t packsize<unsigned char>( const unsigned char& rhs )
{
    return sizeof(unsigned char);
}
template<>
void pack<unsigned char>( const unsigned char& rhs, char *buffer )
{
    memcpy(buffer,&rhs,sizeof(unsigned char));
}
template<>
void unpack<unsigned char>( unsigned char& data, const char *buffer )
{
    memcpy(&data,buffer,sizeof(unsigned char));
}
// std::string
template<>
size_t packsize<std::string>( const std::string& rhs )
{
    return rhs.size()+1;
}
template<>
void pack<std::string>( const std::string& rhs, char *buffer )
{
    memcpy(buffer,rhs.c_str(),rhs.size()+1);
}
template<>
void unpack<std::string>( std::string& data, const char *buffer )
{
    data = std::string(buffer);
}
// MeshType
template<>
size_t packsize<MeshType>( const MeshType& rhs )
{
    return sizeof(MeshType);
}
template<>
void pack<MeshType>( const MeshType& rhs, char *buffer )
{
    memcpy(buffer,&rhs,sizeof(MeshType));
}
template<>
void unpack<MeshType>( MeshType& data, const char *buffer )
{
    memcpy(&data,buffer,sizeof(MeshType));
}
// DatabaseEntry
template<>
size_t packsize<DatabaseEntry>( const DatabaseEntry& rhs )
{
    return packsize(rhs.name)+packsize(rhs.file)+packsize(rhs.offset);
}
template<>
void pack<DatabaseEntry>( const DatabaseEntry& rhs, char *buffer )
{
    size_t i=0;
    pack(rhs.name,&buffer[i]);      i+=packsize(rhs.name);
    pack(rhs.file,&buffer[i]);      i+=packsize(rhs.file);
    pack(rhs.offset,&buffer[i]);    i+=packsize(rhs.offset);
}
template<>
void unpack<DatabaseEntry>( DatabaseEntry& data, const char *buffer )
{
    size_t i=0;
    unpack(data.name,&buffer[i]);       i+=packsize(data.name);
    unpack(data.file,&buffer[i]);       i+=packsize(data.file);
    unpack(data.offset,&buffer[i]);     i+=packsize(data.offset);
}
// MeshDatabase
template<>
size_t packsize<MeshDatabase>( const MeshDatabase& data )
{
    return packsize(data.name) 
         + packsize(data.type)
         + packsize(data.meshClass)
         + packsize(data.format)
         + packsize(data.domains)
         + packsize(data.variables)
         + packsize(data.variable_data);
}
template<>
void pack<MeshDatabase>( const MeshDatabase& rhs, char *buffer )
{
    size_t i = 0;
    pack(rhs.name,&buffer[i]);          i+=packsize(rhs.name);
    pack(rhs.type,&buffer[i]);          i+=packsize(rhs.type);
    pack(rhs.meshClass,&buffer[i]);     i+=packsize(rhs.meshClass);
    pack(rhs.format,&buffer[i]);        i+=packsize(rhs.format);
    pack(rhs.domains,&buffer[i]);       i+=packsize(rhs.domains);
    pack(rhs.variables,&buffer[i]);     i+=packsize(rhs.variables);
    pack(rhs.variable_data,&buffer[i]); i+=packsize(rhs.variable_data);
}
template<>
void unpack<MeshDatabase>( MeshDatabase& data, const char *buffer )
{
    size_t i=0;
    unpack(data.name,&buffer[i]);       i+=packsize(data.name);
    unpack(data.type,&buffer[i]);       i+=packsize(data.type);
    unpack(data.meshClass,&buffer[i]);  i+=packsize(data.meshClass);
    unpack(data.format,&buffer[i]);     i+=packsize(data.format);
    unpack(data.domains,&buffer[i]);    i+=packsize(data.domains);
    unpack(data.variables,&buffer[i]);  i+=packsize(data.variables);
    unpack(data.variable_data,&buffer[i]); i+=packsize(data.variable_data);
}


/****************************************************
* MeshDatabase                                      *
****************************************************/
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
    meshClass = rhs.meshClass;
    format = rhs.format;
    domains = rhs.domains;
    variables = rhs.variables;
    variable_data = rhs.variable_data;
}
MeshDatabase& MeshDatabase::operator=(const MeshDatabase& rhs)
{
    this->name = rhs.name;
    this->type = rhs.type;
    this->meshClass = rhs.meshClass;
    this->format = rhs.format;
    this->domains = rhs.domains;
    this->variables = rhs.variables;
    this->variable_data = rhs.variable_data;
    return *this;
}


/****************************************************
* Split a line into pieces                          *
****************************************************/
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


/****************************************************
* DatabaseEntry                                     *
****************************************************/
std::string DatabaseEntry::write( ) const
{
    char tmp[1000];
    sprintf(tmp,"%s; %s; %lu",name.c_str(),file.c_str(),offset);
    return std::string(tmp);
}
DatabaseEntry::DatabaseEntry( const char* line )
{
    std::vector<std::string> list = splitList(line);
    name = list[0];
    file = list[1];
    offset = atol(list[2].c_str());
}
void DatabaseEntry::read( const char* line )
{
    std::vector<std::string> list = splitList(line);
    name = list[0];
    file = list[1];
    offset = atol(list[2].c_str());
}
void DatabaseEntry::read( const std::string& line )
{
    std::vector<std::string> list = splitList(line.c_str());
    name = list[0];
    file = list[1];
    offset = atol(list[2].c_str());
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
        localsize += packsize(meshes[i]);
    char *localbuf = new char[localsize];
    size_t pos = 0;
    for (size_t i=0; i<meshes.size(); i++) {
        pack( meshes[i], &localbuf[pos] );
        pos += packsize(meshes[i]);
    }
    // Get the number of bytes each processor will be sending/recieving
    int sendsize = static_cast<int>(localsize);
    int *recvsize = new int[size];
    MPI_Allgather(&sendsize,1,MPI_INT,recvsize,1,MPI_INT,comm);
    size_t globalsize = recvsize[0];
    int *disp = new int[size];
    disp[0] = 0;
    for (int i=1; i<size; i++) {
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
        MeshDatabase tmp;
        unpack(tmp,&globalbuf[pos]);
        pos += packsize(tmp);
        std::map<std::string,MeshDatabase>::iterator it = data.find(tmp.name);
        if ( it==data.end() ) {
            data[tmp.name] = tmp;
        } else {
            for (size_t i=0; i<tmp.domains.size(); i++)
                it->second.domains.push_back(tmp.domains[i]);
            for (size_t i=0; i<tmp.variables.size(); i++)
                it->second.variables.push_back(tmp.variables[i]);
            it->second.variable_data.insert(tmp.variable_data.begin(),tmp.variable_data.end());
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
        fprintf(fid,"   type: %i\n",static_cast<int>(meshes[i].type));
        fprintf(fid,"   meshClass: %s\n",meshes[i].meshClass.c_str());
        fprintf(fid,"   format: %i\n",static_cast<int>(meshes[i].format));
        for (size_t j=0; j<meshes[i].domains.size(); j++)
            fprintf(fid,"   domain: %s\n",meshes[i].domains[j].write().c_str());
        fprintf(fid,"   variables: ");
        for (size_t j=0; j<meshes[i].variables.size(); j++) {
            fprintf(fid,"%s; ",meshes[i].variables[j].c_str());
        }
        fprintf(fid,"\n");
        std::map<std::pair<std::string,std::string>,DatabaseEntry>::const_iterator it;
        for (it=meshes[i].variable_data.begin(); it!=meshes[i].variable_data.end(); ++it) {
            const char* domain = it->first.first.c_str();
            const char* variable = it->first.second.c_str();
            fprintf(fid,"   variable(%s,%s): %s\n",domain,variable,it->second.write().c_str());
        }
    }
    fclose(fid);
    PROFILE_STOP("write");
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
        } else if ( strncmp(line,"   meshClass:",13)==0 ) {
            meshes.back().meshClass = deblank(std::string(&line[13]));
        } else if ( strncmp(line,"   domain:",10)==0 ) {
            DatabaseEntry data(&line[10]);
            meshes.back().domains.push_back(data);
        } else if ( strncmp(line,"   variables:",13)==0 ) {
            meshes.back().variables = splitList(&line[13]);
        } else if ( strncmp(line,"   variable(",12)==0 ) {
            size_t i1 = find(line,',');
            size_t i2 = find(line,':');
            std::string domain = deblank(std::string(line,13,i1-13));
            std::string variable = deblank(std::string(line,i1+1,i2-i1));
            std::pair<std::string,std::string> key(domain,variable);
            DatabaseEntry data(&line[i2+1]);
            meshes.back().variable_data.insert( 
                std::pair<std::pair<std::string,std::string>,DatabaseEntry>(key,data) );
        } else {
            ERROR("Error reading line");
        }
    }
    fclose(fid);
    delete [] line;
    PROFILE_STOP("read");
    return meshes;
}


// Return the mesh type
IO::MeshType meshType( std::shared_ptr<IO::Mesh> mesh )
{
    IO::MeshType type = IO::MeshType::Unknown;
    if ( std::dynamic_pointer_cast<IO::PointList>(mesh)!=NULL ) {
        type = IO::MeshType::PointMesh;
    } else if ( std::dynamic_pointer_cast<IO::TriList>(mesh)!=NULL || std::dynamic_pointer_cast<IO::TriMesh>(mesh)!=NULL ) {
        type = IO::MeshType::SurfaceMesh;
    } else {
        ERROR("Unknown mesh");
    }
    return type;
}


} // IO namespace

