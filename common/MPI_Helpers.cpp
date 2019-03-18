#include "common/MPI_Helpers.h"
#include "common/Utilities.h"


/********************************************************
* Return the MPI data type                              *
********************************************************/
template<> MPI_Datatype getMPItype<char>() {
    return MPI_CHAR;
}
template<> MPI_Datatype getMPItype<unsigned char>() {
    return MPI_UNSIGNED_CHAR;
}
template<> MPI_Datatype getMPItype<int>() {
    return MPI_INT;
}
template<> MPI_Datatype getMPItype<long>() {
    return MPI_LONG;
}
template<> MPI_Datatype getMPItype<unsigned long>() {
    return MPI_UNSIGNED_LONG;
}
template<> MPI_Datatype getMPItype<long long>() {
    return MPI_LONG_LONG;
}
template<> MPI_Datatype getMPItype<float>() {
    return MPI_FLOAT;
}
template<> MPI_Datatype getMPItype<double>() {
    return MPI_DOUBLE;
}


/********************************************************
* Concrete implimentations for packing/unpacking        *
********************************************************/
// unsigned char
template<>
size_t packsize<unsigned char>( const unsigned char& )
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
// char
template<>
size_t packsize<char>( const char& )
{
    return sizeof(char);
}
template<>
void pack<char>( const char& rhs, char *buffer )
{
    memcpy(buffer,&rhs,sizeof(char));
}
template<>
void unpack<char>( char& data, const char *buffer )
{
    memcpy(&data,buffer,sizeof(char));
}
// int
template<>
size_t packsize<int>( const int& )
{
    return sizeof(int);
}
template<>
void pack<int>( const int& rhs, char *buffer )
{
    memcpy(buffer,&rhs,sizeof(int));
}
template<>
void unpack<int>( int& data, const char *buffer )
{
    memcpy(&data,buffer,sizeof(int));
}
// unsigned int
template<>
size_t packsize<unsigned int>( const unsigned int& )
{
    return sizeof(unsigned int);
}
template<>
void pack<unsigned int>( const unsigned int& rhs, char *buffer )
{
    memcpy(buffer,&rhs,sizeof(int));
}
template<>
void unpack<unsigned int>( unsigned int& data, const char *buffer )
{
    memcpy(&data,buffer,sizeof(int));
}
// size_t
template<>
size_t packsize<size_t>( const size_t& )
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


/********************************************************
* Fake MPI routines                                     *
********************************************************/
#ifndef USE_MPI
int MPI_Init(int*,char***)
{
    return 0;
}
int MPI_Init_thread(int*,char***, int required, int *provided )
{
    *provided = required;
    return 0;
}
int MPI_Finalize()
{
    return 0;
}
int MPI_Comm_size( MPI_Comm, int *size )
{
    *size = 1;
    return 0;
}
int MPI_Comm_rank( MPI_Comm, int *rank )
{
    *rank = 0;
    return 0;
}
int MPI_Barrier( MPI_Comm )
{
    return 0;
}
int MPI_Waitall( int, MPI_Request[], MPI_Status[] )
{
    return 0;
}
int MPI_Wait( MPI_Request*, MPI_Status* )
{
    return 0;
}
int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
{
    return 0;
}
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
         MPI_Comm comm)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
         MPI_Comm comm, MPI_Status *status)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, const int *recvcounts, const int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                int dest, int sendtag,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int source, int recvtag,
                MPI_Comm comm, MPI_Status *status)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
    ERROR("Not implimented yet");
    return 0;
}
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
    *newcomm = comm;
    return 0;
}
double MPI_Wtime( void )
{
    return 0.0;
}
int MPI_Comm_free(MPI_Comm *group)
{
    return 0;
}
int MPI_Group_free(MPI_Group *group)
{
    return 0;
}
#endif


