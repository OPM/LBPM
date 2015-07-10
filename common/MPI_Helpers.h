// This file contains wrappers for MPI routines and functions to pack/unpack data structures
#ifndef MPI_WRAPPERS_INC
#define MPI_WRAPPERS_INC

#include <string.h>
#include <vector>
#include <set>
#include <map>

#ifdef USE_MPI
    // Inlcude MPI
    #include "mpi.h"
#else
    // Create fake MPI types
    typedef int MPI_Comm;
    typedef int MPI_Request;
    typedef int MPI_Status;
    #define MPI_COMM_WORLD 0
    #define MPI_COMM_SELF 0
    #define MPI_STATUS_IGNORE NULL
    enum MPI_Datatype { MPI_LOGICAL, MPI_CHAR, MPI_UNSIGNED_CHAR, MPI_INT, 
        MPI_UNSIGNED, MPI_LONG, MPI_UNSIGNED_LONG, MPI_LONG_LONG, MPI_FLOAT, MPI_DOUBLE };
    enum MPI_Op { MPI_MIN, MPI_MAX, MPI_SUM };
    enum MPI_Group {  };
    // Fake MPI functions
	int MPI_Init(int*,char***);
	int MPI_Finalize();
    int MPI_Comm_size( MPI_Comm, int *size );
    int MPI_Comm_rank( MPI_Comm, int *rank );
    int MPI_Barrier(MPI_Comm);
    int MPI_Wait(MPI_Request*,MPI_Status*);
    int MPI_Waitall(int,MPI_Request[],MPI_Status[]);
    int MPI_Bcast(void*,int,MPI_Datatype,int,MPI_Comm);
    int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
    int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status);
    int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
    int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request);
    int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
    int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm);
    int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                int dest, int sendtag,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int source, int recvtag,
                MPI_Comm comm, MPI_Status *status);
    double MPI_Wtime( void );
    int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
    int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
#endif


//! Get the size of MPI_COMM_WORLD
inline int MPI_WORLD_SIZE( ) {
    int size = 1;
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    return size;
}

//! Get the size of MPI_COMM_WORLD
inline int MPI_WORLD_RANK( ) {
    int rank = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    return rank;
}

//! Return the appropriate MPI datatype for a class
template<class TYPE>
MPI_Datatype getMPItype();


//! Template function to return the buffer size required to pack a class
template<class TYPE>
size_t packsize( const TYPE& rhs );

//! Template function to pack a class to a buffer
template<class TYPE>
void pack( const TYPE& rhs, char *buffer );

//! Template function to unpack a class from a buffer
template<class TYPE>
void unpack( TYPE& data, const char *buffer );


//! Template function to return the buffer size required to pack a std::vector
template<class TYPE>
size_t packsize( const std::vector<TYPE>& rhs );

//! Template function to pack a class to a buffer
template<class TYPE>
void pack( const std::vector<TYPE>& rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE>
void unpack( std::vector<TYPE>& data, const char *buffer );


//! Template function to return the buffer size required to pack a std::pair
template<class TYPE1, class TYPE2>
size_t packsize( const std::pair<TYPE1,TYPE2>& rhs );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void pack( const std::pair<TYPE1,TYPE2>& rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void unpack( std::pair<TYPE1,TYPE2>& data, const char *buffer );


//! Template function to return the buffer size required to pack a std::map
template<class TYPE1, class TYPE2>
size_t packsize( const std::map<TYPE1,TYPE2>& rhs );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void pack( const std::map<TYPE1,TYPE2>& rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void unpack( std::map<TYPE1,TYPE2>& data, const char *buffer );


//! Template function to return the buffer size required to pack a std::set
template<class TYPE>
size_t packsize( const std::set<TYPE>& rhs );

//! Template function to pack a class to a buffer
template<class TYPE>
void pack( const std::set<TYPE>& rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE>
void unpack( std::set<TYPE>& data, const char *buffer );



#endif


#include "common/MPI_Helpers.hpp"


