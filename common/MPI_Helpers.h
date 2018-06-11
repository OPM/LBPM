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
    #define MPI_COMM_NULL -1
    #define MPI_GROUP_NULL -2
    #define MPI_STATUS_IGNORE NULL
    enum MPI_Datatype { MPI_LOGICAL, MPI_CHAR, MPI_UNSIGNED_CHAR, MPI_INT, 
        MPI_UNSIGNED, MPI_LONG, MPI_UNSIGNED_LONG, MPI_LONG_LONG, MPI_FLOAT, MPI_DOUBLE };
    enum MPI_Op { MPI_MIN, MPI_MAX, MPI_SUM };
    typedef int MPI_Group;
    #define MPI_THREAD_SINGLE 0
    #define MPI_THREAD_FUNNELED 1
    #define MPI_THREAD_SERIALIZED 2
    #define MPI_THREAD_MULTIPLE 3
    // Fake MPI functions
	int MPI_Init(int*,char***);
    int MPI_Init_thread( int *argc, char ***argv, int required, int *provided );
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
    int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, const int *recvcounts, const int *displs,
                  MPI_Datatype recvtype, MPI_Comm comm);
    int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                int dest, int sendtag,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int source, int recvtag,
                MPI_Comm comm, MPI_Status *status);
    int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm);
    double MPI_Wtime( void );
    int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
    int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
    int MPI_Comm_free(MPI_Comm *group);
    int MPI_Group_free(MPI_Group *group);
    int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
#endif


//! Get the size of the MPI_Comm
//  Note: this is a thread and interrupt safe function
inline int comm_size( MPI_Comm comm ) {
    int size = 1;
    MPI_Comm_size( comm, &size );
    return size;
}
    

//! Get the rank of the MPI_Comm
//  Note: this is a thread and interrupt safe function
inline int comm_rank( MPI_Comm comm ) {
    int rank = 1;
    MPI_Comm_rank( comm, &rank );
    return rank;
}
    

//! Get the size of MPI_COMM_WORLD
inline int MPI_WORLD_SIZE( ) {
    return comm_size( MPI_COMM_WORLD );
}

//! Get the size of MPI_COMM_WORLD
inline int MPI_WORLD_RANK( ) {
    return comm_rank( MPI_COMM_WORLD );
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



// Helper functions
inline double sumReduce( MPI_Comm comm, double x )
{
    double y = 0;
	MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_SUM,comm);
    return y;
}
inline float sumReduce( MPI_Comm comm, float x )
{
    float y = 0;
	MPI_Allreduce(&x,&y,1,MPI_FLOAT,MPI_SUM,comm);
    return y;
}
inline int sumReduce( MPI_Comm comm, int x )
{
    int y = 0;
	MPI_Allreduce(&x,&y,1,MPI_INT,MPI_SUM,comm);
    return y;
}
inline bool sumReduce( MPI_Comm comm, bool x )
{
    int y = sumReduce( comm, x?1:0 );
    return y>0;
}
inline std::vector<float> sumReduce( MPI_Comm comm, const std::vector<float>& x )
{
    auto y = x;
	MPI_Allreduce(x.data(),y.data(),x.size(),MPI_FLOAT,MPI_SUM,comm);
    return y;
}
inline std::vector<int> sumReduce( MPI_Comm comm, const std::vector<int>& x )
{
    auto y = x;
	MPI_Allreduce(x.data(),y.data(),x.size(),MPI_INT,MPI_SUM,comm);
    return y;
}
inline double maxReduce( MPI_Comm comm, double x )
{
    double y = 0;
	MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_MAX,comm);
    return y;
}
inline float maxReduce( MPI_Comm comm, float x )
{
    float y = 0;
	MPI_Allreduce(&x,&y,1,MPI_FLOAT,MPI_MAX,comm);
    return y;
}
inline int maxReduce( MPI_Comm comm, int x )
{
    int y = 0;
	MPI_Allreduce(&x,&y,1,MPI_INT,MPI_MAX,comm);
    return y;
}


#endif


#include "common/MPI_Helpers.hpp"


