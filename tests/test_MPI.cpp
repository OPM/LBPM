#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "common/MPI.h"
#include "common/UnitTest.h"
#include "common/Utilities.h"
#include "common/Utilities.hpp"
#include "common/ScaLBL.h"
#include "ProfilerApp.h"


#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
#include <windows.h>
#define sched_yield() Sleep( 0 )
#else
#include <sched.h>
#endif

#undef MPI_CLASS
#define MPI_CLASS Utilities::MPI
#define MPI_ASSERT ASSERT


// Return the time elapsed in seconds
static inline double time() { return MPI_CLASS::time(); }


struct mytype {
    int a;
    double b;
    mytype()
    {
        a = -1;
        b = -1.0;
    }
    mytype( int i )
    {
        a = i;
        b = -1.0;
    }
    mytype( int i, double d )
    {
        a = i;
        b = d;
    }
    bool operator==( const mytype &other )
    {
        if ( a == other.a && b == other.b )
            return true;
        return false;
    }
    bool operator!=( const mytype &other )
    {
        if ( a != other.a || b != other.b )
            return true;
        return false;
    }
};


// Routines to test Reduce with known data types
// flag - 0: all tests should pass
//        1: basic reduce should pass, reduce with rank should fail with error message
template<class type>
int testReduce( MPI_CLASS comm, UnitTest *ut, int flag );
template<>
int testReduce<std::complex<double>>( MPI_CLASS comm, UnitTest *ut, int )
{
    PROFILE_START( "testReduce<complex double>" );
    char message[128];
    std::complex<double> rank = comm.getRank() + 1;
    std::complex<double> N    = ( ( comm.getSize() * ( comm.getSize() + 1 ) ) / 2 );
    // Test sumReduce
    sprintf( message, "sumReduce (%s)", typeid( std::complex<double> ).name() );
    if ( comm.sumReduce<std::complex<double>>( rank ) == N )
        ut->passes( message );
    else
        ut->failure( message );
    sprintf( message, "sumReduce (%s) (x,y)", typeid( std::complex<double> ).name() );
    std::complex<double> y;
    comm.sumReduce<std::complex<double>>( &rank, &y, 1 );
    if ( y == N )
        ut->passes( message );
    else
        ut->failure( message );
    PROFILE_STOP( "testReduce<complex double>" );
    return 2; // Return the number of tests
}
template<class type>
int testReduce( MPI_CLASS comm, UnitTest *ut, int flag )
{
    PROFILE_START( "testReduce" );
    char message[128];
    auto rank = (type) comm.getRank();
    auto size = (type) comm.getSize();
    if ( (int) ( size ) != comm.getSize() ) {
        sprintf( message,
                 "Reduce (%s) cannot represent the number of processors",
                 typeid( type ).name() );
        ut->expected_failure( message );
        PROFILE_STOP2( "testReduce<class type>" );
        return 0;
    }
    type x, y;
    int N = ( ( comm.getSize() * ( comm.getSize() + 1 ) ) / 2 );
    // Test sumReduce
    sprintf( message, "sumReduce (%s)", typeid( type ).name() );
    if ( ( (int) ( (type) N ) ) != N )
        ut->expected_failure( message ); // type cannot represent N
    else if ( comm.sumReduce<type>( rank + 1 ) == (type) N )
        ut->passes( message );
    else
        ut->failure( message );
    sprintf( message, "sumReduce (%s) (x,y)", typeid( type ).name() );
    x = rank + 1;
    comm.sumReduce<type>( &x, &y, 1 );
    if ( ( (int) ( (type) N ) ) != N )
        ut->expected_failure( message );
    else if ( y == (type) N )
        ut->passes( message );
    else
        ut->failure( message );
    // Test minReduce
    sprintf( message, "minReduce (%s)", typeid( type ).name() );
    if ( comm.minReduce<type>( rank + 1 ) == 1 )
        ut->passes( message );
    else
        ut->failure( message );
    sprintf( message, "minReduce (%s) (x,y)", typeid( type ).name() );
    comm.minReduce<type>( &x, &y, 1, nullptr );
    if ( y == 1 )
        ut->passes( message );
    else
        ut->failure( message );
    // Test maxReduce
    sprintf( message, "maxReduce (%s)", typeid( type ).name() );
    if ( comm.maxReduce<type>( rank + 1 ) == size )
        ut->passes( message );
    else
        ut->failure( message );
    sprintf( message, "maxReduce (%s) (x,y)", typeid( type ).name() );
    comm.maxReduce<type>( &x, &y, 1, nullptr );
    if ( y == size )
        ut->passes( message );
    else
        ut->failure( message );
    // Test minReduce with rank
    int rank_of_min = -1;
    int rank_of_max = -1;
    type rank_min   = rank + 1;
    type rank_max   = rank + 1;
    sprintf( message, "minReduce-rank (%s)", typeid( type ).name() );
    try {
        comm.minReduce<type>( &rank_min, 1, &rank_of_min );
        if ( rank_min == 1 && rank_of_min == 0 )
            ut->passes( message );
        else
            ut->failure( message );
        if ( flag == 1 && comm.getSize() > 1 )
            ut->failure( message );
    } catch ( ... ) {
        if ( flag == 1 && comm.getSize() > 1 )
            ut->expected_failure( message );
        else
            ut->failure( message );
    }
    sprintf( message, "minReduce-rank (%s) (x,y)", typeid( type ).name() );
    try {
        comm.minReduce<type>( &x, &rank_min, 1, &rank_of_min );
        if ( rank_min == 1 && rank_of_min == 0 )
            ut->passes( message );
        else
            ut->failure( message );
        if ( flag == 1 && comm.getSize() > 1 )
            ut->failure( message );
    } catch ( ... ) {
        if ( flag == 1 && comm.getSize() > 1 )
            ut->expected_failure( message );
        else
            ut->failure( message );
    }
    // Test maxReduce with rank
    sprintf( message, "maxReduce-rank (%s)", typeid( type ).name() );
    try {
        comm.maxReduce<type>( &rank_max, 1, &rank_of_max );
        if ( rank_max == size && rank_of_max == comm.getSize() - 1 )
            ut->passes( message );
        else
            ut->failure( message );
        if ( flag == 1 && comm.getSize() > 1 )
            ut->failure( message );
    } catch ( ... ) {
        if ( flag == 1 && comm.getSize() > 1 )
            ut->expected_failure( message );
        else
            ut->failure( message );
    }
    sprintf( message, "maxReduce-rank (%s) (x,y)", typeid( type ).name() );
    try {
        comm.maxReduce<type>( &x, &rank_max, 1, &rank_of_max );
        if ( rank_max == size && rank_of_max == comm.getSize() - 1 )
            ut->passes( message );
        else
            ut->failure( message );
        if ( flag == 1 && comm.getSize() > 1 )
            ut->failure( message );
    } catch ( ... ) {
        if ( flag == 1 && comm.getSize() > 1 )
            ut->expected_failure( message );
        else
            ut->failure( message );
    }
    PROFILE_STOP( "testReduce" );
    return 10; // Return the number of tests
}


// Routine to test Scan with known data types
// flag - 0: all tests should pass
//        1: only sumScan is valid (complex<double>)
template<class type>
int testScan( MPI_CLASS comm, UnitTest *ut, int flag = 0 )
{
    PROFILE_START( "testScan" );
    char message[500];
    auto x = ( type )( comm.getRank() + 1 );
    type y;
    sprintf( message, "sumScan (%s)", typeid( type ).name() );
    comm.sumScan<type>( &x, &y, 1 );
    auto N = ( type )( ( ( comm.getRank() + 1 ) * ( comm.getRank() + 2 ) ) / 2 );
    if ( y == N )
        ut->passes( message );
    else
        ut->failure( message );
    if ( flag == 1 ) {
        PROFILE_STOP2( "testScan" );
        return 1;
    }
    sprintf( message, "minScan (%s)", typeid( type ).name() );
    comm.minScan<type>( &x, &y, 1 );
    if ( y == (type) 1 )
        ut->passes( message );
    else
        ut->failure( message );
    sprintf( message, "maxScan (%s)", typeid( type ).name() );
    comm.maxScan<type>( &x, &y, 1 );
    if ( y == x )
        ut->passes( message );
    else
        ut->failure( message );
    PROFILE_STOP( "testScan" );
    return 3; // Return the number of tests
}


// Routine to test bcast
template<class type>
int testBcast( MPI_CLASS comm, UnitTest *ut, type default_val, type new_val )
{
    PROFILE_START( "testBcast" );
    char message[128];
    for ( int i = 0; i < comm.getSize(); i++ ) {
        type tmp1 = default_val;
        if ( comm.getRank() == i )
            tmp1 = new_val;
        sprintf( message, "bcast scalar (%s) from rank %i", typeid( type ).name(), i );
        if ( comm.bcast( tmp1, i ) == new_val )
            ut->passes( message );
        else
            ut->failure( message );
        type tmp2[2];
        tmp2[0] = default_val;
        tmp2[1] = default_val;
        if ( comm.getRank() == i ) {
            tmp2[0] = new_val;
            tmp2[1] = new_val;
        }
        sprintf( message, "bcast vector (%s) from rank %i", typeid( type ).name(), i );
        comm.bcast( tmp2, 2, i );
        if ( tmp2[0] == new_val && tmp2[1] == new_val )
            ut->passes( message );
        else
            ut->failure( message );
    }
    PROFILE_STOP( "testBcast" );
    return 2 * comm.getSize(); // Return the number of tests
}


// Routine to test allGather
template<class type>
int testAllGather( MPI_CLASS comm, UnitTest *ut )
{
    PROFILE_START( "testAllGather" );
    char message[128];
    // Test scalar allGather
    auto x1  = (type) comm.getRank();
    auto *x2 = new type[comm.getSize()];
    comm.allGather( x1, x2 );
    bool pass = true;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        type test = i;
        if ( x2[i] != test )
            pass = false;
    }
    sprintf( message, "allGather scalar (%s)", typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    // Test vector allGather
    int N     = ( comm.getSize() * ( comm.getSize() + 1 ) ) / 2;
    auto *x3  = new type[comm.getRank() + 1];
    auto *x4  = new type[N];
    auto *x5  = new type[N];
    int *size = new int[comm.getSize()];
    for ( int i = 0; i <= comm.getRank(); i++ )
        x3[i] = (type) comm.getRank();
    int tot1 = comm.allGather( x3, comm.getRank() + 1, x4 );
    int tot2 = comm.allGather( x3, comm.getRank() + 1, x5, size );
    pass     = true;
    if ( tot1 != N || tot2 != N )
        pass = false;
    int k = 0;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        if ( size[i] != i + 1 )
            pass = false;
        if ( !pass )
            break;
        for ( int j = 0; j <= i; j++ ) {
            type test = i;
            if ( x4[k] != test || x5[k] != test )
                pass = false;
            k++;
        }
    }
    sprintf( message, "allGather vector (%s)", typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    delete[] x2;
    delete[] x3;
    delete[] x4;
    delete[] x5;
    delete[] size;
    // Test vector allGather with know recive sizes and non-zero displacements
    auto *send      = new type[comm.getRank() + 1];
    auto *recv      = new type[comm.getSize() * comm.getSize() + 1];
    auto *recv_size = new int[comm.getSize()];
    auto *recv_disp = new int[comm.getSize()];
    for ( int i = 0; i <= comm.getRank(); i++ )
        send[i] = i;
    for ( int i = 0; i < comm.getSize(); i++ )
        recv_size[i] = i + 1;
    for ( int i = 0; i < comm.getSize(); i++ )
        recv_disp[i] = 1 + i * comm.getSize() + comm.getSize() - i - 1;
    for ( int i = 0; i <= comm.getSize() * comm.getSize(); i++ )
        recv[i] = (type) -1;
    int tot = comm.allGather( send, comm.getRank() + 1, recv, recv_size, recv_disp, true );
    pass    = true;
    if ( tot != N )
        pass = false;
    auto test = (type) -1;
    if ( recv[0] != test )
        pass = false;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        for ( int j = 0; j < comm.getSize(); j++ ) {
            int l = j + i * comm.getSize() + 1 - recv_disp[i];
            if ( l >= 0 )
                test = l;
            else
                test = (type) -1;
            if ( recv[j + i * comm.getSize() + 1] != test )
                pass = false;
        }
    }
    sprintf( message,
             "allGather vector with known recv and non-zero displacements (%s)",
             typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    delete[] send;
    delete[] recv;
    delete[] recv_size;
    delete[] recv_disp;
    // Test vector allGather with no elements
    size = new int[comm.getSize()];
    sprintf( message, "allGather scalar (%s)", typeid( type ).name() );
    try {
        comm.allGather( &x1, 0, (type *) nullptr, size );
        ut->passes( message );
    } catch ( ... ) {
        ut->failure( message );
    }
    delete[] size;
    PROFILE_STOP( "testAllGather" );
    return 4; // Return the number of tests
}


// Routine to test setGather
template<class type>
int testSetGather( MPI_CLASS comm, UnitTest *ut )
{
    PROFILE_START( "testSetGather" );
    char message[500];
    auto x1 = (type) comm.getRank();
    std::set<type> set;
    set.insert( x1 );
    comm.setGather( set );
    bool pass = true;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        type x2 = i;
        if ( set.find( x2 ) == set.end() )
            pass = false;
    }
    sprintf( message, "setGather (%s)", typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    PROFILE_STOP( "testSetGather" );
    return 1; // Return the number of tests
}


// Routine to test mapGather
template<class type>
int testMapGather( MPI_CLASS comm, UnitTest *ut )
{
    PROFILE_START( "testMapGather" );
    char message[128];
    auto x1 = (type) comm.getRank();
    std::map<int, type> map;
    map.insert( std::pair<int, type>( comm.getRank(), x1 ) );
    comm.mapGather( map );
    bool pass = true;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        type x2 = i;
        auto it = map.find( i );
        if ( it == map.end() )
            pass = false;
        else if ( it->second != x2 )
            pass = false;
    }
    sprintf( message, "mapGather (%s)", typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    PROFILE_STOP( "testMapGather" );
    return 1; // Return the number of tests
}


// Routine to test allToAll
template<class type>
int testAllToAll( MPI_CLASS comm, UnitTest *ut )
{
    PROFILE_START( "testAllToAll" );
    bool pass;
    char message[128];
    int size = 0;
    type *send_data, *recv_data;
    auto *send_cnt  = new int[comm.getSize()];
    auto *recv_cnt  = new int[comm.getSize()];
    auto *send_disp = new int[comm.getSize()];
    auto *recv_disp = new int[comm.getSize()];
    // Test allToAll with a scalar value to each processor
    send_data = new type[comm.getSize()];
    recv_data = new type[comm.getSize()];
    for ( int i = 0; i < comm.getSize(); i++ )
        send_data[i] = comm.getSize();
    comm.allToAll( 1, send_data, recv_data );
    pass = true;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        type test = comm.getSize();
        if ( recv_data[i] != test )
            pass = false;
    }
    delete[] send_data;
    delete[] recv_data;
    sprintf( message, "allToAll with scalar (%s)", typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    // Test allToAll vector with a scalar value to each processor
    send_data = new type[comm.getSize()];
    recv_data = new type[comm.getSize()];
    for ( int i = 0; i < comm.getSize(); i++ ) {
        send_cnt[i]  = 1;
        recv_cnt[i]  = 1;
        send_disp[i] = i;
        recv_disp[i] = i;
        send_data[i] = comm.getSize();
        recv_data[i] = 0;
    }
    size = comm.allToAll( send_data, send_cnt, send_disp, recv_data, recv_cnt, recv_disp, true );
    pass = true;
    if ( size != comm.getSize() )
        pass = false;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        type test = comm.getSize();
        if ( recv_data[i] != test )
            pass = false;
    }
    delete[] send_data;
    delete[] recv_data;
    sprintf( message, "allToAll vector with scalar (%s)", typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    // Test allToAll with a variable number of values per processor and spacing
    send_data = new type[comm.getSize() * comm.getSize()];
    recv_data = new type[2 * comm.getRank() * comm.getSize()];
    for ( int i = 0; i < comm.getSize(); i++ ) {
        send_cnt[i]  = i;
        recv_cnt[i]  = comm.getRank();
        send_disp[i] = i * comm.getSize();
        recv_disp[i] = 2 * i * comm.getRank();
        for ( int j = 0; j < comm.getSize(); j++ ) {
            if ( j < i )
                send_data[j + send_disp[i]] = i;
            else
                send_data[j + send_disp[i]] = (type) -1;
        }
    }
    for ( int i = 0; i < 2 * comm.getRank() * comm.getSize(); i++ )
        recv_data[i] = (type) -2;
    size = comm.allToAll( send_data, send_cnt, send_disp, recv_data, recv_cnt, recv_disp, true );
    pass = true;
    if ( size != comm.getRank() * comm.getSize() )
        pass = false;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        for ( int j = 0; j < 2 * comm.getRank(); j++ ) {
            if ( j < comm.getRank() ) {
                type test = comm.getRank();
                if ( recv_data[j + recv_disp[i]] != test )
                    pass = false;
            } else {
                auto test = (type) -2;
                if ( recv_data[j + recv_disp[i]] != test )
                    pass = false;
            }
        }
    }
    delete[] send_data;
    delete[] recv_data;
    sprintf( message,
             "allToAll with vector of known size and displacements (%s)",
             typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    // Test allToAll with a unknown recieve length
    send_data        = new type[comm.getSize() * comm.getSize()];
    auto *recv_data1 = new type[comm.getSize() * comm.getSize()];
    auto *recv_data2 = new type[comm.getSize() * comm.getSize()];
    for ( int i = 0; i < comm.getSize(); i++ ) {
        send_cnt[i]  = i;
        recv_cnt[i]  = -1;
        send_disp[i] = i * comm.getSize();
        recv_disp[i] = -1;
        for ( int j = 0; j < comm.getSize(); j++ ) {
            if ( j < i )
                send_data[j + send_disp[i]] = i;
            else
                send_data[j + send_disp[i]] = (type) -1;
        }
    }
    for ( int i = 0; i < comm.getSize() * comm.getSize(); i++ ) {
        recv_data1[i] = (type) -2;
        recv_data2[i] = (type) -2;
    }
    int size1 =
        comm.allToAll( send_data, send_cnt, send_disp, recv_data1, recv_cnt, recv_disp, false );
    int size2  = comm.allToAll( send_data, send_cnt, send_disp, recv_data2 );
    bool pass1 = true;
    bool pass2 = true;
    if ( size1 != comm.getRank() * comm.getSize() )
        pass1 = false;
    if ( size2 != comm.getRank() * comm.getSize() )
        pass2 = false;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        if ( recv_cnt[i] != comm.getRank() || recv_disp[i] != i * comm.getRank() )
            pass1 = false;
    }
    for ( int i = 0; i < comm.getRank() * comm.getSize(); i++ ) {
        type test = comm.getRank();
        if ( recv_data1[i] != test )
            pass1 = false;
        if ( recv_data2[i] != test )
            pass2 = false;
    }
    delete[] send_data;
    delete[] recv_data1;
    delete[] recv_data2;
    sprintf( message, "allToAll with vector of unknown size (%s)", typeid( type ).name() );
    if ( pass1 )
        ut->passes( message );
    else
        ut->failure( message );
    sprintf(
        message, "allToAll with vector of unknown size with NULL recv(%s)", typeid( type ).name() );
    if ( pass2 )
        ut->passes( message );
    else
        ut->failure( message );
    // Free temporary variables
    delete[] send_cnt;
    delete[] recv_cnt;
    delete[] send_disp;
    delete[] recv_disp;
    PROFILE_STOP( "testAllToAll" );
    return 5; // Return the number of tests
}


// Routine to test send/recv
template<class type>
int testSendRecv( MPI_CLASS comm, UnitTest *ut, type v1, type v2 )
{
    PROFILE_START( "testSendRecv" );
    char message[128];
    // Test send-recv with a known length
    for ( int i = 0; i < comm.getSize(); i++ ) {
        for ( int j = 0; j < comm.getSize(); j++ ) {
            type x  = v1;
            int tag = i + j * comm.getSize();
            sprintf( message, "send-recv %i-%i known length (%s)", i, j, typeid( type ).name() );
            if ( i == j ) {
                // We are not allowed to send/recieve from the same processor
                continue;
            } else if ( i == comm.getRank() ) {
                // We are sending
                x = v2;
                comm.send( &x, 1, j, tag );
            } else if ( j == comm.getRank() ) {
                // We are recieving
                int size = 1;
                comm.recv( &x, size, i, false, tag );
                if ( size == 1 && x == v2 )
                    ut->passes( message );
                else
                    ut->failure( message );
            }
        }
    }
    // Test send-recv with an unknown length
    for ( int i = 0; i < comm.getSize(); i++ ) {
        for ( int j = 0; j < comm.getSize(); j++ ) {
            type x  = v1;
            int tag = i + j * comm.getSize();
            sprintf( message, "send-recv %i-%i unknown length (%s)", i, j, typeid( type ).name() );
            if ( i == j ) {
                // We are not allowed to send/recieve from the same processor
                continue;
            } else if ( i == comm.getRank() ) {
                // We are sending
                x = v2;
                comm.send( &x, 1, j, tag );
            } else if ( j == comm.getRank() ) {
                // We are recieving
                int size = 1;
                comm.recv( &x, size, i, true, tag );
                if ( size == 1 && x == v2 )
                    ut->passes( message );
                else
                    ut->failure( message );
            }
        }
    }
    // Test send-recv with an empty length
    for ( int i = 0; i < comm.getSize(); i++ ) {
        for ( int j = 0; j < comm.getSize(); j++ ) {
            type x  = v1;
            int tag = i + j * comm.getSize();
            sprintf( message, "send-recv %i-%i empty length (%s)", i, j, typeid( type ).name() );
            if ( i == j ) {
                // We are not allowed to send/recieve from the same processor
                continue;
            } else if ( i == comm.getRank() ) {
                // We are sending
                x = v2;
                comm.send( &x, 0, j, tag );
            } else if ( j == comm.getRank() ) {
                // We are recieving
                int size = comm.probe( i, tag );
                comm.recv( &x, size, i, false, tag );
                if ( size == 0 )
                    ut->passes( message );
                else
                    ut->failure( message );
            }
        }
    }
    PROFILE_STOP( "testSendRecv" );
    return 3 * comm.getSize() * comm.getSize(); // Return the number of tests
}


// Routine to test Isend/Irecv
template<class type>
int testIsendIrecv( MPI_CLASS comm, UnitTest *ut, type v1, type v2 )
{
    PROFILE_START( "testIsendIrecv" );
    char message[128];
    std::vector<MPI_Request> sendRequest;
    std::vector<MPI_Request> recvRequest;
    // Send all messages
    for ( int i = 0; i < comm.getSize(); i++ ) {
        // Check if the current rank is sending
        if ( i != comm.getRank() )
            continue;
        for ( int j = 0; j < comm.getSize(); j++ ) {
            // Start a non-blocking send
            int tag             = i + j * comm.getSize();
            MPI_Request request = comm.Isend( &v1, 1, j, tag );
            sendRequest.insert( sendRequest.begin(), request );
        }
    }
    // Recv all messages
    auto *recv_buffer = new type[comm.getSize()];
    for ( int i = 0; i < comm.getSize(); i++ )
        recv_buffer[i] = v2;
    recv_buffer[comm.getRank()] = v1;
    for ( int j = 0; j < comm.getSize(); j++ ) {
        // Check if the current rank is recieving
        if ( j != comm.getRank() )
            continue;
        for ( int i = 0; i < comm.getSize(); i++ ) {
            // Start a non-blocking recv
            int tag             = i + j * comm.getSize();
            MPI_Request request = comm.Irecv( &recv_buffer[i], 1, i, tag );
            recvRequest.insert( recvRequest.begin(), request );
        }
    }
    // Wait for all communications to finish
    MPI_CLASS::wait( sendRequest[0] );
    sendRequest.erase( sendRequest.begin() + 0 );
    while ( !sendRequest.empty() ) {
        int index = comm.waitAny( sendRequest.size(), &( sendRequest[0] ) );
        sendRequest.erase( sendRequest.begin() + index );
    }
    auto finished = MPI_CLASS::waitSome( recvRequest.size(), recvRequest.data() );
    if ( !recvRequest.empty() ) {
        MPI_ASSERT( !finished.empty() );
        for ( auto it = finished.rbegin(); it != finished.rend(); ++it )
            recvRequest.erase( recvRequest.begin() + ( *it ) );
    }
    if ( !recvRequest.empty() )
        MPI_CLASS::waitAll( recvRequest.size(), &( recvRequest[0] ) );
    Utilities::unique( finished );
    // Check the recieved values
    bool pass = true;
    for ( int i = 0; i < comm.getSize(); i++ ) {
        if ( recv_buffer[i] != v1 )
            pass = false;
    }
    sprintf( message, "Isend-Irecv (%s)", typeid( type ).name() );
    if ( pass )
        ut->passes( message );
    else
        ut->failure( message );
    delete[] recv_buffer;
    PROFILE_STOP( "testIsendIrecv" );
    return comm.getSize() * comm.getSize(); // Return the number of tests
}


// Routine to test CommRanks
int testCommRanks( MPI_CLASS comm, UnitTest *ut )
{
    std::vector<int> neighbors;
    for ( int i = 0; i < comm.getSize(); i++ )
        if ( ( i % 2 ) == 0 )
            neighbors.push_back( i );
    std::vector<int> neighbors2 = comm.commRanks( neighbors );
    bool pass                   = true;
    if ( comm.getRank() % 2 == 0 ) {
        pass = static_cast<int>( neighbors2.size() ) == comm.getSize();
        if ( pass ) {
            for ( int i = 0; i < comm.getSize(); i++ )
                pass = pass && neighbors2[i] == i;
        }
    } else {
        pass = neighbors2.empty();
    }
    auto ranks = comm.globalRanks();
    pass       = pass && (int) ranks.size() == comm.getSize();
    for ( int rank : ranks )
        pass = pass && rank >= 0;
    auto ranks2 = ranks;
    Utilities::unique( ranks2 );
    pass = pass && ranks.size() == ranks2.size();
    comm.barrier();
    if ( pass )
        ut->passes( "commRanks" );
    else
        ut->failure( "commRanks" );
    return 1; // Return the number of tests
}


// Structure to contain timer results
struct testCommTimerResults {
    int N_reduce;
    int N_scan;
    int N_bcast;
    int N_allGather;
    int N_setGather;
    int N_mapGather;
    int N_allToAll;
    int N_sendRecv;
    int N_IsendIrecv;
    double t_reduce;
    double t_scan;
    double t_bcast;
    double t_allGather;
    double t_setGather;
    double t_mapGather;
    double t_allToAll;
    double t_sendRecv;
    double t_IsendIrecv;
    // Constructor
    testCommTimerResults()
    {
        N_reduce     = 0;
        N_scan       = 0;
        N_bcast      = 0;
        N_allGather  = 0;
        N_setGather  = 0;
        N_mapGather  = 0;
        N_allToAll   = 0;
        N_sendRecv   = 0;
        N_IsendIrecv = 0;
        t_reduce     = 0.0;
        t_scan       = 0.0;
        t_bcast      = 0.0;
        t_allGather  = 0.0;
        t_setGather  = 0.0;
        t_mapGather  = 0.0;
        t_allToAll   = 0.0;
        t_sendRecv   = 0.0;
        t_IsendIrecv = 0.0;
    }
    // Print the results
    void print()
    {
        printf( "   Reduce:      N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_reduce,
                t_reduce,
                1e6 * t_reduce / N_reduce );
        printf( "   Scan:        N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_scan,
                t_scan,
                1e6 * t_scan / N_scan );
        printf( "   Bcast:       N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_bcast,
                t_bcast,
                1e6 * t_bcast / N_bcast );
        printf( "   allGather:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_allGather,
                t_allGather,
                1e6 * t_allGather / N_allGather );
        printf( "   allToAll:    N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_allToAll,
                t_allToAll,
                1e6 * t_allToAll / N_allToAll );
        printf( "   send-recv:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_sendRecv,
                t_sendRecv,
                1e6 * t_sendRecv / N_sendRecv );
        printf( "   Isend-Irecv: N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_IsendIrecv,
                t_IsendIrecv,
                1e6 * t_IsendIrecv / N_IsendIrecv );
        printf( "   setGather:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_setGather,
                t_setGather,
                1e6 * t_setGather / N_setGather );
        printf( "   mapGather:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",
                N_mapGather,
                t_mapGather,
                1e6 * t_mapGather / N_mapGather );
    }
};


// This routine will test a single MPI communicator
testCommTimerResults testComm( MPI_CLASS comm, UnitTest *ut )
{
    PROFILE_START( "testComm" );
    testCommTimerResults timer;
    double start_time;
    // Test the tag
    int tag0        = comm.newTag();
    MPI_CLASS comm2 = comm;
    MPI_CLASS comm3( comm );
    bool pass = tag0 > 0 && tag0 < comm.maxTag();
    for ( int i = 1; i < 64; i++ ) {
        if ( comm.newTag() != tag0 + i )
            pass = false;
    }
    for ( int i = 1; i <= 64; i++ ) {
        if ( comm2.newTag() != tag0 + 63 + i )
            pass = false;
    }
    for ( int i = 1; i <= 128; i++ ) {
        if ( comm3.newTag() != tag0 + 127 + i )
            pass = false;
    }
    if ( pass )
        ut->passes( "newTag" );
    else
        ut->failure( "newTag" );
    // Test all and any reduce
    bool test1 = !comm.allReduce( comm.getRank() != 0 );
    bool test2 = comm.allReduce( true );
    if ( test1 && test2 )
        ut->passes( "allReduce" );
    else
        ut->failure( "allReduce" );
    test1 = comm.anyReduce( comm.getRank() == 0 );
    test2 = !comm.anyReduce( false );
    if ( test1 && test2 )
        ut->passes( "anyReduce" );
    else
        ut->failure( "anyReduce" );
    // Test min, max, and sum reduce
    start_time = time();
    timer.N_reduce += testReduce<unsigned char>( comm, ut, 0 );
    timer.N_reduce += testReduce<char>( comm, ut, 0 );
    timer.N_reduce += testReduce<unsigned int>( comm, ut, 0 );
    timer.N_reduce += testReduce<int>( comm, ut, 0 );
    timer.N_reduce += testReduce<unsigned long int>( comm, ut, 0 );
    timer.N_reduce += testReduce<long int>( comm, ut, 0 );
    timer.N_reduce += testReduce<size_t>( comm, ut, 0 );
    timer.N_reduce += testReduce<float>( comm, ut, 0 );
    timer.N_reduce += testReduce<double>( comm, ut, 0 );
    timer.N_reduce += testReduce<std::complex<double>>(
        comm, ut, 2 ); // only sumreduce is valid for complex numbers
    mytype tmp1( 1, -1.0 );
    mytype tmp2;
    if ( comm.getSize() > 1 ) {
        // We can't perform a reduce on an unknown data type (this should throw an error)
        try {
            // This should fail
            tmp2 = comm.sumReduce<mytype>( tmp1 );
            ut->failure( "sumReduce should give an error with an unknown type" );
        } catch ( ... ) {
            ut->passes( "sumReduce should give an error with an unknown type" );
        }
        try {
            // This should fail
            tmp2 = comm.minReduce<mytype>( tmp1 );
            ut->failure( "minReduce should give an error with an unknown type" );
        } catch ( ... ) {
            ut->passes( "minReduce should give an error with an unknown type" );
        }
        try {
            // This should fail
            tmp2 = comm.maxReduce<mytype>( tmp1 );
            ut->failure( "maxReduce should give an error with an unknown type" );
        } catch ( ... ) {
            ut->passes( "maxReduce should give an error with an unknown type" );
        }
        timer.N_reduce += 3;
    }
    timer.t_reduce = time() - start_time;
    // Test min, max, and sum scan
    start_time = time();
    timer.N_scan += testScan<unsigned char>( comm, ut );
    timer.N_scan += testScan<char>( comm, ut );
    timer.N_scan += testScan<unsigned int>( comm, ut );
    timer.N_scan += testScan<int>( comm, ut );
    timer.N_scan += testScan<unsigned long int>( comm, ut );
    timer.N_scan += testScan<long int>( comm, ut );
    timer.N_scan += testScan<size_t>( comm, ut );
    timer.N_scan += testScan<float>( comm, ut );
    timer.N_scan += testScan<double>( comm, ut );
    timer.N_scan +=
        testScan<std::complex<double>>( comm, ut, 1 ); // Only sumScan is valid with complex data
    if ( comm.getSize() > 1 ) {
        // We can't perform a reduce on an unknown data type (this should throw an error)
        try {
            // This should fail
            comm.sumScan<mytype>( &tmp1, &tmp2, 1 );
            ut->failure( "sumReduce should give an error with an unknown type" );
        } catch ( ... ) {
            ut->passes( "sumReduce should give an error with an unknown type" );
        }
        try {
            // This should fail
            comm.minScan<mytype>( &tmp1, &tmp2, 1 );
            ut->failure( "minReduce should give an error with an unknown type" );
        } catch ( ... ) {
            ut->passes( "minReduce should give an error with an unknown type" );
        }
        try {
            // This should fail
            comm.maxScan<mytype>( &tmp1, &tmp2, 1 );
            ut->failure( "maxReduce should give an error with an unknown type" );
        } catch ( ... ) {
            ut->passes( "maxReduce should give an error with an unknown type" );
        }
        timer.N_scan += 3;
    }
    timer.t_scan = time() - start_time;
    // Test bcast
    start_time = time();
    timer.N_bcast += testBcast<unsigned char>( comm, ut, 0, 1 );
    timer.N_bcast += testBcast<char>( comm, ut, -1, 1 );
    timer.N_bcast += testBcast<unsigned int>( comm, ut, 0, 1 );
    timer.N_bcast += testBcast<int>( comm, ut, -1, 1 );
    timer.N_bcast += testBcast<unsigned long int>( comm, ut, 0, 1 );
    timer.N_bcast += testBcast<long int>( comm, ut, -1, 1 );
    timer.N_bcast += testBcast<size_t>( comm, ut, 0, 1 );
    timer.N_bcast += testBcast<float>( comm, ut, -1.0, 1.0 );
    timer.N_bcast += testBcast<double>( comm, ut, -1.0, 1.0 );
    mytype tmp3( -1, -1.0 );
    mytype tmp4( 1, 1.0 );
    timer.N_bcast += testBcast<mytype>( comm, ut, tmp3, tmp4 );
    timer.t_bcast = time() - start_time;
    // Test barrier
    comm.barrier();
    // Test gather
    start_time = time();
    timer.N_allGather += testAllGather<unsigned char>( comm, ut );
    timer.N_allGather += testAllGather<char>( comm, ut );
    timer.N_allGather += testAllGather<unsigned int>( comm, ut );
    timer.N_allGather += testAllGather<int>( comm, ut );
    timer.N_allGather += testAllGather<unsigned long int>( comm, ut );
    timer.N_allGather += testAllGather<long int>( comm, ut );
    timer.N_allGather += testAllGather<size_t>( comm, ut );
    timer.N_allGather += testAllGather<float>( comm, ut );
    timer.N_allGather += testAllGather<double>( comm, ut );
    timer.N_allGather += testAllGather<std::complex<double>>( comm, ut );
    timer.N_allGather += testAllGather<mytype>( comm, ut );
    timer.t_allGather = time() - start_time;
    // Test std::set gather
    start_time = time();
    timer.N_setGather += testSetGather<unsigned char>( comm, ut );
    timer.N_setGather += testSetGather<char>( comm, ut );
    timer.N_setGather += testSetGather<unsigned int>( comm, ut );
    timer.N_setGather += testSetGather<int>( comm, ut );
    timer.N_setGather += testSetGather<unsigned long int>( comm, ut );
    timer.N_setGather += testSetGather<long int>( comm, ut );
    timer.N_setGather += testSetGather<size_t>( comm, ut );
    timer.N_setGather += testSetGather<float>( comm, ut );
    timer.N_setGather += testSetGather<double>( comm, ut );
    timer.t_setGather = time() - start_time;
    // Test std::map gather
    start_time = time();
    timer.N_mapGather += testMapGather<unsigned char>( comm, ut );
    timer.N_mapGather += testMapGather<char>( comm, ut );
    timer.N_mapGather += testMapGather<unsigned int>( comm, ut );
    timer.N_mapGather += testMapGather<int>( comm, ut );
    timer.N_mapGather += testMapGather<unsigned long int>( comm, ut );
    timer.N_mapGather += testMapGather<long int>( comm, ut );
    timer.N_mapGather += testMapGather<size_t>( comm, ut );
    timer.N_mapGather += testMapGather<float>( comm, ut );
    timer.N_mapGather += testMapGather<double>( comm, ut );
    timer.t_mapGather = time() - start_time;
    // Test allToAlll
    start_time = time();
    timer.N_allToAll += testAllToAll<unsigned char>( comm, ut );
    timer.N_allToAll += testAllToAll<char>( comm, ut );
    timer.N_allToAll += testAllToAll<unsigned int>( comm, ut );
    timer.N_allToAll += testAllToAll<int>( comm, ut );
    timer.N_allToAll += testAllToAll<unsigned long int>( comm, ut );
    timer.N_allToAll += testAllToAll<long int>( comm, ut );
    timer.N_allToAll += testAllToAll<size_t>( comm, ut );
    timer.N_allToAll += testAllToAll<float>( comm, ut );
    timer.N_allToAll += testAllToAll<double>( comm, ut );
    timer.N_allToAll += testAllToAll<std::complex<double>>( comm, ut );
    timer.N_allToAll += testAllToAll<mytype>( comm, ut );
    timer.t_allToAll = time() - start_time;
    // Test send/recv
    start_time = time();
    timer.N_sendRecv += testSendRecv<unsigned char>( comm, ut, 0, 1 );
    timer.N_sendRecv += testSendRecv<char>( comm, ut, -1, 1 );
    timer.N_sendRecv += testSendRecv<unsigned int>( comm, ut, 0, 1 );
    timer.N_sendRecv += testSendRecv<int>( comm, ut, -1, 1 );
    timer.N_sendRecv += testSendRecv<unsigned long int>( comm, ut, 0, 1 );
    timer.N_sendRecv += testSendRecv<long int>( comm, ut, -1, 1 );
    timer.N_sendRecv += testSendRecv<size_t>( comm, ut, 0, 1 );
    timer.N_sendRecv += testSendRecv<float>( comm, ut, -1.0, 1.0 );
    timer.N_sendRecv += testSendRecv<double>( comm, ut, -1.0, 1.0 );
    timer.N_sendRecv += testSendRecv<mytype>( comm, ut, tmp3, tmp4 );
    timer.t_sendRecv = time() - start_time;
    // Test Isend/Irecv
    start_time = time();
    timer.N_IsendIrecv += testIsendIrecv<unsigned char>( comm, ut, 0, 1 );
    timer.N_IsendIrecv += testIsendIrecv<char>( comm, ut, -1, 1 );
    timer.N_IsendIrecv += testIsendIrecv<unsigned int>( comm, ut, 0, 1 );
    timer.N_IsendIrecv += testIsendIrecv<int>( comm, ut, -1, 1 );
    timer.N_IsendIrecv += testIsendIrecv<unsigned long int>( comm, ut, 0, 1 );
    timer.N_IsendIrecv += testIsendIrecv<long int>( comm, ut, -1, 1 );
    timer.N_IsendIrecv += testIsendIrecv<size_t>( comm, ut, 0, 1 );
    timer.N_IsendIrecv += testIsendIrecv<float>( comm, ut, -1.0, 1.0 );
    timer.N_IsendIrecv += testIsendIrecv<double>( comm, ut, -1.0, 1.0 );
    timer.N_IsendIrecv += testIsendIrecv<mytype>( comm, ut, tmp3, tmp4 );
    timer.t_IsendIrecv = time() - start_time;
    // Test commRanks
    testCommRanks( comm, ut );
    PROFILE_STOP( "testComm" );
    return timer;
}


// Test comm dup and the number of communicators that can be created
void testCommDup( UnitTest *ut )
{
#if defined( USING_CLANG ) && defined( __APPLE__ )
    // The MPI error handler crashes so this test fails
    // This seems to be a MAC? + Clang + MPICH? issue only
    ut->expected_failure( "testCommDup skipped for this architecture/compiler" );
#else
    MPI_CLASS globalComm( MPI_COMM_WORLD );
    MPI_CLASS dupComm = globalComm.dup();
    if ( globalComm.getCommunicator() != dupComm.getCommunicator() &&
         dupComm.getSize() == globalComm.getSize() && dupComm.getRank() == globalComm.getRank() ) {
        ut->passes( "dup comm" );
    } else {
        ut->failure( "dup comm" );
        return;
    }
#if defined( USE_PETSC ) && !defined( USE_MPI )
    ut->expected_failure( "Skipping dup tests, PETSc (no-mpi) has a limit of 128 unique comms" );
    return;
#endif
    int N_comm_try = 2000; // Maximum number of comms to try and create
    std::vector<MPI_CLASS> comms;
    comms.reserve( N_comm_try );
    try {
        for ( int i = 0; i < N_comm_try; i++ ) {
            MPI_CLASS tmp_comm = globalComm.dup();
            comms.push_back( tmp_comm );
            MPI_ASSERT( globalComm.getCommunicator() != comms[i].getCommunicator() );
            MPI_ASSERT( comms.back().sumReduce<int>( 1 ) ==
                        globalComm.getSize() ); // We need to communicate as part of the test
        }
        ut->passes( Utilities::stringf( "Created %i comms", N_comm_try ) );
    } catch ( ... ) {
        if ( comms.size() < 252 ) {
            ut->failure( "Could not create 252 different communicators" );
        } else {
            int N = comms.size();
            ut->expected_failure( Utilities::stringf(
                "Failed to create an unlimited number of comms (%i)", N ) );
        }
        std::cout << "Maximum number of concurrent communicators: " << comms.size() << std::endl;
    }
    comms.clear();
    size_t N_dup = 0;
    globalComm.barrier();
    try {
        double start = MPI_CLASS::time();
        for ( int i = 0; i < N_comm_try; i++ ) {
            MPI_CLASS tmp_comm1 = globalComm.dup();
            MPI_CLASS tmp_comm2 = globalComm.dup();
            MPI_ASSERT( globalComm.getCommunicator() != tmp_comm1.getCommunicator() );
            MPI_ASSERT( globalComm.getCommunicator() != tmp_comm2.getCommunicator() );
            MPI_ASSERT( tmp_comm1.getCommunicator() != tmp_comm2.getCommunicator() );
            MPI_ASSERT( tmp_comm1.sumReduce<int>( 1 ) ==
                        globalComm.getSize() ); // We need to communicate as part of the test
            MPI_ASSERT( tmp_comm2.sumReduce<int>( 1 ) ==
                        globalComm.getSize() ); // We need to communicate as part of the test
            N_dup += 2;
        }
        double stop = MPI_CLASS::time();
        ut->passes( "Created/Destroyed an unlimited number of comms" );
        char message[128];
        sprintf( message,
                 "Time to create/destroy comm using MPI_CLASS::dup() is: %0.1f us",
                 1e6 * ( stop - start ) / N_dup );
        std::cout << message << std::endl;
    } catch ( ... ) {
        ut->failure( "Failed to create/destroy an unlimited number of comms" );
        std::cout << "Maximum number of communicators created with destruction: " << N_dup
                  << std::endl;
    }
#endif
}

class gpuWrapper{
public:
	gpuWrapper(MPI_CLASS MPI_COMM, int MSG_SIZE);
	~gpuWrapper();
	void Send(double *values);
	void Recv(double *values);
	double *sendbuf, *recvbuf;
private:
	MPI_Request req1[1],req2[1];
	MPI_CLASS comm;
	int sendCount;
	int recvCount;
	int rank, rank_x, rank_X, nprocs;
	int sendtag, recvtag;
};

gpuWrapper::gpuWrapper(MPI_CLASS MPI_COMM, int MSG_SIZE){
    comm = MPI_COMM.dup();
    rank = comm.getRank();
    nprocs = comm.getSize();
	sendCount=MSG_SIZE;
	recvCount=MSG_SIZE;
	ScaLBL_AllocateZeroCopy((void **) &sendbuf, sendCount*sizeof(double));	// Allocate device memory
	ScaLBL_AllocateZeroCopy((void **) &recvbuf, sendCount*sizeof(double));	// Allocate device memory
	rank_X = rank+1;
	rank_x = rank-1;
	if (rank_x < 0)          rank_x = nprocs-1;
	if (!(rank_X < nprocs))  rank_X = 0;
}

gpuWrapper::~gpuWrapper(){
	ScaLBL_FreeDeviceMemory(sendbuf);
	ScaLBL_FreeDeviceMemory(recvbuf);
}

void gpuWrapper::Send(double *values){
	sendtag = recvtag = 130;
    ScaLBL_CopyToDevice(sendbuf,values,sendCount*sizeof(double));
	req1[0] = comm.Isend(sendbuf,sendCount,rank_x,sendtag+0);
	req2[0] = comm.Irecv(recvbuf,recvCount,rank_X,recvtag+0);
}

void gpuWrapper::Recv(double *values){
	comm.waitAll(1,req1);
	comm.waitAll(1,req2);
	ScaLBL_DeviceBarrier();
    ScaLBL_CopyToHost(values,recvbuf,recvCount*sizeof(double));
}

// Test GPU aware MPI
void test_GPU_aware( UnitTest *ut )
{
    constexpr size_t N = 1024*1024;
    constexpr size_t N_msg = 64;
    bool test = true;
    // Get the comm to use
    MPI_CLASS comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
    int size = comm.getSize();
    try {
        // Initialize the device
        int device = ScaLBL_SetDevice(rank); 
        NULL_USE( device );
        // create wrapper for communications
        gpuWrapper gpuComm(comm, N); 
        // Allocate and initialize the buffers
        size_t bytes = N*sizeof(double);
        double *device_send[N_msg] = { nullptr };
        double *device_recv[N_msg] = { nullptr };
        double *host_send[N_msg] = { nullptr };
        double *host_recv[N_msg] = { nullptr };
        for ( size_t k=0; k<N_msg; k++ ) {
            ScaLBL_AllocateDeviceMemory((void**)&device_send[k],bytes);
            ScaLBL_AllocateDeviceMemory((void**)&device_recv[k],bytes);
            host_send[k] = new double[N];
            host_recv[k] = new double[N];
            // Initialize the data
            for ( size_t i=0; i<N; i++ ) {
                host_send[k][i] = 1000 * k * rank + i;
                host_recv[k][i] = 0;
            }
            ScaLBL_CopyToDevice(device_send[k],host_send[k],bytes);
            ScaLBL_CopyToDevice(device_recv[k],host_recv[k],bytes);
        }
        ScaLBL_DeviceBarrier();
        // Send/recieve the data
        int rank_send = ( rank + 1 ) % size;
        int rank_recv = ( rank - 1 + size ) % size;
        MPI_Request req1[N_msg];
        MPI_Request req2[N_msg];
        for ( size_t k=0; k<N_msg; k++ ) {
            req1[k] = comm.Isend( device_send[k], N, rank_send, k );
            req2[k] = comm.Irecv( device_recv[k], N, rank_recv, k );
        }
        comm.waitAll(N_msg,req1);
        comm.waitAll(N_msg,req2);
        // Copy
        for ( size_t k=0; k<N_msg; k++ ) {
            ScaLBL_CopyToHost(host_send[k],device_send[k],bytes);
            ScaLBL_CopyToHost(host_recv[k],device_recv[k],bytes);
        }
        ScaLBL_DeviceBarrier();
        // Check the data
        for ( size_t k=0; k<N_msg; k++ ) {
            for ( size_t i=0; i<N; i++ )
                test = test && host_recv[k][i] == 1000 * k * rank_recv + i;
        }
        // Check the gpu wrapper communications the same way
        for ( size_t k=0; k<N_msg; k++ ) {
        	gpuComm.Send(host_send[k]);
        	gpuComm.Recv(host_recv[k]);
        }
        // Check the data
        for ( size_t k=0; k<N_msg; k++ ) {
            for ( size_t i=0; i<N; i++ )
                test = test && host_recv[k][i] == 1000 * k * rank_recv + i;
        }
        
        // Free buffers
        for ( size_t k=0; k<N_msg; k++ ) {
            ScaLBL_FreeDeviceMemory(device_send[k]);
            ScaLBL_FreeDeviceMemory(device_recv[k]);
            delete [] host_send[k];
            delete [] host_recv[k];
        }
    } catch ( ... ) {
        test = false;
    }
    comm.barrier();
    if ( test ) {
        std::cout << "MPI is GPU aware" << std::endl;
        ut->passes("GPU aware MPI" );
    } else {
        std::cout << "MPI is NOT GPU aware" << std::endl;
        ut->failure("GPU aware MPI" );
    }
}


//  This test will test the MPI class
int main( int argc, char *argv[] )
{
    // Start MPI
    Utilities::MPI::start_MPI( argc, argv );

    // Create the unit test
    UnitTest ut;
    PROFILE_ENABLE( 0 );
    PROFILE_START( "Main" );


    // Limit the scope so objects are destroyed
    {

        // Get the start time for the tests
        double start_time = time();

        // Print the global size (if we are using MPI)
        int global_size = 0;
#ifdef USE_MPI
        MPI_Comm_size( MPI_COMM_WORLD, &global_size );
#else
        global_size = 1;
#endif

        // Test the global communicator (MPI_COMM_WORLD)
        MPI_CLASS globalComm = MPI_CLASS( MPI_COMM_WORLD );
        if ( !globalComm.isNull() )
            ut.passes( "Global communicator created" );
        else
            ut.failure( "Global communicator created" );
        if ( globalComm.getSize() == global_size )
            ut.passes( "Global communicator size" );
        else
            ut.failure( "Global communicator size" );
        if ( globalComm.getRank() == 0 ) {
            std::cout << "MPI_COMM_WORLD = " << global_size << " processors" << std::endl;
            std::cout << "   Largest tag value = " << globalComm.maxTag() << std::endl << std::endl;
        }
#ifdef USE_MPI
        if ( globalComm.getCommunicator() == MPI_COMM_WORLD )
            ut.passes( "Communicator == MPI_COMM_WORLD" );
        else
            ut.failure( "Communicator == MPI_COMM_WORLD" );
#endif
        testCommTimerResults commTimer = testComm( globalComm, &ut );
        if ( globalComm.getRank() == 0 ) {
            std::cout << "Results for global timer (rank 0)" << std::endl;
            commTimer.print();
            std::cout << std::endl;
        }

        // Test bcast with std::string
        std::string rank_string;
        if ( globalComm.getRank() == 0 )
            rank_string = "Rank 0";
        rank_string = globalComm.bcast( rank_string, 0 );
        if ( rank_string == "Rank 0" )
            ut.passes( "Bcast std::string" );
        else
            ut.failure( "Bcast std::string" );

        // Test MPI_COMM_SELF
        MPI_CLASS selfComm = MPI_CLASS( MPI_COMM_SELF );
        if ( !selfComm.isNull() )
            ut.passes( "Self communicator created" );
        else
            ut.failure( "Self communicator created" );
#ifdef USE_MPI
        if ( selfComm.getCommunicator() == MPI_COMM_SELF )
            ut.passes( "Communicator == MPI_COMM_SELF" );
        else
            ut.failure( "Communicator == MPI_COMM_SELF" );
#endif
        testComm( selfComm, &ut );

        // Test == and !=
        if ( globalComm == globalComm && !( selfComm == globalComm ) )
            ut.passes( "==" );
        else
            ut.failure( "==" );
        if ( selfComm != globalComm && !( globalComm != globalComm ) )
            ut.passes( "!=" );
        else
            ut.failure( "!=" );

        // Test MPI_COMM_NULL
        MPI_CLASS nullComm = MPI_CLASS( MPI_COMM_NULL );
        if ( nullComm.isNull() )
            ut.passes( "Null communicator created" );
        else
            ut.failure( "Null communicator created" );
        if ( nullComm.getSize() == 0 )
            ut.passes( "Null communicator has zero size" );
        else
            ut.failure( "Null communicator has zero size" );
#ifdef USE_MPI
        if ( nullComm.getCommunicator() == MPI_COMM_NULL )
            ut.passes( "Communicator == MPI_COMM_NULL" );
        else
            ut.failure( "Communicator == MPI_COMM_NULL" );
#endif

        // Test dup
        MPI_CLASS dupComm = globalComm.dup();
        if ( nullComm.dup().isNull() )
            ut.passes( "Null communicator duplicates a Null communicator" );
        else
            ut.failure( "Null communicator duplicates a Null communicator" );
        testCommDup( &ut );

        // Test compare
        if ( globalComm.compare( globalComm ) == 1 )
            ut.passes( "compare comm global==global" );
        else
            ut.failure( "compare comm global==global" );
        if ( globalComm.compare( dupComm ) == 3 )
            ut.passes( "compare comm global~=dup" );
        else
            ut.failure( "compare comm global~=dup" );
        if ( global_size == 1 ) {
            if ( globalComm.compare( selfComm ) == 3 )
                ut.passes( "compare comm global~=self (global size=1)" );
            else
                ut.failure( "compare comm global~=self (global size=1)" );
        } else {
            if ( globalComm.compare( selfComm ) == 0 )
                ut.passes( "compare comm global!=self" );
            else
                ut.failure( "compare comm global!=self" );
        }

        // Split the global comm and test
        PROFILE_START( "Split" );
        int color;
        if ( globalComm.getRank() == 0 )
            color = 0;
        else if ( globalComm.getRank() < 3 )
            color = 1;
        else
            color = 2 + ( globalComm.getRank() - 2 ) / 4;
        std::vector<MPI_CLASS> splitComms( 4 );
        splitComms[0] = globalComm.split( color );
        splitComms[1] = globalComm.split( color, globalComm.getRank() );
        if ( splitComms[0].getCommunicator() != globalComm.getCommunicator() &&
             splitComms[1].getCommunicator() != globalComm.getCommunicator() &&
             splitComms[0].getCommunicator() != splitComms[1].getCommunicator() )
            ut.passes( "split comm has different communicator" );
        else
            ut.failure( "split comm has different communicator" );
        if ( globalComm.getSize() > 1 ) {
            if ( splitComms[0].getSize() < globalComm.getSize() )
                ut.passes( "split comm is smaller" );
            else
                ut.failure( "split comm is smaller" );
        }
        if ( splitComms[0].getRank() == splitComms[1].getRank() )
            ut.passes( "split sort by rank" );
        else
            ut.failure( "split sort by rank" );
        testComm( splitComms[0], &ut );
        splitComms[2] = globalComm.split( -1 );
        if ( splitComms[2].isNull() )
            ut.passes( "split with color=-1 returns NULL communicator" );
        else
            ut.failure( "split with color=-1 returns NULL communicator" );
        splitComms[3] = splitComms[0]; // Make a copy to ensure there are no memory leaks
        splitComms[3] = splitComms[2]; // Perform assignement to check memory leaks
        MPI_ASSERT( splitComms[3] == splitComms[2] );
        PROFILE_STOP( "Split" );

        // Test  <  <=  >  >=
        if ( globalComm.getSize() > 1 ) {
            if ( splitComms[0] < globalComm && splitComms[1] < globalComm &&
                 !( globalComm < globalComm ) && !( globalComm < splitComms[0] ) )
                ut.passes( " < comm" );
            else
                ut.failure( " < comm" );
            if ( splitComms[0] <= globalComm && splitComms[1] <= globalComm &&
                 globalComm <= globalComm && !( globalComm <= splitComms[0] ) )
                ut.passes( " <= comm" );
            else
                ut.failure( " <= comm" );
            if ( globalComm > splitComms[0] && globalComm > splitComms[1] &&
                 !( globalComm > globalComm ) && !( splitComms[0] > globalComm ) )
                ut.passes( " > comm" );
            else
                ut.failure( " > comm" );
            if ( globalComm >= splitComms[0] && globalComm >= splitComms[1] &&
                 globalComm >= globalComm && !( splitComms[0] >= globalComm ) )
                ut.passes( " >= comm" );
            else
                ut.failure( " >= comm" );
        }

        // Test intersection
        // Test globalComm with selfComm
        if ( globalComm.getSize() > 1 ) {
            MPI_CLASS comm1 = MPI_CLASS::intersect( globalComm, selfComm );
            MPI_CLASS comm2 = MPI_CLASS::intersect( selfComm, globalComm );
            MPI_CLASS comm3 = MPI_CLASS::intersect( globalComm, globalComm );
            if ( comm1.compare( globalComm ) == 0 && comm1.compare( selfComm ) != 0 &&
                 comm2.compare( globalComm ) == 0 && comm2.compare( selfComm ) != 0 &&
                 comm3.compare( globalComm ) != 0 && comm3.compare( selfComm ) == 0 )
                ut.passes( "intersection of globalComm and selfComm" );
            else
                ut.failure( "intersection of globalComm and selfComm" );
        }

        // Test case where we have disjoint sets (this can only happen of one of the comms is null)
        {
            MPI_CLASS intersection = MPI_CLASS::intersect( globalComm, nullComm );
            if ( intersection.isNull() )
                ut.passes( "intersection of non-overlapping comms" );
            else
                ut.failure( "intersection of non-overlapping comms" );
        }

        // Test case where the comms partially overlap
        if ( globalComm.getSize() > 2 ) {
            int n = globalComm.getSize() - 1;
            // Intersect 2 comms (all other ranks will be null)
            MPI_CLASS split1       = globalComm.split( globalComm.getRank() == 0 ? -1 : 0 );
            MPI_CLASS split2       = globalComm.split( globalComm.getRank() == n ? -1 : 0 );
            MPI_CLASS intersection = MPI_CLASS::intersect( split1, split2 );
            bool pass              = true;
            if ( globalComm.getRank() == 0 || globalComm.getRank() == n ) {
                if ( !intersection.isNull() )
                    pass = false;
            } else {
                if ( intersection.compare( split1 ) != 0 || intersection.compare( split2 ) != 0 ||
                     intersection.getSize() != globalComm.getSize() - 2 )
                    pass = false;
            }
            // Intersect 2 sets for ranks (3 groups should result)
            // split1 = globalComm.split(globalComm.getRank()==0?1:2);
            // split2 = globalComm.split(globalComm.getRank()==n?1:2);
            // intersection = MPI_CLASS::intersect( split1, split2 );
            // bool pass = true;
            // if ( globalComm.getRank()==0 || globalComm.getRank()==n ) {
            //    if ( intersection.compare(selfComm)==0 )
            //        pass = false;
            //} else {
            //    if ( intersection.compare(split1)!=0 || intersection.compare(split2)!=0 ||
            //         intersection.getSize()!=globalComm.getSize()-2 )
            //        pass = false;
            //}
            if ( pass )
                ut.passes( "intersection of partially overlapping comms" );
            else
                ut.failure( "intersection of partially overlapping comms" );
        }

        // Test splitByNode
        MPI_CLASS nodeComm = globalComm.splitByNode();
        std::string localName = MPI_CLASS::getNodeName();
        std::vector<std::string> globalStrings( globalComm.getSize() );
        std::vector<std::string> nodeStrings( nodeComm.getSize() );
        globalComm.allGather<std::string>( localName, &globalStrings[0] );
        nodeComm.allGather<std::string>( localName, &nodeStrings[0] );
        int N_local = 0;
        for ( auto &nodeString : nodeStrings ) {
            if ( nodeString == localName )
                N_local++;
        }
        int N_global = 0;
        for ( auto &globalString : globalStrings ) {
            if ( globalString == localName )
                N_global++;
        }
        if ( !nodeComm.isNull() && N_local == nodeComm.getSize() && N_local == N_global )
            ut.passes( "splitByNode" );
        else
            ut.failure( "splitByNode" );

        // Test the call to load balance the processes
        MPI_CLASS::balanceProcesses( globalComm, 1 );
        std::vector<int> cpus = MPI_CLASS::getProcessAffinity();
        size_t maxProcNode    = nodeComm.maxReduce( cpus.size() );
        bool pass_balance     = cpus.size() == maxProcNode && !cpus.empty();
        MPI_CLASS::balanceProcesses( globalComm, 2 );
        cpus = MPI_CLASS::getProcessAffinity();
        if ( cpus.size() < 1 || cpus.size() > maxProcNode / nodeComm.getSize() )
            pass_balance = false;
        if ( pass_balance ) {
            ut.passes( "balanceProcesses" );
        } else {
#ifdef __APPLE__
            ut.expected_failure( "balanceProcesses" );
#else
            ut.failure( "balanceProcesses" );
#endif
        }

        // Test the performance of sched_yield (used internally by MPI wait routines)
        globalComm.barrier();
        double start_yield = time();
        for ( int i = 0; i < 10000; i++ )
            sched_yield();
        double time_yield = ( time() - start_yield ) / 10000;
        if ( globalComm.getRank() == 0 )
            std::cout << "Time to yield: " << time_yield * 1e6 << " us" << std::endl;

        // Test time and tick
        double end_time = MPI_CLASS::time();
        double time_res = MPI_CLASS::tick();
        if ( globalComm.getRank() == 0 ) {
            std::cout << "Time to run tests: " << end_time - start_time << std::endl;
            std::cout << "Timer resolution: " << time_res << std::endl;
            if ( time_res > 0 && time_res < 1 && ( end_time - start_time ) >= time_res )
                ut.passes( "time and tick" );
            else
                ut.failure( "time and tick" );
            std::cout << std::endl;
        }

        // Test GPU aware MPI
        test_GPU_aware( &ut );

    } // Limit the scope so objects are destroyed

    // Finished testing, report the results
    PROFILE_START( "Report" );
    double start_time = time();
    ut.report();
    int num_failed  = ut.NumFailGlobal();
    double end_time = time();
    if ( MPI_CLASS( MPI_COMM_WORLD ).getRank() == 0 )
        std::cout << "Time to report: " << end_time - start_time << std::endl << std::endl;
    PROFILE_STOP( "Report" );
    PROFILE_STOP( "Main" );

    // Shutdown
    PROFILE_SAVE( "test_MPI" );
    ut.reset();
    Utilities::MPI::stop_MPI();
    return num_failed;
}
