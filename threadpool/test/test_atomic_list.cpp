#include "threadpool/atomic_list.h"
#include "common/UnitTest.h"
#include "common/Utilities.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <string>
#include <thread>
#include <vector>



static void modify_list( AtomicList<int, 1024> &list )
{
    const int N_count = 50000;
    for ( int i = 0; i < N_count; i++ ) {
        auto v1 = list.remove_first();
        auto v2 = list.remove( []( int ) { return true; } );
        auto v3 = list.remove( []( int v ) { return v >= ( rand() / 8 ); } );
        auto v4 = list.remove( []( int v ) { return v >= ( rand() / 4 ); } );
        auto v5 = list.remove( []( int v ) { return v >= ( rand() / 2 ); } );
        if ( v1 != -1 ) {
            list.insert( v1 );
        }
        if ( v2 != -1 ) {
            list.insert( v2 );
        }
        if ( v3 != -1 ) {
            list.insert( v3 );
        }
        if ( v4 != -1 ) {
            list.insert( v4 );
        }
        if ( v5 != -1 ) {
            list.insert( v5 );
        }
    }
}


static bool check_list( const std::vector<int> &x, AtomicList<int, 1024> &list )
{
    bool pass = list.check();
    pass      = pass && (int) x.size() == list.size();
    if ( pass ) {
        for ( int i : x )
            pass = pass && i == list.remove( []( int ) { return true; } );
    }
    // Restore the list
    for ( int i = 0; i < list.size(); i++ )
        list.remove_first();
    for ( int i : x )
        list.insert( i );
    return pass;
}


static inline void clear_list( AtomicList<int, 1024> &list )
{
    for ( int i = 0; i < list.size(); i++ )
        list.remove_first();
}


/******************************************************************
 * The main program                                                *
 ******************************************************************/
int main( int, char *[] )
{
    UnitTest ut;

    int N_threads = 8; // Number of threads

    // Create the list
    AtomicList<int, 1024> list( -1 );
    if ( list.size() == 0 && list.check() )
        ut.passes( "Initialize" );
    else
        ut.failure( "Initialize" );

    // Initialize the list with some empty values
    for ( int i = 0; i < 80; i++ )
        list.insert( rand() );
    list.insert( 2 );
    list.insert( 1 );
    list.insert( rand() );

    // Try to pull off a couple of values
    int v1 = list.remove( []( int a ) { return a == 1; } ); // Find the entry with 1
    int v2 = list.remove( []( int ) { return true; } );     // Get the first entry
    int v3 = list.remove( []( int ) { return false; } );    // Fail to get an entry
    if ( v1 == 1 && v2 == 2 && v3 == -1 && list.size() == 81 && list.check() )
        ut.passes( "Basic sanity test" );
    else
        ut.failure( "Basic sanity test" );

    // Clear the list
    while ( list.remove( []( int ) { return true; } ) != -1 ) {
    }

    // Create a list of known values
    // std::vector<int> data0(512);
    std::vector<int> data0( 5 * N_threads );
    for ( int &i : data0 )
        i = rand();
    auto data = data0;
    std::sort( data.begin(), data.end() );

    // Test the cost to insert
    int N_it = 20;
    for ( int i = 0; i < list.size(); i++ )
        list.remove( []( int ) { return true; } );
    std::chrono::duration<double> time;
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        clear_list( list );
        start = std::chrono::high_resolution_clock::now();
        for ( int i : data0 )
            list.insert( i );
        stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    printf( "insert time/item = %0.0f ns\n", 1e9 * time.count() / ( N_it * data0.size() ) );

    // Test the cost to remove (first)
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        check_list( data, list );
        start = std::chrono::high_resolution_clock::now();
        for ( size_t i = 0; i < data0.size(); i++ )
            list.remove_first();
        stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    printf( "remove (first) time/item = %0.0f ns\n", 1e9 * time.count() / ( N_it * data0.size() ) );

    // Test the cost to remove (in order)
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        check_list( data, list );
        start = std::chrono::high_resolution_clock::now();
        for ( size_t i = 0; i < data0.size(); i++ )
            list.remove( []( int ) { return true; } );
        stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    printf(
        "remove (ordered) time/item = %0.0f ns\n", 1e9 * time.count() / ( N_it * data0.size() ) );

    // Test the cost to remove (out order)
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        check_list( data, list );
        start = std::chrono::high_resolution_clock::now();
        for ( int tmp : data0 ) {
            list.remove( [tmp]( int v ) { return v == tmp; } );
        }
        stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    printf(
        "remove (unordered) time/item = %0.0f ns\n", 1e9 * time.count() / ( N_it * data0.size() ) );

    // Read/write to the list and check the results
    int64_t N0 = list.N_remove();
    check_list( data, list );
    start = std::chrono::high_resolution_clock::now();
    modify_list( list );
    stop               = std::chrono::high_resolution_clock::now();
    double time_serial = std::chrono::duration<double>( stop - start ).count();
    int64_t N1         = list.N_remove();
    bool pass          = check_list( data, list );
    if ( pass )
        ut.passes( "Serial get/insert" );
    else
        ut.failure( "Serial get/insert" );
    printf( "serial time = %0.5f s\n", time_serial );
    printf( "serial time/item = %0.0f ns\n", 1e9 * time_serial / ( N1 - N0 ) );

    // Have multiple threads reading/writing to the list simultaneously
    std::vector<std::thread> threads( N_threads );
    start = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < N_threads; i++ )
        threads[i] = std::thread( modify_list, std::ref( list ) );
    for ( int i = 0; i < N_threads; i++ )
        threads[i].join();
    stop                 = std::chrono::high_resolution_clock::now();
    double time_parallel = std::chrono::duration<double>( stop - start ).count();
    int64_t N2           = list.N_remove();
    pass                 = check_list( data, list );
    if ( pass )
        ut.passes( "Parallel get/insert" );
    else
        ut.failure( "Parallel get/insert" );
    printf( "parallel time = %0.5f s\n", time_parallel );
    printf( "parallel time/item = %0.0f ns\n", 1e9 * time_parallel / ( N2 - N1 ) );

    // Try to over-fill the list
    while ( !list.empty() )
        list.remove_first();
    for ( int i = 1; i <= list.capacity(); i++ )
        list.insert( i );
    try {
        list.insert( list.capacity() + 1 );
        ut.failure( "List overflow" );
    } catch ( const std::exception &e ) {
        ut.passes( "List overflow" );
    } catch ( ... ) {
        ut.failure( "List overflow (unknown exception)" );
    }

    // Finished
    ut.report();
    auto N_errors = static_cast<int>( ut.NumFailGlobal() );
    return N_errors;
}
