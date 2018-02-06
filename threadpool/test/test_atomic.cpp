#include "threadpool/atomic_helpers.h"
#include "common/UnitTest.h"
#include "common/Utilities.h"
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <string>
#include <thread>
#include <vector>


#define perr std::cerr
#define pout std::cout
#define printp printf


// Function to increment/decrement a counter N times
static void modify_counter( int N, AtomicOperations::counter_t &counter )
{
    if ( N > 0 ) {
        for ( int i = 0; i < N; i++ )
            counter.increment();
    } else if ( N < 0 ) {
        for ( int i = 0; i < -N; i++ )
            counter.decrement();
    }
}


/******************************************************************
 * The main program                                                *
 ******************************************************************/
#ifdef USE_WINDOWS
int __cdecl main( int, char ** )
{
#elif defined( USE_LINUX ) || defined( USE_MAC )
int main( int, char *[] )
{
#else
#error Unknown OS
#endif
    UnitTest ut;

    int N_threads = 64;      // Number of threads
    int N_count   = 1000000; // Number of work items

// Ensure we are using all processors
#ifdef __USE_GNU
    int N_procs = sysconf( _SC_NPROCESSORS_ONLN );
    cpu_set_t mask;
    CPU_ZERO( &mask );
    for ( int i = 0; i < N_procs; i++ )
        CPU_SET( i, &mask );
    sched_setaffinity( getpid(), sizeof( cpu_set_t ), &mask );
#endif

    // Create the counter we want to test
    AtomicOperations::counter_t count;
    if ( count.increment() == 1 )
        ut.passes( "increment count" );
    else
        ut.failure( "increment count" );
    if ( count.decrement() == 0 )
        ut.passes( "decrement count" );
    else
        ut.failure( "decrement count" );
    count.setCount( 3 );
    if ( count.getCount() == 3 )
        ut.passes( "set count" );
    else
        ut.failure( "set count" );
    count.setCount( 0 );

    // Increment the counter in serial
    auto start = std::chrono::high_resolution_clock::now();
    modify_counter( N_count, count );
    auto stop              = std::chrono::high_resolution_clock::now();
    double time_inc_serial = std::chrono::duration<double>( stop - start ).count() / N_count;
    int val                = count.getCount();
    if ( val != N_count ) {
        char tmp[100];
        sprintf( tmp, "Count of %i did not match expected count of %i", val, N_count );
        ut.failure( tmp );
    }
    printp( "Time to increment (serial) = %0.1f ns\n", 1e9 * time_inc_serial );

    // Decrement the counter in serial
    start = std::chrono::high_resolution_clock::now();
    modify_counter( -N_count, count );
    stop                   = std::chrono::high_resolution_clock::now();
    double time_dec_serial = std::chrono::duration<double>( stop - start ).count() / N_count;
    val                    = count.getCount();
    if ( val != 0 ) {
        char tmp[100];
        sprintf( tmp, "Count of %i did not match expected count of %i", val, 0 );
        ut.failure( tmp );
    }
    printp( "Time to decrement (serial) = %0.1f ns\n", 1e9 * time_dec_serial );

    // Increment the counter in parallel
    std::vector<std::thread> threads( N_threads );
    start = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < N_threads; i++ )
        threads[i] = std::thread( modify_counter, N_count, std::ref( count ) );
    for ( int i = 0; i < N_threads; i++ )
        threads[i].join();
    stop = std::chrono::high_resolution_clock::now();
    double time_inc_parallel =
        std::chrono::duration<double>( stop - start ).count() / ( N_count * N_threads );
    val = count.getCount();
    if ( val != N_count * N_threads ) {
        char tmp[100];
        sprintf( tmp, "Count of %i did not match expected count of %i", val, N_count * N_threads );
        ut.failure( tmp );
    }
    printp( "Time to increment (parallel) = %0.1f ns\n", 1e9 * time_inc_parallel );

    // Decrement the counter in parallel
    start = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < N_threads; i++ )
        threads[i] = std::thread( modify_counter, -N_count, std::ref( count ) );
    for ( int i = 0; i < N_threads; i++ )
        threads[i].join();
    stop = std::chrono::high_resolution_clock::now();
    double time_dec_parallel =
        std::chrono::duration<double>( stop - start ).count() / ( N_count * N_threads );
    val = count.getCount();
    if ( val != 0 ) {
        char tmp[100];
        sprintf( tmp, "Count of %i did not match expected count of %i", val, 0 );
        ut.failure( tmp );
    }
    printp( "Time to decrement (parallel) = %0.1f ns\n", 1e9 * time_dec_parallel );

    // Check the time to increment/decrement
    if ( time_inc_serial > 100e-9 || time_dec_serial > 100e-9 || time_inc_parallel > 100e-9 ||
         time_dec_serial > 100e-9 ) {
#if USE_GCOV
        ut.expected_failure( "Time to increment/decrement count is too expensive" );
#else
        ut.failure( "Time to increment/decrement count is too expensive" );
#endif
    } else {
        ut.passes( "Time to increment/decrement passed" );
    }

    // Finished
    ut.report();
    auto N_errors = static_cast<int>( ut.NumFailGlobal() );
    return N_errors;
}
