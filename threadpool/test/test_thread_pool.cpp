#include "ProfilerApp.h"
#ifdef USE_TIMER
#include "MemoryApp.h"
#endif
#include "threadpool/thread_pool.h"
#include "common/UnitTest.h"
#include "common/Utilities.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>


#define MAX( x, y ) ( ( x ) > ( y ) ? ( x ) : ( y ) )


#define perr std::cerr
#define pout std::cout
#define printp printf


#ifdef USE_MPI
#include "mpi.h"
#endif

#define to_ns( x ) std::chrono::duration_cast<std::chrono::nanoseconds>( x ).count()
#define to_ms( x ) std::chrono::duration_cast<std::chrono::milliseconds>( x ).count()


// Wrapper functions for mpi
static inline void barrier()
{
#ifdef USE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
}
static inline int getRank()
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif
    return rank;
}
static inline int getSize()
{
    int size = 0;
#ifdef USE_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif
    return size;
}


// Function to waste CPU cycles
void waste_cpu( int N )
{
    if ( N > 10000 ) {
        PROFILE_START( "waste_cpu", 2 );
    }
    double pi = 3.141592653589793;
    double x  = 1.0;
    N         = std::max( 10, N );
    {
        for ( int i = 0; i < N; i++ )
            x = sqrt( x * exp( pi / x ) );
    } // style to limit gcov hits
    if ( fabs( x - 2.926064057273157 ) > 1e-12 ) {
        abort();
    }
    if ( N > 10000 ) {
        PROFILE_STOP( "waste_cpu", 2 );
    }
}


// Sleep for the given time
// Note: since we may encounter interrupts, we may not sleep for the desired time
//   so we need to perform the sleep in a loop
void sleep_ms( int64_t N )
{
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    while ( to_ms( t2 - t1 ) < N ) {
        int N2 = N - to_ms( t2 - t1 );
        std::this_thread::sleep_for( std::chrono::milliseconds( N2 ) );
        t2 = std::chrono::high_resolution_clock::now();
    }
}
void sleep_s( int N ) { sleep_ms( 1000 * N ); }


// Function to sleep for N seconds then increment a global count
static volatile int global_sleep_count = 0;
void sleep_inc( int N )
{
    PROFILE_START( "sleep_inc" );
    sleep_s( N );
    ++global_sleep_count;
    PROFILE_STOP( "sleep_inc" );
}
void sleep_inc2( double x )
{
    sleep_ms( static_cast<int>( round( x * 1000 ) ) );
    ++global_sleep_count;
}
void sleep_msg( double x, std::string msg )
{
    PROFILE_START( msg );
    sleep_ms( static_cast<int>( round( x * 1000 ) ) );
    NULL_USE( msg );
    PROFILE_STOP( msg );
}
bool check_inc( int N ) { return global_sleep_count == N; }


// Function to return the processor for the given thread
std::mutex print_processor_mutex;

void print_processor( ThreadPool *tpool )
{
    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif
    int thread    = tpool->getThreadNumber();
    int processor = ThreadPool::getCurrentProcessor();
    char tmp[100];
    sprintf( tmp, "%i:  Thread,proc = %i,%i\n", rank, thread, processor );
    sleep_ms( 10 * rank );
    print_processor_mutex.lock();
    pout << tmp;
    print_processor_mutex.unlock();
    sleep_ms( 100 );
}


// Function to test how a member thread interacts with the thread pool
int test_member_thread( ThreadPool *tpool )
{
    int N_errors = 0;
    // Member threads are not allowed to wait for the pool to finish
    try {
        tpool->wait_pool_finished();
        N_errors++;
    } catch ( ... ) {
    }
    // Member threads are not allowed to change the size of the pool
    try {
        tpool->wait_pool_finished();
        N_errors++;
    } catch ( ... ) {
    }
    return N_errors;
}


/******************************************************************
 * Test the TPOOL_ADD_WORK macro with variable number of arguments *
 ******************************************************************/
static int myfun0() { return 0; }
static int myfun1( int ) { return 1; }
static int myfun2( int, float ) { return 2; }
static int myfun3( int, float, double ) { return 3; }
static int myfun4( int, float, double, char ) { return 4; }
static int myfun5( int, float, double, char, std::string ) { return 5; }
static int myfun6( int, float, double, char, std::string, int ) { return 6; }
static int myfun7( int, float, double, char, std::string, int, int ) { return 7; }
static int test_function_arguements( ThreadPool *tpool )
{
    int N_errors = 0;
    // Test some basic types of instantiations
    ThreadPool::thread_id_t id0 = TPOOL_ADD_WORK( tpool, myfun0, ( nullptr ) );
    ThreadPool::thread_id_t id1 = TPOOL_ADD_WORK( tpool, myfun1, ( (int) 1 ) );
    ThreadPool::thread_id_t id2 = TPOOL_ADD_WORK( tpool, myfun2, ( (int) 1, (float) 2 ) );
    ThreadPool::thread_id_t id3 =
        TPOOL_ADD_WORK( tpool, myfun3, ( (int) 1, (float) 2, (double) 3 ) );
    ThreadPool::thread_id_t id4 =
        TPOOL_ADD_WORK( tpool, myfun4, ( (int) 1, (float) 2, (double) 3, (char) 4 ) );
    ThreadPool::thread_id_t id5 = TPOOL_ADD_WORK(
        tpool, myfun5, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ) ) );
    ThreadPool::thread_id_t id52 = TPOOL_ADD_WORK(
        tpool, myfun5, ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ) ), -1 );
    ThreadPool::thread_id_t id6 = TPOOL_ADD_WORK( tpool, myfun6,
        ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ), (int) 1 ) );
    ThreadPool::thread_id_t id7 = TPOOL_ADD_WORK( tpool, myfun7,
        ( (int) 1, (float) 2, (double) 3, (char) 4, std::string( "test" ), (int) 1, (int) 1 ) );
    tpool->wait_pool_finished();
    if ( !tpool->isFinished( id0 ) ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id0 ) != 0 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id1 ) != 1 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id2 ) != 2 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id3 ) != 3 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id4 ) != 4 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id5 ) != 5 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id52 ) != 5 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id6 ) != 6 ) {
        N_errors++;
    }
    if ( tpool->getFunctionRet<int>( id7 ) != 7 ) {
        N_errors++;
    }
    return N_errors;
}


/******************************************************************
 * Examples to derive a user work item                             *
 ******************************************************************/
class UserWorkItemVoid : public ThreadPool::WorkItem
{
public:
    // User defined constructor (does not need to match any interfaces)
    explicit UserWorkItemVoid( int dummy )
    {
        // User initialized variables
        NULL_USE( dummy );
    }
    // User defined run (can do anything)
    void run() override
    {
        // Perform the tasks
        printf( "Hello work from UserWorkItem (void)" );
    }
    // Will the routine return a result
    bool has_result() const override { return false; }
    // User defined destructor
    ~UserWorkItemVoid() override = default;
};
class UserWorkItemInt : public ThreadPool::WorkItemRet<int>
{
public:
    // User defined constructor (does not need to match any interfaces)
    explicit UserWorkItemInt( int dummy )
    {
        // User initialized variables
        NULL_USE( dummy );
    }
    // User defined run (can do anything)
    void run() override
    {
        // Perform the tasks
        printf( "Hello work from UserWorkItem (int)" );
        // Store the results (it's type will match the template)
        ThreadPool::WorkItemRet<int>::d_result = 1;
    }
    // User defined destructor
    ~UserWorkItemInt() override = default;
};


/******************************************************************
 * test the time to run N tasks in parallel                        *
 ******************************************************************/
template<class Ret, class... Args>
inline double launchAndTime( ThreadPool &tpool, int N, Ret ( *routine )( Args... ), Args... args )
{
    tpool.wait_pool_finished();
    auto start = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < N; i++ )
        ThreadPool_add_work( &tpool, 0, routine, args... );
    tpool.wait_pool_finished();
    auto stop = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>( stop - start ).count();
}


// Move constructor function
volatile ThreadPool::thread_id_t f1( volatile ThreadPool::thread_id_t a ) { return a; }
ThreadPool::thread_id_t f2( ThreadPool::thread_id_t a ) { return a; }


/******************************************************************
 * Test the basic functionallity of the atomics                    *
 ******************************************************************/
int test_atomics()
{
    using namespace AtomicOperations;
    int N_errors = 0;
    volatile int32_atomic i32;
    volatile int64_atomic i64;
    i32 = 32;
    i64 = 64;
    if ( atomic_increment( &i32 ) != 33 || atomic_increment( &i64 ) != 65 )
        N_errors++;
    if ( atomic_decrement( &i32 ) != 32 || atomic_decrement( &i64 ) != 64 )
        N_errors++;
    if ( atomic_add( &i32, 2 ) != 34 || atomic_add( &i64, 4 ) != 68 )
        N_errors++;
    if ( atomic_compare_and_swap( &i32, 0, 0 ) || atomic_compare_and_swap( &i64, 0, 0 ) )
        N_errors++;
    if ( !atomic_compare_and_swap( &i32, 34, 32 ) || !atomic_compare_and_swap( &i64, 68, 64 ) )
        N_errors++;
    if ( i32 != 32 || i64 != 64 )
        N_errors++;
    return N_errors;
}


/******************************************************************
 * Test FIFO behavior                                              *
 ******************************************************************/
void test_FIFO( UnitTest &ut, ThreadPool &tpool )
{
    int rank    = getRank();
    int size    = getSize();
    const int N = 4000;
    for ( int r = 0; r < size; r++ ) {
        barrier();
        if ( r != rank )
            continue;
        std::vector<ThreadPool::thread_id_t> ids;
        ids.reserve( N );
        for ( size_t i = 0; i < N; i++ )
            ids.emplace_back( TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.001 ) ) );
        bool pass = true;
        while ( tpool.N_queued() > 0 ) {
            int i1 = -1, i2 = ids.size();
            for ( int i = N - 1; i >= 0; i-- ) {
                bool started = ids[i].started();
                if ( started )
                    i1 = std::max<int>( i1, i ); // Last index to processing item
                else
                    i2 = std::min<int>( i2, i ); // First index to queued item
            }
            int diff = i1 == -1 ? 0 : ( i2 - i1 - 1 );
            if ( abs( diff ) > 4 ) {
                printf( "%i %i %i\n", i1, i2, diff );
                pass = pass && abs( i2 - i1 - 1 ) <= 2;
            }
        }
        ids.clear();
        tpool.wait_pool_finished();
        if ( pass )
            ut.passes( "Thread pool behaves as FIFO" );
        else
            ut.failure( "Thread pool does not behave as FIFO" );
    }
}


/******************************************************************
 * The main program                                                *
 ******************************************************************/
#ifdef USE_WINDOWS
int __cdecl main( int argc, char **argv )
{
#elif defined( USE_LINUX ) || defined( USE_MAC )
int main( int argc, char *argv[] )
{
#else
#error Unknown OS
#endif

    int N_threads = 4;    // Number of threads
    int N_work    = 2000; // Number of work items
    int N_it      = 10;   // Number of cycles to run
    int N_problem = 5;    // Problem size
    PROFILE_ENABLE( 3 );
    PROFILE_ENABLE_TRACE();
    PROFILE_DISABLE_MEMORY();
    UnitTest ut;


    // Initialize MPI and set the error handlers
#ifdef USE_MPI
    int provided_thread_support = -1;
    MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided_thread_support );
    Utilities::setErrorHandlers();
    // Disable OS specific warnings for all non-root ranks
#endif
    int rank = getRank();
    int size = getSize();
    if ( rank > 0 )
        ThreadPool::set_OS_warnings( 1 );
    NULL_USE( size );
    NULL_USE( argc );
    NULL_USE( argv );


    // Test the atomics
    if ( test_atomics() == 0 )
        ut.passes( "Atomics passed" );
    else
        ut.failure( "Atomics failed" );

    // Initialize the data
    std::vector<int> data1( N_work, 0 );
    std::vector<int> priority( N_work, 0 );
    for ( int i = 0; i < N_work; i++ ) {
        data1[i]    = N_problem;
        priority[i] = i % 128;
    }


    // Print the size of the thread pool class
    printp( "Size of ThreadPool = %i\n", (int) sizeof( ThreadPool ) );


    // Get the number of processors availible
    barrier();
    int N_procs = ThreadPool::getNumberOfProcessors();
    if ( N_procs > 0 )
        ut.passes( "getNumberOfProcessors" );
    else
        ut.failure( "getNumberOfProcessors" );
    printp( "%i processors availible\n", N_procs );


    // Get the processor affinities for the process
    barrier();
    std::vector<int> cpus = ThreadPool::getProcessAffinity();
    printp( "%i cpus for current process: ", (int) cpus.size() );
    for ( int cpu : cpus )
        printp( "%i ", cpu );
    printp( "\n" );
    if ( !cpus.empty() ) {
        ut.passes( "getProcessAffinity" );
    } else {
#ifdef __APPLE__
        ut.expected_failure( "getProcessAffinity" );
#else
        ut.failure( "getProcessAffinity" );
#endif
    }


    // Test setting the process affinities
    barrier();
    bool pass = false;
    if ( !cpus.empty() && N_procs > 0 ) {
        if ( cpus.size() == 1 ) {
            cpus.resize( N_procs );
            for ( int i = 0; i < N_procs; i++ )
                cpus.push_back( i );
            try {
                ThreadPool::setProcessAffinity( cpus );
            } catch ( ... ) {
            }
            cpus                  = ThreadPool::getProcessAffinity();
            std::vector<int> cpus = ThreadPool::getProcessAffinity();
            printp( "%i cpus for current process (updated): ", (int) cpus.size() );
            for ( int cpu : cpus )
                printp( "%i ", cpu );
            printp( "\n" );
            pass = cpus.size() > 1;
        } else {
            std::vector<int> cpus_orig = cpus;
            std::vector<int> cpus_tmp( 1, cpus[0] );
            try {
                ThreadPool::setProcessAffinity( cpus_tmp );
            } catch ( ... ) {
            }
            cpus = ThreadPool::getProcessAffinity();
            if ( cpus.size() == 1 )
                pass = true;
            try {
                ThreadPool::setProcessAffinity( cpus_orig );
            } catch ( ... ) {
            }
            cpus = ThreadPool::getProcessAffinity();
            if ( cpus.size() != cpus_orig.size() )
                pass = false;
        }
    }
    if ( pass ) {
        ut.passes( "setProcessAffinity" );
    } else {
#ifdef __APPLE__
        ut.expected_failure( "setProcessAffinity" );
#else
        ut.failure( "setProcessAffinity" );
#endif
    }
    int N_procs_used = std::min<int>( N_procs, N_threads );
    printp( "%i processors used\n", N_procs_used );


    // Create the thread pool
    barrier();
    printp( "Creating thread pool\n" );
    ThreadPool tpool0;
    ThreadPool tpool;
    ThreadPool::thread_id_t id;
    id = TPOOL_ADD_WORK( &tpool, waste_cpu, ( data1[0] ) );
    if ( id == ThreadPool::thread_id_t() || !tpool.isValid( id ) )
        ut.failure( "Errors with id" );
    tpool.setNumThreads( N_threads );
    if ( tpool.getNumThreads() == N_threads )
        ut.passes( "Created thread pool" );
    else
        ut.failure( "Failed to create tpool with desired number of threads" );


    // Test setting the thread affinities
    barrier();
    if ( cpus.size() > 1 ) {
        sleep_ms( 50 );
        // First make sure we can get the thread affinities
        std::vector<int> procs = ThreadPool::getThreadAffinity();
        if ( procs == cpus ) {
            ut.passes( "getThreadAffinity() matches procs" );
        } else {
            char msg[100];
            sprintf( msg, "getThreadAffinity() does not match procs (%i,%i)",
                static_cast<int>( procs.size() ), static_cast<int>( cpus.size() ) );
            ut.failure( msg );
        }
        pass = true;
        for ( int i = 0; i < N_threads; i++ ) {
            std::vector<int> procs_thread = tpool.getThreadAffinity( i );
            if ( procs_thread != procs ) {
                printp( "%i: Initial thread affinity: ", rank );
                for ( int i : procs_thread )
                    printp( "%i ", i );
                printp( "\n" );
                pass = false;
            }
        }
        if ( pass )
            ut.passes( "getThreadAffinity(thread) matches procs" );
        else
            ut.failure( "getThreadAffinity(thread) does not match procs" );
        // Try to set the thread affinities
        pass = true;
        if ( !procs.empty() ) {
            int N_procs_thread = std::max<int>( (int) cpus.size() / N_threads, 1 );
            for ( int i = 0; i < N_threads; i++ ) {
                std::vector<int> procs_thread( N_procs_thread, -1 );
                for ( int j = 0; j < N_procs_thread; j++ )
                    procs_thread[j] = procs[( i * N_procs_thread + j ) % procs.size()];
                tpool.setThreadAffinity( i, procs_thread );
                sleep_ms( 10 ); // Give time for OS to update thread affinities
                std::vector<int> procs_thread2 = tpool.getThreadAffinity( i );
                if ( procs_thread2 != procs_thread ) {
                    printp( "%i: Final thread affinity: ", rank );
                    for ( int i : procs_thread )
                        printp( "%i ", i );
                    printp( "\n" );
                    pass = false;
                }
            }
        }
        if ( pass )
            ut.passes( "setThreadAffinity passes" );
        else
            ut.failure( "setThreadAffinity failed to change affinity" );
    }


    // Reset the thread affinities
    barrier();
    tpool.setNumThreads( tpool.getNumThreads(), "none" );
    // tpool.setNumThreads(tpool.getNumThreads(),"independent");
    for ( int i = 0; i < N_threads; i++ ) {
        std::vector<int> procs_thread = tpool.getThreadAffinity( i );
        printp( "Thread affinity: " );
        for ( int i : procs_thread )
            printp( "%i ", i );
        printp( "\n" );
    }

    // Print the current processors by thread id
    barrier();
    ThreadPool::set_OS_warnings( 1 );
    print_processor( &tpool );
    launchAndTime( tpool, N_threads, print_processor, &tpool );

    // Run some basic tests
    barrier();
    auto start = std::chrono::high_resolution_clock::now();
    for ( int n = 0; n < N_it; n++ ) {
        for ( int i = 0; i < N_work; i++ )
            waste_cpu( data1[i] );
    }
    auto stop   = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>( stop - start ).count();
    printp( "Time for serial cycle = %0.0f us\n", 1e6 * time / N_it );
    printp( "Time for serial item = %0.0f ns\n", 1e9 * time / ( N_it * N_work ) );
    id = TPOOL_ADD_WORK( &tpool, waste_cpu, ( data1[0] ) );
    tpool.wait( id );
    std::vector<ThreadPool::thread_id_t> ids2;
    ids2.push_back( TPOOL_ADD_WORK( &tpool, waste_cpu, ( data1[0] ) ) );
    tpool.wait( ids2[0] );

    // Test the move operator for thread_id
    ThreadPool::thread_id_t id1          = f1( id );         // move-construct from rvalue temporary
    ThreadPool::thread_id_t id2          = std::move( id1 ); // move-construct from xvalue
    volatile ThreadPool::thread_id_t id3 = f2( id );         // move-construct from rvalue temporary
    volatile ThreadPool::thread_id_t id4 = std::move( id3 ); // move-construct from xvalue
    id2.reset();
    id4.reset();

    // Test calling functions with different number of arguments
    barrier();
    printp( "Testing arguments:\n" );
    int N_errors_args = test_function_arguements( &tpool );
    if ( N_errors_args == 0 )
        ut.passes( "Calling function with default arguments" );
    else
        ut.failure( "Error calling function with default arguments" );


    // Check that the threads can sleep in parallel (this does not depend on the number of
    // processors)
    barrier();
    tpool.wait_pool_finished();
    start = std::chrono::high_resolution_clock::now();
    sleep_inc( 1 );
    stop                  = std::chrono::high_resolution_clock::now();
    double sleep_serial   = std::chrono::duration<double>( stop - start ).count();
    double sleep_parallel = launchAndTime( tpool, N_threads, sleep_inc, 1 );
    double sleep_speedup  = N_procs_used * sleep_serial / sleep_parallel;
    printf( "%i:  Speedup on %i sleeping threads: %0.3f\n", rank, N_procs_used, sleep_speedup );
    printf( "%i:    ts = %0.3f, tp = %0.3f\n", rank, sleep_serial, sleep_parallel );
    if ( fabs( sleep_serial - 1.0 ) < 0.05 && fabs( sleep_parallel - 1.0 ) < 0.25 &&
         sleep_speedup > 3 )
        ut.passes( "Passed thread sleep" );
    else
        ut.failure( "Failed thread sleep" );


    // Check that the threads are actually working in parallel
    barrier();
    if ( N_procs_used > 1 ) {
#ifdef USE_MPI
        // Use a non-blocking serialization of the MPI processes
        // if we do not have a sufficient number of processors
        bool serialize_mpi = N_procs < N_threads * size;
        int buf;
        MPI_Request request;
        MPI_Status status;
        if ( serialize_mpi && rank > 0 ) {
            MPI_Irecv( &buf, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &request );
            int flag = false;
            while ( !flag ) {
                MPI_Test( &request, &flag, &status );
                sleep_s( 1 );
            }
        }
#endif
        int N = 20000000; // Enough work to keep the processor busy for ~ 1 s
        // Run in serial
        start = std::chrono::high_resolution_clock::now();
        waste_cpu( N );
        stop               = std::chrono::high_resolution_clock::now();
        double time_serial = std::chrono::duration<double>( stop - start ).count();
        // Run in parallel
        double time_parallel  = launchAndTime( tpool, N_procs_used, waste_cpu, N );
        double time_parallel2 = launchAndTime( tpool, N_procs_used, waste_cpu, N / 1000 );
        double speedup        = N_procs_used * time_serial / time_parallel;
        printf( "%i:  Speedup on %i procs: %0.3f\n", rank, N_procs_used, speedup );
        printf( "%i:    ts = %0.3f, tp = %0.3f, tp2 = %0.3f\n", rank, time_serial, time_parallel,
            time_parallel2 );
        if ( speedup > 1.4 ) {
            ut.passes( "Passed speedup test" );
        } else {
#ifdef USE_GCOV
            ut.expected_failure( "Times do not indicate tests are running in parallel (gcov)" );
#else
            ut.failure( "Times do not indicate tests are running in parallel" );
#endif
        }
#ifdef USE_MPI
        if ( serialize_mpi ) {
            if ( rank < size - 1 )
                MPI_Send( &N, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD );
            if ( rank == size - 1 ) {
                for ( int i = 0; i < size - 1; i++ )
                    MPI_Send( &N, 1, MPI_INT, i, 1, MPI_COMM_WORLD );
            } else {
                MPI_Irecv( &buf, 1, MPI_INT, size - 1, 1, MPI_COMM_WORLD, &request );
                int flag = false;
                MPI_Status status;
                while ( !flag ) {
                    MPI_Test( &request, &flag, &status );
                    sleep_s( 1 );
                }
            }
        }
#endif
    } else {
        ut.expected_failure( "Testing thread performance with less than 1 processor" );
    }


    // Test first-in-first-out scheduler (also ensures priorities)
    test_FIFO( ut, tpool );


    // Test adding a work item with a dependency
    barrier();
    {
        // Test that we sucessfully wait on the work items
        std::vector<ThreadPool::thread_id_t> ids;
        ids.reserve( 5 );
        global_sleep_count = 0; // Reset the count before this test
        ThreadPool::thread_id_t id0;
        auto id1    = TPOOL_ADD_WORK( &tpool, sleep_inc, ( 1 ) );
        auto id2    = TPOOL_ADD_WORK( &tpool, sleep_inc, ( 2 ) );
        auto *wait1 = new WorkItemFull<bool, int>( check_inc, 1 );
        auto *wait2 = new WorkItemFull<bool, int>( check_inc, 2 );
        wait1->add_dependency( id0 );
        wait1->add_dependency( id1 );
        wait2->add_dependency( id1 );
        wait2->add_dependency( id2 );
        ids.clear();
        ids.push_back( tpool.add_work( wait1 ) );
        ids.push_back( tpool.add_work( wait2 ) );
        tpool.wait_all( ids.size(), &ids[0] );
        if ( !tpool.getFunctionRet<bool>( ids[0] ) || !tpool.getFunctionRet<bool>( ids[1] ) )
            ut.failure( "Failed to wait on required dependency" );
        else
            ut.passes( "Dependencies" );
        tpool.wait_pool_finished();
        // Test waiting on more dependencies than in the thread pool (changing priorities)
        ids.clear();
        for ( size_t i = 0; i < 20; i++ )
            ids.push_back( TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.1 ) ) );
        auto *wait3 = new WorkItemFull<void, double>( sleep_inc2, 0 );
        wait3->add_dependencies( ids );
        id = tpool.add_work( wait3, 50 );
        tpool.wait( id );
        bool pass = true;
        for ( auto &id : ids )
            pass = pass && id.finished();
        ids.clear();
        if ( pass )
            ut.passes( "Dependencies2" );
        else
            ut.failure( "Dependencies2" );
        // Check that we can handle more complex dependencies
        id1 = TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.5 ) );
        for ( int i = 0; i < 10; i++ ) {
            wait1 = new WorkItemFull<bool, int>( check_inc, 1 );
            wait1->add_dependency( id1 );
            tpool.add_work( wait1 );
        }
        tpool.wait_pool_finished();
        ids.clear();
        for ( int i = 0; i < 5; i++ )
            ids.push_back( TPOOL_ADD_WORK( &tpool, sleep_inc2, ( 0.5 ) ) );
        sleep_inc2( 0.002 );
        ThreadPool::WorkItem *work = new WorkItemFull<void, int>( waste_cpu, 100 );
        work->add_dependencies( ids );
        id = tpool.add_work( work, 10 );
        tpool.wait( id );
    }

    // Test the timing creating and running a work item
    barrier();
    {
        printp( "Testing timmings (creating/running work item):\n" );
        std::string timer_name = "Create/Run work item";
        PROFILE_START( timer_name );
        int64_t time_create = 0;
        int64_t time_run    = 0;
        int64_t time_delete = 0;
        std::vector<ThreadPool::WorkItem *> work( N_work );
        start = std::chrono::high_resolution_clock::now();
        for ( int n = 0; n < N_it; n++ ) {
            auto t1 = std::chrono::high_resolution_clock::now();
            for ( int i = 0; i < N_work; i++ )
                work[i] = ThreadPool::createWork<void, int>( waste_cpu, data1[i] );
            auto t2 = std::chrono::high_resolution_clock::now();
            for ( int i = 0; i < N_work; i++ )
                work[i]->run();
            auto t3 = std::chrono::high_resolution_clock::now();
            for ( int i = 0; i < N_work; i++ )
                delete work[i];
            auto t4 = std::chrono::high_resolution_clock::now();
            time_create += to_ns( t2 - t1 );
            time_run += to_ns( t3 - t2 );
            time_delete += to_ns( t4 - t3 );
            if ( ( n + 1 ) % 100 == 0 )
                printp( "Cycle %i of %i finished\n", n + 1, N_it );
        }
        stop = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration<double>( stop - start ).count();
        PROFILE_STOP( timer_name );
        printp( "   time = %0.0f ms\n", 1e3 * time );
        printp( "   time / cycle = %0.0f us\n", 1e6 * time / N_it );
        printp( "   average time / item = %0.0f ns\n", 1e9 * time / ( N_it * N_work ) );
        printp( "      create = %i ns\n", time_create / ( N_it * N_work ) );
        printp( "      run    = %i ns\n", time_run / ( N_it * N_work ) );
        printp( "      delete = %i us\n", time_delete / ( N_it * N_work ) );
    }

    // Test the timing adding a single item
    barrier();
    for ( int it = 0; it < 2; it++ ) {
        ThreadPool *tpool_ptr = nullptr;
        std::string timer_name;
        if ( it == 0 ) {
            printp( "Testing timmings (adding a single item to empty tpool):\n" );
            timer_name = "Add single item to empty pool";
            tpool_ptr  = &tpool0;
        } else if ( it == 1 ) {
            printp( "Testing timmings (adding a single item):\n" );
            timer_name = "Add single item to tpool";
            tpool_ptr  = &tpool;
        }
        PROFILE_START( timer_name );
        std::vector<ThreadPool::thread_id_t> ids( N_work );
        int64_t time_add  = 0;
        int64_t time_wait = 0;
        start             = std::chrono::high_resolution_clock::now();
        for ( int n = 0; n < N_it; n++ ) {
            auto t1 = std::chrono::high_resolution_clock::now();
            for ( int i = 0; i < N_work; i++ )
                ids[i] = TPOOL_ADD_WORK( tpool_ptr, waste_cpu, ( data1[i] ), priority[i] );
            auto t2 = std::chrono::high_resolution_clock::now();
            tpool_ptr->wait_all( N_work, &ids[0] );
            auto t3 = std::chrono::high_resolution_clock::now();
            time_add += to_ns( t2 - t1 );
            time_wait += to_ns( t3 - t2 );
            if ( ( n + 1 ) % 100 == 0 )
                printp( "Cycle %i of %i finished\n", n + 1, N_it );
        }
        stop = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration<double>( stop - start ).count();
        PROFILE_STOP( timer_name );
        printp( "   time = %0.0f ms\n", 1e3 * time );
        printp( "   time / cycle = %0.0f us\n", 1e6 * time / N_it );
        printp( "   average time / item = %0.0f ns\n", 1e9 * time / ( N_it * N_work ) );
        printp( "      create and add = %i ns\n", time_add / ( N_it * N_work ) );
        printp( "      wait = %i us\n", time_wait / ( N_it * N_work ) );
    }

    // Test the timing pre-creating the work items and adding multiple at a time
    barrier();
    for ( int it = 0; it < 2; it++ ) {
        ThreadPool *tpool_ptr = nullptr;
        std::string timer_name;
        if ( it == 0 ) {
            printp( "Testing timmings (adding a block of items to empty tpool):\n" );
            timer_name = "Add multiple items to empty pool";
            tpool_ptr  = &tpool0;
        } else if ( it == 1 ) {
            printp( "Testing timmings (adding a block of items):\n" );
            timer_name = "Add multiple items to tpool";
            tpool_ptr  = &tpool;
        }
        PROFILE_START( timer_name );
        int64_t time_create_work = 0;
        int64_t time_add_work    = 0;
        int64_t time_wait_work   = 0;
        std::vector<ThreadPool::WorkItem *> work( N_work );
        start = std::chrono::high_resolution_clock::now();
        for ( int n = 0; n < N_it; n++ ) {
            auto t1 = std::chrono::high_resolution_clock::now();
            for ( int i = 0; i < N_work; i++ )
                work[i] = ThreadPool::createWork<void, int>( waste_cpu, data1[i] );
            auto t2  = std::chrono::high_resolution_clock::now();
            auto ids = tpool_ptr->add_work( work, priority );
            auto t3  = std::chrono::high_resolution_clock::now();
            tpool_ptr->wait_all( ids );
            auto t4 = std::chrono::high_resolution_clock::now();
            time_create_work += to_ns( t2 - t1 );
            time_add_work += to_ns( t3 - t2 );
            time_wait_work += to_ns( t4 - t3 );
            if ( ( n + 1 ) % 100 == 0 )
                printp( "Cycle %i of %i finished\n", n + 1, N_it );
        }
        stop = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration<double>( stop - start ).count();
        PROFILE_STOP( timer_name );
        printp( "   time = %0.0f ms\n", 1e3 * time );
        printp( "   time / cycle = %0.0f us\n", 1e6 * time / N_it );
        printp( "   average time / item = %0.0f ns\n", 1e9 * time / ( N_it * N_work ) );
        printp( "      create = %i ns\n", time_create_work / ( N_it * N_work ) );
        printp( "      add = %i ns\n", time_add_work / ( N_it * N_work ) );
        printp( "      wait = %i ns\n", time_wait_work / ( N_it * N_work ) );
    }

    // Run a dependency test that tests a simple case that should keep the thread pool busy
    // Note: Checking the results requires looking at the trace data
    tpool.wait_pool_finished();
    PROFILE_START( "Dependency test" );
    for ( int i = 0; i < 10; i++ ) {
        char msg[3][100];
        sprintf( msg[0], "Item %i-%i", i, 0 );
        sprintf( msg[1], "Item %i-%i", i, 1 );
        sprintf( msg[2], "Item %i-%i", i, 2 );
        ThreadPool::WorkItem *work =
            new WorkItemFull<void, double, std::string>( sleep_msg, 0.5, msg[0] );
        ThreadPool::WorkItem *work1 =
            new WorkItemFull<void, double, std::string>( sleep_msg, 0.1, msg[1] );
        ThreadPool::WorkItem *work2 =
            new WorkItemFull<void, double, std::string>( sleep_msg, 0.1, msg[2] );
        ThreadPool::thread_id_t id = tpool.add_work( work );
        work1->add_dependency( id );
        work2->add_dependency( id );
        tpool.add_work( work1 );
        tpool.add_work( work2 );
    }
    tpool.wait_pool_finished();
    PROFILE_STOP( "Dependency test" );

    // Close the thread pool
    tpool.setNumThreads( 0 );

    // Save the profiling results
    PROFILE_SAVE( "test_thread_pool" );
    PROFILE_DISABLE();

    // Test creating/destroying a thread pool using new
    barrier();
    pass = true;
    try {
        ThreadPool *tpool = new ThreadPool( ThreadPool::MAX_NUM_THREADS - 1 );
        if ( tpool->getNumThreads() != ThreadPool::MAX_NUM_THREADS - 1 )
            pass = false;
        if ( !ThreadPool::is_valid( tpool ) )
            pass = false;
        delete tpool;
        // Check that tpool is invalid
        // Note: valgrind will report this as an invalid memory read, but we want to keep the test)
        if ( ThreadPool::is_valid( tpool ) )
            pass = false;
    } catch ( ... ) {
        pass = false;
    }
    if ( pass )
        ut.passes( "Created/destroyed thread pool with new" );
    else
        ut.failure( "Created/destroyed thread pool with new" );

    // Print the test results
    barrier();
    ut.report();
    auto N_errors = static_cast<int>( ut.NumFailGlobal() );

    // Shudown MPI
    pout << "Shutting down\n";
    barrier();
#ifdef USE_TIMER
    if ( rank == 0 )
        MemoryApp::print( pout );
#endif
#ifdef USE_MPI
    MPI_Finalize();
    sleep_ms( 10 );
#endif
    return N_errors;
}
