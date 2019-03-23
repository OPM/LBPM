#include "threadpool/thread_pool.h"
#include "common/Utilities.h"
#include "StackTrace/StackTrace.h"
#include "StackTrace/Utilities.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <bitset>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>
#include <typeinfo>


// Add profile timers or performance counters to the threadpool
#define PROFILE_THREADPOOL_PERFORMANCE 0
#define MONITOR_THREADPOOL_PERFORMANCE 0


#define perr std::cerr
#define pout std::cout
#define printp printf


// OS specific includes / definitions
// clang-format off
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
    #define USE_WINDOWS
#elif defined( __APPLE__ )
    #define USE_MAC
#elif defined( __linux ) || defined( __linux__ ) || defined( __unix ) || defined( __posix )
    #define USE_LINUX
#else
    #error Unknown OS
#endif
#if defined( USE_WINDOWS )
    #include <process.h>
    #include <windows.h>
    // Disable warning: the inline specifier cannot be used when a friend
    // declaration refers to a specialization of a function template
    #pragma warning( disable : 4396 )
#endif
#if defined(USE_LINUX) || defined(USE_MAC)
    #include <pthread.h>
    #include <unistd.h>
#endif
#ifdef USE_MAC
    // https://developer.apple.com/library/mac/#releasenotes/Performance/RN-AffinityAPI
    // http://plugins.svn.wordpress.org/wp-xhprof-profiler/trunk/facebook-xhprof/extension/xhprof..c
    #include <mach/mach_init.h>
    #include <mach/thread_policy.h>
    #define cpu_set_t thread_affinity_policy_data_t
    #define CPU_SET( cpu_id, new_mask ) *new_mask.affinity_tag = ( cpu_id + 1 )
    #define CPU_ZERO( new_mask ) ( *( new_mask ) ).affinity_tag = THREAD_AFFINITY_TAG_NULL
    #define sched_setaffinity( pid, size, mask ) \
        thread_policy_set(                       \
            mach_thread_self(), THREAD_AFFINITY_POLICY, mask, THREAD_AFFINITY_POLICY_COUNT )
    #define sched_getaffinity( pid, size, mask ) \
        thread_policy_get(                       \
            mach_thread_self(), THREAD_AFFINITY_POLICY, mask, THREAD_AFFINITY_POLICY_COUNT )
#endif
// clang-format on


// Set some macros
// clang-format off
#if PROFILE_THREADPOOL_PERFORMANCE == 1
#define PROFILE_THREADPOOL_START(X)  PROFILE_START(X,3)
#define PROFILE_THREADPOOL_START2(X) PROFILE_START2(X,3)
#define PROFILE_THREADPOOL_STOP(X)   PROFILE_STOP(X,3)
#define PROFILE_THREADPOOL_STOP2(X)  PROFILE_STOP2(X,3)
#else
#define PROFILE_THREADPOOL_START(X)  do {} while ( 0 )
#define PROFILE_THREADPOOL_START2(X) do {} while ( 0 )
#define PROFILE_THREADPOOL_STOP(X)   do {} while ( 0 )
#define PROFILE_THREADPOOL_STOP2(X)  do {} while ( 0 )
#endif
#if MONITOR_THREADPOOL_PERFORMANCE == 1
#define accumulate(x,t1,t2)  AtomicOperations::atomic_add( \
    &x, std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() );
#endif
// clang-format on


#if MONITOR_THREADPOOL_PERFORMANCE == 1
static AtomicOperations::int64_atomic total_add_work_time[5] = { 0, 0, 0, 0, 0 };
#endif


// Set env
static std::mutex Utilities_mutex;
void setenv( const std::string &name, const std::string &value )
{
    Utilities_mutex.lock();
#if defined( USE_LINUX ) || defined( USE_MAC )
    bool pass = false;
    if ( value.empty() )
        pass = ::setenv( name.data(), value.data(), 1 ) == 0;
    else
        pass = ::unsetenv( name.data() ) == 0;
#elif defined( USE_WINDOWS )
    bool pass = SetEnvironmentVariable( name.data(), value.data() ) != 0;
#else
#error Unknown OS
#endif
    Utilities_mutex.unlock();
    if ( !pass ) {
        char msg[1024];
        if ( !value.empty() )
            sprintf(
                msg, "Error setting enviornmental variable: %s=%s\n", name.data(), value.data() );
        else
            sprintf( msg, "Error clearing enviornmental variable: %s\n", name.data() );
        ERROR( msg );
    }
}


// Helper functions
template<class T>
void quicksort( int N, T *data );
template<class T>
inline void quicksort( std::vector<T> &x )
{
    quicksort( x.size(), x.data() );
}
static inline int find_id( int, const ThreadPool::thread_id_t *, const ThreadPool::thread_id_t & );


// Function to generate a random number for checking if tpool is valid
static inline bool validHeadTail( uint32_t key )
{
    return ( key > 10 ) && ( ~key > 10 ) && ( key % 2 != 0 ) && ( key % 3 == 2 );
}
static inline uint32_t generateHeadTail()
{
    uint32_t key = 0;
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_int_distribution<> dis( 1, 0xFFFFFF );
    while ( !validHeadTail( key ) )
        key = static_cast<uint32_t>( dis( gen ) ) * 0x9E3779B9; // 2^32*0.5*(sqrt(5)-1)
    return key;
}


/******************************************************************
 * Run some basic compile-time checks                              *
 ******************************************************************/
static_assert( ThreadPool::MAX_THREADS % 64 == 0, "MAX_THREADS must be a multiple of 64" );
static_assert( ThreadPool::MAX_THREADS < 65535, "MAX_THREADS must < 65535" );
static_assert( sizeof( AtomicOperations::int32_atomic ) == 4, "atomic32 must be a 32-bit integer" );
static_assert( sizeof( AtomicOperations::int64_atomic ) == 8, "atomic64 must be a 64-bit integer" );


/******************************************************************
 * Get/Set a bit                                                   *
 * Note: these functions are thread-safe                           *
 ******************************************************************/
static inline void set_bit( volatile AtomicOperations::int64_atomic *x, size_t index )
{
    uint64_t mask = 0x01;
    mask <<= index % 64;
    size_t i  = index / 64;
    bool test = false;
    while ( !test ) {
        AtomicOperations::int64_atomic y = x[i];
        test = AtomicOperations::atomic_compare_and_swap( &x[i], y, ( y | mask ) );
    }
}
static inline void unset_bit( volatile AtomicOperations::int64_atomic *x, size_t index )
{
    uint64_t mask = 0x01;
    mask <<= index % 64;
    mask      = ~mask;
    size_t i  = index / 64;
    bool test = false;
    while ( !test ) {
        AtomicOperations::int64_atomic y = x[i];
        test = AtomicOperations::atomic_compare_and_swap( &x[i], y, ( y & mask ) );
    }
}
static inline bool get_bit( const volatile AtomicOperations::int64_atomic *x, size_t index )
{
    uint64_t mask = 0x01;
    mask <<= index % 64;
    // This is thread-safe since we only care about a single bit
    AtomicOperations::int64_atomic y = x[index / 64];
    return ( y & mask ) != 0;
}


/******************************************************************
 * Simple function to check if the parity is odd (true) or even    *
 ******************************************************************/
static inline bool is_odd8( size_t x )
{ // This only works for 64-bit integers
    x ^= ( x >> 1 );
    x ^= ( x >> 2 );
    x ^= ( x >> 4 );
    x ^= ( x >> 8 );
    x ^= ( x >> 16 );
    x ^= ( x >> 32 );
    return ( x & 0x01 ) > 0;
}
template<class int_type>
static inline int count_bits( int_type x )
{
    int count = 0;
    for ( size_t i = 0; i < 8 * sizeof( int_type ); i++ ) {
        if ( ( x >> i ) & 0x01 )
            ++count;
    }
    return count;
}


/******************************************************************
 * Set the global constants                                        *
 ******************************************************************/
constexpr uint16_t ThreadPool::MAX_THREADS;
constexpr uint16_t ThreadPool::MAX_WAIT;


/******************************************************************
 * Set the behavior of OS warnings                                 *
 ******************************************************************/
static int global_OS_behavior = 0;
static std::mutex OS_warning_mutex;
void ThreadPool::set_OS_warnings( int behavior )
{
    ASSERT( behavior >= 0 && behavior <= 2 );
    global_OS_behavior = behavior;
}
static void OS_warning( const std::string &message )
{
    OS_warning_mutex.lock();
    if ( global_OS_behavior == 0 ) {
        pout << "Warning: " << message << std::endl;
    } else if ( global_OS_behavior == 2 ) {
        perr << "Error: " << message << std::endl;
    }
    OS_warning_mutex.unlock();
}
void ThreadPool::setErrorHandler( std::function<void( const std::string & )> fun )
{
    d_errorHandler = fun;
}

/******************************************************************
 * Function to return the number of processors availible           *
 ******************************************************************/
int ThreadPool::getNumberOfProcessors() { return std::thread::hardware_concurrency(); }


/******************************************************************
 * Function to return the processor number of the current thread   *
 ******************************************************************/
int ThreadPool::getCurrentProcessor()
{
#if defined( USE_LINUX )
    return sched_getcpu() + 1;
#elif defined( USE_MAC )
    OS_warning( "MAC does not support getCurrentProcessor" );
    return 0;
#elif defined( USE_WINDOWS )
    return GetCurrentProcessorNumber() + 1;
#else
#error Unknown OS
#endif
}


/******************************************************************
 * Function to get/set the affinity of the current process         *
 ******************************************************************/
std::vector<int> ThreadPool::getProcessAffinity()
{
    std::vector<int> procs;
#ifdef USE_LINUX
#ifdef _GNU_SOURCE
    cpu_set_t mask;
    int error = sched_getaffinity( getpid(), sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        throw std::logic_error( "Error getting process affinity" );
    for ( size_t i = 0; i < sizeof( cpu_set_t ) * CHAR_BIT; i++ ) {
        if ( CPU_ISSET( i, &mask ) )
            procs.push_back( i );
    }
#else
#warning sched_getaffinity is not supported for this compiler/OS
    OS_warning( "sched_getaffinity is not supported for this compiler/OS" );
#endif
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    OS_warning( "MAC does not support getting the process affinity" );
#elif defined( USE_WINDOWS )
    HANDLE hProc = GetCurrentProcess();
    size_t procMask;
    size_t sysMask;
    PDWORD_PTR procMaskPtr = reinterpret_cast<PDWORD_PTR>( &procMask );
    PDWORD_PTR sysMaskPtr  = reinterpret_cast<PDWORD_PTR>( &sysMask );
    GetProcessAffinityMask( hProc, procMaskPtr, sysMaskPtr );
    for ( size_t i = 0; i < sizeof( size_t ) * CHAR_BIT; i++ ) {
        if ( ( procMask & 0x1 ) != 0 )
            procs.push_back( i );
        procMask >>= 1;
    }
#else
#error Unknown OS
#endif
    return procs;
}
void ThreadPool::setProcessAffinity( const std::vector<int> &procs )
{
#ifdef USE_LINUX
#ifdef _GNU_SOURCE
    cpu_set_t mask;
    CPU_ZERO( &mask );
    for ( size_t i = 0; i < procs.size(); i++ )
        CPU_SET( procs[i], &mask );
    int error = sched_setaffinity( getpid(), sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        throw std::logic_error( "Error setting process affinity" );
#else
#warning sched_setaffinity is not supported for this compiler/OS
    OS_warning( "sched_setaffinity is not supported for this compiler/OS" );
#endif
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    OS_warning( "MAC does not support setting the process affinity" );
#elif defined( USE_WINDOWS )
    DWORD mask = 0;
    for ( size_t i = 0; i < procs.size(); i++ )
        mask |= ( (DWORD) 1 ) << procs[i];
    HANDLE hProc = GetCurrentProcess();
    SetProcessAffinityMask( hProc, mask );
#else
#error Unknown OS
#endif
}


/******************************************************************
 * Function to get the thread affinities                           *
 ******************************************************************/
#ifdef USE_WINDOWS
DWORD GetThreadAffinityMask( HANDLE thread )
{
    DWORD mask = 1;
    DWORD old  = 0;
    // try every CPU one by one until one works or none are left
    while ( mask ) {
        old = SetThreadAffinityMask( thread, mask );
        if ( old ) {                              // this one worked
            SetThreadAffinityMask( thread, old ); // restore original
            return old;
        } else {
            if ( GetLastError() != ERROR_INVALID_PARAMETER )
                return 0; // fatal error, might as well throw an exception
        }
        mask <<= 1;
    }
    return 0;
}
#endif
std::vector<int> ThreadPool::getThreadAffinity()
{
    std::vector<int> procs;
#ifdef USE_LINUX
#ifdef _GNU_SOURCE
    cpu_set_t mask;
    int error = pthread_getaffinity_np( pthread_self(), sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        throw std::logic_error( "Error getting thread affinity" );
    for ( size_t i = 0; i < sizeof( cpu_set_t ) * CHAR_BIT; i++ ) {
        if ( CPU_ISSET( i, &mask ) )
            procs.push_back( i );
    }
#else
#warning pthread_getaffinity_np is not supported
    OS_warning( "pthread does not support pthread_getaffinity_np" );
#endif
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    OS_warning( "MAC does not support getting the thread affinity" );
#elif defined( USE_WINDOWS )
    size_t procMask = GetThreadAffinityMask( GetCurrentThread() );
    for ( size_t i = 0; i < sizeof( size_t ) * CHAR_BIT; i++ ) {
        if ( ( procMask & 0x1 ) != 0 )
            procs.push_back( i );
        procMask >>= 1;
    }
#else
#error Unknown OS
#endif
    return procs;
}
std::vector<int> ThreadPool::getThreadAffinity( int thread ) const
{
    if ( thread >= getNumThreads() )
        std::logic_error( "Invalid thread number" );
    std::vector<int> procs;
#ifdef USE_LINUX
#ifdef _GNU_SOURCE
    auto handle = const_cast<std::thread &>( d_thread[thread] ).native_handle();
    cpu_set_t mask;
    int error = pthread_getaffinity_np( handle, sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        throw std::logic_error( "Error getting thread affinity" );
    for ( size_t i = 0; i < sizeof( cpu_set_t ) * CHAR_BIT; i++ ) {
        if ( CPU_ISSET( i, &mask ) )
            procs.push_back( i );
    }
#else
#warning pthread_getaffinity_np is not supported
    OS_warning( "pthread does not support pthread_getaffinity_np" );
#endif
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    OS_warning( "MAC does not support getting the thread affinity" );
#elif defined( USE_WINDOWS )
    auto handle     = const_cast<std::thread &>( d_thread[thread] ).native_handle();
    size_t procMask = GetThreadAffinityMask( handle );
    for ( size_t i = 0; i < sizeof( size_t ) * CHAR_BIT; i++ ) {
        if ( ( procMask & 0x1 ) != 0 )
            procs.push_back( i );
        procMask >>= 1;
    }
#else
#error Unknown OS
#endif
    return procs;
}


/******************************************************************
 * Function to set the thread affinity                             *
 ******************************************************************/
void ThreadPool::setThreadAffinity( const std::vector<int> &procs )
{
#ifdef USE_LINUX
#ifdef _GNU_SOURCE
    cpu_set_t mask;
    CPU_ZERO( &mask );
    for ( size_t i = 0; i < procs.size(); i++ )
        CPU_SET( procs[i], &mask );
    int error = pthread_setaffinity_np( pthread_self(), sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        throw std::logic_error( "Error setting thread affinity" );
#else
#warning pthread_getaffinity_np is not supported
    OS_warning( "pthread does not support pthread_setaffinity_np" );
#endif
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    NULL_USE( procs );
    OS_warning( "MAC does not support setting the thread affinity" );
#elif defined( USE_WINDOWS )
    DWORD mask = 0;
    for ( size_t i = 0; i < procs.size(); i++ )
        mask |= ( (DWORD) 1 ) << procs[i];
    SetThreadAffinityMask( GetCurrentThread(), mask );
#else
#error Unknown OS
#endif
}
void ThreadPool::setThreadAffinity( int thread, const std::vector<int> &procs ) const
{
    if ( thread >= getNumThreads() )
        std::logic_error( "Invalid thread number" );
#ifdef USE_LINUX
#ifdef __USE_GNU
    cpu_set_t mask;
    CPU_ZERO( &mask );
    for ( size_t i = 0; i < procs.size(); i++ )
        CPU_SET( procs[i], &mask );
    auto handle = const_cast<std::thread &>( d_thread[thread] ).native_handle();
    int error   = pthread_setaffinity_np( handle, sizeof( cpu_set_t ), &mask );
    if ( error != 0 )
        throw std::logic_error( "Error setting thread affinity" );
#else
#warning pthread_getaffinity_np is not supported
    OS_warning( "pthread does not support pthread_setaffinity_np" );
#endif
#elif defined( USE_MAC )
    // MAC does not support getting or setting the affinity
    NULL_USE( procs );
    OS_warning( "MAC does not support getting the process affinity" );
#elif defined( USE_WINDOWS )
    DWORD mask = 0;
    for ( size_t i = 0; i < procs.size(); i++ )
        mask |= ( (DWORD) 1 ) << procs[i];
    auto handle = const_cast<std::thread &>( d_thread[thread] ).native_handle();
    SetThreadAffinityMask( handle, mask );
#else
#error Unknown OS
#endif
}


/******************************************************************
 * Function to perform some basic checks before we start           *
 ******************************************************************/
void ThreadPool::check_startup()
{
    // Check getting/setting a bit
    AtomicOperations::int64_atomic x[2] = { 0x0, 0x7 };
    set_bit( x, 2 );
    unset_bit( x, 66 );
    if ( x[0] != 4 || x[1] != 3 || !get_bit( x, 2 ) || get_bit( x, 66 ) )
        throw std::logic_error( "Getting/setting a bit failed" );
    // Check the thread id
    bool pass = true;
    ThreadPool::thread_id_t id;
    if ( id.getPriority() != -128 )
        pass = false;
    id.reset( 3, 564, nullptr );
    if ( id.getPriority() != 3 || id.getLocalID() != 564 )
        pass = false;
    if ( count_bits( 0x0 ) != 0 || count_bits( 0x03 ) != 2 )
        pass = false;
    if ( count_bits( ~( (size_t) 0 ) ) != 8 * sizeof( size_t ) )
        pass = false;
    if ( sizeof( size_t ) == 8 ) {
        if ( is_odd8( 0x0 ) || !is_odd8( 0x02 ) || is_odd8( 0x03 ) )
            pass = false;
        if ( is_odd8( ~( (size_t) 0 ) ) || !is_odd8( thread_id_t::maxThreadID ) )
            pass = false;
        for ( size_t i = 0; i < 1024; i++ ) {
            if ( ( count_bits( thread_id_t::maxThreadID - i ) % 2 == 1 ) !=
                 is_odd8( thread_id_t::maxThreadID - i ) ) {
                printp( "%i %i %s\n", count_bits( thread_id_t::maxThreadID - i ),
                    is_odd8( thread_id_t::maxThreadID - i ) ? 1 : 0,
                    std::bitset<64>( thread_id_t::maxThreadID - i ).to_string().c_str() );
                pass = false;
            }
        }
    }
    d_id_assign = thread_id_t::maxThreadID;
    AtomicOperations::atomic_decrement( &d_id_assign ); // Advance the id
    AtomicOperations::atomic_decrement( &d_id_assign ); // Advance the id
    ThreadPool::thread_id_t id2;
    id2.reset( 3, d_id_assign, nullptr );
    if ( isValid( id ) || !isValid( id2 ) )
        pass = false;
    if ( !pass )
        throw std::logic_error( "thread id test failed" );
}


/******************************************************************
 * Constructors/destructor                                         *
 ******************************************************************/
ThreadPool::ThreadPool(
    const int N, const std::string &affinity, const std::vector<int> &procs, int queueSize )
    : d_queue_list( queueSize )
{
    // Run some basic tests on startup
    check_startup();
    // Initialize the header/tail
    d_NULL_HEAD = generateHeadTail();
    d_NULL_TAIL = d_NULL_HEAD;
    // Initialize the variables to NULL values
    d_id_assign     = 0;
    d_signal_empty  = false;
    d_signal_count  = 0;
    d_N_threads     = 0;
    d_num_active    = 0;
    d_N_added       = 0;
    d_N_started     = 0;
    d_N_finished    = 0;
    d_max_wait_time = 600;
    memset( (void *) d_active, 0, MAX_THREADS / 8 );
    memset( (void *) d_cancel, 0, MAX_THREADS / 8 );
    d_wait_last = nullptr;
    for ( auto &i : d_wait )
        i = nullptr;
    // Initialize the id
    d_id_assign = thread_id_t::maxThreadID;
    // Create the threads
    setNumThreads( N, affinity, procs );
    // Verify that the threadpool is valid
    if ( !is_valid( this ) )
        throw std::logic_error( "Thread pool is not valid" );
}
ThreadPool::~ThreadPool()
{
    DISABLE_WARNINGS
    if ( !is_valid( this ) ) {
        std::cerr << "Thread pool is not valid, error calling destructor\n";
        return;
    }
    ENABLE_WARNINGS
    // Destroy the threads
    setNumThreads( 0 );
    // Delete all remaining data
    d_N_threads = ~0;
    d_NULL_HEAD = 0;
    d_NULL_TAIL = 0;
    delete d_wait_last;
#if MONITOR_THREADPOOL_PERFORMANCE == 1
    // Print the performance metrics
    printp( "ThreadPool Performance:\n" );
    printp( "add_work:  %lu us,  %lu us,  %lu us,  %lu us,  %lu us\n",
        total_add_work_time[0] / 1000, total_add_work_time[1] / 1000, total_add_work_time[2] / 1000,
        total_add_work_time[3] / 1000, total_add_work_time[4] / 1000 );
#endif
}


/******************************************************************
 * Check if the pointer points to a valid thread pool object       *
 ******************************************************************/
bool ThreadPool::is_valid( const ThreadPool *tpool )
{
    if ( tpool == nullptr )
        return false;
    if ( tpool->d_N_threads > MAX_THREADS )
        return false;
    if ( !validHeadTail( tpool->d_NULL_HEAD ) || tpool->d_NULL_HEAD != tpool->d_NULL_TAIL )
        return false;
    return true;
}


/******************************************************************
 * This function creates the threads in the thread pool            *
 ******************************************************************/
void ThreadPool::setNumThreads(
    int num_worker_threads, const std::string &affinity, const std::vector<int> &procs )
{
    // Check if we are a member thread
    if ( isMemberThread() )
        throw std::logic_error(
            "Member threads are not allowed to change the number of threads in the pool" );
    // Determing the number of threads we need to create or destroy
    if ( num_worker_threads > MAX_THREADS ) {
        printp( "Warning: Maximum Number of Threads is %i\n", MAX_THREADS );
        printp( "         Only that number will be created\n" );
        num_worker_threads = MAX_THREADS;
    } else if ( num_worker_threads < 0 ) {
        printp( "Error: cannot have a negitive number of threads\n" );
        printp( "       Setting the number of threads to 0\n" );
        num_worker_threads = 0;
    }
    int d_N_threads_diff = num_worker_threads - d_N_threads;
    if ( d_N_threads_diff > 0 ) {
        // Check that no threads are in the process of being deleted
        for ( long i : d_cancel ) {
            if ( i != 0 )
                throw std::logic_error(
                    "Threads are being created and destroyed at the same time" );
        }
        // Create the threads
        int j = d_N_threads;
        for ( int i = 0; i < d_N_threads_diff; i++ ) {
            d_N_threads++;
            set_bit( d_cancel, j );
            d_thread[j] = std::thread( create_new_thread, this, j );
            j++;
        }
        // Wait for all of the threads to finish initialization
        while ( true ) {
            std::this_thread::sleep_for( std::chrono::milliseconds( 25 ) );
            bool wait = false;
            for ( long i : d_cancel ) {
                if ( i != 0 )
                    wait = true;
            }
            if ( !wait )
                break;
        }
        std::this_thread::sleep_for( std::chrono::milliseconds( 25 ) );
    } else if ( d_N_threads_diff < 0 ) {
        // Reduce the number of threads
        if ( num_worker_threads == 0 ) {
            // Special case if we want to delete all of the threads
            wait_pool_finished();
        }
        // Tell the threads to shutdown
        for ( int i = 0; i > d_N_threads_diff; i-- )
            set_bit( d_cancel, d_N_threads - 1 + i );
        // Wake all threads to process the shutdown
        d_wait_work.notify_all();
        std::this_thread::sleep_for( std::chrono::milliseconds( 25 ) );
        // Wait for the threads to close
        for ( int i = 0; i > d_N_threads_diff; i-- ) {
            d_thread[d_N_threads - 1 + i].join();
            d_thread[d_N_threads - 1 + i] = std::thread();
            unset_bit( d_cancel, d_N_threads - 1 + i );
            d_threadId[d_N_threads - 1 + i] = std::thread::id();
        }
        d_N_threads += d_N_threads_diff;
    }
    if ( d_N_threads == 0 )
        return;
    // Get the default thread affinity to use
    std::vector<int> cpus;
    int tmp            = global_OS_behavior;
    global_OS_behavior = 1;
    OS_warning( "Dummy message (should not print)" );
    try {
        cpus = ThreadPool::getProcessAffinity();
    } catch ( ... ) {
        pout << "Warning: Unable to get default cpus for thread affinities\n";
    }
    if ( !cpus.empty() && !procs.empty() ) {
        cpus.resize( procs.size() );
        for ( size_t i = 0; i < procs.size(); i++ )
            cpus[i] = procs[i];
    }
    // Set the affinity model and the associated thread affinities
    // Note: not all OS's support setting the thread affinities
    std::vector<std::vector<int>> t_procs( d_N_threads );
    if ( cpus.empty() ) {
        // We do not have a list of cpus to use, do nothing (OS not supported)
    } else if ( affinity == "none" ) {
        // We are using the default thread affinities (all threads get all procs of the program)
        for ( int i = 0; i < d_N_threads; i++ )
            t_procs[i] = cpus;
    } else if ( affinity == "independent" ) {
        // We want to use an independent set of processors for each thread
        if ( cpus.size() == d_N_threads ) {
            // The number of cpus matches the number of threads
            for ( int i = 0; i < d_N_threads; i++ )
                t_procs[i] = std::vector<int>( 1, cpus[i] );
        } else if ( cpus.size() > d_N_threads ) {
            // There are more cpus than threads, threads will use more the one processor
            int N_procs_thread = ( cpus.size() + d_N_threads - 1 ) / d_N_threads;
            size_t k           = 0;
            for ( int i = 0; i < d_N_threads; i++ ) {
                for ( int j = 0; j < N_procs_thread && k < cpus.size(); j++ ) {
                    t_procs[i].push_back( cpus[k] );
                    k++;
                }
            }
        } else {
            // There are fewer cpus than threads, threads will share a processor
            auto N_threads_proc = ( cpus.size() + d_N_threads - 1 ) / cpus.size();
            for ( int i = 0; i < d_N_threads; i++ )
                t_procs[i].push_back( cpus[i / N_threads_proc] );
        }
    } else {
        global_OS_behavior = tmp;
        throw std::logic_error( "Unknown affinity model" );
    }
    try {
        for ( int i = 0; i < d_N_threads; i++ ) {
            ThreadPool::setThreadAffinity( i, t_procs[i] );
            auto cpus2 = getThreadAffinity( i );
            if ( cpus2 != t_procs[i] )
                pout << "Warning: error setting affinities (failed to set)\n";
        }
    } catch ( ... ) {
        pout << "Warning: error setting affinities (exception)\n";
    }
    global_OS_behavior = tmp;
}


/******************************************************************
 * This is the function that controls the individual thread and    *
 * allows it to do work.                                           *
 * Note: this function is lock free                                *
 ******************************************************************/
void ThreadPool::tpool_thread( int thread_id )
{
    bool shutdown         = false;
    bool printInfo        = false;
    d_threadId[thread_id] = std::this_thread::get_id();
    if ( get_bit( d_active, thread_id ) )
        throw std::logic_error( "Thread cannot already be active" );
    AtomicOperations::atomic_increment( &d_num_active );
    set_bit( d_active, thread_id );
    unset_bit( d_cancel, thread_id );
    setenv( "OMP_NUM_THREADS", "1" );
    setenv( "MKL_NUM_THREADS", "1" );
    if ( printInfo ) {
        // Print the pid
        printp( "pid = %i\n", (int) getpid() );
        // Print the processor affinities for the process
        try {
            auto cpus = ThreadPool::getProcessAffinity();
            printp( "%i cpus for current thread: ", (int) cpus.size() );
            for ( int cpu : cpus )
                printp( "%i ", cpu );
            printp( "\n" );
        } catch ( ... ) {
            printp( "Unable to get process affinity\n" );
        }
    }
    // Check for shutdown
    PROFILE_THREADPOOL_START( "thread active" );
    shutdown = false;
    while ( !shutdown ) {
        // Check if there is work to do
        if ( !d_queue_list.empty() ) {
            // Get next work item to process
            auto work_id =
                d_queue_list.remove( []( const thread_id_t &id ) { return id.ready(); } );
            if ( work_id.isNull() ) {
                std::this_thread::yield();
                continue;
            }
            WorkItem *work = work_id.work();
            AtomicOperations::atomic_increment( &d_N_started );
            // Start work here
            PROFILE_THREADPOOL_START( "thread working" );
            work->d_state = 2;
            if ( d_errorHandler ) {
                try {
                    work->run();
                } catch ( std::exception &e ) {
                    auto msg = Utilities::stringf(
                        "Error, caught exception in thread %i:\n  %s\n", thread_id, e.what() );
                    d_errorHandler( msg );
                } catch ( ... ) {
                    auto msg = Utilities::stringf(
                        "Error, caught unknown exception in thread %i\n", thread_id );
                    d_errorHandler( msg );
                }
            } else {
                work->run();
            }
            work->d_state = 3;
            PROFILE_THREADPOOL_STOP( "thread working" );
            AtomicOperations::atomic_increment( &d_N_finished );
            // Check if any threads are waiting on the current work item
            // This can be done without blocking
            for ( auto &i : d_wait ) {
                auto wait = AtomicOperations::atomic_get( &i );
                if ( wait != nullptr )
                    wait->id_finished( work_id );
            }
            // Check the signal count and signal if desired
            // This can be done without blocking
            if ( d_signal_count > 0 ) {
                int count = AtomicOperations::atomic_decrement( &d_signal_count );
                if ( count == 0 )
                    d_wait_finished.notify_all();
            }
        } else {
            int N_active = AtomicOperations::atomic_decrement( &d_num_active );
            unset_bit( d_active, thread_id );
            // Yield to give the main thread a chance to update
            std::this_thread::yield();
            // Alert main thread that a thread finished processing
            if ( ( N_active == 0 ) && d_signal_empty ) {
                d_wait_finished.notify_all();
                d_signal_empty = false;
            }
            // Wait for work
            PROFILE_THREADPOOL_STOP2( "thread active" );
            double wait_time = thread_id <= 2 ? 0.01 : 0.1;
            if ( d_queue_list.empty() )
                d_wait_work.wait_for( wait_time );
            PROFILE_THREADPOOL_START2( "thread active" );
            AtomicOperations::atomic_increment( &d_num_active );
            set_bit( d_active, thread_id );
        }
        // Check if there is a shutdown requested
        shutdown = get_bit( d_cancel, thread_id );
    }
    PROFILE_THREADPOOL_STOP( "thread active" );
    AtomicOperations::atomic_decrement( &d_num_active );
    unset_bit( d_active, thread_id );
    return;
}


/******************************************************************
 * This is the function that adds work to the thread pool          *
 * Note: this version uses a last in - first out work scheduling.  *
 ******************************************************************/
inline void ThreadPool::add_work( const ThreadPool::thread_id_t &id )
{
    auto work     = id.work();
    work->d_state = 1;
    // Check and change priorities of dependency ids
    int priority = id.getPriority();
    auto compare = []( const thread_id_t &a, const thread_id_t &b ) { return a == b; };
    for ( int i = 0; i < work->d_N_ids; i++ ) {
        const auto &id1 = work->d_ids[i];
        if ( !id1.started() && id1 < id ) {
            // Remove and add the id back with a higher priority
            auto id2 = d_queue_list.remove( compare, id1 );
            id2.setPriority( std::max( priority, id2.getPriority() ) );
            d_queue_list.insert( id2 );
        }
    }
    d_queue_list.insert( id );
    AtomicOperations::atomic_increment( &d_N_added );
}
void ThreadPool::add_work(
    size_t N, ThreadPool::WorkItem *work[], const int *priority, ThreadPool::thread_id_t *ids )
{
    // If we have a very long list, break it up into smaller pieces to keep the threads busy
    constexpr size_t block_size = 256;
    if ( N > block_size ) {
        size_t i = 0;
        while ( i < N ) {
            add_work( std::min( N - i, block_size ), &work[i], &priority[i], &ids[i] );
            i += block_size;
        }
        return;
    }
    PROFILE_THREADPOOL_START( "add_work" );
#if MONITOR_THREADPOOL_PERFORMANCE == 1
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    // Create the thread ids (can be done without blocking)
    for ( size_t i = 0; i < N; i++ )
        ids[i].reset( priority[i], AtomicOperations::atomic_decrement( &d_id_assign ), work[i] );
#if MONITOR_THREADPOOL_PERFORMANCE == 1
    auto t2 = std::chrono::high_resolution_clock::now();
    accumulate( total_add_work_time[0], t1, t2 );
#endif
    // If there are no threads, perform the work immediately
    if ( d_N_threads < 1 ) {
        for ( size_t i = 0; i < N; i++ ) {
            work[i]->d_state = 2;
            work[i]->run();
            work[i]->d_state = 3;
        }
#if MONITOR_THREADPOOL_PERFORMANCE == 1
        auto t5 = std::chrono::high_resolution_clock::now();
        accumulate( total_add_work_time[4], t2, t5 );
#endif
        PROFILE_THREADPOOL_STOP2( "add_work" );
        return;
    }
    // Wait for enough room in the queue (doesn't need blocking since it isn't that precise)
    if ( N > d_queue_list.capacity() - d_queue_list.size() ) {
        int N_wait = N - ( d_queue_list.capacity() - d_queue_list.size() );
        while ( N_wait > 0 ) {
            d_signal_count = std::min<int>( N_wait, 255 );
            d_wait_finished.wait_for( 1e-4 );
            N_wait = N - ( d_queue_list.capacity() - d_queue_list.size() );
        }
    }
#if MONITOR_THREADPOOL_PERFORMANCE == 1
    auto t3 = std::chrono::high_resolution_clock::now();
    accumulate( total_add_work_time[1], t2, t3 );
#endif
    // Get add the work items to the queue
    for ( size_t i = 0; i < N; i++ )
        add_work( ids[i] );
#if MONITOR_THREADPOOL_PERFORMANCE == 1
    auto t4 = std::chrono::high_resolution_clock::now();
    accumulate( total_add_work_time[2], t3, t4 );
#endif
    // Activate sleeping threads
    if ( d_num_active == d_N_threads ) {
        // All threads are active, no need to wake anybody
    } else if ( d_queue_list.empty() ) {
        // Queue is empty, no need to activate
    } else if ( N == 1 ) {
        // Added 1 item to the queue, wake 1 worker
        d_wait_work.notify_one();
    } else {
        // Added multple items in the queue, wake all workers
        d_wait_work.notify_all();
    }
#if MONITOR_THREADPOOL_PERFORMANCE == 1
    auto t5 = std::chrono::high_resolution_clock::now();
    accumulate( total_add_work_time[3], t4, t5 );
#endif
    PROFILE_THREADPOOL_STOP( "add_work" );
}


/******************************************************************
 * This function waits for a some of the work items to finish      *
 ******************************************************************/
static inline void check_finished(
    size_t N_work, const ThreadPool::thread_id_t *ids, size_t &N_finished, bool *finished )
{
    for ( size_t k = 0; k < N_work; k++ ) {
        if ( !finished[k] && ids[k].finished() ) {
            N_finished++;
            finished[k] = true;
        }
    }
}
int ThreadPool::wait_some( size_t N_work, const ThreadPool::thread_id_t *ids, size_t N_wait,
    bool *finished, int max_wait ) const
{
    // Check the inputs
    if ( N_wait > N_work )
        throw std::logic_error( "Invalid arguments in thread pool wait" );
    size_t N_finished = 0;
    memset( finished, 0, N_work * sizeof( bool ) );
    // Check that all the ids are valid
    size_t next_id = d_id_assign - 1;
    for ( size_t k = 0; k < N_work; k++ ) {
        if ( !ids[k].initialized() ) {
            finished[k] = true;
            N_finished++;
        }
        size_t local_id = ids[k].getLocalID();
        bool test = local_id == 0 || local_id > thread_id_t::maxThreadID || local_id <= next_id;
        test      = test && !finished[k];
        if ( test )
            throw std::logic_error( "Invalid ids for wait" );
    }
    // Check which ids have finished
    check_finished( N_work, ids, N_finished, finished );
    // If enough ids have finished return
    if ( N_finished >= N_wait )
        return N_finished;
    // Create the wait event struct
    auto tmp = new wait_ids_struct( N_work, ids, N_wait, d_cond_pool, MAX_WAIT, d_wait );
    // Wait for the ids
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = t1;
    int dt1 = 0;
    while ( dt1 < max_wait ) {
        if ( tmp->wait_for( std::min( max_wait, d_max_wait_time ), 0.01 ) )
            break;
        auto t3 = std::chrono::high_resolution_clock::now();
        dt1     = std::chrono::duration_cast<std::chrono::seconds>( t3 - t1 ).count();
        int dt2 = std::chrono::duration_cast<std::chrono::seconds>( t3 - t2 ).count();
        if ( dt2 >= d_max_wait_time ) {
            print_wait_warning();
            t2 = t3;
        }
    }
    // Update the ids that have finished
    check_finished( N_work, ids, N_finished, finished );
    // Delete the wait event struct
    // Note: we want to maintain the reference in case a thread is still using it
    // Note: technically this should be atomic, but it really isn't necessary here
    std::swap( d_wait_last, tmp );
    delete tmp;
    return N_finished;
}


/******************************************************************
 * This function waits for all of the threads to finish their work *
 ******************************************************************/
void ThreadPool::print_wait_warning() const
{
    pout << "Warning: Maximum wait time in ThreadPool exceeded, threads may be hung\n";
    pout << "N_active: " << d_num_active << std::endl;
    pout << "N_queued: " << d_queue_list.size() << std::endl;
    pout << "N_added: " << d_N_added << std::endl;
    pout << "N_started: " << d_N_started << std::endl;
    pout << "N_finished: " << d_N_finished << std::endl;
    pout << "queue.insert(): " << d_queue_list.N_insert() << std::endl;
    pout << "queue.remove(): " << d_queue_list.N_remove() << std::endl;
    pout << "Stack Trace:\n";
    auto call_stack = StackTrace::getAllCallStacks();
    StackTrace::cleanupStackTrace( call_stack );
    auto text = call_stack.print( "  " );
    for ( auto &line : text )
        pout << line << std::endl;
}
void ThreadPool::wait_pool_finished() const
{
    // First check that we are not one of the threads
    if ( isMemberThread() )
        throw std::logic_error( "Member thread attempted to call wait_pool_finished" );
    // Wait for all threads to finish their work
    auto t1 = std::chrono::high_resolution_clock::now();
    while ( d_num_active > 0 || !d_queue_list.empty() ) {
        // Wait for signal from last thread
        d_signal_empty = true;
        d_wait_finished.wait_for( 5e-4 );
        if ( d_num_active == 0 && d_queue_list.empty() )
            break;
        // Check that we have not exceeded the maximum time
        auto t2     = std::chrono::high_resolution_clock::now();
        int seconds = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
        if ( seconds > d_max_wait_time ) {
            print_wait_warning();
            t1 = t2;
        }
    }
    d_signal_empty = false;
}


/******************************************************************
 * Member functions of wait_ids_struct                             *
 ******************************************************************/
ThreadPool::wait_ids_struct::wait_ids_struct( size_t N, const ThreadPool::thread_id_t *ids,
    size_t N_wait, AtomicOperations::pool<condition_variable, 128> &cv_pool, int N_wait_list,
    volatile wait_ids_struct **list )
    : d_wait( N_wait ), d_N( 0 ), d_cv_pool( cv_pool ), d_wait_event( cv_pool.get() )
{
    d_ids = new ThreadPool::thread_id_t[N];
    for ( size_t i = 0; i < N; i++ ) {
        if ( ids[i].finished() )
            d_wait = std::max( d_wait - 1, 0 );
        else
            d_ids[d_N++] = ids[i];
    }
    quicksort( d_N, d_ids );
    d_finished = new bool[d_N];
    memset( (void *) d_finished, 0, d_N );
    int i = 0;
    while (
        !AtomicOperations::atomic_compare_and_swap( (void *volatile *) &list[i], nullptr, this ) ) {
        i = ( i + 1 ) % N_wait_list;
    }
    d_ptr = &list[i];
}
ThreadPool::wait_ids_struct::~wait_ids_struct()
{
    d_cv_pool.put( d_wait_event );
    delete[] d_finished;
    delete[] d_ids;
}
void ThreadPool::wait_ids_struct::id_finished( const ThreadPool::thread_id_t &id ) const
{
    int index = find_id( d_N, d_ids, id );
    if ( index >= 0 ) {
        d_finished[index] = true;
        int N_finished    = 0;
        for ( int i = 0; i < d_N; i++ )
            N_finished += d_finished[i] ? 1 : 0;
        if ( N_finished >= d_wait ) {
            d_N    = 0;
            d_wait = 0;
            AtomicOperations::atomic_compare_and_swap(
                (void *volatile *) d_ptr, (void *) *d_ptr, nullptr );
            d_wait_event->notify_all();
        }
    }
}
inline bool ThreadPool::wait_ids_struct::check()
{
    int N_finished = 0;
    for ( int i = 0; i < d_N; i++ )
        N_finished += d_finished[i] ? 1 : 0;
    if ( N_finished >= d_wait || d_N == 0 ) {
        *d_ptr = nullptr;
        d_wait = 0;
        d_N    = 0;
        return true;
    }
    return false;
}
bool ThreadPool::wait_ids_struct::wait_for( double total_time, double recheck_time )
{
    int total   = 1e6 * total_time;
    int recheck = 1e6 * recheck_time;
    auto t1     = std::chrono::high_resolution_clock::now();
    auto t2     = t1;
    int us1     = 0;
    while ( us1 < total ) {
        for ( int i = 0; i < d_N; i++ ) {
            if ( d_ids[i].finished() )
                d_finished[i] = true;
        }
        if ( check() )
            return true;
        int us2 = 0;
        while ( us2 < recheck ) {
            double dt = 1e-6 * std::max( 10, recheck - us2 );
            d_wait_event->wait_for( dt );
            if ( check() )
                return true;
            auto t3 = std::chrono::high_resolution_clock::now();
            us2     = std::chrono::duration_cast<std::chrono::microseconds>( t3 - t2 ).count();
            t2      = t3;
        }
        us1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    }
    return false;
}


/******************************************************************
 * templated quicksort routine                                     *
 ******************************************************************/
template<class T>
void quicksort( int n, T *arr )
{
    if ( n <= 1 )
        return;
    bool test;
    int i, ir, j, jstack, k, l, istack[100];
    T a, tmp_a;
    jstack = 0;
    l      = 0;
    ir     = n - 1;
    while ( true ) {
        if ( ir - l < 7 ) { // Insertion sort when subarray small enough.
            for ( j = l + 1; j <= ir; j++ ) {
                a    = arr[j];
                test = true;
                for ( i = j - 1; i >= 0; i-- ) {
                    if ( arr[i] < a ) {
                        arr[i + 1] = a;
                        test       = false;
                        break;
                    }
                    arr[i + 1] = arr[i];
                }
                if ( test ) {
                    i          = l - 1;
                    arr[i + 1] = a;
                }
            }
            if ( jstack == 0 )
                return;
            ir = istack[jstack]; // Pop stack and begin a new round of partitioning.
            l  = istack[jstack - 1];
            jstack -= 2;
        } else {
            k = ( l + ir ) / 2; // Choose median of left, center and right elements as partitioning
                                // element a. Also rearrange so that a(l) < a(l+1) < a(ir).
            tmp_a      = arr[k];
            arr[k]     = arr[l + 1];
            arr[l + 1] = tmp_a;
            if ( arr[l] > arr[ir] ) {
                tmp_a   = arr[l];
                arr[l]  = arr[ir];
                arr[ir] = tmp_a;
            }
            if ( arr[l + 1] > arr[ir] ) {
                tmp_a      = arr[l + 1];
                arr[l + 1] = arr[ir];
                arr[ir]    = tmp_a;
            }
            if ( arr[l] > arr[l + 1] ) {
                tmp_a      = arr[l];
                arr[l]     = arr[l + 1];
                arr[l + 1] = tmp_a;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l + 1]; // Partitioning element.
            for ( i = l + 2; i <= ir; i++ ) {
                if ( arr[i] < a )
                    continue;
                while ( arr[j] > a ) // Scan down to find element < a.
                    j--;
                if ( j < i )
                    break;       // Pointers crossed. Exit with partitioning complete.
                tmp_a  = arr[i]; // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
            }
            arr[l + 1] = arr[j]; // Insert partitioning element in both arrays.
            arr[j]     = a;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if ( ir - i + 1 >= j - l ) {
                istack[jstack]     = ir;
                istack[jstack - 1] = i;
                ir                 = j - 1;
            } else {
                istack[jstack]     = j - 1;
                istack[jstack - 1] = l;
                l                  = i;
            }
        }
    }
}


/************************************************************************
 * Function to find the id in a sorted vector                            *
 ************************************************************************/
inline int find_id( int n, const ThreadPool::thread_id_t *x, const ThreadPool::thread_id_t &id )
{
    if ( n == 0 )
        return -1;
    // Check if value is within the range of x
    if ( id == x[0] )
        return 0;
    if ( id < x[0] )
        return -1;
    if ( id == x[n - 1] )
        return n - 1;
    if ( id > x[n - 1] )
        return -1;
    // Perform the search
    size_t lower = 0;
    size_t upper = n - 1;
    while ( ( upper - lower ) != 1 ) {
        size_t index = ( upper + lower ) / 2;
        if ( x[index] == id )
            return index;
        if ( x[index] >= id )
            upper = index;
        else
            lower = index;
    }
    return -1;
}


/************************************************************************
 * Function to add dependencies to the work item                         *
 * Note: when expanding the size of d_ids, we need to allocate space for *
 * one extra entry for a spinlock.                                       *
 ************************************************************************/
void ThreadPool::WorkItem::add_dependencies( size_t N, const ThreadPool::thread_id_t *ids )
{
    if ( d_state != 0 ) {
        // The item has already been added to the threadpool,
        // we are not allowed to add dependencies
        throw std::logic_error(
            "Cannot add dependency to work item once it has been added the the threadpool" );
    }
    if ( d_N_ids + N > 0xFFFF )
        throw std::logic_error( "Cannot add more than 65000 dependencies" );
    if ( d_N_ids + N + 1 > d_size ) {
        thread_id_t *tmp = d_ids;
        unsigned int N2  = d_size;
        if ( N2 == 0 ) {
            N2 = 8;
        }
        while ( N2 < d_N_ids + N + 1 )
            N2 *= 2;
        d_ids = new thread_id_t[N2];
        for ( size_t i = 0; i < d_N_ids; i++ )
            const_cast<thread_id_t &>( ids[i] ).swap( tmp[i] );
        delete[] tmp;
        d_size     = N2;
        auto *lock = reinterpret_cast<int *>( &d_ids[d_size - 1] );
        *lock      = 0;
    }
    const ThreadPool::thread_id_t id0;
    for ( size_t i = 0; i < N; i++ ) {
        if ( ids[i] != id0 ) {
            if ( !ids[i].finished() ) {
                d_ids[d_N_ids] = ids[i];
                d_N_ids++;
            }
        }
    }
}

