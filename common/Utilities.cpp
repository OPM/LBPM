#include "common/Utilities.h"
#include "StackTrace/StackTrace.h"
#include "StackTrace/ErrorHandlers.h"

#ifdef USE_TIMER
#include "MemoryApp.h"
#include "ProfilerApp.h"
#endif

#ifdef USE_MPI
#include "common/MPI.h"
#endif

#include <algorithm>
#include <math.h>
#include <mutex>


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
// clang-format on


// Mutex for Utility functions
static std::mutex Utilities_mutex;


/****************************************************************************
 *  Function to perform the default startup/shutdown sequences               *
 ****************************************************************************/
void Utilities::startup( int argc, char **argv )
{
    NULL_USE( argc );
    NULL_USE( argv );
    // Disable OpenMP
    Utilities::setenv( "OMP_NUM_THREADS", "1" );
    Utilities::setenv( "MKL_NUM_THREADS", "1" );
    // Start MPI
#ifdef USE_MPI
    int provided;
    MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
    if ( provided < MPI_THREAD_MULTIPLE ) {
        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        if ( rank == 0 )
            std::cerr << "Warning: Failed to start MPI with necessary thread support, thread support will be disabled" << std::endl;
    }
    StackTrace::globalCallStackInitialize( MPI_COMM_WORLD );
#endif
    // Set the error handlers
    Utilities::setAbortBehavior( true, 3 );
    Utilities::setErrorHandlers();
}
void Utilities::shutdown()
{
    // Clear the error handlers
    Utilities::clearErrorHandlers();
    StackTrace::clearSignals();
    StackTrace::clearSymbols();
    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    StackTrace::globalCallStackFinalize();
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
#endif
#ifdef USE_TIMER
    PROFILE_DISABLE();
    auto memory = MemoryApp::getMemoryStats();
    if ( rank == 0 && memory.N_new > memory.N_delete )
        MemoryApp::print( std::cout );
#endif
}


/****************************************************************************
 *  Function to set an environemental variable                               *
 ****************************************************************************/
void Utilities::setenv( const std::string &name, const std::string &value )
{
    Utilities_mutex.lock();
#if defined( USE_LINUX ) || defined( USE_MAC )
    bool pass = false;
    if ( !value.empty() )
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
std::string Utilities::getenv( const std::string &name )
{
    std::string var;
    Utilities_mutex.lock();
    auto tmp = std::getenv( name.data() );
    if ( tmp )
        var = std::string( tmp );
    Utilities_mutex.unlock();
    return var;
}


/****************************************************************************
 *  Factor a number into it's prime factors                                  *
 ****************************************************************************/
std::vector<int> Utilities::factor(size_t number)
{
    if ( number<=3 ) 
        return std::vector<int>(1,(int)number);
    size_t i, n, n_max;
    bool factor_found;
    // Compute the maximum number of factors
    int N_primes_max = 1;
    n = number;
    while (n >>= 1) ++N_primes_max;
    // Initialize n, factors 
    n = number;
    std::vector<int> factors;
    factors.reserve(N_primes_max);
    while ( 1 ) {
        // Check if n is a trivial prime number
        if ( n==2 || n==3 || n==5 ) {
            factors.push_back( (int) n );
            break;
        } 
        // Check if n is divisible by 2
        if ( n%2 == 0 ) {
            factors.push_back( 2 );
            n/=2;
            continue;
        } 
        // Check each odd number until a factor is reached
        n_max = (size_t) floor(sqrt((double) n));
        factor_found = false;
        for (i=3; i<=n_max; i+=2) {
            if ( n%i == 0 ) {
                factors.push_back( i );
                n/=i;
                factor_found = true;
                break;
            } 
        }
        if ( factor_found )
            continue;
        // No factors were found, the number must be prime
        factors.push_back( (int) n );
        break;
    }
    // Sort the factors
    std::sort( factors.begin(), factors.end() );
    return factors;
}


/****************************************************************************
 *  Dummy function to prevent compiler from optimizing away variable         *
 ****************************************************************************/
void Utilities::nullUse( void* data )
{
    NULL_USE(data);
}



