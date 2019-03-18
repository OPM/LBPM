#include "threadpool/atomic_helpers.h"
#include <stdexcept>

#ifdef USE_PTHREAD_ATOMIC_LOCK
// Print a warning if we defaulted to use pthreads for atomic operations
// This can decrease the performance of atomic operations
// We print the message here so it is only printed once
#warning using pthreads for atomic operations, this may affect performance
#endif


namespace AtomicOperations {

#ifdef USE_PTHREAD_ATOMIC_LOCK
pthread_mutex_t atomic_pthread_lock;
static pthread_mutexattr_t threadpool_global_attr;
static int create_atomic_pthread_lock()
{
    pthread_mutexattr_init( &threadpool_global_attr );
    int error = pthread_mutex_init( &atomic_pthread_lock, &threadpool_global_attr );
    if ( error != 0 )
        throw std::logic_error( "Error initializing mutex:" );
    return error;
}
int atomic_pthread_lock_initialized = create_atomic_pthread_lock();
#endif


// Atomic operations for floating types
double atomic_add( double volatile *x, double y )
{
    static_assert( sizeof( double ) == sizeof( int64_atomic ), "Unexpected size" );
    union U {
        double d;
        int64_atomic i;
    };
    U a, b;
    bool swap = false;
    auto x2   = reinterpret_cast<int64_atomic volatile *>( x );
    while ( !swap ) {
        a.i  = atomic_add( x2, 0 );
        b.d  = a.d + y;
        swap = atomic_compare_and_swap( x2, a.i, b.i );
    }
    return b.d;
}
float atomic_add( float volatile *x, float y )
{
    static_assert( sizeof( float ) == sizeof( int32_atomic ), "Unexpected size" );
    union U {
        float d;
        int32_atomic i;
    };
    U a, b;
    bool swap = false;
    auto x2   = reinterpret_cast<int32_atomic volatile *>( x );
    while ( !swap ) {
        a.i  = atomic_add( x2, 0 );
        b.d  = a.d + y;
        swap = atomic_compare_and_swap( x2, a.i, b.i );
    }
    return b.d;
}


} // AtomicOperations namespace

