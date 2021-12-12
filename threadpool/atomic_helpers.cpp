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

