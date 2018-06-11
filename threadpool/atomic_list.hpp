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
#ifndef included_AtomicList_hpp
#define included_AtomicList_hpp


#include <iostream>
#include <stdexcept>
#include <thread>


/******************************************************************
 * Constructor                                                     *
 ******************************************************************/
template<class TYPE, int MAX_SIZE, class COMPARE>
AtomicList<TYPE, MAX_SIZE, COMPARE>::AtomicList( const TYPE &default_value, const COMPARE &comp )
    : d_compare( comp ), d_default( default_value )
{
    d_N        = 0;
    d_next[0]  = -1;
    d_unused   = 1;
    d_N_insert = 0;
    d_N_remove = 0;
    for ( int i = 0; i < MAX_SIZE; i++ ) {
        d_next[i + 1] = -5 - i;
        d_objects[i]  = d_default;
    }
}


/******************************************************************
 * Remove an item                                                  *
 ******************************************************************/
template<class TYPE, int MAX_SIZE, class COMPARE>
template<class Compare, class... Args>
inline TYPE AtomicList<TYPE, MAX_SIZE, COMPARE>::remove( Compare compare, Args... args )
{
    // Acquiring temporary ownership
    int pos   = 0;
    auto next = lock( 0 );
    while ( true ) {
        if ( next == -1 ) {
            // We have no more entires to search
            unlock( pos, -1 );
            pos = -1;
            break;
        }
        if ( next < 0 )
            throw std::logic_error( "Internal error" );
        // Acquire ownership of the next item
        int next2 = lock( next );
        // Test to see if the object passes compare
        bool test = compare( const_cast<TYPE &>( d_objects[next - 1] ), args... );
        if ( test ) {
            // We want to return this object, update next to point to another entry and remove the
            // entry
            unlock( next, -3 );
            unlock( pos, next2 );
            pos = next;
            break;
        }
        // Release the ownership and move on
        unlock( pos, next );
        pos  = next;
        next = next2;
    }
    TYPE rtn( d_default );
    if ( pos != -1 ) {
        std::swap( rtn, const_cast<TYPE &>( d_objects[pos - 1] ) );
        put_unused( pos );
        AtomicOperations::atomic_decrement( &d_N );
        AtomicOperations::atomic_increment( &d_N_remove );
    }
    return rtn;
}
template<class TYPE, int MAX_SIZE, class COMPARE>
inline TYPE AtomicList<TYPE, MAX_SIZE, COMPARE>::remove_first()
{
    TYPE rtn( d_default );
    auto next = lock( 0 );
    if ( next != -1 ) {
        int next2 = lock( next );
        unlock( next, -3 );
        unlock( 0, next2 );
        std::swap( rtn, const_cast<TYPE &>( d_objects[next - 1] ) );
        put_unused( next );
        AtomicOperations::atomic_decrement( &d_N );
        AtomicOperations::atomic_increment( &d_N_remove );
    } else {
        unlock( 0, next );
    }
    return rtn;
}


/******************************************************************
 * Insert an item                                                  *
 ******************************************************************/
template<class TYPE, int MAX_SIZE, class COMPARE>
inline void AtomicList<TYPE, MAX_SIZE, COMPARE>::insert( TYPE x )
{
    int N_used = AtomicOperations::atomic_increment( &d_N );
    if ( N_used > MAX_SIZE ) {
        AtomicOperations::atomic_decrement( &d_N );
        throw std::logic_error( "No room in list" );
    }
    // Get an index to store the entry
    auto index = get_unused();
    if ( index < 1 )
        throw std::logic_error( "Internal error" );
    // Store the object in d_objects
    AtomicOperations::atomic_increment( &d_N_insert );
    d_objects[index - 1] = x;
    d_next[index]        = -1;
    // Find the position to store and update the next entires
    int pos   = 0;
    auto next = lock( pos );
    while ( true ) {
        // Get the next item in the list (acquiring temporary ownership)
        if ( next == -1 ) {
            // We have no more entires to search, store here
            unlock( pos, index );
            break;
        }
        // Test to see if the object is < the value being compared
        bool test = d_compare.operator()( x, const_cast<TYPE &>( d_objects[next - 1] ) );
        if ( test ) {
            // We want to store this object before next
            d_next[index] = next;
            unlock( pos, index );
            break;
        }
        // Release the ownership and move on
        int last = pos;
        pos      = next;
        next     = lock( next );
        unlock( last, pos );
    }
}


/******************************************************************
 * Check the internal structures of the list                       *
 * This is mostly thread-safe, but blocks all threads              *
 ******************************************************************/
template<class TYPE, int MAX_SIZE, class COMPARE>
inline bool AtomicList<TYPE, MAX_SIZE, COMPARE>::check()
{
    // Get the lock and check for any other threads modifying the list
    auto start = lock( 0 );
    std::this_thread::sleep_for( std::chrono::microseconds( 100 ) );
    // Perform the checks on the list
    bool pass    = true;
    int N1       = 0;
    int N2       = 0;
    int N_unused = 0;
    int N_tail   = 0;
    for ( int i = 0; i < MAX_SIZE; i++ ) {
        if ( d_objects[i] != d_default )
            N1++;
    }
    for ( int i = 0; i < MAX_SIZE + 1; i++ ) {
        int next = i == 0 ? start : d_next[i];
        if ( next > 0 ) {
            N2++;
        } else if ( next < -3 ) {
            N_unused++;
        } else if ( next == -1 ) {
            N_tail++;
        } else {
            pass = false;
        }
    }
    pass    = pass && N_tail == 1 && N1 == d_N && N2 == d_N && N_unused + d_N == MAX_SIZE;
    int it  = 0;
    int pos = 0;
    while ( true ) {
        int next = pos == 0 ? start : d_next[pos];
        if ( next == -1 )
            break;
        pos = next;
        it++;
    }
    pass = pass && it == d_N;
    // Unlock the list and return the results
    unlock( 0, start );
    return pass;
}


/******************************************************************
 * MemoryPool                                                      *
 ******************************************************************/
template<class TYPE, class INT_TYPE>
MemoryPool<TYPE, INT_TYPE>::MemoryPool( size_t size )
{
    static_assert( sizeof( TYPE ) >= sizeof( int ),
        "sizeof(TYPE) must be >= sizeof(int) to ensure proper operation" );
    static_assert( sizeof( TYPE ) >= sizeof( INT_TYPE ),
        "sizeof(TYPE) must be >= sizeof(INT_TYPE) to ensure proper operation" );
    d_objects = reinterpret_cast<TYPE *>( malloc( sizeof( TYPE ) * size ) );
    d_next    = 1;
    for ( size_t i = 0; i < size; i++ )
        reinterpret_cast<volatile INT_TYPE &>( d_objects[i] ) = i + 1;
    reinterpret_cast<volatile INT_TYPE &>( d_objects[size - 1] ) = -1;
}
template<class TYPE, class INT_TYPE>
MemoryPool<TYPE, INT_TYPE>::~MemoryPool()
{
    free( const_cast<TYPE *>( d_objects ) );
    d_objects = nullptr;
}
template<class TYPE, class INT_TYPE>
inline TYPE *MemoryPool<TYPE, INT_TYPE>::allocate()
{
    AtomicOperations::int32_atomic i = 0;
    while ( i == 0 )
        AtomicOperations::atomic_swap( &d_next, &i );
    TYPE *ptr = nullptr;
    if ( i != -1 ) {
        INT_TYPE j = reinterpret_cast<volatile INT_TYPE &>( d_objects[i - 1] );
        ptr        = const_cast<TYPE *>( &d_objects[i - 1] );
        new ( ptr ) TYPE();
        i = j + 1;
    }
    AtomicOperations::atomic_fetch_and_or( &d_next, i );
    return ptr;
}
template<class TYPE, class INT_TYPE>
inline void MemoryPool<TYPE, INT_TYPE>::free( TYPE *ptr )
{
    ptr->~TYPE();
    AtomicOperations::int32_atomic i = 0;
    while ( i == 0 )
        AtomicOperations::atomic_swap( &d_next, &i );
    reinterpret_cast<INT_TYPE &>( *ptr ) = i - 1;
    i                                    = ptr - d_objects + 1;
    AtomicOperations::atomic_fetch_and_or( &d_next, i );
}


#endif
