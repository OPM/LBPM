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
#ifndef included_AtomicModelAtomicList
#define included_AtomicModelAtomicList

#include <atomic>
#include <csignal>
#include <functional>

#include "threadpool/atomic_helpers.h"


/** \class AtomicList
 *
 * \brief Maintain a sorted list of entries
 * \details This class implements a basic sorted list that is thread-safe and lock-free.
 *    Entries are stored smallest to largest according to the compare operator
 */
template<class TYPE, int MAX_SIZE, class COMPARE = std::less<TYPE>>
class AtomicList final
{
public:
    //! Default constructor
    AtomicList( const TYPE &default_value = TYPE(), const COMPARE &comp = COMPARE() );

    /*!
     * \brief   Remove an item from the list
     * \details Find and remove first entry that meets the given criteria
     * @return          Return the item that matches the criteria,
     *                  or the default item if no item matches
     * @param compare   Comparison function object (i.e. an object that satisfies
     *                  the requirements of Compare) which returns ​true if the
     *                  given value meets the selection criteria.
     *                  The signature of the comparison function should be equivalent to:
     *                      bool cmp( const TYPE& value, ... );
     * @param args      Additional arguments for the comparison
     */
    template<class Compare, class... Args>
    inline TYPE remove( Compare compare, Args... args );

    //! Remove the first from the list
    inline TYPE remove_first();

    /*!
     * \brief   Insert an item
     * \details Insert an item into the list
     * @param x         Item to insert
     */
    inline void insert( TYPE x );

    /*!
     * \brief   Return the size of the list
     * \details Return the number of items in the list
     */
    inline int size() const { return AtomicOperations::atomic_get( &d_N ); }

    /*!
     * \brief   Check if the list is empty
     * \details Return true if the list is empty
     */
    inline bool empty() const { return AtomicOperations::atomic_get( &d_N ) == 0; }

    /*!
     * \brief   Return the capacity of the list
     * \details Return the maximum number of items the list can hold
     */
    inline int capacity() const { return MAX_SIZE; }

    /*!
     * \brief   Check the list
     * \details Perform a series of checks to verify the list is in a stable state.
     *    Note: This function is only partially thread-safe: it will block all other
     *    operations on the list, but check may fail if we caught a thread modifing the list.
     *    It is intended for debugging purposes only!
     * @return          This function returns true if the list is in a good working state
     */
    inline bool check();


    //! Return the total number of inserts since object creation
    inline int64_t N_insert() const { return AtomicOperations::atomic_get( &d_N_insert ); }


    //! Return the total number of removals since object creation
    inline int64_t N_remove() const { return AtomicOperations::atomic_get( &d_N_remove ); }

private:
    // Data members
    COMPARE d_compare;
    volatile TYPE d_default;
    volatile TYPE d_objects[MAX_SIZE];
    volatile AtomicOperations::int32_atomic d_N;
    volatile AtomicOperations::int32_atomic d_next[MAX_SIZE + 1];
    volatile AtomicOperations::int32_atomic d_unused;
    volatile AtomicOperations::int64_atomic d_N_insert;
    volatile AtomicOperations::int64_atomic d_N_remove;

private:
    inline int lock( int i )
    {
        if ( i == -1 )
            return -1;
        int tmp = 0;
        while ( tmp == 0 )
            tmp = AtomicOperations::atomic_fetch_and_and( &d_next[i], 0 );
        return tmp;
    }
    inline void unlock( int i, int value )
    {
        if ( i != -1 )
            AtomicOperations::atomic_fetch_and_or( &d_next[i], value );
    }
    inline int get_unused()
    {
        int i = 0;
        while ( i == 0 )
            i = AtomicOperations::atomic_fetch_and_and( &d_unused, 0 );
        AtomicOperations::atomic_fetch_and_or( &d_unused, -( d_next[i] + 4 ) + 1 );
        d_next[i] = -3;
        return i;
    }
    inline void put_unused( int i )
    {
        int j = 0;
        while ( j == 0 )
            AtomicOperations::atomic_swap( &d_unused, &j );
        d_next[i] = -3 - j;
        AtomicOperations::atomic_fetch_and_or( &d_unused, i );
    }


private:
    AtomicList( const AtomicList & );
    AtomicList &operator=( const AtomicList & );
};


/** \class MemoryPool
 *
 * \brief Pool allocator
 * \details This class implements a basic fast pool allocator that is thread-safe.
 */
template<class TYPE, class INT_TYPE = int>
class MemoryPool final
{
public:
    //! Default constructor
    explicit MemoryPool( size_t size );

    //! destructor
    ~MemoryPool();

    /*!
     * \brief   Allocate an object
     * \details Allocates a new object from the pool
     * @return          Return the new pointer, or nullptr if there is no more room in the pool
     */
    inline TYPE *allocate();

    /*!
     * \brief   Insert an item
     * \details Insert an item into the list
     * @param ptr       The pointer to free
     */
    inline void free( TYPE *ptr );

private:
    // Data members
    volatile TYPE *d_objects;
    volatile AtomicOperations::int32_atomic d_next;

private:
    MemoryPool( const MemoryPool & );
    MemoryPool &operator=( const MemoryPool & );
};


#include "threadpool/atomic_list.hpp"

#endif

