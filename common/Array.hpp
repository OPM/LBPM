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
#ifndef included_ArrayClass_hpp
#define included_ArrayClass_hpp

#include "common/Array.h"
#include "common/FunctionTable.h"
#include "common/FunctionTable.hpp"
#include "common/Utilities.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <memory>


/********************************************************
 *  External instantiations                              *
 ********************************************************/
extern template class Array<char, FunctionTable>;
extern template class Array<uint8_t, FunctionTable>;
extern template class Array<uint16_t, FunctionTable>;
extern template class Array<uint32_t, FunctionTable>;
extern template class Array<uint64_t, FunctionTable>;
extern template class Array<int8_t, FunctionTable>;
extern template class Array<int16_t, FunctionTable>;
extern template class Array<int32_t, FunctionTable>;
extern template class Array<int64_t, FunctionTable>;
extern template class Array<double, FunctionTable>;
extern template class Array<float, FunctionTable>;


/********************************************************
 *  Helper functions                                     *
 ********************************************************/
template<class TYPE>
inline typename std::enable_if<std::is_integral<TYPE>::value, size_t>::type getRangeSize(
    const Range<TYPE> &range )
{
    return ( static_cast<int64_t>( range.j ) - static_cast<int64_t>( range.i ) ) /
           static_cast<int64_t>( range.k );
}
template<class TYPE>
inline typename std::enable_if<std::is_floating_point<TYPE>::value, size_t>::type getRangeSize(
    const Range<TYPE> &range )
{
    double tmp = static_cast<double>( ( range.j - range.i ) ) / static_cast<double>( range.k );
    return static_cast<size_t>( floor( tmp + 1e-12 ) + 1 );
}
template<class TYPE>
inline typename std::enable_if<std::is_same<TYPE, std::complex<float>>::value ||
                                   std::is_same<TYPE, std::complex<double>>::value,
    size_t>::type
getRangeSize( const Range<TYPE> &range )
{
    double tmp = std::real( ( range.j - range.i ) / ( range.k ) );
    return static_cast<size_t>( floor( tmp + 1e-12 ) + 1 );
}
template<class TYPE>
inline typename std::enable_if<std::is_integral<TYPE>::value, TYPE>::type getRangeValue(
    const Range<TYPE> &range, size_t index )
{
    return range.i + index * range.k;
}
template<class TYPE>
inline typename std::enable_if<std::is_floating_point<TYPE>::value, TYPE>::type getRangeValue(
    const Range<TYPE> &range, size_t index )
{
    return range.k * ( range.i / range.k + index );
}
template<class TYPE>
inline typename std::enable_if<std::is_same<TYPE, std::complex<float>>::value ||
                                   std::is_same<TYPE, std::complex<double>>::value,
    TYPE>::type
getRangeValue( const Range<TYPE> &range, size_t index )
{
    return range.k * ( range.i / range.k + static_cast<TYPE>( index ) );
}


/********************************************************
 *  Constructors                                         *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array()
    : d_isCopyable( true ), d_isFixedSize( false ), d_data( nullptr )
{
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( const ArraySize &N )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( N );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( size_t N ) : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( ArraySize( N ) );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( size_t N_rows, size_t N_cols )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( ArraySize( N_rows, N_cols ) );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( size_t N1, size_t N2, size_t N3 )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( ArraySize( N1, N2, N3 ) );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( size_t N1, size_t N2, size_t N3, size_t N4 )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( ArraySize( N1, N2, N3, N4 ) );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( size_t N1, size_t N2, size_t N3, size_t N4, size_t N5 )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( ArraySize( N1, N2, N3, N4, N5 ) );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( const std::vector<size_t> &N, const TYPE *data )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( N );
    if ( data ) {
        for ( size_t i = 0; i < d_size.length(); i++ )
            d_data[i] = data[i];
    }
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( const Range<TYPE> &range )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    size_t N = getRangeSize( range );
    allocate( { N } );
    for ( size_t i = 0; i < N; i++ )
        d_data[i] = getRangeValue( range, i );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( std::string str ) : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( 0 );
    if ( (int) std::count( str.begin(), str.end(), ' ' ) == (int) str.length() ) {
        // Empty string
        return;
    }
    // Remove unnecessary whitespace
    while ( str.front() == ' ' )
        str.erase( 0, 1 );
    while ( str.back() == ' ' )
        str.resize( str.length() - 1 );
    while ( str.find( ',' ) != std::string::npos )
        str[str.find( ',' )] = ' ';
    while ( str.find( "  " ) != std::string::npos )
        str.replace( str.find( "  " ), 2, " " );
    // Check if the string is of the format [...]
    if ( str.front() == '[' && str.back() == ']' ) {
        str.erase( 0, 1 );
        str.resize( str.length() - 1 );
        *this = Array( str );
        return;
    }
    // Check if we are dealing with a 2D array
    if ( str.find( ';' ) != std::string::npos ) {
        size_t i1 = 0;
        size_t i2 = str.find( ';' );
        std::vector<Array> x;
        while ( i2 > i1 ) {
            auto tmp = str.substr( i1, i2 - i1 );
            x.emplace_back( Array( tmp ) );
            i1 = i2 + 1;
            i2 = str.find( ';', i1 + 1 );
            if ( i2 == std::string::npos )
                i2 = str.length();
        }
        for ( auto &y : x )
            y.reshape( { 1, y.length() } );
        *this = cat( x, 0 );
        return;
    }
    // Begin parsing the array constructor
    size_t i1 = 0;
    size_t i2 = str.find( ' ' );
    std::vector<TYPE> data;
    while ( i2 > i1 ) {
        auto tmp = str.substr( i1, i2 - i1 );
        int type = std::count( tmp.begin(), tmp.end(), ':' );
        if ( type == 0 ) {
            data.push_back( std::stod( tmp ) );
        } else if ( type == 1 ) {
            size_t k   = tmp.find( ':' );
            TYPE x1    = std::stod( tmp.substr( 0, k ) );
            TYPE x2    = std::stod( tmp.substr( k + 1 ) );
            double tol = 1e-8 * ( x2 - x1 );
            for ( TYPE x = x1; x <= x2 + tol; x += 1 )
                data.push_back( x );
        } else if ( type == 2 ) {
            size_t k1  = tmp.find( ':' );
            size_t k2  = tmp.find( ':', k1 + 1 );
            TYPE x1    = std::stod( tmp.substr( 0, k1 ) );
            TYPE dx    = std::stod( tmp.substr( k1 + 1, k2 - k1 - 1 ) );
            TYPE x2    = std::stod( tmp.substr( k2 + 1 ) );
            double tol = 1e-8 * ( x2 - x1 );
            for ( TYPE x = x1; x <= x2 + tol; x += dx )
                data.push_back( x );
        } else {
            throw std::logic_error( "Failed to parse string constructor: " + str );
        }
        i1 = i2;
        i2 = str.find( ' ', i1 + 1 );
        if ( i2 == std::string::npos )
            i2 = str.length();
    }
    allocate( data.size() );
    for ( size_t i = 0; i < data.size(); i++ )
        d_data[i] = data[i];
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( std::initializer_list<TYPE> x )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    allocate( { x.size() } );
    auto it = x.begin();
    for ( size_t i = 0; i < x.size(); ++i, ++it )
        d_data[i] = *it;
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::allocate( const ArraySize &N )
{
    if ( d_isFixedSize )
        throw std::logic_error( "Array cannot be resized" );
    d_size      = N;
    auto length = d_size.length();
    if ( length == 0 )
        d_ptr.reset();
    else
        d_ptr.reset( new ( std::nothrow ) TYPE[length], []( TYPE *p ) { delete[] p; } );
    d_data = d_ptr.get();
    if ( length > 0 && d_data == nullptr )
        throw std::logic_error( "Failed to allocate array" );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( const Array &rhs )
    : d_isCopyable( true ), d_isFixedSize( false )
{
    if ( !rhs.d_isCopyable )
        throw std::logic_error( "Array cannot be copied" );
    allocate( rhs.size() );
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] = rhs.d_data[i];
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::Array( Array &&rhs )
    : d_isCopyable( rhs.d_isCopyable ),
      d_isFixedSize( rhs.d_isFixedSize ),
      d_size( rhs.d_size ),
      d_data( rhs.d_data )
{
    rhs.d_data = nullptr;
    d_ptr      = std::move( rhs.d_ptr );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> &Array<TYPE, FUN, Allocator>::operator=( const Array &rhs )
{
    if ( this == &rhs )
        return *this;
    if ( !rhs.d_isCopyable )
        throw std::logic_error( "Array cannot be copied" );
    this->allocate( rhs.size() );
    for ( size_t i = 0; i < d_size.length(); i++ )
        this->d_data[i] = rhs.d_data[i];
    return *this;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> &Array<TYPE, FUN, Allocator>::operator=( Array &&rhs )
{
    if ( this == &rhs )
        return *this;
    d_isCopyable  = rhs.d_isCopyable;
    d_isFixedSize = rhs.d_isFixedSize;
    d_size        = rhs.d_size;
    d_data        = rhs.d_data;
    d_ptr         = std::move( rhs.d_ptr );
    rhs.d_data    = nullptr;
    return *this;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> &Array<TYPE, FUN, Allocator>::operator=( const std::vector<TYPE> &rhs )
{
    this->allocate( ArraySize( rhs.size() ) );
    for ( size_t i = 0; i < rhs.size(); i++ )
        this->d_data[i] = rhs[i];
    return *this;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator>::~Array()
{
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::clear()
{
    d_isCopyable  = true;
    d_isFixedSize = false;
    d_size        = ArraySize();
    d_ptr.reset();
    d_data = nullptr;
}


/********************************************************
 *  Copy/move values from one array to another (resize)  *
 ********************************************************/
template<class TYPE>
static inline void moveValues( const ArraySize &N1, const ArraySize &N2, TYPE *data1, TYPE *data2 )
{
    for ( size_t i5 = 0; i5 < std::min( N1[4], N2[4] ); i5++ ) {
        for ( size_t i4 = 0; i4 < std::min( N1[3], N2[3] ); i4++ ) {
            for ( size_t i3 = 0; i3 < std::min( N1[2], N2[2] ); i3++ ) {
                for ( size_t i2 = 0; i2 < std::min( N1[1], N2[1] ); i2++ ) {
                    for ( size_t i1 = 0; i1 < std::min( N1[0], N2[0] ); i1++ ) {
                        size_t index1 = N1.index( i1, i2, i3, i4, i5 );
                        size_t index2 = N2.index( i1, i2, i3, i4, i5 );
                        data2[index2] = std::move( data1[index1] );
                    }
                }
            }
        }
    }
}
template<bool test, class TYPE>
static inline typename std::enable_if<test, void>::type copyValues(
    const ArraySize &N1, const ArraySize &N2, const TYPE *data1, TYPE *data2 )
{
    for ( size_t i5 = 0; i5 < std::min( N1[4], N2[4] ); i5++ ) {
        for ( size_t i4 = 0; i4 < std::min( N1[3], N2[3] ); i4++ ) {
            for ( size_t i3 = 0; i3 < std::min( N1[2], N2[2] ); i3++ ) {
                for ( size_t i2 = 0; i2 < std::min( N1[1], N2[1] ); i2++ ) {
                    for ( size_t i1 = 0; i1 < std::min( N1[0], N2[0] ); i1++ ) {
                        size_t index1 = N1.index( i1, i2, i3, i4, i5 );
                        size_t index2 = N2.index( i1, i2, i3, i4, i5 );
                        data2[index2] = data1[index1];
                    }
                }
            }
        }
    }
}
template<bool test, class TYPE>
static inline typename std::enable_if<!test, void>::type copyValues(
    const ArraySize &, const ArraySize &, const TYPE *, TYPE * )
{
    throw std::logic_error( "No copy constructor" );
}


/********************************************************
 *  Resize the array                                     *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::resize( const ArraySize &N )
{
    // Check if the array actually changed size
    bool equal = true;
    for ( size_t i = 0; i < ArraySize::maxDim(); i++ )
        equal = equal && N[i] == d_size[i];
    if ( equal ) {
        d_size = N;
        return;
    }
    // Store the old data
    auto N0    = d_size;
    auto data0 = d_ptr;
    // Allocate new data
    allocate( N );
    // Copy the old values
    if ( N.length() > 0 && d_size.length() > 0 ) {
        if ( data0.use_count() <= 1 ) {
            // We own the data, use std:move
            moveValues( N0, N, data0.get(), d_data );
        } else {
            // We do not own the data, copy
            copyValues<std::is_copy_constructible<TYPE>::value, TYPE>( N0, N, data0.get(), d_data );
        }
    }
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::resizeDim( int dim, size_t N, const TYPE &value )
{
    if ( dim < 0 || dim > d_size.ndim() )
        throw std::out_of_range( "Invalid dimension" );
    size_t N0 = d_size[dim];
    auto size = d_size;
    size.resize( dim, N );
    resize( size );
    size_t n1 = 1, n2 = 1;
    for ( int d = 0; d < dim; d++ )
        n1 *= size[d];
    for ( size_t d = dim + 1; d < size.ndim(); d++ )
        n2 *= size[d];
    for ( size_t k = 0; k < n2; k++ ) {
        for ( size_t j = N0; j < N; j++ ) {
            for ( size_t i = 0; i < n1; i++ ) {
                d_data[i + j * n1 + k * n1 * N] = value;
            }
        }
    }
}


/********************************************************
 *  Reshape the array                                     *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::reshape( const ArraySize &N )
{
    if ( N.length() != d_size.length() )
        throw std::logic_error( "reshape is not allowed to change the array size" );
    d_size = N;
}


/********************************************************
 *  Subset the array                                     *
 ********************************************************/
// Helper function to check subset indices
template<class TYPE, class FUN, class Allocator>
inline void Array<TYPE, FUN, Allocator>::checkSubsetIndex(
    const std::vector<Range<size_t>> &range ) const
{
    bool test = (int) range.size() == d_size.ndim();
    for ( size_t d = 0; d < range.size(); d++ )
        test = test && range[d].j <= d_size[d];
    if ( !test )
        throw std::logic_error( "indices for subset are invalid" );
}
template<class TYPE, class FUN, class Allocator>
std::vector<Range<size_t>> Array<TYPE, FUN, Allocator>::convert(
    const std::vector<size_t> &index ) const
{
    std::vector<Range<size_t>> range( d_size.ndim() );
    if ( index.size() % 2 != 0 || static_cast<int>( index.size() / 2 ) < d_size.ndim() )
        throw std::logic_error( "indices for subset are invalid" );
    for ( int d = 0; d < d_size.ndim(); d++ )
        range[d] = Range<size_t>( index[2 * d + 0], index[2 * d + 1] );
    return range;
}
// Helper function to return dimensions for the subset array
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::getSubsetArrays( const std::vector<Range<size_t>> &index,
    std::array<size_t, 5> &first, std::array<size_t, 5> &last, std::array<size_t, 5> &inc,
    std::array<size_t, 5> &N )
{
    first.fill( 0 );
    last.fill( 0 );
    inc.fill( 1 );
    N.fill( 1 );
    size_t ndim = index.size();
    for ( size_t d = 0; d < ndim; d++ ) {
        first[d] = index[d].i;
        last[d]  = index[d].j;
        inc[d]   = index[d].k;
        N[d]     = ( last[d] - first[d] + inc[d] ) / inc[d];
    }
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::subset(
    const std::vector<Range<size_t>> &index ) const
{
    // Get the subset indicies
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( index, first, last, inc, N1 );
    ArraySize S1( d_size.ndim(), N1.data() );
    // Create the new array
    Array<TYPE> subset_array( S1 );
    // Fill the new array
    static_assert( ArraySize::maxDim() == 5, "Not programmed for more than 5 dimensions" );
    TYPE *subset_data = subset_array.data();
    for ( size_t i4 = first[4], k1 = 0; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0], k1++ ) {
                        size_t k2       = d_size.index( i0, i1, i2, i3, i4 );
                        subset_data[k1] = d_data[k2];
                    }
                }
            }
        }
    }
    return subset_array;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::subset(
    const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return subset( range );
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::copySubset(
    const std::vector<Range<size_t>> &index, const Array<TYPE, FUN, Allocator> &subset )
{
    // Get the subset indices
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( index, first, last, inc, N1 );
    // Copy the sub-array
    static_assert( ArraySize::maxDim() == 5, "Not programmed for more than 5 dimensions" );
    const TYPE *src_data = subset.data();
    for ( size_t i4 = first[4], k1 = 0; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0], k1++ ) {
                        size_t k2  = d_size.index( i0, i1, i2, i3, i4 );
                        d_data[k2] = src_data[k1];
                    }
                }
            }
        }
    }
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::addSubset(
    const std::vector<Range<size_t>> &index, const Array<TYPE, FUN, Allocator> &subset )
{
    // Get the subset indices
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( index, first, last, inc, N1 );
    // add the sub-array
    static_assert( ArraySize::maxDim() == 5, "Not programmed for more than 5 dimensions" );
    for ( size_t i4 = first[4], k1 = 0; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0], k1++ ) {
                        size_t k2 = d_size.index( i0, i1, i2, i3, i4 );
                        d_data[k2] += subset.d_data[k1];
                    }
                }
            }
        }
    }
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::copySubset(
    const std::vector<size_t> &index, const Array<TYPE, FUN, Allocator> &subset )
{
    auto range = convert( index );
    copySubset( range, subset );
}

template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::addSubset(
    const std::vector<size_t> &index, const Array<TYPE, FUN, Allocator> &subset )
{
    auto range = convert( index );
    addSubset( range, subset );
}


/********************************************************
 *  Operator overloading                                 *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
bool Array<TYPE, FUN, Allocator>::operator==( const Array &rhs ) const
{
    if ( this == &rhs )
        return true;
    if ( d_size != rhs.d_size )
        return false;
    bool match = true;
    for ( size_t i = 0; i < d_size.length(); i++ )
        match = match && d_data[i] == rhs.d_data[i];
    return match;
}


/********************************************************
 *  Get a view of an C array                             *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
std::unique_ptr<Array<TYPE, FUN, Allocator>> Array<TYPE, FUN, Allocator>::view(
    const ArraySize &N, std::shared_ptr<TYPE> &data )
{
    auto array    = std::make_unique<Array<TYPE, FUN, Allocator>>();
    array->d_size = N;
    array->d_ptr  = data;
    array->d_data = array->d_ptr.get();
    return array;
}
template<class TYPE, class FUN, class Allocator>
std::unique_ptr<const Array<TYPE, FUN, Allocator>> Array<TYPE, FUN, Allocator>::constView(
    const ArraySize &N, std::shared_ptr<const TYPE> const &data )
{
    auto array    = std::make_unique<Array<TYPE, FUN, Allocator>>();
    array->d_size = N;
    array->d_ptr  = std::const_pointer_cast<TYPE>( data );
    array->d_data = array->d_ptr.get();
    return array;
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::view2( Array<TYPE, FUN, Allocator> &src )
{
    view2( src.size(), src.getPtr() );
    d_data = src.d_data;
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::view2( const ArraySize &N, std::shared_ptr<TYPE> const &data )
{
    d_size = N;
    d_ptr  = data;
    d_data = d_ptr.get();
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::viewRaw(
    const ArraySize &N, TYPE *data, bool isCopyable, bool isFixedSize )
{
    d_isCopyable  = isCopyable;
    d_isFixedSize = isFixedSize;
    d_ptr.reset();
    d_size = N;
    d_data = data;
}


/********************************************************
 *  Basic functions                                      *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::swap( Array &other )
{
    // check that dimensions match
    if ( d_size != other.d_size )
        throw std::logic_error( "length of arrays do not match" );
    // swap the data
    std::swap( d_data, other.d_data );
    std::swap( d_ptr, other.d_ptr );
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::fill( const TYPE &value )
{
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] = value;
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::scale( const TYPE &value )
{
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] *= value;
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::pow(
    const Array<TYPE, FUN, Allocator> &baseArray, const TYPE &exp )
{
    // not insisting on the shapes being the same
    // but insisting on the total size being the same
    if ( d_size.length() != baseArray.length() )
        throw std::logic_error( "length of arrays do not match" );

    const auto base_data = baseArray.data();
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] = std::pow( base_data[i], exp );
}


/********************************************************
 *  Replicate the array                                  *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::repmat(
    const std::vector<size_t> &N_rep ) const
{
    std::vector<size_t> N2( d_size.begin(), d_size.end() );
    if ( N2.size() < N_rep.size() )
        N2.resize( N_rep.size(), 1 );
    std::array<size_t, 5> N1, Nr;
    N1.fill( 1 );
    Nr.fill( 1 );
    for ( size_t d = 0; d < N_rep.size(); d++ ) {
        N1[d] = d_size[d];
        Nr[d] = N_rep[d];
        N2[d] *= N_rep[d];
    }
    Array<TYPE, FUN, Allocator> y( N2 );
    static_assert( ArraySize::maxDim() <= 5, "Not programmed for dimensions > 5" );
    TYPE *y2 = y.data();
    for ( size_t i4 = 0, index = 0; i4 < N1[4]; i4++ ) {
        for ( size_t j4 = 0; j4 < Nr[4]; j4++ ) {
            for ( size_t i3 = 0; i3 < N1[3]; i3++ ) {
                for ( size_t j4 = 0; j4 < Nr[3]; j4++ ) {
                    for ( size_t i2 = 0; i2 < N1[2]; i2++ ) {
                        for ( size_t j4 = 0; j4 < Nr[2]; j4++ ) {
                            for ( size_t i1 = 0; i1 < N1[1]; i1++ ) {
                                for ( size_t j4 = 0; j4 < Nr[1]; j4++ ) {
                                    for ( size_t i0 = 0; i0 < N1[0]; i0++ ) {
                                        size_t k = d_size.index( i0, i1, i2, i3, i4 );
                                        TYPE x   = d_data[k];
                                        for ( size_t j4 = 0; j4 < Nr[0]; j4++, index++ )
                                            y2[index] = x;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return y;
}


/********************************************************
 *  Simple math operations                               *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
bool Array<TYPE, FUN, Allocator>::NaNs() const
{
    bool test = false;
    for ( size_t i = 0; i < d_size.length(); i++ )
        test = test || d_data[i] != d_data[i];
    return test;
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::mean( void ) const
{
    TYPE x = this->sum() / d_size.length();
    return x;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::min( int dir ) const
{
    auto size_ans = d_size;
    size_ans.resize( dir, 1 );
    Array<TYPE, FUN, Allocator> ans( size_ans );
    size_t N1 = 1, N2 = 1, N3 = 1;
    for ( int d = 0; d < std::min<int>( dir, d_size.ndim() ); d++ )
        N1 *= d_size[d];
    N2 = d_size[dir];
    for ( size_t d = dir + 1; d < d_size.ndim(); d++ )
        N3 *= d_size[d];
    TYPE *data2 = ans.d_data;
    for ( size_t i3 = 0; i3 < N3; i3++ ) {
        for ( size_t i1 = 0; i1 < N1; i1++ ) {
            TYPE x = d_data[i1 + i3 * N1 * N2];
            for ( size_t i2 = 0; i2 < N2; i2++ )
                x = std::min( x, d_data[i1 + i2 * N1 + i3 * N1 * N2] );
            data2[i1 + i3 * N1] = x;
        }
    }
    return ans;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::max( int dir ) const
{
    auto size_ans = d_size;
    size_ans.resize( dir, 1 );
    Array<TYPE, FUN, Allocator> ans( size_ans );
    size_t N1 = 1, N2 = 1, N3 = 1;
    for ( int d = 0; d < std::min<int>( dir, d_size.ndim() ); d++ )
        N1 *= d_size[d];
    N2 = d_size[dir];
    DISABLE_WARNINGS // Suppress false array subscript is above array bounds
        for ( size_t d = dir + 1; d < d_size.ndim(); d++ ) N3 *= d_size[d];
    ENABLE_WARNINGS // Enable warnings
        TYPE *data2 = ans.d_data;
    for ( size_t i3 = 0; i3 < N3; i3++ ) {
        for ( size_t i1 = 0; i1 < N1; i1++ ) {
            TYPE x = d_data[i1 + i3 * N1 * N2];
            for ( size_t i2 = 0; i2 < N2; i2++ )
                x = std::max( x, d_data[i1 + i2 * N1 + i3 * N1 * N2] );
            data2[i1 + i3 * N1] = x;
        }
    }
    return ans;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::sum( int dir ) const
{
    auto size_ans = d_size;
    size_ans.resize( dir, 1 );
    Array<TYPE, FUN, Allocator> ans( size_ans );
    size_t N1 = 1, N2 = 1, N3 = 1;
    for ( int d = 0; d < std::min<int>( dir, d_size.ndim() ); d++ )
        N1 *= d_size[d];
    N2 = d_size[dir];
    DISABLE_WARNINGS
    for ( size_t d = dir + 1; d < d_size.ndim(); d++ )
        N3 *= d_size[d];
    ENABLE_WARNINGS
    TYPE *data2 = ans.d_data;
    for ( size_t i3 = 0; i3 < N3; i3++ ) {
        for ( size_t i1 = 0; i1 < N1; i1++ ) {
            TYPE x = 0;
            for ( size_t i2 = 0; i2 < N2; i2++ )
                x += d_data[i1 + i2 * N1 + i3 * N1 * N2];
            data2[i1 + i3 * N1] = x;
        }
    }
    return ans;
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::min( const std::vector<Range<size_t>> &range ) const
{
    // Get the subset indicies
    checkSubsetIndex( range );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( range, first, last, inc, N1 );
    static_assert( ArraySize::maxDim() <= 5, "Function programmed for more than 5 dimensions" );
    TYPE x = std::numeric_limits<TYPE>::max();
    for ( size_t i4 = first[4]; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0] ) {
                        size_t k1 = d_size.index( i0, i1, i2, i3, i4 );
                        x         = std::min( x, d_data[k1] );
                    }
                }
            }
        }
    }
    return x;
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::max( const std::vector<Range<size_t>> &range ) const
{
    // Get the subset indicies
    checkSubsetIndex( range );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( range, first, last, inc, N1 );
    static_assert( ArraySize::maxDim() <= 5, "Function programmed for more than 5 dimensions" );
    TYPE x = std::numeric_limits<TYPE>::min();
    for ( size_t i4 = first[4]; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0] ) {
                        size_t k1 = d_size.index( i0, i1, i2, i3, i4 );
                        x         = std::max( x, d_data[k1] );
                    }
                }
            }
        }
    }
    return x;
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::sum( const std::vector<Range<size_t>> &range ) const
{
    // Get the subset indicies
    checkSubsetIndex( range );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( range, first, last, inc, N1 );
    static_assert( ArraySize::maxDim() <= 5, "Function programmed for more than 5 dimensions" );
    TYPE x = 0;
    for ( size_t i4 = first[4]; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0] ) {
                        size_t k1 = d_size.index( i0, i1, i2, i3, i4 );
                        x += d_data[k1];
                    }
                }
            }
        }
    }
    return x;
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::mean( const std::vector<Range<size_t>> &range ) const
{
    // Get the subset indicies
    checkSubsetIndex( range );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( range, first, last, inc, N1 );
    static_assert( ArraySize::maxDim() <= 5, "Function programmed for more than 5 dimensions" );
    size_t n = 1;
    for ( auto &d : N1 )
        n *= d;
    TYPE x = sum( range ) / n;
    return x;
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::min( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return min( range );
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::max( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return max( range );
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::sum( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return sum( range );
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::mean( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return mean( range );
}


/********************************************************
 *  Find all elements that match the given operation     *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
std::vector<size_t> Array<TYPE, FUN, Allocator>::find(
    const TYPE &value, std::function<bool( const TYPE &, const TYPE & )> compare ) const
{
    std::vector<size_t> result;
    result.reserve( d_size.length() );
    for ( size_t i = 0; i < d_size.length(); i++ ) {
        if ( compare( d_data[i], value ) )
            result.push_back( i );
    }
    return result;
}


/********************************************************
 *  Print an array to an output stream                   *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::print(
    std::ostream &os, const std::string &name, const std::string &prefix ) const
{
    if ( d_size.ndim() == 1 ) {
        for ( size_t i = 0; i < d_size[0]; i++ )
            os << prefix << name << "[" << i << "] = " << d_data[i] << std::endl;
    } else if ( d_size.ndim() == 2 ) {
        os << prefix << name << ":" << std::endl;
        for ( size_t i = 0; i < d_size[0]; i++ ) {
            for ( size_t j = 0; j < d_size[1]; j++ )
                os << prefix << "  " << operator()( i, j );
            os << std::endl;
        }
    } else {
        throw std::logic_error( "Not programmed for this dimension" );
    }
}


/********************************************************
 *  Reverse dimensions (transpose)                       *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::reverseDim() const
{
    size_t N2[ArraySize::maxDim()];
    for ( int d = 0; d < ArraySize::maxDim(); d++ )
        N2[d] = d_size[ArraySize::maxDim() - d - 1];
    ArraySize S2( ArraySize::maxDim(), N2 );
    Array<TYPE, FUN, Allocator> y( S2 );
    static_assert( ArraySize::maxDim() == 5, "Not programmed for dimensions other than 5" );
    TYPE *y2 = y.data();
    for ( size_t i0 = 0; i0 < d_size[0]; i0++ ) {
        for ( size_t i1 = 0; i1 < d_size[1]; i1++ ) {
            for ( size_t i2 = 0; i2 < d_size[2]; i2++ ) {
                for ( size_t i3 = 0; i3 < d_size[3]; i3++ ) {
                    for ( size_t i4 = 0; i4 < d_size[4]; i4++ ) {
                        y2[S2.index( i4, i3, i2, i1, i0 )] =
                            d_data[d_size.index( i0, i1, i2, i3, i4 )];
                    }
                }
            }
        }
    }
    for ( int d = 0; d < d_size.ndim(); d++ )
        N2[d] = d_size[d_size.ndim() - d - 1];
    y.reshape( ArraySize( d_size.ndim(), N2 ) );
    return y;
}


/********************************************************
 *  Coarsen the array                                    *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::coarsen(
    const Array<TYPE, FUN, Allocator> &filter ) const
{
    auto S2 = size();
    for ( size_t i = 0; i < S2.size(); i++ ) {
        size_t s = S2[i] / filter.size( i );
        S2.resize( i, s );
        if ( S2[i] * filter.size( i ) != size( i ) )
            throw std::invalid_argument( "Array must be multiple of filter size" );
    }
    Array<TYPE, FUN, Allocator> y( S2 );
    if ( d_size.ndim() > 3 )
        throw std::logic_error( "Function not programmed for more than 3 dimensions" );
    const auto &Nh = filter.d_size;
    for ( size_t k1 = 0; k1 < y.d_size[2]; k1++ ) {
        for ( size_t j1 = 0; j1 < y.d_size[1]; j1++ ) {
            for ( size_t i1 = 0; i1 < y.d_size[0]; i1++ ) {
                TYPE tmp = 0;
                for ( size_t k2 = 0; k2 < Nh[2]; k2++ ) {
                    for ( size_t j2 = 0; j2 < Nh[1]; j2++ ) {
                        for ( size_t i2 = 0; i2 < Nh[0]; i2++ ) {
                            tmp += filter( i2, j2, k2 ) * this->operator()( i1 *Nh[0] + i2,
                                                              j1 * Nh[1] + j2, k1 * Nh[2] + k2 );
                        }
                    }
                }
                y( i1, j1, k1 ) = tmp;
            }
        }
    }
    return y;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::coarsen( const std::vector<size_t> &ratio,
    std::function<TYPE( const Array<TYPE, FUN, Allocator> & )> filter ) const
{
    if ( ratio.size() != d_size.ndim() )
        throw std::logic_error( "ratio size does not match ndim" );
    auto S2 = size();
    for ( size_t i = 0; i < S2.size(); i++ ) {
        S2.resize( i, S2[i] / ratio[i] );
        if ( S2[i] * ratio[i] != size( i ) )
            throw std::invalid_argument( "Array must be multiple of filter size" );
    }
    Array<TYPE, FUN, Allocator> tmp( ratio );
    Array<TYPE, FUN, Allocator> y( S2 );
    if ( d_size.ndim() > 3 )
        throw std::logic_error( "Function not programmed for more than 3 dimensions" );
    for ( size_t k1 = 0; k1 < y.d_size[2]; k1++ ) {
        for ( size_t j1 = 0; j1 < y.d_size[1]; j1++ ) {
            for ( size_t i1 = 0; i1 < y.d_size[0]; i1++ ) {
                for ( size_t k2 = 0; k2 < ratio[2]; k2++ ) {
                    for ( size_t j2 = 0; j2 < ratio[1]; j2++ ) {
                        for ( size_t i2 = 0; i2 < ratio[0]; i2++ ) {
                            tmp( i2, j2, k2 ) = this->operator()(
                                i1 *ratio[0] + i2, j1 * ratio[1] + j2, k1 * ratio[2] + k2 );
                        }
                    }
                }
                y( i1, j1, k1 ) = filter( tmp );
            }
        }
    }
    return y;
}


/********************************************************
 *  Concatenates the arrays                              *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::cat( const Array<TYPE, FUN, Allocator> &x, int dim )
{
    std::vector<Array<TYPE, FUN, Allocator>> tmp( 2 );
    tmp[0].view2( *this );
    tmp[1].view2( const_cast<Array<TYPE, FUN, Allocator> &>( x ) );
    *this = cat( tmp, dim );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::cat( const std::vector<Array> &x, int dim )
{
    if ( x.empty() )
        return Array<TYPE, FUN, Allocator>();
    // Check that the dimensions match
    bool check = true;
    for ( size_t i = 1; i < x.size(); i++ ) {
        check = check && x[i].ndim() == x[0].ndim();
        for ( int d = 0; d < x[0].ndim(); d++ )
            if ( d != dim )
                check = check && x[i].size( d ) == x[0].size( d );
    }
    if ( !check )
        throw std::logic_error( "Array dimensions do not match for concatenation" );
    // Create the output array
    auto size = x[0].d_size;
    for ( size_t i = 1; i < x.size(); i++ )
        size.resize( dim, size[dim] + x[i].size( dim ) );
    Array<TYPE, FUN, Allocator> out( size );
    size_t N1 = 1;
    size_t N2 = size[dim];
    size_t N3 = 1;
    for ( int d = 0; d < dim; d++ )
        N1 *= size[d];
    for ( size_t d = dim + 1; d < size.ndim(); d++ )
        N3 *= size[d];
    TYPE *data = out.data();
    for ( size_t i = 0, i0 = 0; i < x.size(); i++ ) {
        const TYPE *src = x[i].data();
        size_t N22      = x[i].size( dim );
        for ( size_t j2 = 0; j2 < N3; j2++ ) {
            for ( size_t i1 = 0; i1 < N22; i1++ ) {
                for ( size_t j1 = 0; j1 < N1; j1++ ) {
                    data[j1 + ( i1 + i0 ) * N1 + j2 * N1 * N2] = src[j1 + i1 * N1 + j2 * N1 * N22];
                }
            }
        }
        i0 += N22;
    }
    return out;
}


/********************************************************
 *  Interpolate                                          *
 ********************************************************/
template<class T>
struct is_compatible_double : std::integral_constant<bool,
                                  std::is_floating_point<typename std::remove_cv<T>::type>::value ||
                                      std::is_integral<typename std::remove_cv<T>::type>::value> {
};
template<class TYPE>
inline typename std::enable_if<is_compatible_double<TYPE>::value, TYPE>::type Array_interp_1D(
    double x, int N, const TYPE *data )
{
    int i = floor( x );
    i     = std::max( i, 0 );
    i     = std::min( i, N - 2 );
    return ( i + 1 - x ) * data[i] + ( x - i ) * data[i + 1];
}
template<class TYPE>
inline typename std::enable_if<is_compatible_double<TYPE>::value, TYPE>::type Array_interp_2D(
    double x, double y, int Nx, int Ny, const TYPE *data )
{
    int i       = floor( x );
    i           = std::max( i, 0 );
    i           = std::min( i, Nx - 2 );
    double dx   = x - i;
    double dx2  = 1.0 - dx;
    int j       = floor( y );
    j           = std::max( j, 0 );
    j           = std::min( j, Ny - 2 );
    double dy   = y - j;
    double dy2  = 1.0 - dy;
    double f[4] = { (double) data[i + j * Nx], (double) data[i + 1 + j * Nx],
        (double) data[i + ( j + 1 ) * Nx], (double) data[i + 1 + ( j + 1 ) * Nx] };
    return ( dx * f[1] + dx2 * f[0] ) * dy2 + ( dx * f[3] + dx2 * f[2] ) * dy;
}
template<class TYPE>
inline typename std::enable_if<is_compatible_double<TYPE>::value, TYPE>::type Array_interp_3D(
    double x, double y, double z, int Nx, int Ny, int Nz, const TYPE *data )
{
    int i       = floor( x );
    i           = std::max( i, 0 );
    i           = std::min( i, Nx - 2 );
    double dx   = x - i;
    double dx2  = 1.0 - dx;
    int j       = floor( y );
    j           = std::max( j, 0 );
    j           = std::min( j, Ny - 2 );
    double dy   = y - j;
    double dy2  = 1.0 - dy;
    int k       = floor( z );
    k           = std::max( k, 0 );
    k           = std::min( k, Nz - 2 );
    double dz   = z - k;
    double dz2  = 1.0 - dz;
    double f[8] = { (double) data[i + j * Nx + k * Nx * Ny],
        (double) data[i + 1 + j * Nx + k * Nx * Ny],
        (double) data[i + ( j + 1 ) * Nx + k * Nx * Ny],
        (double) data[i + 1 + ( j + 1 ) * Nx + k * Nx * Ny],
        (double) data[i + j * Nx + ( k + 1 ) * Nx * Ny],
        (double) data[i + 1 + j * Nx + ( k + 1 ) * Nx * Ny],
        (double) data[i + ( j + 1 ) * Nx + ( k + 1 ) * Nx * Ny],
        (double) data[i + 1 + ( j + 1 ) * Nx + ( k + 1 ) * Nx * Ny] };
    double h0   = ( dx * f[1] + dx2 * f[0] ) * dy2 + ( dx * f[3] + dx2 * f[2] ) * dy;
    double h1   = ( dx * f[5] + dx2 * f[4] ) * dy2 + ( dx * f[7] + dx2 * f[6] ) * dy;
    return h0 * dz2 + h1 * dz;
}
template<class TYPE>
inline typename std::enable_if<!is_compatible_double<TYPE>::value, TYPE>::type Array_interp_1D(
    double, int, const TYPE * )
{
    throw std::logic_error( "Invalid conversion" );
}
template<class TYPE>
inline typename std::enable_if<!is_compatible_double<TYPE>::value, TYPE>::type Array_interp_2D(
    double, double, int, int, const TYPE * )
{
    throw std::logic_error( "Invalid conversion" );
}
template<class TYPE>
inline typename std::enable_if<!is_compatible_double<TYPE>::value, TYPE>::type Array_interp_3D(
    double, double, double, int, int, int, const TYPE * )
{
    throw std::logic_error( "Invalid conversion" );
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::interp( const std::vector<double> &x ) const
{
    int ndim = 0, dim[5];
    double x2[5];
    for ( int d = 0; d < d_size.ndim(); d++ ) {
        if ( d_size[d] > 1 ) {
            x2[ndim]  = x[d];
            dim[ndim] = d_size[d];
            ndim++;
        }
    }
    TYPE f = 0;
    if ( ndim == 0 ) {
        // No data, do nothing
    } else if ( ndim == 1 ) {
        f = Array_interp_1D( x2[0], dim[0], d_data );
    } else if ( ndim == 2 ) {
        f = Array_interp_2D( x2[0], x2[1], dim[0], dim[1], d_data );
    } else if ( ndim == 3 ) {
        f = Array_interp_3D( x2[0], x2[1], x2[2], dim[0], dim[1], dim[2], d_data );
    } else {
        throw std::logic_error( "Not finished" );
    }
    return f;
}


/********************************************************
 *  Math operations (should call the Math class)         *
 ********************************************************/
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::rand()
{
    FUN::rand( *this );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> &Array<TYPE, FUN, Allocator>::operator+=(
    const Array<TYPE, FUN, Allocator> &rhs )
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a + b; };
    FUN::transform( fun, *this, rhs, *this );
    return *this;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> &Array<TYPE, FUN, Allocator>::operator-=(
    const Array<TYPE, FUN, Allocator> &rhs )
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a - b; };
    FUN::transform( fun, *this, rhs, *this );
    return *this;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> &Array<TYPE, FUN, Allocator>::operator+=( const TYPE &rhs )
{
    const auto &fun = [rhs]( const TYPE &x ) { return x + rhs; };
    FUN::transform( fun, *this, *this );
    return *this;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> &Array<TYPE, FUN, Allocator>::operator-=( const TYPE &rhs )
{
    const auto &fun = [rhs]( const TYPE &x ) { return x - rhs; };
    FUN::transform( fun, *this, *this );
    return *this;
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::min() const
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a < b ? a : b; };
    return FUN::reduce( fun, *this );
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::max() const
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a > b ? a : b; };
    return FUN::reduce( fun, *this );
}
template<class TYPE, class FUN, class Allocator>
TYPE Array<TYPE, FUN, Allocator>::sum() const
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a + b; };
    return FUN::reduce( fun, *this );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::multiply(
    const Array<TYPE, FUN, Allocator> &a, const Array<TYPE, FUN, Allocator> &b )
{
    Array<TYPE, FUN, Allocator> c;
    FUN::multiply( a, b, c );
    return c;
}
template<class TYPE, class FUN, class Allocator>
void Array<TYPE, FUN, Allocator>::axpby(
    const TYPE &alpha, const Array<TYPE, FUN, Allocator> &x, const TYPE &beta )
{
    const auto &fun = [alpha, beta](
                          const TYPE &x, const TYPE &y ) { return alpha * x + beta * y; };
    return FUN::transform( fun, x, *this, *this );
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::transform(
    std::function<TYPE( const TYPE & )> fun, const Array<TYPE, FUN, Allocator> &x )
{
    Array<TYPE, FUN, Allocator> y;
    FUN::transform( fun, x, y );
    return y;
}
template<class TYPE, class FUN, class Allocator>
Array<TYPE, FUN, Allocator> Array<TYPE, FUN, Allocator>::transform(
    std::function<TYPE( const TYPE &, const TYPE & )> fun, const Array<TYPE, FUN, Allocator> &x,
    const Array<TYPE, FUN, Allocator> &y )
{
    Array<TYPE, FUN, Allocator> z;
    FUN::transform( fun, x, y, z );
    return z;
}
template<class TYPE, class FUN, class Allocator>
bool Array<TYPE, FUN, Allocator>::equals( const Array &rhs, TYPE tol ) const
{
    return FUN::equals( *this, rhs, tol );
}

#endif
