#ifndef included_ArrayClass_hpp
#define included_ArrayClass_hpp

#include "common/Array.h"
#include "common/FunctionTable.h"
#include "common/Utilities.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>


/********************************************************
 *  ArraySize                                            *
 ********************************************************/
inline ArraySize::ArraySize()
{
    d_ndim   = 1;
    d_N[0]   = 0;
    d_N[1]   = 1;
    d_N[2]   = 1;
    d_N[3]   = 1;
    d_N[4]   = 1;
    d_length = 0;
}
inline ArraySize::ArraySize( size_t N1 )
{
    d_ndim   = 1;
    d_N[0]   = N1;
    d_N[1]   = 1;
    d_N[2]   = 1;
    d_N[3]   = 1;
    d_N[4]   = 1;
    d_length = N1;
}
inline ArraySize::ArraySize( size_t N1, size_t N2 )
{
    d_ndim   = 2;
    d_N[0]   = N1;
    d_N[1]   = N2;
    d_N[2]   = 1;
    d_N[3]   = 1;
    d_N[4]   = 1;
    d_length = N1 * N2;
}
inline ArraySize::ArraySize( size_t N1, size_t N2, size_t N3 )
{
    d_ndim   = 3;
    d_N[0]   = N1;
    d_N[1]   = N2;
    d_N[2]   = N3;
    d_N[3]   = 1;
    d_N[4]   = 1;
    d_length = N1 * N2 * N3;
}
inline ArraySize::ArraySize( size_t N1, size_t N2, size_t N3, size_t N4 )
{
    d_ndim   = 4;
    d_N[0]   = N1;
    d_N[1]   = N2;
    d_N[2]   = N3;
    d_N[3]   = N4;
    d_N[4]   = 1;
    d_length = N1 * N2 * N3 * N4;
}
inline ArraySize::ArraySize( size_t N1, size_t N2, size_t N3, size_t N4, size_t N5 )
{
    d_ndim   = 5;
    d_N[0]   = N1;
    d_N[1]   = N2;
    d_N[2]   = N3;
    d_N[3]   = N4;
    d_N[4]   = N5;
    d_length = N1 * N2 * N3 * N4 * N5;
}
inline ArraySize::ArraySize( std::initializer_list<size_t> N )
{
    d_ndim  = N.size();
    d_N[0]  = 0;
    d_N[1]  = 1;
    d_N[2]  = 1;
    d_N[3]  = 1;
    d_N[4]  = 1;
    auto it = N.begin();
    for ( size_t i = 0; i < d_ndim; i++, ++it )
        d_N[i] = *it;
    d_length = 1;
    for ( size_t i = 0; i < maxDim(); i++ )
        d_length *= d_N[i];
    if ( d_ndim == 0 )
        d_length = 0;
}
inline ArraySize::ArraySize( size_t ndim, const size_t *dims )
{
    d_ndim = ndim;
    d_N[0] = 0;
    d_N[1] = 1;
    d_N[2] = 1;
    d_N[3] = 1;
    d_N[4] = 1;
    for ( size_t i = 0; i < ndim; i++ )
        d_N[i] = dims[i];
    d_length = 1;
    for ( size_t i = 0; i < maxDim(); i++ )
        d_length *= d_N[i];
    if ( d_ndim == 0 )
        d_length = 0;
}
inline ArraySize::ArraySize( const std::vector<size_t> &N )
{
    d_ndim = N.size();
    d_N[0] = 0;
    d_N[1] = 1;
    d_N[2] = 1;
    d_N[3] = 1;
    d_N[4] = 1;
    for ( size_t i = 0; i < d_ndim; i++ )
        d_N[i] = N[i];
    d_length = 1;
    for ( size_t i = 0; i < maxDim(); i++ )
        d_length *= d_N[i];
    if ( d_ndim == 0 )
        d_length = 0;
}
inline ArraySize::ArraySize( const ArraySize &rhs ) { memcpy( this, &rhs, sizeof( *this ) ); }
inline ArraySize::ArraySize( ArraySize &&rhs ) { memcpy( this, &rhs, sizeof( *this ) ); }
inline ArraySize &ArraySize::operator=( const ArraySize &rhs )
{
    if ( this != &rhs )
        memcpy( this, &rhs, sizeof( *this ) );
    return *this;
}
inline ArraySize &ArraySize::operator=( ArraySize &&rhs )
{
    if ( this != &rhs )
        memcpy( this, &rhs, sizeof( *this ) );
    return *this;
}
inline void ArraySize::resize( uint8_t dim, size_t N )
{
    if ( dim >= d_ndim )
        throw std::out_of_range( "Invalid dimension" );
    d_N[dim] = N;
    d_length = 1;
    for ( size_t i = 0; i < maxDim(); i++ )
        d_length *= d_N[i];
}


/********************************************************
 *  Constructors                                         *
 ********************************************************/
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array()
{
    d_data = nullptr;
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( const ArraySize &N )
{
    allocate( N );
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( size_t N )
{
    allocate( ArraySize( N ) );
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( size_t N_rows, size_t N_cols )
{
    allocate( ArraySize( N_rows, N_cols ) );
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( size_t N1, size_t N2, size_t N3 )
{
    allocate( ArraySize( N1, N2, N3 ) );
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( size_t N1, size_t N2, size_t N3, size_t N4 )
{
    allocate( ArraySize( N1, N2, N3, N4 ) );
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( size_t N1, size_t N2, size_t N3, size_t N4, size_t N5 )
{
    allocate( ArraySize( N1, N2, N3, N4, N5 ) );
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( const std::vector<size_t> &N, const TYPE *data )
{
    allocate( N );
    if ( data ) {
        for ( size_t i = 0; i < d_size.length(); i++ )
            d_data[i] = data[i];
    }
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( const Range<TYPE> &range )
{
    double tmp = static_cast<double>( ( range.j - range.i ) ) / static_cast<double>( range.k );
    size_t N   = static_cast<size_t>( floor( tmp + 1e-12 ) + 1 );
    allocate( { N } );
    for ( size_t i = 0; i < N; i++ )
        d_data[i] = range.k * ( range.i / range.k + i );
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( std::initializer_list<TYPE> x )
{
    allocate( { x.size() } );
    auto it = x.begin();
    for ( size_t i = 0; i < x.size(); ++i, ++it )
        d_data[i] = *it;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::allocate( const ArraySize &N )
{
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
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( const Array &rhs ) : d_size( rhs.d_size ), d_data( nullptr )
{
    allocate( rhs.size() );
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] = rhs.d_data[i];
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::Array( Array &&rhs ) : d_size( rhs.d_size ), d_data( rhs.d_data )
{
    rhs.d_size = ArraySize();
    rhs.d_data = nullptr;
    d_ptr      = std::move( rhs.d_ptr );
}
template<class TYPE, class FUN>
Array<TYPE, FUN> &Array<TYPE, FUN>::operator=( const Array &rhs )
{
    if ( this == &rhs )
        return *this;
    this->allocate( rhs.size() );
    for ( size_t i = 0; i < d_size.length(); i++ )
        this->d_data[i] = rhs.d_data[i];
    return *this;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> &Array<TYPE, FUN>::operator=( Array &&rhs )
{
    if ( this == &rhs )
        return *this;
    d_size     = rhs.d_size;
    rhs.d_size = ArraySize();
    d_data     = rhs.d_data;
    rhs.d_data = nullptr;
    d_ptr      = std::move( rhs.d_ptr );
    return *this;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> &Array<TYPE, FUN>::operator=( const std::vector<TYPE> &rhs )
{
    this->allocate( ArraySize( rhs.size() ) );
    for ( size_t i = 0; i < rhs.size(); i++ )
        this->d_data[i] = rhs[i];
    return *this;
}
template<class TYPE, class FUN>
Array<TYPE, FUN>::~Array()
{
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::clear()
{
    d_size = ArraySize();
    d_ptr.reset();
    d_data = nullptr;
}


/********************************************************
 *  Access elements                                      *
 ********************************************************/


/********************************************************
 *  Copy/move values from one array to another (resize)  *
 ********************************************************/
template<class TYPE>
inline void moveValues( const ArraySize &N1, const ArraySize &N2, TYPE *data1, TYPE *data2 )
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
inline typename std::enable_if<test, void>::type copyValues(
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
inline typename std::enable_if<!test, void>::type copyValues(
    const ArraySize &, const ArraySize &, const TYPE *, TYPE * )
{
    throw std::logic_error( "No copy constructor" );
}


/********************************************************
 *  Resize the array                                     *
 ********************************************************/
template<class TYPE, class FUN>
void Array<TYPE, FUN>::resize( size_t N )
{
    resize( ArraySize( N ) );
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::resize( size_t N1, size_t N2 )
{
    resize( ArraySize( N1, N2 ) );
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::resize( size_t N1, size_t N2, size_t N3 )
{
    resize( ArraySize( N1, N2, N3 ) );
}

template<class TYPE, class FUN>
void Array<TYPE, FUN>::resize( const ArraySize &N )
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
template<class TYPE, class FUN>
void Array<TYPE, FUN>::resizeDim( int dim, size_t N, const TYPE &value )
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
template<class TYPE, class FUN>
void Array<TYPE, FUN>::reshape( const ArraySize &N )
{
    if ( N.length() != d_size.length() )
        throw std::logic_error( "reshape is not allowed to change the array size" );
    d_size = N;
}


/********************************************************
 *  Subset the array                                     *
 ********************************************************/
// Helper function to check subset indices
template<class TYPE, class FUN>
inline void Array<TYPE, FUN>::checkSubsetIndex( const std::vector<Range<size_t>> &range ) const
{
    bool test = (int) range.size() == d_size.ndim();
    for ( size_t d = 0; d < range.size(); d++ )
        test = test && range[d].i >= 0 && range[d].j <= d_size[d];
    if ( !test )
        throw std::logic_error( "indices for subset are invalid" );
}
template<class TYPE, class FUN>
inline std::vector<Range<size_t>> Array<TYPE, FUN>::convert(
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
template<class TYPE, class FUN>
inline void Array<TYPE, FUN>::getSubsetArrays( const std::vector<Range<size_t>> &index,
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
template<class TYPE, class FUN>
template<class TYPE2>
Array<TYPE2, FUN> Array<TYPE, FUN>::subset( const std::vector<Range<size_t>> &index ) const
{
    // Get the subset indicies
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( index, first, last, inc, N1 );
    ArraySize S1( d_size.ndim(), N1.data() );
    // Create the new array
    Array<TYPE2> subset_array( S1 );
    // Fill the new array
    static_assert( ArraySize::maxDim() == 5, "Not programmed for more than 5 dimensions" );
    TYPE2 *subset_data = subset_array.data();
    for ( size_t i4 = first[4], k1 = 0; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0], k1++ ) {
                        size_t k2       = d_size.index( i0, i1, i2, i3, i4 );
                        subset_data[k1] = static_cast<TYPE2>( d_data[k2] );
                    }
                }
            }
        }
    }
    return subset_array;
}
template<class TYPE, class FUN>
template<class TYPE2>
Array<TYPE2, FUN> Array<TYPE, FUN>::subset( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return subset( range );
}
template<class TYPE, class FUN>
template<class TYPE2>
void Array<TYPE, FUN>::copySubset(
    const std::vector<Range<size_t>> &index, const Array<TYPE2, FUN> &subset )
{
    // Get the subset indices
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, inc, N1;
    getSubsetArrays( index, first, last, inc, N1 );
    // Copy the sub-array
    static_assert( ArraySize::maxDim() == 5, "Not programmed for more than 5 dimensions" );
    const TYPE2 *src_data = subset.data();
    for ( size_t i4 = first[4], k1 = 0; i4 <= last[4]; i4 += inc[4] ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3 += inc[3] ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2 += inc[2] ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1 += inc[1] ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0 += inc[0], k1++ ) {
                        size_t k2  = d_size.index( i0, i1, i2, i3, i4 );
                        d_data[k2] = static_cast<TYPE>( src_data[k1] );
                    }
                }
            }
        }
    }
}

template<class TYPE, class FUN>
void Array<TYPE, FUN>::addSubset(
    const std::vector<Range<size_t>> &index, const Array<TYPE, FUN> &subset )
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
template<class TYPE, class FUN>
template<class TYPE2>
void Array<TYPE, FUN>::copySubset(
    const std::vector<size_t> &index, const Array<TYPE2, FUN> &subset )
{
    auto range = convert( index );
    copySubset( range, subset );
}

template<class TYPE, class FUN>
void Array<TYPE, FUN>::addSubset( const std::vector<size_t> &index, const Array<TYPE, FUN> &subset )
{
    auto range = convert( index );
    addSubset( range, subset );
}


/********************************************************
 *  Operator overloading                                 *
 ********************************************************/
template<class TYPE, class FUN>
bool Array<TYPE, FUN>::operator==( const Array &rhs ) const
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
template<class TYPE, class FUN>
std::shared_ptr<Array<TYPE, FUN>> Array<TYPE, FUN>::view(
    size_t N, std::shared_ptr<TYPE> const &data )
{
    return view( ArraySize( N ), data );
}
template<class TYPE, class FUN>
std::shared_ptr<Array<TYPE, FUN>> Array<TYPE, FUN>::view(
    size_t N1, size_t N2, std::shared_ptr<TYPE> const &data )
{
    return view( ArraySize( N1, N2 ), data );
}
template<class TYPE, class FUN>
std::shared_ptr<Array<TYPE, FUN>> Array<TYPE, FUN>::view(
    size_t N1, size_t N2, size_t N3, std::shared_ptr<TYPE> const &data )
{
    return view( ArraySize( N1, N2, N3 ), data );
}
template<class TYPE, class FUN>
std::shared_ptr<const Array<TYPE, FUN>> Array<TYPE, FUN>::constView(
    size_t N, std::shared_ptr<const TYPE> const &data )
{
    return constView( ArraySize( N ), data );
}
template<class TYPE, class FUN>
std::shared_ptr<const Array<TYPE, FUN>> Array<TYPE, FUN>::constView(
    size_t N1, size_t N2, std::shared_ptr<const TYPE> const &data )
{
    return constView( ArraySize( N1, N2 ), data );
}
template<class TYPE, class FUN>
std::shared_ptr<const Array<TYPE, FUN>> Array<TYPE, FUN>::constView(
    size_t N1, size_t N2, size_t N3, std::shared_ptr<const TYPE> const &data )
{
    return constView( ArraySize( N1, N2, N3 ), data );
}
template<class TYPE, class FUN>
std::shared_ptr<Array<TYPE, FUN>> Array<TYPE, FUN>::view(
    const ArraySize &N, std::shared_ptr<TYPE> const &data )
{
    std::shared_ptr<Array<TYPE, FUN>> array( new Array<TYPE, FUN>() );
    array->d_size = N;
    array->d_ptr  = data;
    array->d_data = array->d_ptr.get();
    return array;
}
template<class TYPE, class FUN>
std::shared_ptr<const Array<TYPE, FUN>> Array<TYPE, FUN>::constView(
    const ArraySize &N, std::shared_ptr<const TYPE> const &data )
{
    std::shared_ptr<Array<TYPE, FUN>> array( new Array<TYPE, FUN>() );
    array->d_size = N;
    array->d_ptr  = data;
    array->d_data = array->d_ptr.get();
    return array;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::view2( Array<TYPE, FUN> &src )
{
    view2( src.size(), src.getPtr() );
    d_data = src.d_data;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::view2( const ArraySize &N, std::shared_ptr<TYPE> const &data )
{
    d_size = N;
    d_ptr  = data;
    d_data = d_ptr.get();
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::viewRaw( int ndim, const size_t *dims, TYPE *data )
{
    d_size = ArraySize( ndim, dims );
    d_ptr.reset();
    d_data = data;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::viewRaw( const ArraySize &N, TYPE *data )
{
    d_size = N;
    d_ptr.reset();
    d_data = data;
}


/********************************************************
 *  Convert array types                                  *
 ********************************************************/
template<class TYPE, class FUN>
template<class TYPE2>
std::shared_ptr<Array<TYPE2>> Array<TYPE, FUN>::convert( std::shared_ptr<Array<TYPE, FUN>> array )
{
    if ( std::is_same<TYPE, TYPE2>() )
        return array;
    std::shared_ptr<Array<TYPE2>> array2( new Array<TYPE2>( array->size() ) );
    array2.copy( *array );
    return array2;
}
template<class TYPE, class FUN>
template<class TYPE2>
std::shared_ptr<const Array<TYPE2>> Array<TYPE, FUN>::convert(
    std::shared_ptr<const Array<TYPE, FUN>> array )
{
    return Array<TYPE, FUN>::convert( std::const_pointer_cast<Array<TYPE2>>( array ) );
}
template<class TYPE, class FUN>
template<class TYPE2>
void Array<TYPE, FUN>::copy( const Array<TYPE2> &array )
{
    resize( array.size() );
    const TYPE2 *src = array.data();
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] = static_cast<TYPE>( src[i] );
}
template<class TYPE, class FUN>
template<class TYPE2>
void Array<TYPE, FUN>::copy( const TYPE2 *src )
{
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] = static_cast<TYPE>( src[i] );
}
template<class TYPE, class FUN>
template<class TYPE2>
void Array<TYPE, FUN>::copyTo( TYPE2 *dst ) const
{
    for ( size_t i = 0; i < d_size.length(); i++ )
        dst[i] = static_cast<TYPE2>( d_data[i] );
}
template<class TYPE, class FUN>
template<class TYPE2>
Array<TYPE2, FUN> Array<TYPE, FUN>::cloneTo() const
{
    Array<TYPE2, FUN> dst( this->size() );
    auto dst_data = dst.data();
    for ( size_t i = 0; i < d_size.length(); i++ )
        dst_data[i] = static_cast<TYPE2>( d_data[i] );
    return dst;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::fill( const TYPE &value )
{
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] = value;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::scale( const TYPE &value )
{
    for ( size_t i = 0; i < d_size.length(); i++ )
        d_data[i] *= value;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::pow( const Array<TYPE, FUN> &baseArray, const TYPE &exp )
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
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::repmat( const std::vector<size_t> &N_rep ) const
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
    Array<TYPE, FUN> y( N2 );
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
template<class TYPE, class FUN>
bool Array<TYPE, FUN>::NaNs() const
{
    bool test = false;
    for ( size_t i = 0; i < d_size.length(); i++ )
        test = test || d_data[i] != d_data[i];
    return test;
}

template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::mean( void ) const
{
    TYPE x = this->sum() / d_size.length();
    return x;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::min( int dir ) const
{
    auto size_ans = d_size;
    size_ans.resize( dir, 1 );
    Array<TYPE, FUN> ans( size_ans );
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
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::max( int dir ) const
{
    auto size_ans = d_size;
    size_ans.resize( dir, 1 );
    Array<TYPE, FUN> ans( size_ans );
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
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::sum( int dir ) const
{
    auto size_ans = d_size;
    size_ans.resize( dir, 1 );
    Array<TYPE, FUN> ans( size_ans );
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
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::min( const std::vector<Range<size_t>> &range ) const
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
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::max( const std::vector<Range<size_t>> &range ) const
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
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::sum( const std::vector<Range<size_t>> &range ) const
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
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::mean( const std::vector<Range<size_t>> &range ) const
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
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::min( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return min( range );
}
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::max( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return max( range );
}
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::sum( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return sum( range );
}
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::mean( const std::vector<size_t> &index ) const
{
    auto range = convert( index );
    return mean( range );
}


/********************************************************
 *  Find all elements that match the given operation     *
 ********************************************************/
template<class TYPE, class FUN>
std::vector<size_t> Array<TYPE, FUN>::find(
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
template<class TYPE, class FUN>
void Array<TYPE, FUN>::print(
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
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::reverseDim() const
{
    size_t N2[ArraySize::maxDim()];
    for ( int d = 0; d < ArraySize::maxDim(); d++ )
        N2[d] = d_size[ArraySize::maxDim() - d - 1];
    ArraySize S2( ArraySize::maxDim(), N2 );
    Array<TYPE, FUN> y( S2 );
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
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::coarsen( const Array<TYPE, FUN> &filter ) const
{
  auto S2 = size();
  for ( size_t i = 0; i < S2.size(); i++ ) {
    S2.resize( i, S2[i] / filter.size(i) );
    if ( S2[i] * filter.size( i ) != size( i ) )
      throw std::invalid_argument( "Array must be multiple of filter size" );
  }
  Array<TYPE, FUN> y( S2 );
  if ( d_size.ndim() <= 3 )
    throw std::logic_error( "Function programmed for more than 3 dimensions" );
    const auto& Nh = filter.d_size;
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
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::coarsen(
    const std::vector<size_t> &ratio, std::function<TYPE( const Array<TYPE, FUN> & )> filter ) const
{
  //if ( ratio.size() != d_size.ndim() )
  //     throw std::logic_error( "ratio size does not match ndim" );
    auto S2 = size();
    for ( size_t i = 0; i < S2.size(); i++ ) {
        S2.resize( i, S2[i] / ratio[i] );
        if ( S2[i] * ratio[i] != size( i ) )
            throw std::invalid_argument( "Array must be multiple of filter size" );
    }
    Array<TYPE, FUN> tmp( ratio );
    Array<TYPE, FUN> y( S2 );
    if ( d_size.ndim() <= 3 )
        throw std::logic_error( "Function programmed for more than 3 dimensions" );
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
template<class TYPE, class FUN>
void Array<TYPE, FUN>::cat( const Array<TYPE, FUN> &x, int dim )
{
    std::vector<Array<TYPE, FUN>> tmp( 2 );
    tmp[0].view2( *this );
    tmp[1].view2( const_cast<Array<TYPE, FUN> &>( x ) );
    *this = cat( tmp, dim );
}
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::cat( const std::vector<Array> &x, int dim )
{
    if ( x.empty() )
        return Array<TYPE, FUN>();
    // Check that the dimensions match
    bool check = true;
    for ( size_t i = 1; i < x.size(); i++ ) {
        check = check && x[i].ndim() == x[0].ndim();
        for ( int d = 0; d < x[0].ndim(); d++ )
            check = check && d == dim;
    }
    if ( !check )
        throw std::logic_error( "Array dimensions do not match for concatenation" );
    // Create the output array
    auto size = x[0].d_size;
    for ( size_t i = 1; i < x.size(); i++ )
        size.resize( dim, size[dim] + x[i].size( dim ) );
    Array<TYPE, FUN> out( size );
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
 *  Math operations (should call the Math class)         *
 ********************************************************/
template<class TYPE, class FUN>
void Array<TYPE, FUN>::rand()
{
    FUN::rand( *this );
}
template<class TYPE, class FUN>
Array<TYPE, FUN> &Array<TYPE, FUN>::operator+=( const Array<TYPE, FUN> &rhs )
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a + b; };
    FUN::transform( fun, *this, rhs, *this );
    return *this;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> &Array<TYPE, FUN>::operator-=( const Array<TYPE, FUN> &rhs )
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a - b; };
    FUN::transform( fun, *this, rhs, *this );
    return *this;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> &Array<TYPE, FUN>::operator+=( const TYPE &rhs )
{
    const auto &fun = [rhs]( const TYPE &x ) { return x + rhs; };
    FUN::transform( fun, *this, *this );
    return *this;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> &Array<TYPE, FUN>::operator-=( const TYPE &rhs )
{
    const auto &fun = [rhs]( const TYPE &x ) { return x - rhs; };
    FUN::transform( fun, *this, *this );
    return *this;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> operator+( const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b )
{
    Array<TYPE, FUN> c;
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a + b; };
    FUN::transform( fun, a, b, c );
    return c;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> operator-( const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b )
{
    Array<TYPE, FUN> c;
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a - b; };
    FUN::transform( fun, a, b, c );
    return c;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> operator*( const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b )
{
    return Array<TYPE, FUN>::multiply( a, b );
}
template<class TYPE, class FUN>
inline Array<TYPE, FUN> operator*( const Array<TYPE, FUN> &a, const std::vector<TYPE> &b )
{
    Array<TYPE, FUN> b2;
    b2.viewRaw( { b.size() }, const_cast<TYPE *>( b.data() ) );
    return Array<TYPE, FUN>::multiply( a, b2 );
}
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::min() const
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a < b ? a : b; };
    return FUN::reduce( fun, *this );
}
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::max() const
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a > b ? a : b; };
    return FUN::reduce( fun, *this );
}
template<class TYPE, class FUN>
TYPE Array<TYPE, FUN>::sum() const
{
    const auto &fun = []( const TYPE &a, const TYPE &b ) { return a + b; };
    return FUN::reduce( fun, *this );
}
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::multiply( const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b )
{
    Array<TYPE, FUN> c;
    FUN::multiply( a, b, c );
    return c;
}
template<class TYPE, class FUN>
void Array<TYPE, FUN>::axpby( const TYPE &alpha, const Array<TYPE, FUN> &x, const TYPE &beta )
{
    const auto &fun = [alpha, beta](
                          const TYPE &x, const TYPE &y ) { return alpha * x + beta * y; };
    return FUN::transform( fun, x, *this );
}
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::transform(
    std::function<TYPE( const TYPE & )> fun, const Array<TYPE, FUN> &x )
{
    Array<TYPE, FUN> y;
    FUN::transform( fun, x, y );
    return y;
}
template<class TYPE, class FUN>
Array<TYPE, FUN> Array<TYPE, FUN>::transform( std::function<TYPE( const TYPE &, const TYPE & )> fun,
    const Array<TYPE, FUN> &x, const Array<TYPE, FUN> &y )
{
    Array<TYPE, FUN> z;
    FUN::transform( fun, x, y, z );
    return z;
}


#endif
