#ifndef included_ArrayClass_hpp
#define included_ArrayClass_hpp

#include "common/Array.h"
#include "common/Utilities.h"
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>



/********************************************************
*  Constructors                                         *
********************************************************/
template <class TYPE>
Array<TYPE>::Array()
{
    d_ndim   = 1;
    d_length = 0;
    for ( size_t i = 0; i < ARRAY_NDIM_MAX; i++ )
        d_N[i]     = 1;
    d_N[0]         = 0;
    d_data         = nullptr;
}
template <class TYPE>
Array<TYPE>::Array( size_t N )
{
    allocate( std::vector<size_t>( 1, N ) );
}
template <class TYPE>
Array<TYPE>::Array( size_t N_rows, size_t N_columns )
{
    std::vector<size_t> N( 2 );
    N[0] = N_rows;
    N[1] = N_columns;
    allocate( N );
}
template <class TYPE>
Array<TYPE>::Array( size_t N1, size_t N2, size_t N3 )
{
    std::vector<size_t> N( 3 );
    N[0] = N1;
    N[1] = N2;
    N[2] = N3;
    allocate( N );
}
template <class TYPE>
Array<TYPE>::Array( const std::vector<size_t> &N, const TYPE *data )
{
    allocate( N );
    if ( data != NULL ) {
        for ( size_t i = 0; i < d_length; i++ )
            d_data[i]  = data[i];
    }
}
template <class TYPE>
void Array<TYPE>::allocate( const std::vector<size_t> &N )
{
    d_ndim   = static_cast<int>( N.size() );
    d_length = 1;
    for ( size_t i = 0; i < ARRAY_NDIM_MAX; i++ )
        d_N[i]     = 1;
    for ( size_t i = 0; i < N.size(); i++ ) {
        d_N[i] = N[i];
        d_length *= N[i];
    }
    if ( N.empty() ) {
        d_N[0]   = 0;
        d_length = 0;
    }
    if ( d_length == 0 )
        d_ptr.reset();
    else
        d_ptr.reset( new ( std::nothrow ) TYPE[d_length], []( TYPE *p ) { delete[] p; } );
    d_data = d_ptr.get();
    if ( d_length > 0 && d_data == nullptr )
        throw std::logic_error( "Failed to allocate array" );
}
template <class TYPE>
Array<TYPE>::Array( const Array &rhs )
    : d_ndim( rhs.d_ndim ), d_length( rhs.d_length ), d_data( nullptr )
{
    allocate( rhs.size() );
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i]  = rhs.d_data[i];
}
template <class TYPE>
Array<TYPE>::Array( Array &&rhs )
    : d_ndim( rhs.d_ndim ), d_length( rhs.d_length ), d_data( rhs.d_data )
{
    rhs.d_ndim = 0;
    memcpy( d_N, rhs.d_N, sizeof( rhs.d_N ) );
    memset( rhs.d_N, 0, sizeof( rhs.d_N ) );
    rhs.d_length = 0;
    rhs.d_data   = nullptr;
    d_ptr        = std::move( rhs.d_ptr );
}
template <class TYPE>
Array<TYPE> &Array<TYPE>::operator=( const Array &rhs )
{
    if ( this == &rhs )
        return *this;
    this->allocate( rhs.size() );
    for ( size_t i      = 0; i < d_length; i++ )
        this->d_data[i] = rhs.d_data[i];
    return *this;
}
template <class TYPE>
Array<TYPE> &Array<TYPE>::operator=( Array &&rhs )
{
    if ( this == &rhs )
        return *this;
    d_ndim     = rhs.d_ndim;
    rhs.d_ndim = 0;
    memcpy( d_N, rhs.d_N, sizeof( rhs.d_N ) );
    memset( rhs.d_N, 0, sizeof( rhs.d_N ) );
    d_length     = rhs.d_length;
    rhs.d_length = 0;
    d_data       = rhs.d_data;
    rhs.d_data   = nullptr;
    d_ptr        = std::move( rhs.d_ptr );
    return *this;
}
template <class TYPE>
Array<TYPE> &Array<TYPE>::operator=( const std::vector<TYPE> &rhs )
{
    this->allocate( std::vector<size_t>( 1, rhs.size() ) );
    for ( size_t i      = 0; i < rhs.size(); i++ )
        this->d_data[i] = rhs[i];
    return *this;
}
template <class TYPE>
Array<TYPE>::~Array()
{
}
template <class TYPE>
void Array<TYPE>::clear()
{
    d_ndim   = 0;
    d_length = 0;
    for ( size_t i = 0; i < ARRAY_NDIM_MAX; i++ )
        d_N[i]     = 1;
    d_N[0]         = 0;
    d_ptr.reset();
    d_data = nullptr;
}


/********************************************************
*  Check if the size of the array matches rhs           *
********************************************************/
template <class TYPE>
template <class TYPE2>
bool Array<TYPE>::sizeMatch( const Array<TYPE2>& rhs ) const
{
    bool test = d_ndim == rhs.d_ndim;
    for ( int d = 0; d < d_ndim; d++ )
        test = test && d_N[d] == rhs.d_N[d];
    return test;
}


/********************************************************
*  Resize the array                                     *
********************************************************/
template <class TYPE>
void Array<TYPE>::resize( size_t N )
{
    resize( std::vector<size_t>{N} );
}
template <class TYPE>
void Array<TYPE>::resize( size_t N1, size_t N2 )
{
    resize( std::vector<size_t>{N1,N2} );
}
template <class TYPE>
void Array<TYPE>::resize( size_t N1, size_t N2, size_t N3 )
{
    resize( std::vector<size_t>{N1,N2,N3} );
}
template <class TYPE>
void Array<TYPE>::resize( const std::vector<size_t> &N )
{
    // Check if the array actually changed size
    size_t new_length = 1;
    for ( size_t i = 0; i < N.size(); i++ )
        new_length *= N[i];
    if ( N.empty() )
        new_length = 0;
    bool changed = new_length != d_length || (int) N.size() != d_ndim;
    for ( size_t i = 0; i < N.size(); i++ )
        changed = changed || N[i] != d_N[i];
    if ( !changed )
        return;
// Store the old data
#if ARRAY_NDIM_MAX > 5
#error Function programmed for more than 5 dimensions
#endif
    std::array<size_t, 5> N1{ { 1, 1, 1, 1, 1 } };
    std::array<size_t, 5> N2{ { 1, 1, 1, 1, 1 } };
    for ( int d = 0; d < d_ndim; d++ )
        N1[d]   = d_N[d];
    for ( size_t d = 0; d < N.size(); d++ )
        N2[d]      = N[d];
    if ( d_ndim == 0 ) {
        N1[0] = 0;
    }
    if ( N.empty() ) {
        N2[0] = 0;
    }
    std::shared_ptr<TYPE> old_data = d_ptr;
    // Allocate new data
    allocate( N );
    // Copy the old values
    if ( d_length > 0 ) {
        TYPE *data1 = old_data.get();
        TYPE *data2 = d_data;
        if ( old_data.unique() ) {
            // We own the data, use std:move
            for ( size_t i5 = 0; i5 < std::min( N1[4], N2[4] ); i5++ ) {
                for ( size_t i4 = 0; i4 < std::min( N1[3], N2[3] ); i4++ ) {
                    for ( size_t i3 = 0; i3 < std::min( N1[2], N2[2] ); i3++ ) {
                        for ( size_t i2 = 0; i2 < std::min( N1[1], N2[1] ); i2++ ) {
                            for ( size_t i1 = 0; i1 < std::min( N1[0], N2[0] ); i1++ ) {
                                size_t index1 = GET_ARRAY_INDEX5D( N1, i1, i2, i3, i4, i5 );
                                size_t index2 = GET_ARRAY_INDEX5D( N2, i1, i2, i3, i4, i5 );
                                data2[index2] = std::move( data1[index1] );
                            }
                        }
                    }
                }
            }
        } else {
            // We do not own the data, copy
            for ( size_t i5 = 0; i5 < std::min( N1[4], N2[4] ); i5++ ) {
                for ( size_t i4 = 0; i4 < std::min( N1[3], N2[3] ); i4++ ) {
                    for ( size_t i3 = 0; i3 < std::min( N1[2], N2[2] ); i3++ ) {
                        for ( size_t i2 = 0; i2 < std::min( N1[1], N2[1] ); i2++ ) {
                            for ( size_t i1 = 0; i1 < std::min( N1[0], N2[0] ); i1++ ) {
                                size_t index1 = GET_ARRAY_INDEX5D( N1, i1, i2, i3, i4, i5 );
                                size_t index2 = GET_ARRAY_INDEX5D( N2, i1, i2, i3, i4, i5 );
                                data2[index2] = data1[index1];
                            }
                        }
                    }
                }
            }
        }
    }
}
template <class TYPE>
void Array<TYPE>::resizeDim( int dim, size_t N, const TYPE &value )
{
    if ( dim >= d_ndim )
        throw std::logic_error( "Invalid dimension" );
    std::vector<size_t> N2 = size();
    size_t N0              = N2[dim];
    N2[dim]                = N;
    resize( N2 );
    size_t n1 = 1, n2 = 1;
    for ( int d = 0; d < dim; d++ )
        n1 *= N2[d];
    for ( size_t d = dim + 1; d < N2.size(); d++ )
        n2 *= N2[d];
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
template <class TYPE>
void Array<TYPE>::reshape( const std::vector<size_t> &N )
{
    size_t new_length = 1;
    for ( size_t i = 0; i < N.size(); i++ )
        new_length *= N[i];
    if ( new_length != d_length )
        throw std::logic_error( "reshape is not allowed to change the array size" );
    d_ndim = N.size();
    for ( size_t i = 0; i < ARRAY_NDIM_MAX; i++ )
        d_N[i]     = 1;
    for ( size_t i = 0; i < N.size(); i++ )
        d_N[i]     = N[i];
}


/********************************************************
*  Subset the array                                     *
********************************************************/
// clang-format off
// Helper function to check subset indices
template <class TYPE>
inline void Array<TYPE>::checkSubsetIndex( const std::vector<size_t> &index ) const
{
    bool test = index.size() % 2 == 0 && (int) index.size() / 2 <= d_ndim;
    for ( size_t d = 0; d < index.size() / 2; d++ )
        test = test && index[2 * d + 0] < d_N[d] && index[2 * d + 1] < d_N[d];
    if ( !test )
        throw std::logic_error( "indices for subset are invalid" );
}
// Helper function to return dimensions as a std::array for hard coded loops
template <class TYPE>
inline std::array<size_t,5> Array<TYPE>::getDimArray() const
{
    #if ARRAY_NDIM_MAX > 5
        #error Function programmed for more than 5 dimensions
    #endif
    std::array<size_t,5> N{ { 1, 1, 1, 1, 1 } };
    for ( int d = 0; d < d_ndim; d++ )
        N[d] = d_N[d];
    return N;
}
// Helper function to return dimensions for the subset array
template <class TYPE>
inline void Array<TYPE>::getSubsetArrays( const std::vector<size_t> &index,
                                          std::array<size_t, 5> &first,
                                          std::array<size_t, 5> &last,
                                          std::array<size_t, 5> &N )
{
    #if ARRAY_NDIM_MAX > 5
        #error Function programmed for more than 5 dimensions
    #endif
    size_t ndim = index.size() / 2;
    for ( size_t d = 0; d < ndim; d++ ) {
        first[d] = index[2 * d + 0];
        last[d]  = index[2 * d + 1];
        N[d]     = last[d] - first[d] + 1;
    }
    for ( size_t d = ndim; d < 5; d++ ) {
        first[d] = 0;
        last[d]  = 0;
        N[d]     = 1;
    }
}
template <class TYPE>
template <class TYPE2>
Array<TYPE2> Array<TYPE>::subset( const std::vector<size_t> &index ) const
{
    // Get the subset indicies
    checkSubsetIndex( index );
    std::array<size_t,5> first, last, N1;
    getSubsetArrays( index, first, last, N1 );
    std::array<size_t,5> N2 = getDimArray();
    // Create the new array
    std::vector<size_t> dim( d_ndim );
    for ( int d = 0; d < d_ndim; d++ )
        dim[d]  = last[d] - first[d] + 1;
    Array<TYPE2> subset( dim );
    // Fill the new array
    #if ARRAY_NDIM_MAX > 5
        #error Function programmed for more than 5 dimensions
    #endif
    TYPE2 *subset_data = subset.data();
    for (size_t i4=first[4]; i4<=last[4]; i4++) {
        for (size_t i3=first[3]; i3<=last[3]; i3++) {
            for (size_t i2=first[2]; i2<=last[2]; i2++) {
                for (size_t i1=first[1]; i1<=last[1]; i1++) {
                    for (size_t i0=first[0]; i0<=last[0]; i0++) {
                        size_t k1 = GET_ARRAY_INDEX5D( N1, i0-first[0], 
                            i1-first[1], i2-first[2], i3-first[3], i4-first[4] );
                        size_t k2       = GET_ARRAY_INDEX5D( N2, i0, i1, i2, i3, i4 );
                        subset_data[k1] = static_cast<TYPE2>( d_data[k2] );
                    }
                }
            }
        }
    }
    return subset;
}
template <class TYPE>
template <class TYPE2>
void Array<TYPE>::copySubset( const std::vector<size_t> &index, const Array<TYPE2> &subset )
{
    // Get the subset indices
    checkSubsetIndex( index );
    std::array<size_t,5> first, last, N1;
    getSubsetArrays( index, first, last, N1 );
    std::array<size_t,5> N2 = getDimArray();
    // Copy the sub-array
    #if ARRAY_NDIM_MAX > 5
        #error Function programmed for more than 5 dimensions
    #endif
    const TYPE2 *src_data = subset.data();
    for (size_t i4=first[4]; i4<=last[4]; i4++) {
        for (size_t i3=first[3]; i3<=last[3]; i3++) {
            for (size_t i2=first[2]; i2<=last[2]; i2++) {
                for (size_t i1=first[1]; i1<=last[1]; i1++) {
                    for (size_t i0=first[0]; i0<=last[0]; i0++) {
                        size_t k1 = GET_ARRAY_INDEX5D( N1, i0-first[0], 
                            i1-first[1], i2-first[2], i3-first[3], i4-first[4] );
                        size_t k2  = GET_ARRAY_INDEX5D( N2, i0, i1, i2, i3, i4 );
                        d_data[k2] = static_cast<TYPE>( src_data[k1] );
                    }
                }
            }
        }
    }
}

template <class TYPE>
void Array<TYPE>::addSubset( const std::vector<size_t> &index, const Array<TYPE> &subset )
{
    // Get the subset indices
    checkSubsetIndex( index );
    std::array<size_t,5> first, last, N1;
    getSubsetArrays( index, first, last, N1 );
    std::array<size_t,5> N2 = getDimArray();
    // add the sub-array
    #if ARRAY_NDIM_MAX > 5
        #error Function programmed for more than 5 dimensions
    #endif
    for (size_t i4=first[4]; i4<=last[4]; i4++) {
        for (size_t i3=first[3]; i3<=last[3]; i3++) {
            for (size_t i2=first[2]; i2<=last[2]; i2++) {
                for (size_t i1=first[1]; i1<=last[1]; i1++) {
                    for (size_t i0=first[0]; i0<=last[0]; i0++) {
                        size_t k1 = GET_ARRAY_INDEX5D( N1, i0-first[0], 
                            i1-first[1], i2-first[2], i3-first[3], i4-first[4] );
                        size_t k2  = GET_ARRAY_INDEX5D( N2, i0, i1, i2, i3, i4 );
                        d_data[k2] += subset.d_data[k1];
                    }
                }
            }
        }
    }
}
// clang-format on


/********************************************************
*  Operator overloading                                 *
********************************************************/
template <class TYPE>
bool Array<TYPE>::operator==( const Array &rhs ) const
{
    if ( this == &rhs )
        return true;
    if ( d_length != rhs.d_length )
        return false;
    if ( d_ndim != rhs.d_ndim )
        return false;
    for ( int d = 0; d < d_ndim; d++ ) {
        if ( d_N[d] != rhs.d_N[d] )
            return false;
    }
    bool match = true;
    for ( size_t i = 0; i < d_length; i++ )
        match = match && d_data[i] == rhs.d_data[i];
    return match;
}


/********************************************************
*  Get a view of an C array                             *
********************************************************/
template <class TYPE>
std::shared_ptr<Array<TYPE>> Array<TYPE>::view( size_t N, std::shared_ptr<TYPE> const &data )
{
    return view( std::vector<size_t>{N}, data );
}
template <class TYPE>
std::shared_ptr<Array<TYPE>> Array<TYPE>::view(
    size_t N1, size_t N2, std::shared_ptr<TYPE> const &data )
{
    return view( std::vector<size_t>{N1,N2}, data );
}
template <class TYPE>
std::shared_ptr<Array<TYPE>> Array<TYPE>::view(
    size_t N1, size_t N2, size_t N3, std::shared_ptr<TYPE> const &data )
{
    return view( std::vector<size_t>{N1,N2,N3}, data );
}
template <class TYPE>
std::shared_ptr<const Array<TYPE>> Array<TYPE>::constView(
    size_t N, std::shared_ptr<const TYPE> const &data )
{
    return constView( std::vector<size_t>{N}, data );
}
template <class TYPE>
std::shared_ptr<const Array<TYPE>> Array<TYPE>::constView(
    size_t N1, size_t N2, std::shared_ptr<const TYPE> const &data )
{
    return constView( std::vector<size_t>{N1,N2}, data );
}
template <class TYPE>
std::shared_ptr<const Array<TYPE>> Array<TYPE>::constView(
    size_t N1, size_t N2, size_t N3, std::shared_ptr<const TYPE> const &data )
{
    return constView( std::vector<size_t>{N1,N2,N3}, data );
}
template <class TYPE>
std::shared_ptr<Array<TYPE>> Array<TYPE>::view(
    const std::vector<size_t> &N, std::shared_ptr<TYPE> const &data )
{
    std::shared_ptr<Array<TYPE>> array( new Array<TYPE>() );
    array->d_ndim   = N.size();
    array->d_length = 1;
    for ( size_t i = 0; i < N.size(); i++ ) {
        array->d_N[i] = N[i];
        array->d_length *= N[i];
    }
    if ( array->d_ndim == 0 )
        array->d_length = 0;
    array->d_ptr        = data;
    array->d_data       = array->d_ptr.get();
    return array;
}
template <class TYPE>
std::shared_ptr<const Array<TYPE>> Array<TYPE>::constView(
    const std::vector<size_t> &N, std::shared_ptr<const TYPE> const &data )
{
    return view( N, std::const_pointer_cast<TYPE>( data ) );
}
template <class TYPE>
void Array<TYPE>::view2( Array<TYPE> &src )
{
    view2( src.size(), src.getPtr() );
    d_data = src.d_data;
}
template <class TYPE>
void Array<TYPE>::view2( const std::vector<size_t> &N, std::shared_ptr<TYPE> const &data )
{
    d_ndim = static_cast<int>( N.size() );
    for ( size_t i = 0; i < ARRAY_NDIM_MAX; i++ ) {
        d_N[i] = 1;
    }
    d_length = d_ndim == 0 ? 0 : 1;
    for ( size_t i = 0; i < N.size(); i++ ) {
        d_N[i] = N[i];
        d_length *= d_N[i];
    }
    d_ptr  = data;
    d_data = d_ptr.get();
}

template <class TYPE>
void Array<TYPE>::viewRaw( const std::initializer_list<size_t> &N, TYPE *data )
{
    d_ndim = static_cast<int>( N.size() );
    for ( size_t i = 0; i < ARRAY_NDIM_MAX; i++ ) {
        d_N[i] = 1;
    }
    d_length = d_ndim == 0 ? 0 : 1;
    size_t i = 0;
    for ( auto it = N.begin(); it != N.end(); ++it, ++i ) {
        d_N[i] = *it;
        d_length *= *it;
    }
    d_ptr.reset();
    d_data = data;
}
template <class TYPE>
void Array<TYPE>::viewRaw( const std::vector<size_t> &N, TYPE *data )
{
    d_ndim = static_cast<int>( N.size() );
    for ( size_t i = 0; i < ARRAY_NDIM_MAX; i++ ) {
        d_N[i] = 1;
    }
    d_length = d_ndim == 0 ? 0 : 1;
    size_t i = 0;
    for ( auto it = N.begin(); it != N.end(); ++it, ++i ) {
        d_N[i] = *it;
        d_length *= *it;
    }
    d_ptr.reset();
    d_data = data;
}


/********************************************************
*  Convert array types                                  *
********************************************************/
template <class TYPE>
template <class TYPE2>
std::shared_ptr<Array<TYPE2>> Array<TYPE>::convert( std::shared_ptr<Array<TYPE>> array )
{
    if ( std::is_same<TYPE, TYPE2>() )
        return array;
    std::shared_ptr<Array<TYPE2>> array2( new Array<TYPE2>( array->size() ) );
    array2.copy( *array );
    return array2;
}
template <class TYPE>
template <class TYPE2>
std::shared_ptr<const Array<TYPE2>> Array<TYPE>::convert( std::shared_ptr<const Array<TYPE>> array )
{
    return Array<TYPE>::convert( std::const_pointer_cast<Array<TYPE2>>( array ) );
}
template <class TYPE>
template <class TYPE2>
void Array<TYPE>::copy( const Array<TYPE2> &array )
{
    resize( array.size() );
    const TYPE2 *src = array.data();
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i]  = static_cast<TYPE>( src[i] );
}
template <class TYPE>
template <class TYPE2>
void Array<TYPE>::copy( const TYPE2 *src )
{
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i]  = static_cast<TYPE>( src[i] );
}
template <class TYPE>
template <class TYPE2>
void Array<TYPE>::copyTo( TYPE2 *dst ) const
{
    for ( size_t i = 0; i < d_length; i++ )
        dst[i]     = static_cast<TYPE2>( d_data[i] );
}
template <class TYPE>
void Array<TYPE>::fill( const TYPE &value )
{
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i]  = value;
}
template <class TYPE>
void Array<TYPE>::scale( const TYPE &value )
{
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i] *= value;
}
template <class TYPE>
    void Array<TYPE>::pow(const Array<TYPE> &baseArray, const TYPE &exp )
{
    // not insisting on the shapes being the same
    // but insisting on the total size being the same
    AMP_ASSERT(d_length==baseArray.length());

    const auto base_data = baseArray.data();
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i]  = std::pow(base_data[i], exp);
}

/********************************************************
*  Simple math operations                               *
********************************************************/
template <class TYPE>
bool Array<TYPE>::NaNs() const
{
    bool test = false;
    for ( size_t i = 0; i < d_length; i++ )
        test = test || d_data[i] != d_data[i];
    return test;
}
template <class TYPE>
TYPE Array<TYPE>::min() const
{
    TYPE x = std::numeric_limits<TYPE>::max();
    for ( size_t i = 0; i < d_length; i++ )
        x = std::min( x, d_data[i] );
    return x;
}
template <class TYPE>
TYPE Array<TYPE>::max() const
{
    TYPE x = std::numeric_limits<TYPE>::min();
    for ( size_t i = 0; i < d_length; i++ )
        x = std::max( x, d_data[i] );
    return x;
}
template <class TYPE>
TYPE Array<TYPE>::sum() const
{
    TYPE x = 0;
    for ( size_t i = 0; i < d_length; i++ )
        x += d_data[i];
    return x;
}
template <class TYPE>
TYPE Array<TYPE>::mean( void ) const
{
    TYPE x = sum() / d_length;
    return x;
}
template <class TYPE>
Array<TYPE> Array<TYPE>::min( int dir ) const
{
    std::vector<size_t> size_ans = size();
    size_ans[dir]                = 1;
    Array<TYPE> ans( size_ans );
    size_t N1 = 1, N2 = 1, N3 = 1;
    for ( int d = 0; d < std::min( dir, d_ndim ); d++ )
        N1 *= d_N[d];
    N2 = d_N[dir];
    for ( int d = dir + 1; d < std::min( d_ndim, ARRAY_NDIM_MAX ); d++ )
        N3 *= d_N[d];
    TYPE *data2 = ans.d_data;
    for ( size_t i3 = 0; i3 < N3; i3++ ) {
        for ( size_t i1 = 0; i1 < N1; i1++ ) {
            TYPE x = d_data[i1 + i3 * N1 * N2];
            for ( size_t i2 = 0; i2 < N2; i2++ )
                x               = std::min( x, d_data[i1 + i2 * N1 + i3 * N1 * N2] );
            data2[i1 + i3 * N1] = x;
        }
    }
    return ans;
}
template <class TYPE>
Array<TYPE> Array<TYPE>::max( int dir ) const
{
    std::vector<size_t> size_ans = size();
    size_ans[dir]                = 1;
    Array<TYPE> ans( size_ans );
    size_t N1 = 1, N2 = 1, N3 = 1;
    for ( int d = 0; d < std::min( dir, d_ndim ); d++ )
        N1 *= d_N[d];
    N2 = d_N[dir];
    for ( int d = dir + 1; d < std::min( d_ndim, ARRAY_NDIM_MAX ); d++ )
        N3 *= d_N[d];
    TYPE *data2 = ans.d_data;
    for ( size_t i3 = 0; i3 < N3; i3++ ) {
        for ( size_t i1 = 0; i1 < N1; i1++ ) {
            TYPE x = d_data[i1 + i3 * N1 * N2];
            for ( size_t i2 = 0; i2 < N2; i2++ )
                x               = std::max( x, d_data[i1 + i2 * N1 + i3 * N1 * N2] );
            data2[i1 + i3 * N1] = x;
        }
    }
    return ans;
}
template <class TYPE>
Array<TYPE> Array<TYPE>::sum( int dir ) const
{
    std::vector<size_t> size_ans = size();
    size_ans[dir]                = 1;
    Array<TYPE> ans( size_ans );
    size_t N1 = 1, N2 = 1, N3 = 1;
    for ( int d = 0; d < std::min( dir, d_ndim ); d++ )
        N1 *= d_N[d];
    N2 = d_N[dir];
    for ( int d = dir + 1; d < std::min( d_ndim, ARRAY_NDIM_MAX ); d++ )
        N3 *= d_N[d];
    TYPE *data2 = ans.d_data;
    for ( int i3 = 0; i3 < N3; i3++ ) {
        for ( int i1 = 0; i1 < N1; i1++ ) {
            TYPE x = 0;
            for ( size_t i2 = 0; i2 < N2; i2++ )
                x += d_data[i1 + i2 * N1 + i3 * N1 * N2];
            data2[i1 + i3 * N1] = x;
        }
    }
    return ans;
}
template <class TYPE>
TYPE Array<TYPE>::min( const std::vector<size_t> &index ) const
{
    // Get the subset indicies
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, N1;
    getSubsetArrays( index, first, last, N1 );
    std::array<size_t, 5> N2 = getDimArray();
#if ARRAY_NDIM_MAX > 5
#error Function programmed for more than 5 dimensions
#endif
    TYPE x = std::numeric_limits<TYPE>::max();
    for ( size_t i4 = first[4]; i4 <= last[4]; i4++ ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3++ ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2++ ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1++ ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0++ ) {
                        size_t k1 = GET_ARRAY_INDEX5D( N2, i0, i1, i2, i3, i4 );
                        x         = std::min( x, d_data[k1] );
                    }
                }
            }
        }
    }

    return x;
}
template <class TYPE>
TYPE Array<TYPE>::max( const std::vector<size_t> &index ) const
{
    // Get the subset indicies
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, N1;
    getSubsetArrays( index, first, last, N1 );
    std::array<size_t, 5> N2 = getDimArray();
#if ARRAY_NDIM_MAX > 5
#error Function programmed for more than 5 dimensions
#endif
    TYPE x = std::numeric_limits<TYPE>::min();
    for ( size_t i4 = first[4]; i4 <= last[4]; i4++ ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3++ ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2++ ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1++ ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0++ ) {
                        size_t k1 = GET_ARRAY_INDEX5D( N2, i0, i1, i2, i3, i4 );
                        x         = std::max( x, d_data[k1] );
                    }
                }
            }
        }
    }
    return x;
}
template <class TYPE>
TYPE Array<TYPE>::sum( const std::vector<size_t> &index ) const
{
    // Get the subset indicies
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, N1;
    getSubsetArrays( index, first, last, N1 );
    std::array<size_t, 5> N2 = getDimArray();
#if ARRAY_NDIM_MAX > 5
#error Function programmed for more than 5 dimensions
#endif
    TYPE x = 0;
    for ( size_t i4 = first[4]; i4 <= last[4]; i4++ ) {
        for ( size_t i3 = first[3]; i3 <= last[3]; i3++ ) {
            for ( size_t i2 = first[2]; i2 <= last[2]; i2++ ) {
                for ( size_t i1 = first[1]; i1 <= last[1]; i1++ ) {
                    for ( size_t i0 = first[0]; i0 <= last[0]; i0++ ) {
                        size_t k1 = GET_ARRAY_INDEX5D( N2, i0, i1, i2, i3, i4 );
                        x += d_data[k1];
                    }
                }
            }
        }
    }
    return x;
}
template <class TYPE>
TYPE Array<TYPE>::mean( const std::vector<size_t> &index ) const
{
    // Get the subset indicies
    checkSubsetIndex( index );
    std::array<size_t, 5> first, last, N1;
    getSubsetArrays( index, first, last, N1 );
#if ARRAY_NDIM_MAX > 5
#error Function programmed for more than 5 dimensions
#endif
    size_t n = 1;
    for ( auto &d : N1 )
        n *= d;
    TYPE x = sum( index ) / n;
    return x;
}

template <class TYPE>
Array<TYPE> &Array<TYPE>::operator+=( const Array<TYPE> &rhs )
{
    if ( !sizeMatch(rhs) )
        throw std::logic_error( "Array don't match" );
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i] += rhs.d_data[i];
    return *this;
}
template <class TYPE>
Array<TYPE> &Array<TYPE>::operator-=( const Array<TYPE> &rhs )
{
    if ( !sizeMatch(rhs) )
        throw std::logic_error( "Array don't match" );
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i] -= rhs.d_data[i];
    return *this;
}
template <class TYPE>
Array<TYPE> &Array<TYPE>::operator+=( const TYPE &rhs )
{
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i] += rhs;
    return *this;
}
template <class TYPE>
Array<TYPE> &Array<TYPE>::operator-=( const TYPE &rhs )
{
    for ( size_t i = 0; i < d_length; i++ )
        d_data[i] -= rhs;
    return *this;
}
template <class TYPE>
Array<TYPE> operator+( const Array<TYPE>& a, const Array<TYPE>& b )
{
    Array<TYPE> c = a;
    c += b;
    return c;
}
template <class TYPE>
Array<TYPE> operator-( const Array<TYPE>& a, const Array<TYPE>& b )
{
    Array<TYPE> c = a;
    c -= b;
    return c;
}
template <class TYPE>
Array<TYPE> operator*( const Array<TYPE>& a, const Array<TYPE>& b )
{
    return Array<TYPE>::multiply(a,b);
}
template <class TYPE>
Array<TYPE> Array<TYPE>::multiply( const Array<TYPE>& a, const Array<TYPE>& b )
{
    Array<TYPE> c;
    if ( a.d_ndim==2 && b.d_ndim==2 ) {
        c.resize( a.size(0), b.size(1) );
        c.fill(0);
        for (size_t k=0; k<b.size(1); k++) {
            for (size_t j=0; j<a.size(1); j++) {
                for (size_t i=0; i<a.size(0); i++) {
                    c(i,k) += a(i,j) * b(j,k);
                }
            }
        }
    } else {
        throw std::logic_error("Not finished yet");
    }
    return c;
}


/********************************************************
*  Find all elements that match the given operation     *
********************************************************/
template <class TYPE>
std::vector<size_t> Array<TYPE>::find(
    const TYPE &value, std::function<bool( const TYPE &, const TYPE & )> compare ) const
{
    std::vector<size_t> result;
    result.reserve( d_length );
    for ( size_t i = 0; i < d_length; i++ ) {
        if ( compare( d_data[i], value ) )
            result.push_back( i );
    }
    return result;
}


/********************************************************
*  Print an array to an output stream                   *
********************************************************/
template <class TYPE>
void Array<TYPE>::print( std::ostream& os, const std::string& name, const std::string& prefix ) const
{
    if ( d_ndim==1 ) {
        for (size_t i=0; i<d_N[0]; i++)
            os << prefix << name << "[" << i << "] = " << operator()(i) << std::endl;
    } else if ( d_ndim==2 ) {
        os << prefix << name << ":" << std::endl;
        for (size_t i=0; i<d_N[0]; i++) {
            for (size_t j=0; j<d_N[1]; j++)
                os << prefix << "  " << operator()(i,j);
            os << std::endl;
        }
    } else {
        throw std::logic_error("Not programmed for this dimension");
    }
}


/********************************************************
*  Reverse dimensions (transpose)                       *
********************************************************/
template <class TYPE>
Array<TYPE> Array<TYPE>::reverseDim( ) const
{
    std::vector<size_t> N2(ARRAY_NDIM_MAX);
    for ( int d=0; d<ARRAY_NDIM_MAX; d++)
        N2[d] = d_N[ARRAY_NDIM_MAX-d-1];
    Array<TYPE> y( N2 );
#if ARRAY_NDIM_MAX != 5
    #error Function programmed for dimensions other than 5
#endif
    TYPE* y2 = y.data();
    for (size_t i0=0; i0<d_N[0]; i0++) {
        for (size_t i1=0; i1<d_N[1]; i1++) {
            for (size_t i2=0; i2<d_N[2]; i2++) {
                for (size_t i3=0; i3<d_N[3]; i3++) {
                    for (size_t i4=0; i4<d_N[4]; i4++) {
                        y2[GET_ARRAY_INDEX5D(N2,i4,i3,i2,i1,i0)] = d_data[GET_ARRAY_INDEX5D(d_N,i0,i1,i2,i3,i4)];
                    }
                }
            }
        }
    }
    auto S2 = size();
    for ( int d=0; d<d_ndim; d++)
        S2[d] = size(d_ndim-d-1);
    y.reshape( S2 );
    return y;
}


/********************************************************
*  Coarsen the array                                    *
********************************************************/
template <class TYPE>
Array<TYPE> Array<TYPE>::coarsen( const Array<TYPE>& filter ) const
{
    auto S2 = size();
    for (size_t i=0; i<S2.size(); i++) {
        S2[i] /= filter.size(i);
        INSIST(S2[i]*filter.size(i)==size(i),"Array must be multiple of filter size");
    }
    Array<TYPE> y( S2 );
    INSIST(d_ndim<=3,"Function programmed for more than 5 dimensions");
    const size_t *Nh = filter.d_N;
    for (size_t k1=0; k1<y.d_N[2]; k1++) {
        for (size_t j1=0; j1<y.d_N[1]; j1++) {
            for (size_t i1=0; i1<y.d_N[0]; i1++) {
                TYPE tmp = 0;
                for (size_t k2=0; k2<Nh[2]; k2++) {
                    for (size_t j2=0; j2<Nh[1]; j2++) {
                        for (size_t i2=0; i2<Nh[0]; i2++) {
                            tmp += filter(i2,j2,k2) * this->operator()(i1*Nh[0]+i2,j1*Nh[1]+j2,k1*Nh[2]+k2);
                        }
                    }
                }
                y(i1,j1,k1) = tmp;
            }
        }
    }
    return y;
}
template <class TYPE>
Array<TYPE> Array<TYPE>::coarsen( const std::vector<size_t>& ratio, std::function<TYPE(const Array<TYPE>&)> filter ) const
{
    ASSERT((int)ratio.size()==d_ndim);
    auto S2 = size();
    for (size_t i=0; i<S2.size(); i++) {
        S2[i] /= ratio[i];
        INSIST(S2[i]*ratio[i]==size(i),"Array must be multiple of filter size");
    }
    Array<TYPE> tmp(ratio);
    TYPE* tmp2 = tmp.data();
    Array<TYPE> y( S2 );
    INSIST(d_ndim<=3,"Function programmed for more than 3 dimensions");
    for (size_t k1=0; k1<y.d_N[2]; k1++) {
        for (size_t j1=0; j1<y.d_N[1]; j1++) {
            for (size_t i1=0; i1<y.d_N[0]; i1++) {
                for (size_t k2=0; k2<ratio[2]; k2++) {
                    for (size_t j2=0; j2<ratio[1]; j2++) {
                        for (size_t i2=0; i2<ratio[0]; i2++) {
                            tmp2[GET_ARRAY_INDEX3D(tmp.d_N,i2,j2,k2)] = this->operator()(i1*ratio[0]+i2,j1*ratio[1]+j2,k1*ratio[2]+k2);
                        }
                    }
                }
                y(i1,j1,k1) = filter(tmp);
            }
        }
    }
    return y;
}


#endif
