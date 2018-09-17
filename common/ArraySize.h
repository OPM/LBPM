#ifndef included_ArraySizeClass
#define included_ArraySizeClass


#include <array>
#include <cstring>
#include <initializer_list>
#include <vector>


#if defined( __CUDA_ARCH__ )
#include <cuda.h>
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
#if defined( USING_GCC ) || defined( USING_CLANG )
#define ARRAY_ATTRIBUTE HOST_DEVICE __attribute__( ( always_inline ) )
#else
#define ARRAY_ATTRIBUTE HOST_DEVICE
#endif


#if ( defined( DEBUG ) || defined( _DEBUG ) ) && !defined( NDEBUG )
#define CHECK_ARRAY_LENGTH( i )                                      \
    do {                                                             \
        if ( i >= d_length )                                         \
            throw std::out_of_range( "Index exceeds array bounds" ); \
    } while ( 0 )
#else
#define CHECK_ARRAY_LENGTH( i ) \
    do {                        \
    } while ( 0 )
#endif


// Forward declerations
class FunctionTable;
template<class TYPE, class FUN = FunctionTable, class Allocator = std::nullptr_t>
class Array;


//! Simple range class
template<class TYPE = size_t>
class Range final
{
public:
    //! Empty constructor
    constexpr Range() : i( 0 ), j( -1 ), k( 1 ) {}

    /*!
     * Create a range i:k:j (or i:j)
     * @param i_            Starting value
     * @param j_            Ending value
     * @param k_            Increment value
     */
    constexpr Range( TYPE i_, TYPE j_, TYPE k_ = 1 ) : i( i_ ), j( j_ ), k( k_ ) {}

    TYPE i, j, k;
};


//! Simple class to store the array dimensions
class ArraySize final
{
public:
    //! Empty constructor
    constexpr ArraySize() : d_ndim( 1 ), d_length( 0 ), d_N{ 0, 1, 1, 1, 1 } {}

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     */
    constexpr ArraySize( size_t N1 ) : d_ndim( 1 ), d_length( N1 ), d_N{ N1, 1, 1, 1, 1 } {}

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     */
    constexpr ArraySize( size_t N1, size_t N2 )
        : d_ndim( 2 ), d_length( N1 * N2 ), d_N{ N1, N2, 1, 1, 1 }
    {
    }

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     */
    constexpr ArraySize( size_t N1, size_t N2, size_t N3 )
        : d_ndim( 3 ), d_length( N1 * N2 * N3 ), d_N{ N1, N2, N3, 1, 1 }
    {
    }

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     */
    constexpr ArraySize( size_t N1, size_t N2, size_t N3, size_t N4 )
        : d_ndim( 4 ), d_length( N1 * N2 * N3 * N4 ), d_N{ N1, N2, N3, N4, 1 }
    {
    }

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     * @param N5            Number of elements in the fifth dimension
     */
    constexpr ArraySize( size_t N1, size_t N2, size_t N3, size_t N4, size_t N5 )
        : d_ndim( 5 ), d_length( N1 * N2 * N3 * N4 * N5 ), d_N{ N1, N2, N3, N4, N5 }
    {
    }

    /*!
     * Create from initializer list
     * @param N             Size of the array
     */
    constexpr ArraySize( std::initializer_list<size_t> N )
        : d_ndim( N.size() ), d_length( 0 ), d_N{ 0, 1, 1, 1, 1 }
    {
        if ( d_ndim > maxDim() )
            throw std::out_of_range( "Maximum number of dimensions exceeded" );
        auto it = N.begin();
        for ( size_t i = 0; i < d_ndim; i++, ++it )
            d_N[i] = *it;
        d_length = 1;
        for ( unsigned long i : d_N )
            d_length *= i;
        if ( d_ndim == 0 )
            d_length = 0;
    }


    /*!
     * Create from raw pointer
     * @param ndim          Number of dimensions
     * @param dims          Dimensions
     */
    constexpr ArraySize( size_t ndim, const size_t *dims )
        : d_ndim( ndim ), d_length( 0 ), d_N{ 0, 1, 1, 1, 1 }
    {
        if ( d_ndim > maxDim() )
            throw std::out_of_range( "Maximum number of dimensions exceeded" );
        for ( size_t i = 0; i < ndim; i++ )
            d_N[i] = dims[i];
        d_length = 1;
        for ( unsigned long i : d_N )
            d_length *= i;
        if ( d_ndim == 0 )
            d_length = 0;
    }

    /*!
     * Create from std::vector
     * @param N             Size of the array
     */
    ArraySize( const std::vector<size_t> &N );

    // Copy/assignment constructors
    constexpr ArraySize( ArraySize &&rhs )      = default;
    constexpr ArraySize( const ArraySize &rhs ) = default;
    constexpr ArraySize &operator=( ArraySize &&rhs ) = default;
    constexpr ArraySize &operator=( const ArraySize &rhs ) = default;

    /*!
     * Access the ith dimension
     * @param i             Index to access
     */
    constexpr ARRAY_ATTRIBUTE size_t operator[]( size_t i ) const { return d_N[i]; }

    //! Return the number of dimensions
    constexpr ARRAY_ATTRIBUTE uint8_t ndim() const { return d_ndim; }

    //! Return the number of dimensions
    constexpr ARRAY_ATTRIBUTE size_t size() const { return d_ndim; }

    //! Return the total number of elements in the array
    constexpr ARRAY_ATTRIBUTE size_t length() const { return d_length; }

    //! Resize the dimension
    constexpr void resize( uint8_t dim, size_t N )
    {
        if ( dim >= d_ndim )
            throw std::out_of_range( "Invalid dimension" );
        d_N[dim] = N;
        d_length = 1;
        for ( unsigned long i : d_N )
            d_length *= i;
    }

    /*!
     * Reshape the Array so that the number of dimensions is the
     *    max of ndim and the largest dim>1.
     * @param ndim          Desired number of dimensions
     */
    constexpr void setNdim( uint8_t ndim ) { d_ndim = std::max( ndim, d_ndim ); }

    //! Returns an iterator to the beginning
    constexpr const size_t *begin() const { return d_N; }

    //! Returns an iterator to the end
    constexpr const size_t *end() const { return d_N + d_ndim; }

    // Check if two array sizes are equal
    constexpr ARRAY_ATTRIBUTE bool operator==( const ArraySize &rhs ) const
    {
        return d_ndim == rhs.d_ndim && memcmp( d_N, rhs.d_N, sizeof( d_N ) ) == 0;
    }

    // Check if two array sizes are equal (ignoring the dimension)
    constexpr ARRAY_ATTRIBUTE bool approxEqual( const ArraySize &rhs ) const
    {
        return ( length() == 0 && rhs.length() == 0 ) || memcmp( d_N, rhs.d_N, sizeof( d_N ) ) == 0;
    }

    //! Check if two matrices are not equal
    constexpr ARRAY_ATTRIBUTE bool operator!=( const ArraySize &rhs ) const
    {
        return d_ndim != rhs.d_ndim || memcmp( d_N, rhs.d_N, sizeof( d_N ) ) != 0;
    }

    //! Maximum supported dimension
    constexpr ARRAY_ATTRIBUTE static uint8_t maxDim() { return 5u; }

    //! Get the index
    constexpr ARRAY_ATTRIBUTE size_t index( size_t i ) const
    {
        CHECK_ARRAY_LENGTH( i );
        return i;
    }

    //! Get the index
    constexpr ARRAY_ATTRIBUTE size_t index( size_t i1, size_t i2 ) const
    {
        size_t index = i1 + i2 * d_N[0];
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

    //! Get the index
    constexpr ARRAY_ATTRIBUTE size_t index( size_t i1, size_t i2, size_t i3 ) const
    {
        size_t index = i1 + d_N[0] * ( i2 + d_N[1] * i3 );
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

    //! Get the index
    constexpr ARRAY_ATTRIBUTE size_t index( size_t i1, size_t i2, size_t i3, size_t i4 ) const
    {
        size_t index = i1 + d_N[0] * ( i2 + d_N[1] * ( i3 + d_N[2] * i4 ) );
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

    //! Get the index
    constexpr ARRAY_ATTRIBUTE size_t index(
        size_t i1, size_t i2, size_t i3, size_t i4, size_t i5 ) const
    {
        size_t index = i1 + d_N[0] * ( i2 + d_N[1] * ( i3 + d_N[2] * ( i4 + d_N[3] * i5 ) ) );
        CHECK_ARRAY_LENGTH( index );
        return index;
    }

private:
    uint8_t d_ndim;
    size_t d_length;
    size_t d_N[5];
};


// Function to concatenate dimensions of two array sizes
constexpr ArraySize cat( const ArraySize &x, const ArraySize &y )
{
    if ( x.ndim() + y.ndim() > ArraySize::maxDim() )
        throw std::out_of_range( "Maximum number of dimensions exceeded" );
    size_t N[ArraySize::maxDim()] = { 0 };
    for ( int i = 0; i < x.ndim(); i++ )
        N[i] = x[i];
    for ( int i = 0; i < y.ndim(); i++ )
        N[i + x.ndim()] = y[i];
    return ArraySize( x.ndim() + y.ndim(), N );
}


#endif
