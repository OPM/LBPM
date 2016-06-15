#ifndef included_ArrayClass
#define included_ArrayClass

#include <vector>
#include <array>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <iostream>


#define ARRAY_NDIM_MAX 5 // Maximum number of dimensions supported


#define GET_ARRAY_INDEX3D( N, i1, i2, i3 ) i1 + N[0] * ( i2 + N[1] * i3 )
#define GET_ARRAY_INDEX4D( N, i1, i2, i3, i4 ) i1 + N[0] * ( i2 + N[1] * ( i3 + N[2] * i4 ) )
#define GET_ARRAY_INDEX5D( N, i1, i2, i3, i4, i5 ) i1 + N[0] * ( i2 + N[1] * ( i3 + N[2] * ( i4 + N[3] * i5 ) ) )

#if defined( DEBUG ) || defined( _DEBUG )
    #define CHECK_ARRAY_INDEX3D( N, i1, i2, i3 )                  \
        if ( GET_ARRAY_INDEX3D( N, i1, i2, i3 ) < 0 || GET_ARRAY_INDEX3D( N, i1, i2, i3 ) >= d_length ) \
            throw std::logic_error( "Index exceeds array bounds" );
    #define CHECK_ARRAY_INDEX4D( N, i1, i2, i3, i4 )              \
        if ( GET_ARRAY_INDEX4D( N, i1, i2, i3, i4 ) < 0 ||        \
             GET_ARRAY_INDEX4D( N, i1, i2, i3, i4 ) >= d_length ) \
            throw std::logic_error( "Index exceeds array bounds" );
#else
    #define CHECK_ARRAY_INDEX3D( N, i1, i2, i3 )
    #define CHECK_ARRAY_INDEX4D( N, i1, i2, i3, i4 )
#endif


#if defined( __CUDA_ARCH__ )
#include <cuda.h>
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif


/*!
 * Class Array is a multi-dimensional array class written by Mark Berrill
 */
template <class TYPE>
class Array
{
public:
    /*!
     * Create a new empty Array
     */
    Array();

    /*!
     * Create a new 1D Array with the given number of elements
     * @param N             Number of elements in the array
     */
    explicit Array( size_t N );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     */
    explicit Array( size_t N_rows, size_t N_columns );

    /*!
     * Create a new 3D Array with the given number of rows and columns
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     */
    explicit Array( size_t N1, size_t N2, size_t N3 );

    /*!
     * Create a multi-dimensional Array with the given number of elements
     * @param N             Number of elements in each dimension
     * @param data          Optional raw array to copy the src data
     */
    explicit Array( const std::vector<size_t> &N, const TYPE *data = NULL );

    /*!
     * Copy constructor
     * @param rhs           Array to copy
     */
    Array( const Array &rhs );

    /*!
     * Move constructor
     * @param rhs           Array to copy
     */
    Array( Array &&rhs );

    /*!
     * Assignment operator
     * @param rhs           Array to copy
     */
    Array &operator=( const Array &rhs );

    /*!
     * Move assignment operator
     * @param rhs           Array to copy
     */
    Array &operator=( Array &&rhs );

    /*!
     * Assignment operator
     * @param rhs           std::vector to copy
     */
    Array &operator=( const std::vector<TYPE> &rhs );


    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( size_t N, std::shared_ptr<TYPE> const &data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view(
        size_t N_rows, size_t N_columns, std::shared_ptr<TYPE> const &data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view(
        size_t N1, size_t N2, size_t N3, std::shared_ptr<TYPE> const &data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view(
        const std::vector<size_t> &N, std::shared_ptr<TYPE> const &data );


    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView(
        size_t N, std::shared_ptr<const TYPE> const &data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView(
        size_t N_rows, size_t N_columns, std::shared_ptr<const TYPE> const &data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView(
        size_t N1, size_t N2, size_t N3, std::shared_ptr<const TYPE> const &data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView(
        const std::vector<size_t> &N, std::shared_ptr<const TYPE> const &data );


    /*!
     * Make this object a view of the src
     * @param src           Source vector to take the view of
     */
    void view2( Array &src );

    /*!
     * Make this object a view of the data
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void view2( const std::vector<size_t> &N, std::shared_ptr<TYPE> const &data );

    /*!
     * Make this object a view of the raw data (expert use only).
     * Use view2( N, std::shared_ptr(data,[](TYPE*){}) ) instead.
     *   Note: this interface is not recommended as it does not protect from
     *   the src data being deleted while still being used by the Array.
     *   Additionally for maximum performance it does not set the internal shared_ptr
     *   so functions like getPtr and resize will not work correctly.
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void viewRaw( const std::initializer_list<size_t> &N, TYPE *data );

    /*!
     * Make this object a view of the raw data (expert use only).
     * Use view2( N, std::shared_ptr(data,[](TYPE*){}) ) instead.
     *   Note: this interface is not recommended as it does not protect from
     *   the src data being deleted while still being used by the Array.
     *   Additionally for maximum performance it does not set the internal shared_ptr
     *   so functions like getPtr and resize will not work correctly.
     * @param N             Number of elements in each dimension
     * @param data          Pointer to the data
     */
    void viewRaw( const std::vector<size_t> &N, TYPE *data );

    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template <class TYPE2>
    static std::shared_ptr<Array<TYPE2>> convert( std::shared_ptr<Array<TYPE>> array );


    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template <class TYPE2>
    static std::shared_ptr<const Array<TYPE2>> convert( std::shared_ptr<const Array<TYPE>> array );


    /*!
     * Copy and convert data from another array to this array
     * @param array         Source array
     */
    template <class TYPE2>
    void copy( const Array<TYPE2> &array );

    /*!
     * Copy and convert data from a raw vector to this array.
     *    Note: The current array must be allocated to the proper size first.
     * @param array         Source array
     */
    template <class TYPE2>
    void copy( const TYPE2 *array );

    /*!
     * Copy and convert data from this array to a raw vector.
     * @param array         Source array
     */
    template <class TYPE2>
    void copyTo( TYPE2 *array ) const;


    /*!
     * Fill the array with the given value
     * @param value         Value to fill
     */
    void fill( const TYPE &value );

    /*!
     * Scale the array by the given value
     * @param scale         Value to scale by
     */
    void scale( const TYPE &scale );

    /*!
     * Set the values of this array to pow(base, exp)
     * @param base        Base array
     * @param exp         Exponent value
     */
    void pow( const Array<TYPE> &baseArray, const TYPE &exp );

    //! Destructor
    ~Array();


    //! Clear the data in the array
    void clear();


    //! Return the size of the Array
    inline int ndim() const { return d_ndim; }


    //! Return the size of the Array
    inline std::vector<size_t> size() const { return std::vector<size_t>( d_N, d_N + d_ndim ); }


    //! Return the size of the Array
    inline size_t size( int d ) const { return d_N[d]; }


    //! Return the size of the Array
    inline size_t length() const { return d_length; }


    //! Return true if the Array is empty
    inline bool empty() const { return d_length == 0; }


    /*!
     * Resize the Array
     * @param N             NUmber of elements
     */
    void resize( size_t N );

    /*!
     * Resize the Array
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     */
    void resize( size_t N_rows, size_t N_columns );

    /*!
     * Resize the Array
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     */
    void resize( size_t N1, size_t N2, size_t N3 );

    /*!
     * Resize the Array
     * @param N             Number of elements in each dimension
     */
    void resize( const std::vector<size_t> &N );

    /*!
     * Resize the given dimension of the array
     * @param dim           The dimension to resize
     * @param N             Number of elements for the given dimension
     * @param value         Value to initialize new elements
     */
    void resizeDim( int dim, size_t N, const TYPE &value );


    /*!
     * Reshape the Array (total size of array will not change)
     * @param N             Number of elements in each dimension
     */
    void reshape( const std::vector<size_t> &N );


    /*!
     * Subset the Array (total size of array will not change)
     * @param index         Index to subset (imin,imax,jmin,jmax,kmin,kmax,...)
     */
    template<class TYPE2=TYPE>
    Array<TYPE2> subset( const std::vector<size_t> &index ) const;

    /*!
     * Copy data from an array into a subset of this array
     * @param index         Index of the subset (imin,imax,jmin,jmax,kmin,kmax,...)
     * @param subset        The subset array to copy from
     */
    template <class TYPE2>
    void copySubset( const std::vector<size_t> &index, const Array<TYPE2> &subset );

    /*!
     * Add data from an array into a subset of this array
     * @param index         Index of the subset (imin,imax,jmin,jmax,kmin,kmax,...)
     * @param subset        The subset array to add from
     */
    void addSubset( const std::vector<size_t> &index, const Array<TYPE> &subset );


    /*!
     * Access the desired element
     * @param i             The row index
     */
    HOST_DEVICE inline TYPE &operator()( size_t i )
    {
        CHECK_ARRAY_INDEX3D( d_N, i, 0, 0 ) return d_data[i];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     */
    HOST_DEVICE inline const TYPE &operator()( size_t i ) const
    {
        CHECK_ARRAY_INDEX3D( d_N, i, 0, 0 ) return d_data[i];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    HOST_DEVICE inline TYPE &operator()( size_t i, size_t j )
    {
        CHECK_ARRAY_INDEX3D( d_N, i, j, 0 ) return d_data[i + j * d_N[0]];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    HOST_DEVICE inline const TYPE &operator()( size_t i, size_t j ) const
    {
        CHECK_ARRAY_INDEX3D( d_N, i, j, 0 ) return d_data[i + j * d_N[0]];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     */
    HOST_DEVICE inline TYPE &operator()( size_t i, size_t j, size_t k )
    {
        CHECK_ARRAY_INDEX3D( d_N, i, j, k ) return d_data[GET_ARRAY_INDEX3D( d_N, i, j, k )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     */
    HOST_DEVICE inline const TYPE &operator()( size_t i, size_t j, size_t k ) const
    {
        CHECK_ARRAY_INDEX3D( d_N, i, j, k ) return d_data[GET_ARRAY_INDEX3D( d_N, i, j, k )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     * @param l             The fourth index
     */
    HOST_DEVICE inline TYPE &operator()( size_t i, size_t j, size_t k, size_t l )
    {
        CHECK_ARRAY_INDEX4D( d_N, i, j, k, l ) return d_data[GET_ARRAY_INDEX4D( d_N, i, j, k, l )];
    }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     * @param k             The third index
     * @param l             The fourth index
     */
    HOST_DEVICE inline const TYPE &operator()( size_t i, size_t j, size_t k, size_t l ) const
    {
        CHECK_ARRAY_INDEX4D( d_N, i, j, k, l ) return d_data[GET_ARRAY_INDEX4D( d_N, i, j, k, l )];
    }

    //! Check if two matrices are equal
    // Equality means the dimensions and data have to be identical
    bool operator==( const Array &rhs ) const;

    //! Check if two matrices are not equal
    inline bool operator!=( const Array &rhs ) const { return !this->operator==( rhs ); }


    //! Return the pointer to the raw data
    inline std::shared_ptr<TYPE> getPtr() { return d_ptr; }

    //! Return the pointer to the raw data
    inline std::shared_ptr<const TYPE> getPtr() const { return d_ptr; }

    //! Return the pointer to the raw data
    HOST_DEVICE inline TYPE *data() { return d_data; }

    //! Return the pointer to the raw data
    HOST_DEVICE inline const TYPE *data() const { return d_data; }


    //! Return true if NaNs are present
    inline bool NaNs() const;

    //! Return the smallest value
    inline TYPE min() const;

    //! Return the largest value
    inline TYPE max() const;

    //! Return the sum of all elements
    inline TYPE sum() const;

    //! Return the mean of all elements
    inline TYPE mean() const;

    //! Return the min of all elements in a given direction
    Array<TYPE> min( int dir ) const;

    //! Return the max of all elements in a given direction
    Array<TYPE> max( int dir ) const;

    //! Return the sum of all elements in a given direction
    Array<TYPE> sum( int dir ) const;

    //! Return the smallest value
    inline TYPE min( const std::vector<size_t> &index ) const;

    //! Return the largest value
    inline TYPE max( const std::vector<size_t> &index ) const;

    //! Return the sum of all elements
    inline TYPE sum( const std::vector<size_t> &index ) const;

    //! Return the mean of all elements
    inline TYPE mean( const std::vector<size_t> &index ) const;

    //! Find all elements that match the operator
    std::vector<size_t> find(
        const TYPE &value, std::function<bool( const TYPE &, const TYPE & )> compare ) const;

    //! Add another array
    Array &operator+=( const Array &rhs );

    //! Subtract another array
    Array &operator-=( const Array &rhs );

    //! Add a scalar
    Array &operator+=( const TYPE &rhs );

    //! Subtract a scalar
    Array &operator-=( const TYPE &rhs );

    //! Print an array
    void print( std::ostream& os, const std::string& name="A", const std::string& prefix="" ) const;

    //! Multiply two arrays
    static Array multiply( const Array& a, const Array& b );

    //! Transpose an array
    Array<TYPE> reverseDim( ) const;

    //! Coarsen an array using the given filter
    Array<TYPE> coarsen( const Array<TYPE>& filter ) const;

    //! Coarsen an array using the given filter
    Array<TYPE> coarsen( const std::vector<size_t>& ratio, std::function<TYPE(const Array<TYPE>&)> filter ) const;

private:
    int d_ndim;                  // Number of dimensions in array
    size_t d_N[ARRAY_NDIM_MAX];  // Size of each dimension
    size_t d_length;             // Total length of array
    TYPE *d_data;                // Raw pointer to data in array
    std::shared_ptr<TYPE> d_ptr; // Shared pointer to data in array
    void allocate( const std::vector<size_t> &N );

private:
    template<class TYPE2>
    inline bool sizeMatch( const Array<TYPE2>& rhs ) const;
    inline void checkSubsetIndex( const std::vector<size_t> &index ) const;
    inline std::array<size_t, 5> getDimArray() const;
    static inline void getSubsetArrays( const std::vector<size_t> &index,
        std::array<size_t, 5> &first, std::array<size_t, 5> &last, std::array<size_t, 5> &N );
};


typedef Array<double> DoubleArray;
typedef Array<float> FloatArray;
typedef Array<int> IntArray;


#include "common/Array.hpp"

#endif
