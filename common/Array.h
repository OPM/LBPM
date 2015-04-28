#ifndef included_ArrayClass
#define included_ArrayClass

#include <iostream>
#include <memory>
#include <vector>
#include "shared_ptr.h"
#include "common/Utilities.h"


#define GET_ARRAY_INDEX(i1,i2,i3,i4) i1+d_N[0]*(i2+d_N[1]*(i3+d_N[2]*i4))
#ifdef DEBUG
    #define CHECK_ARRAY_INDEX(i1,i2,i3,i4) \
        if ( GET_ARRAY_INDEX(i1,i2,i3,i4)>d_length ) \
            ERROR("Index exceeds array bounds");
#else
    #define CHECK_ARRAY_INDEX(i1,i2,i3,i4) 
#endif


/*!
 * Class Array is a simple array class
 */
template<class TYPE>
class Array
{
public:

    /*!
     * Create a new empty Array
     */
    Array( );

    /*!
     * Create a new 1D Array with the given number of elements
     * @param N             Number of elements in the array
     */
    Array( size_t N );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     */
    Array( size_t N_rows, size_t N_columns );

    /*!
     * Create a new 3D Array with the given number of rows and columns
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     */
    Array( size_t N1, size_t N2, size_t N3 );

    /*!
     * Create a multi-dimensional Array with the given number of elements
     * @param N             Number of elements in each dimension
     */
    Array( const std::vector<size_t>& N );

    /*!
     * Copy constructor
     * @param rhs           Array to copy
     */
    Array( const Array& rhs );

    /*!
     * Assignment operator
     * @param rhs           Array to copy
     */
    Array& operator=( const Array& rhs );
    
    
    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( size_t N, std::shared_ptr<TYPE> data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( size_t N_rows, size_t N_columns, std::shared_ptr<TYPE> data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<Array> view( size_t N1, size_t N2, size_t N3, std::shared_ptr<TYPE> data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     */
    static std::shared_ptr<Array> view( const std::vector<size_t>& N, std::shared_ptr<TYPE> data );

   
    
    /*!
     * Create a 1D Array view to a raw block of data
     * @param N             Number of elements in the array
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView( size_t N, std::shared_ptr<const TYPE> data );

    /*!
     * Create a new 2D Array with the given number of rows and columns
     * @param N_rows        Number of rows
     * @param N_columns     Number of columns
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView( size_t N_rows, size_t N_columns, std::shared_ptr<const TYPE> data );

    /*!
     * Create a new 3D Array view to a raw block of data
     * @param N1            Number of rows
     * @param N2            Number of columns
     * @param N3            Number of elements in the third dimension
     * @param data          Pointer to the data
     */
    static std::shared_ptr<const Array> constView( size_t N1, size_t N2, size_t N3, std::shared_ptr<const TYPE> data );

    /*!
     * Create a multi-dimensional Array view to a raw block of data
     * @param N             Number of elements in each dimension
     */
    static std::shared_ptr<const Array> constView( const std::vector<size_t>& N, std::shared_ptr<const TYPE> data );


    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template<class TYPE2>
    static std::shared_ptr<Array<TYPE2> > convert( std::shared_ptr<Array<TYPE> > array );


    /*!
     * Convert an array of one type to another.  This may or may not allocate new memory.
     * @param array         Input array
     */
    template<class TYPE2>
    static std::shared_ptr<const Array<TYPE2> > convert( std::shared_ptr<const Array<TYPE> > array );


    /*!
     * Copy and convert data from another array to this array
     * @param array         Source array
     */
    template<class TYPE2>
    void copy( const Array<TYPE2>& array );

    /*!
     * Copy and convert data from a raw vector to this array.
     *    Note: The current array must be allocated to the proper size first.
     * @param array         Source array
     */
    template<class TYPE2>
    void copy( const TYPE2* array );

    /*!
     * Fill the array with the given value
     * @param value         Value to fill
     */
    void fill( const TYPE& value );


    //! Destructor
    ~Array( );


    //! Return the size of the Array
    inline int ndim( ) const { return d_ndim; }


    //! Return the size of the Array
    inline std::vector<size_t> size( ) const { return std::vector<size_t>(d_N,d_N+d_ndim); }


    //! Return the size of the Array
    inline size_t size( int d ) const { return d_N[d]; }


    //! Return the size of the Array
    inline size_t length( ) const { return d_length; }


    //! Return true if the Array is empty
    inline bool empty( ) const { return d_length==0; }


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
    void resize( const std::vector<size_t>& N );


    /*!
     * Reshape the Array (total size of array will not change)
     * @param N             Number of elements in each dimension
     */
    void reshape( const std::vector<size_t>& N );


    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline TYPE& operator()( size_t i ) { CHECK_ARRAY_INDEX(i,0,0,0) return d_data[i]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline const TYPE& operator()( size_t i ) const { CHECK_ARRAY_INDEX(i,0,0,0) return d_data[i]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline TYPE& operator()( size_t i, size_t j ) { CHECK_ARRAY_INDEX(i,j,0,0) return d_data[i+j*d_N[0]]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline const TYPE& operator()( size_t i, size_t j ) const { CHECK_ARRAY_INDEX(i,j,0,0) return d_data[i+j*d_N[0]]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline TYPE& operator()( size_t i, size_t j, size_t k ) { CHECK_ARRAY_INDEX(i,j,k,0) return d_data[GET_ARRAY_INDEX(i,j,k,0)]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline const TYPE& operator()( size_t i, size_t j, size_t k ) const { CHECK_ARRAY_INDEX(i,j,k,0) return d_data[GET_ARRAY_INDEX(i,j,k,0)]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline TYPE& operator()( size_t i, size_t j, size_t k, size_t m ) { CHECK_ARRAY_INDEX(i,j,k,m) return d_data[GET_ARRAY_INDEX(i,j,k,m)]; }

    /*!
     * Access the desired element
     * @param i             The row index
     * @param j             The column index
     */
    inline const TYPE& operator()( size_t i, size_t j, size_t k, size_t m ) const { CHECK_ARRAY_INDEX(i,j,k,m) return d_data[GET_ARRAY_INDEX(i,j,k,m)]; }


    //! Check if two matricies are equal
    bool operator==( const Array& rhs ) const;

    //! Check if two matricies are not equal
    inline bool operator!=( const Array& rhs ) const { return !this->operator==(rhs); }


    //! Return the pointer to the raw data
    inline std::shared_ptr<TYPE> getPtr( ) { return d_ptr; }

    //! Return the pointer to the raw data
    inline std::shared_ptr<const TYPE> getPtr( ) const { return d_ptr; }

    //! Return the pointer to the raw data
    inline TYPE* get( ) { return d_data; }

    //! Return the pointer to the raw data
    inline const TYPE* get( ) const { return d_data; }


    //! Return true if NaNs are present
    inline bool NaNs( ) const;

    //! Return the smallest value
    inline TYPE min( ) const;

    //! Return the largest value
    inline TYPE max( ) const;

    //! Return the sum of all elements
    inline TYPE sum( ) const;

    //! Return the sum of all elements in a given direction
    std::shared_ptr<Array<TYPE> > sum( int dir ) const;


private:
    int d_ndim;
    size_t d_N[4];
    size_t d_length;
    TYPE *d_data;
    std::shared_ptr<TYPE> d_ptr;
    void allocate( const std::vector<size_t>& N );
};


typedef Array<int> IntArray;
typedef Array<double> DoubleArray;


#include "common/Array.hpp"

#endif

