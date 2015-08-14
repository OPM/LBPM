#ifndef included_ArrayClass_hpp
#define included_ArrayClass_hpp

#include "common/Array.h"
#include "common/Utilities.h"
#include <algorithm>
#include <limits>


template<class TYPE>
void DeleteArray( TYPE* x )
{
    delete [] x;
}


/********************************************************
*  Constructors                                         *
********************************************************/
template<class TYPE>
Array<TYPE>::Array( )
{
    d_ndim = 0;
    d_length = 0;
    for (size_t i=0; i<sizeof(d_N)/sizeof(size_t); i++)
        d_N[i] = 1;
    d_N[0] = 0;
    d_data = d_ptr.get();
}
template<class TYPE>
Array<TYPE>::Array( size_t N )
{
    allocate(std::vector<size_t>(1,N));
}
template<class TYPE>
Array<TYPE>::Array( size_t N_rows, size_t N_columns )
{
    std::vector<size_t> N(2);
    N[0] = N_rows;
    N[1] = N_columns;
    allocate(N);
}
template<class TYPE>
Array<TYPE>::Array( size_t N1, size_t N2, size_t N3 )
{
    std::vector<size_t> N(3);
    N[0] = N1;
    N[1] = N2;
    N[2] = N3;
    allocate(N);
}
template<class TYPE>
Array<TYPE>::Array( const std::vector<size_t>& N )
{
    allocate(N);
}
template<class TYPE>
void Array<TYPE>::allocate( const std::vector<size_t>& N )
{
    d_ndim = N.size();
    d_length = 1;
    for (size_t i=0; i<sizeof(d_N)/sizeof(size_t); i++)
        d_N[i] = 1;
    for (size_t i=0; i<N.size(); i++) {
        d_N[i] = N[i];
        d_length *= N[i];
    }
    if ( N.empty() ) {
        d_N[0] = 0;
        d_length = 0;
    }
    d_ptr.reset();
    if ( d_length > 0 )
        d_ptr = std::shared_ptr<TYPE>(new TYPE[d_length],DeleteArray<TYPE>);
    d_data = d_ptr.get();
    if ( d_length>0 && d_data==NULL )
        ERROR("Failed to allocate array");
}
template<class TYPE>
Array<TYPE>::Array( const Array& rhs ):
    d_ndim(rhs.d_ndim), d_length(rhs.d_length), d_data(NULL)
{
    allocate( std::vector<size_t>(rhs.d_N,rhs.d_N+rhs.d_ndim) );
    for (size_t i=0; i<d_length; i++)
        d_data[i] = rhs.d_data[i];
}
template<class TYPE>
Array<TYPE>& Array<TYPE>::operator=( const Array& rhs )
{
    if ( this == &rhs ) 
        return *this;
    this->allocate( rhs.size() );
    for (size_t i=0; i<d_length; i++)
        this->d_data[i] = rhs.d_data[i];
    return *this;
}
template<class TYPE>
Array<TYPE>::~Array( )
{
}


/********************************************************
*  Resize the array                                     *
********************************************************/
template<class TYPE>
void Array<TYPE>::resize( size_t N )
{
    resize(std::vector<size_t>(1,N));
}
template<class TYPE>
void Array<TYPE>::resize( size_t N1, size_t N2 )
{
    std::vector<size_t> N(2);
    N[0] = N1;
    N[1] = N2;
    resize(N);
}
template<class TYPE>
void Array<TYPE>::resize( size_t N1, size_t N2, size_t N3 )
{
    std::vector<size_t> N(3);
    N[0] = N1;
    N[1] = N2;
    N[2] = N3;
    resize(N);
}
template<class TYPE>
void Array<TYPE>::resize( const std::vector<size_t>& N )
{
    // Check if the array actually changed size
    size_t new_length = 1;
    for (size_t i=0; i<N.size(); i++)
        new_length *= N[i];
    bool changed = new_length!=d_length;
    for (size_t i=0; i<N.size(); i++)
        changed = changed || N[i]!=d_N[i];
    if ( !changed )
        return;
    // Store the old data
    const size_t ndim_max = sizeof(d_N)/sizeof(size_t);
    std::vector<size_t> N1(ndim_max,1), N2(ndim_max,1);
    for (size_t d=0; d<d_ndim; d++)
        N1[d] = d_N[d];
    for (size_t d=0; d<N.size(); d++)
        N2[d] = N[d];
    if ( d_ndim==0 ) { N1[0] = 0; }
    if ( N.empty() ) { N2[0] = 0; }
    std::shared_ptr<TYPE> old_data = d_ptr;
    // Allocate new data
    allocate(N);
    // Copy the old values
    if ( d_length > 0 ) {
        ASSERT(sizeof(d_N)/sizeof(size_t)==4);
        TYPE *data1 = old_data.get();
        TYPE *data2 = d_data;
        for (size_t m=0; m<std::min(N1[3],N2[3]); m++) {
            for (size_t k=0; k<std::min(N1[2],N2[2]); k++) {
                for (size_t j=0; j<std::min(N1[1],N2[1]); j++) {
                    for (size_t i=0; i<std::min(N1[0],N2[0]); i++) {
                        size_t index1 = i + j*N1[0] + k*N1[0]*N1[1] + m*N1[0]*N1[1]*N1[2];
                        size_t index2 = i + j*N2[0] + k*N2[0]*N2[1] + m*N2[0]*N2[1]*N1[2];
                        data2[index2] = data1[index1];
                    }
                }
            }
        }
    }
}


/********************************************************
*  Rehape the array                                     *
********************************************************/
template<class TYPE>
void Array<TYPE>::reshape( const std::vector<size_t>& N )
{
    size_t new_length = 1;
    for (size_t i=0; i<N.size(); i++)
        new_length *= N[i];
    if ( new_length!=d_length )
        ERROR("reshape is not allowed to change the array size");
    d_ndim = N.size();
    for (size_t i=0; i<sizeof(d_N)/sizeof(size_t); i++)
        d_N[i] = 1;
    for (size_t i=0; i<N.size(); i++)
        d_N[i] = N[i];
}


/********************************************************
*  Operator overloading                                 *
********************************************************/
template<class TYPE>
bool Array<TYPE>::operator==( const Array& rhs ) const
{
    if ( this==&rhs )
        return true;
    if ( d_length!=rhs.d_length )
        return false;
    bool match = true;
    for (size_t i=0; i<d_length; i++)
        match = match && d_data[i]==rhs.d_data[i];
    return match;
}


/********************************************************
*  Get a view of an C array                             *
********************************************************/
template<class TYPE>
std::shared_ptr<Array<TYPE> > Array<TYPE>::view( size_t N, std::shared_ptr<TYPE> data )
{
    view(std::vector<size_t>(1,N),data);
}
template<class TYPE>
std::shared_ptr<Array<TYPE> > Array<TYPE>::view( size_t N1, size_t N2, std::shared_ptr<TYPE> data )
{
    std::vector<size_t> N(2);
    N[0] = N1;
    N[1] = N2;
    view(N,data);
}
template<class TYPE>
std::shared_ptr<Array<TYPE> > Array<TYPE>::view( size_t N1, size_t N2, size_t N3, std::shared_ptr<TYPE> data )
{
    std::vector<size_t> N(3);
    N[0] = N1;
    N[1] = N2;
    N[2] = N3;
    view(N,data);
}
template<class TYPE>
std::shared_ptr<const Array<TYPE> > Array<TYPE>::constView( size_t N, std::shared_ptr<const TYPE> data )
{
    constView(std::vector<size_t>(1,N),data);
}
template<class TYPE>
std::shared_ptr<const Array<TYPE> > Array<TYPE>::constView( size_t N1, size_t N2, std::shared_ptr<const TYPE> data )
{
    std::vector<size_t> N(2);
    N[0] = N1;
    N[1] = N2;
    constView(N,data);
}
template<class TYPE>
std::shared_ptr<const Array<TYPE> > Array<TYPE>::constView( size_t N1, size_t N2, size_t N3, std::shared_ptr<const TYPE> data )
{
    std::vector<size_t> N(3);
    N[0] = N1;
    N[1] = N2;
    N[2] = N3;
    constView(N,data);
}
template<class TYPE>
std::shared_ptr<Array<TYPE> > Array<TYPE>::view( const std::vector<size_t>& N, std::shared_ptr<TYPE> data )
{
    std::shared_ptr<Array<TYPE> > array(new Array<TYPE>());
    array->d_ndim = N.size();
    array->d_length = 1;
    for (size_t i=0; i<N.size(); i++) {
        array->d_N[i] = N[i];
        array->d_length *= N[i];
    }
    array->d_ptr = data;
    array->d_data = array->d_ptr.get();
    return array;
}
template<class TYPE>
std::shared_ptr<const Array<TYPE> > Array<TYPE>::constView( const std::vector<size_t>& N, std::shared_ptr<const TYPE> data )
{
    return view(N,std::const_pointer_cast<TYPE>(data));
}


/********************************************************
*  Convert array types                                  *
********************************************************/
template<class TYPE>
template<class TYPE2>
std::shared_ptr<Array<TYPE2> > Array<TYPE>::convert( std::shared_ptr<Array<TYPE> > array )
{
    std::shared_ptr<Array<TYPE2> > array2( new Array<TYPE2>(array->size()) );
    array2.copy( *array );
    return array2;
}
template<class TYPE>
template<class TYPE2>
std::shared_ptr<const Array<TYPE2> > Array<TYPE>::convert( std::shared_ptr<const Array<TYPE> > array )
{
    return Array<TYPE>::convert( std::const_pointer_cast<Array<TYPE2> >(array) );
}
template<class TYPE>
template<class TYPE2>
void Array<TYPE>::copy( const Array<TYPE2>& array )
{
    resize( array.size() );
    const TYPE2 *src = array.get();
    for (size_t i=0; i<d_length; i++)
        d_data[i] = static_cast<TYPE>(src[i]);
}
template<class TYPE>
template<class TYPE2>
void Array<TYPE>::copy( const TYPE2* src )
{
    for (size_t i=0; i<d_length; i++)
        d_data[i] = static_cast<TYPE>(src[i]);
}
template<class TYPE>
void Array<TYPE>::fill( const TYPE& value )
{
    for (size_t i=0; i<d_length; i++)
        d_data[i] = value;
}

/********************************************************
*  std::swap                                            *
********************************************************/
template<class TYPE>
void std::swap(Array<TYPE>& v1, Array<TYPE>& v2)
{
    std::swap(v1.d_ndim,v2.d_ndim);
    std::swap(v1.d_N,v2.d_N);
    std::swap(v1.d_length,v2.d_length);
    std::swap(v1.d_data,v2.d_data);
    std::swap(v1.d_ptr,v2.d_ptr);
}


/********************************************************
*  Simple math operations                               *
********************************************************/
template<class TYPE>
bool Array<TYPE>::NaNs( ) const
{
    bool test = false;
    for (size_t i=0; i<d_length; i++)
        test = test || d_data[i]!=d_data[i];
    return test;
}
template<class TYPE>
TYPE Array<TYPE>::min( ) const
{
    TYPE x = std::numeric_limits<TYPE>::max();
    for (size_t i=0; i<d_length; i++)
        x = std::min(x,d_data[i]);
    return x;
}
template<class TYPE>
TYPE Array<TYPE>::max( ) const
{
    TYPE x = std::numeric_limits<TYPE>::min();
    for (size_t i=0; i<d_length; i++)
        x = std::max(x,d_data[i]);
    return x;
}
template<class TYPE>
TYPE Array<TYPE>::sum( ) const
{
    TYPE x = 0;
    for (size_t i=0; i<d_length; i++)
        x += d_data[i];
    return x;
}
template<class TYPE>
std::shared_ptr<Array<TYPE> > Array<TYPE>::sum( int dir ) const
{
    std::vector<size_t> size_ans = size();
    std::shared_ptr<Array<TYPE> > ans( new Array<TYPE>(size_ans) );
    size_t N1=1, N2=1, N3=1;
    for (int d=0; d<std::min(dir,d_ndim); d++)
        N1 *= d_N[d];
    N2 = d_N[dir];
    for (int d=dir+1; d<std::min(dir,d_ndim); d++)
        N3 *= d_N[d];
    TYPE* data2 = ans->d_data;
    for (int i3=0; i3<N3; i3++) {
        for (int i1=0; i1<N1; i1++) {
            TYPE x = 0;
            for (size_t i2=0; i2<N2; i2++)
                x += d_data[i1+i2*N1+i3*N1*N2];
            data2[i1+i3*N1] = x;
        }
    }
    return ans;
}



#endif

