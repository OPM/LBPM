// This file contains helper functions and interfaces for reading/writing HDF5
#ifndef included_HDF5_hpp
#define included_HDF5_hpp
#ifdef USE_HDF5

#include "IO/HDF5_IO.h"
#include "common/Array.h"
#include "common/Array.hpp"
#include "common/Utilities.h"

#include <array>
#include <complex>
#include <memory>
#include <type_traits>
#include <vector>


namespace IO {
namespace HDF5 {


/********************************************************
 *  External instantiations (scalar)                     *
 ********************************************************/
// clang-format off
template<> void writeHDF5<char>( hid_t, const std::string &, const char & );
template<> void readHDF5<char>( hid_t, const std::string &, char & );
template<> void writeHDF5<bool>( hid_t, const std::string &, const bool & );
template<> void readHDF5<bool>( hid_t, const std::string &, bool & );
template<> void writeHDF5<int>( hid_t, const std::string &, const int & );
template<> void readHDF5<int>( hid_t, const std::string &, int & );
template<> void writeHDF5<long>( hid_t, const std::string &, const long & );
template<> void readHDF5<long>( hid_t, const std::string &, long & );
template<> void writeHDF5<float>( hid_t, const std::string &, const float & );
template<> void readHDF5<float>( hid_t, const std::string &, float & );
template<> void writeHDF5<double>( hid_t, const std::string &, const double & );
template<> void readHDF5<double>( hid_t, const std::string &, double & );
template<> void writeHDF5<unsigned char>( hid_t, const std::string &, const unsigned char & );
template<> void readHDF5<unsigned char>( hid_t, const std::string &, unsigned char & );
template<> void writeHDF5<unsigned int>( hid_t, const std::string &, const unsigned int & );
template<> void readHDF5<unsigned int>( hid_t, const std::string &, unsigned int & );
template<> void writeHDF5<unsigned long>( hid_t, const std::string &, const unsigned long & );
template<> void readHDF5<unsigned long>( hid_t, const std::string &, unsigned long & );
template<> void writeHDF5<std::string>( hid_t, const std::string &, const std::string & );
template<> void readHDF5<std::string>( hid_t, const std::string &, std::string & );
template<> void writeHDF5<std::complex<double>>( hid_t, const std::string &, const std::complex<double> & );
template<> void readHDF5<std::complex<double>>( hid_t, const std::string &, std::complex<double> & );
template<> void writeHDF5<std::complex<float>>( hid_t, const std::string &, const std::complex<float> & );
template<> void readHDF5<std::complex<float>>( hid_t, const std::string &, std::complex<float> & );
// clang-format on


/********************************************************
 *  External instantiations (Array)                      *
 ********************************************************/
// clang-format off
template<> void writeHDF5<Array<char>>( hid_t, const std::string &, const Array<char> & );
template<> void readHDF5<Array<char>>( hid_t, const std::string &, Array<char> & );
template<> void writeHDF5<Array<bool>>( hid_t, const std::string &, const Array<bool> & );
template<> void readHDF5<Array<bool>>( hid_t, const std::string &, Array<bool> & );
template<> void writeHDF5<Array<int>>( hid_t, const std::string &, const Array<int> & );
template<> void readHDF5<Array<int>>( hid_t, const std::string &, Array<int> & );
template<> void writeHDF5<Array<long>>( hid_t, const std::string &, const Array<long> & );
template<> void readHDF5<Array<long>>( hid_t, const std::string &, Array<long> & );
template<> void writeHDF5<Array<float>>( hid_t, const std::string &, const Array<float> & );
template<> void readHDF5<Array<float>>( hid_t, const std::string &, Array<float> & );
template<> void writeHDF5<Array<double>>( hid_t, const std::string &, const Array<double> & );
template<> void readHDF5<Array<double>>( hid_t, const std::string &, Array<double> & );
template<> void writeHDF5<Array<unsigned char>>( hid_t, const std::string &, const Array<unsigned char> & );
template<> void readHDF5<Array<unsigned char>>( hid_t, const std::string &, Array<unsigned char> & );
template<> void writeHDF5<Array<unsigned int>>( hid_t, const std::string &, const Array<unsigned int> & );
template<> void readHDF5<Array<unsigned int>>( hid_t, const std::string &, Array<unsigned int> & );
template<> void writeHDF5<Array<unsigned long>>( hid_t, const std::string &, const Array<unsigned long> & );
template<> void readHDF5<Array<unsigned long>>( hid_t, const std::string &, Array<unsigned long> & );
template<> void writeHDF5<Array<std::string>>( hid_t, const std::string &, const Array<std::string> & );
template<> void readHDF5<Array<std::string>>( hid_t, const std::string &, Array<std::string> & );
template<> void writeHDF5<Array<std::string>>( hid_t, const std::string &, const Array<std::string> & );
template<> void readHDF5<Array<std::string>>( hid_t, const std::string &, Array<std::string> & );
template<> void writeHDF5<Array<std::complex<double>>>( hid_t, const std::string &, const Array<std::complex<double>> & );
template<> void readHDF5<Array<std::complex<double>>>( hid_t, const std::string &, Array<std::complex<double>> & );
template<> void writeHDF5<Array<std::complex<float>>>( hid_t, const std::string &, const Array<std::complex<float>> & );
template<> void readHDF5<Array<std::complex<float>>>( hid_t, const std::string &, Array<std::complex<float>> & );
// clang-format on


/******************************************************************
 * Default implimentation                                          *
 ******************************************************************/
/*template<class TYPE>
void writeHDF5( hid_t fid, const std::string &name, const TYPE &x )
{
    NULL_USE( fid );
    if constexpr ( is_shared_ptr<TYPE>::value ) {
        // We are dealing with a std::shared_ptr
        writeHDF5( fid, name, *x );
    } else if constexpr ( is_vector<TYPE>::value ) {
        // We are dealing with a std::vector
        typedef decltype( *x.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        Array<TYPE4> y;
        y.viewRaw( { x.size() }, const_cast<TYPE4 *>( x.data() ) );
        writeHDF5( fid, name, y );
    } else if constexpr ( std::is_array<TYPE>::value ) {
        // We are dealing with a std::array
        typedef decltype( *x.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        Array<TYPE4> y;
        y.viewRaw( { x.size() }, const_cast<TYPE4 *>( x.data() ) );
        writeHDF5( fid, name, y );
    } else if constexpr ( is_Array<TYPE>::value ) {
        // We are dealing with an Array
        std::string typeName = Utilities::demangle( typeid( TYPE ).name() );
        throw std::logic_error( "Unsupported type writeHDF5<Array<" + typeName + ">>" );
    } else if constexpr ( std::is_same<TYPE, std::string>::value ) {
        // We are dealing with a std::string (should be handled through specialization)
        throw std::logic_error( "Internal error" );
    } else if constexpr ( std::is_same<TYPE, std::string>::value ||
                          std::is_same<TYPE, char *>::value ||
                          std::is_same<TYPE, const char *>::value ) {
        // We are dealing with a string or char array
        writeHDF5( fid, name, std::string( x ) );
    } else if constexpr ( has_size<TYPE>::value ) {
        // We are dealing with a container
        typedef decltype( *x.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        std::vector<TYPE4> x2( x.begin(), x.end() );
        writeHDF5<std::vector<TYPE4>>( fid, name, x2 );
    } else {
        throw std::logic_error( "Unsupported type" );
    }
}
template<class TYPE>
void readHDF5( hid_t fid, const std::string &name, TYPE &x )
{
    NULL_USE( fid );
    if constexpr ( is_shared_ptr<TYPE>::value ) {
        // We are dealing with a std::shared_ptr
        readHDF5( fid, name, *x );
    } else if constexpr ( is_vector<TYPE>::value ) {
        // We are dealing with a std::vector
        typedef typename std::remove_reference<decltype( *x.begin() )>::type TYPE2;
        Array<TYPE2> y;
        readHDF5( fid, name, y );
        x.resize( y.length() );
        // Swap the elements in the arrays to use the move operator
        for ( size_t i = 0; i < x.size(); i++ )
            std::swap( x[i], y( i ) );
    } else if constexpr ( std::is_array<TYPE>::value ) {
        // We are dealing with a std::array
        typedef typename std::remove_reference<decltype( *x.begin() )>::type TYPE2;
        Array<TYPE2> y;
        readHDF5( fid, name, y );
        ASSERT( y.length() == x.size() );
        // Swap the elements in the arrays to use the move operator
        for ( size_t i = 0; i < x.size(); i++ )
            std::swap( x[i], y( i ) );
    } else if constexpr ( is_Array<TYPE>::value ) {
        // We are dealing with an Array
        std::string typeName = Utilities::demangle( typeid( TYPE ).name() );
        throw std::logic_error( "Unsupported type readHDF5<Array<" + typeName + ">>" );
    } else if constexpr ( std::is_same<TYPE, std::string>::value ) {
        // We are dealing with a std::string (should be handled through specialization)
        throw std::logic_error( "Internal error" );
    } else if constexpr ( std::is_same<TYPE, std::string>::value ||
                          std::is_same<TYPE, char *>::value ||
                          std::is_same<TYPE, const char *>::value ) {
        // We are dealing with a string or char array
        throw std::logic_error(
            "Reading data into a string, char*, const char* is not supported" );
    } else if constexpr ( has_size<TYPE>::value ) {
        // We are dealing with a container
        typedef typename std::remove_reference<decltype( *x.begin() )>::type TYPE2;
        Array<TYPE2> y;
        readHDF5( fid, name, y );
        if ( x.size() == y.length() ) {
            auto it = x.begin();
            for ( size_t i = 0; i < y.length(); i++, ++it )
                *it = y( i );
        } else {
            throw std::logic_error( "Reading data into an arbitrary container is not finished" );
        }
    } else {
        throw std::logic_error( "Unsupported type" );
    }
}*/


/************************************************************************
 * Helper function to get the size of an Array                           *
 * Note that HDF5 uses C ordered arrays so we need to flip the dimensions*
 ************************************************************************/
template<class T>
inline std::vector<hsize_t> arraySize( const Array<T> &x )
{
    int N   = x.ndim();
    auto s1 = x.size();
    std::vector<hsize_t> s2( std::max( N, 1 ), 0 );
    for ( int i = 0; i < N; i++ )
        s2[N - i - 1] = static_cast<hsize_t>( s1[i] );
    return s2;
}
inline std::vector<size_t> convertSize( int N, const hsize_t *dims )
{
    if ( N == 0 )
        return std::vector<size_t>( 1, 1 );
    std::vector<size_t> size( N, 0 );
    for ( int i = 0; i < N; i++ )
        size[N - i - 1] = static_cast<size_t>( dims[i] );
    return size;
}


/************************************************************************
 * readAndConvertHDF5Data                                                *
 ************************************************************************/
template<class T>
typename std::enable_if<std::is_integral<T>::value || std::is_floating_point<T>::value, void>::type
readAndConvertHDF5Data( hid_t dataset, hid_t datatype, Array<T> &data )
{
    if ( H5Tequal( datatype, H5T_NATIVE_CHAR ) ) {
        Array<char> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_UCHAR ) ) {
        Array<unsigned char> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_INT8 ) ) {
        Array<int8_t> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_UINT8 ) ) {
        Array<uint8_t> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_INT ) ) {
        Array<int> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_UINT ) ) {
        Array<unsigned int> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_LONG ) ) {
        Array<long int> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_ULONG ) ) {
        Array<unsigned long int> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_FLOAT ) ) {
        Array<float> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else if ( H5Tequal( datatype, H5T_NATIVE_DOUBLE ) ) {
        Array<double> data2( data.size() );
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2.data() );
        data.copy( data2 );
    } else {
        ERROR( "We need to convert unknown data format" );
    }
}
template<class T>
typename std::enable_if<!std::is_integral<T>::value && !std::is_floating_point<T>::value,
    void>::type
readAndConvertHDF5Data( hid_t, hid_t, Array<T> & )
{
    ERROR( "Unable to convert data" );
}


/************************************************************************
 * Default writeHDF5Array                                                *
 ************************************************************************/
template<class T>
void writeHDF5ArrayDefault( hid_t fid, const std::string &name, const Array<T> &data )
{
    size_t N_bytes = data.length() * sizeof( T );
    auto dim       = arraySize( data );
    hid_t plist    = H5P_DEFAULT;
    if ( N_bytes < 0x7500 ) {
        // Use compact storage (limited to < 30K)
        plist       = H5Pcreate( H5P_DATASET_CREATE );
        auto status = H5Pset_layout( plist, H5D_COMPACT );
        ASSERT( status == 0 );
    } else if ( std::is_same<T, double>::value || std::is_same<T, float>::value ) {
        // Use compression if availible
        plist = createChunk( dim, defaultCompression( fid ) );
    }
    hid_t dataspace = H5Screate_simple( dim.size(), dim.data(), NULL );
    hid_t datatype  = getHDF5datatype<T>();
    hid_t dataset =
        H5Dcreate2( fid, name.data(), datatype, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT );
    const void *ptr = data.data() == NULL ? ( (void *) 1 ) : data.data();
    H5Dwrite( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr );
    H5Dclose( dataset );
    H5Tclose( datatype );
    H5Sclose( dataspace );
    if ( plist != H5P_DEFAULT )
        H5Pclose( plist );
}


/************************************************************************
 * Default readHDF5Array                                                 *
 ************************************************************************/
template<class T>
void readHDF5ArrayDefault( hid_t fid, const std::string &name, Array<T> &data )
{
    if ( !H5Dexists( fid, name ) ) {
        // Dataset does not exist
        data.resize( 0 );
        return;
    }
    hid_t dataset   = H5Dopen2( fid, name.data(), H5P_DEFAULT );
    hid_t datatype  = H5Dget_type( dataset );
    hid_t dataspace = H5Dget_space( dataset );
    hsize_t dims0[10];
    int ndim  = H5Sget_simple_extent_dims( dataspace, dims0, NULL );
    auto dims = convertSize( ndim, dims0 );
    data.resize( dims );
    hid_t datatype2 = getHDF5datatype<T>();
    if ( data.empty() ) {
        // The data is empty
    } else if ( H5Tequal( datatype, datatype2 ) ) {
        // The type of Array and the data in HDF5 match
        H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data() );
    } else {
        // Try to convert the data
        readAndConvertHDF5Data( dataset, datatype, data );
    }
    H5Dclose( dataset );
    H5Tclose( datatype );
    H5Tclose( datatype2 );
    H5Sclose( dataspace );
}


} // namespace HDF5
} // namespace IO


#endif
#endif
