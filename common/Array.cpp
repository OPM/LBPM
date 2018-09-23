#include "common/Array.h"
#include "common/Array.hpp"

#include <complex>


/********************************************************
 *  ArraySize                                            *
 ********************************************************/
ArraySize::ArraySize( const std::vector<size_t>& N )
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
    for ( unsigned long i : d_N )
        d_length *= i;
    if ( d_ndim == 0 )
        d_length = 0;
}


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
template class Array<char, FunctionTable>;
template class Array<uint8_t, FunctionTable>;
template class Array<uint16_t, FunctionTable>;
template class Array<uint32_t, FunctionTable>;
template class Array<uint64_t, FunctionTable>;
template class Array<int8_t, FunctionTable>;
template class Array<int16_t, FunctionTable>;
template class Array<int32_t, FunctionTable>;
template class Array<int64_t, FunctionTable>;
template class Array<float, FunctionTable>;
template class Array<double, FunctionTable>;


/********************************************************
 *  Explicit instantiations of Array<bool>               *
 ********************************************************/
// clang-format off
template Array<bool, FunctionTable>::Array();
template Array<bool, FunctionTable>::~Array();
template Array<bool, FunctionTable>::Array( size_t );
template Array<bool, FunctionTable>::Array( size_t, size_t );
template Array<bool, FunctionTable>::Array( size_t, size_t, size_t );
template Array<bool, FunctionTable>::Array( size_t, size_t, size_t, size_t );
template Array<bool, FunctionTable>::Array( size_t, size_t, size_t, size_t, size_t );
template Array<bool, FunctionTable>::Array( const std::vector<size_t>&, const bool* );
template Array<bool, FunctionTable>::Array( std::string );
template Array<bool, FunctionTable>::Array( std::initializer_list<bool> );
template Array<bool, FunctionTable>::Array( const Array<bool, FunctionTable>& );
template Array<bool, FunctionTable>::Array( Array<bool, FunctionTable>&& );
template Array<bool, FunctionTable>& Array<bool, FunctionTable>::operator=( const Array<bool, FunctionTable>& );
template Array<bool, FunctionTable>& Array<bool, FunctionTable>::operator=( Array<bool, FunctionTable>&& );
template Array<bool, FunctionTable>& Array<bool, FunctionTable>::operator=( const std::vector<bool>& );
template void Array<bool, FunctionTable>::fill(bool const&);
template void Array<bool, FunctionTable>::clear();
template bool Array<bool, FunctionTable>::operator==(Array<bool, FunctionTable> const&) const;
template void Array<bool, FunctionTable>::resize( ArraySize const& );
// clang-format on


/********************************************************
 *  Explicit instantiations of Array<std::complex>       *
 ********************************************************/
// clang-format off
template Array<std::complex<double>, FunctionTable>::Array();
template Array<std::complex<double>, FunctionTable>::~Array();
template Array<std::complex<double>, FunctionTable>::Array( size_t );
template Array<std::complex<double>, FunctionTable>::Array( size_t, size_t );
template Array<std::complex<double>, FunctionTable>::Array( size_t, size_t, size_t );
template Array<std::complex<double>, FunctionTable>::Array( size_t, size_t, size_t, size_t );
template Array<std::complex<double>, FunctionTable>::Array( size_t, size_t, size_t, size_t, size_t );
template Array<std::complex<double>, FunctionTable>::Array( const std::vector<size_t>&, const std::complex<double>* );
template Array<std::complex<double>, FunctionTable>::Array( std::initializer_list<std::complex<double>> );
template Array<std::complex<double>, FunctionTable>::Array( const Range<std::complex<double>>& range );
template Array<std::complex<double>, FunctionTable>::Array( const Array<std::complex<double>, FunctionTable>& );
template Array<std::complex<double>, FunctionTable>::Array( Array<std::complex<double>, FunctionTable>&& );
template Array<std::complex<double>, FunctionTable>& Array<std::complex<double>, FunctionTable>::operator=( const Array<std::complex<double>, FunctionTable>& );
template Array<std::complex<double>, FunctionTable>& Array<std::complex<double>, FunctionTable>::operator=( Array<std::complex<double>, FunctionTable>&& );
template Array<std::complex<double>, FunctionTable>& Array<std::complex<double>, FunctionTable>::operator=( const std::vector<std::complex<double>>& );
template void Array<std::complex<double>, FunctionTable>::resize( ArraySize const& );
template void Array<std::complex<double>, FunctionTable>::viewRaw( ArraySize const&, std::complex<double>*, bool, bool );
template void Array<std::complex<double>, FunctionTable>::fill(std::complex<double> const&);
template void Array<std::complex<double>, FunctionTable>::clear();
template bool Array<std::complex<double>, FunctionTable>::operator==(Array<std::complex<double>, FunctionTable> const&) const;
template Array<std::complex<double>, FunctionTable> Array<std::complex<double>, FunctionTable>::repmat(std::vector<unsigned long, std::allocator<unsigned long> > const&) const;
// clang-format on


/********************************************************
 *  Explicit instantiations of Array<std::string>        *
 ********************************************************/
// clang-format off
template Array<std::string, FunctionTable>::Array();
template Array<std::string, FunctionTable>::~Array();
template Array<std::string, FunctionTable>::Array( size_t );
template Array<std::string, FunctionTable>::Array( size_t, size_t );
template Array<std::string, FunctionTable>::Array( size_t, size_t, size_t );
template Array<std::string, FunctionTable>::Array( size_t, size_t, size_t, size_t );
template Array<std::string, FunctionTable>::Array( size_t, size_t, size_t, size_t, size_t );
template Array<std::string, FunctionTable>::Array( const std::vector<size_t>&, const std::string* );
template Array<std::string, FunctionTable>::Array( std::initializer_list<std::string> );
template Array<std::string, FunctionTable>::Array( const Array<std::string, FunctionTable>& );
template Array<std::string, FunctionTable>::Array( Array<std::string, FunctionTable>&& );
template Array<std::string, FunctionTable>& Array<std::string, FunctionTable>::operator=( const Array<std::string, FunctionTable>& );
template Array<std::string, FunctionTable>& Array<std::string, FunctionTable>::operator=( Array<std::string, FunctionTable>&& );
template Array<std::string, FunctionTable>& Array<std::string, FunctionTable>::operator=( const std::vector<std::string>& );
template void Array<std::string, FunctionTable>::resize( ArraySize const& );
// clang-format on
