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
#ifndef included_FunctionTable_hpp
#define included_FunctionTable_hpp

#include "common/FunctionTable.h"
#include "common/Utilities.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <random>


/********************************************************
 *  Random number initialization                         *
 ********************************************************/
template<class TYPE>
static inline typename std::enable_if<std::is_integral<TYPE>::value>::type genRand(
    size_t N, TYPE* x )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_int_distribution<TYPE> dis;
    for ( size_t i = 0; i < N; i++ )
        x[i] = dis( gen );
}
template<class TYPE>
static inline typename std::enable_if<std::is_floating_point<TYPE>::value>::type genRand(
    size_t N, TYPE* x )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<TYPE> dis( 0, 1 );
    for ( size_t i = 0; i < N; i++ )
        x[i] = dis( gen );
}
template<class TYPE, class FUN>
inline void FunctionTable::rand( Array<TYPE, FUN>& x )
{
    genRand<TYPE>( x.length(), x.data() );
}


/********************************************************
 *  Reduction                                            *
 ********************************************************/
template<class TYPE, class FUN, typename LAMBDA>
inline TYPE FunctionTable::reduce( LAMBDA& op, const Array<TYPE, FUN>& A )
{
    if ( A.length() == 0 )
        return TYPE();
    const TYPE* x  = A.data();
    TYPE y         = x[0];
    const size_t N = A.length();
    for ( size_t i = 1; i < N; i++ )
        y = op( x[i], y );
    return y;
}


/********************************************************
 *  Unary transformation                                 *
 ********************************************************/
template<class TYPE, class FUN, typename LAMBDA>
inline void FunctionTable::transform( LAMBDA& fun, const Array<TYPE, FUN>& x, Array<TYPE, FUN>& y )
{
    y.resize( x.size() );
    const size_t N = x.length();
    for ( size_t i = 0; i < N; i++ )
        y( i ) = fun( x( i ) );
}
template<class TYPE, class FUN, typename LAMBDA>
inline void FunctionTable::transform(
    LAMBDA& fun, const Array<TYPE, FUN>& x, const Array<TYPE, FUN>& y, Array<TYPE, FUN>& z )
{
    if ( x.size() != y.size() )
        throw std::logic_error( "Sizes of x and y do not match" );
    z.resize( x.size() );
    const size_t N = x.length();
    for ( size_t i = 0; i < N; i++ )
        z( i ) = fun( x( i ), y( i ) );
}


/********************************************************
 *  axpy                                                 *
 ********************************************************/
template<class TYPE>
inline void call_axpy( size_t N, const TYPE alpha, const TYPE* x, TYPE* y );
template<>
inline void call_axpy<float>( size_t, const float, const float*, float* )
{
    throw std::logic_error( "LapackWrappers not configured" );
}
template<>
inline void call_axpy<double>( size_t, const double, const double*, double* )
{
    throw std::logic_error( "LapackWrappers not configured" );
}
template<class TYPE>
inline void call_axpy( size_t N, const TYPE alpha, const TYPE* x, TYPE* y )
{
    for ( size_t i = 0; i < N; i++ )
        y[i] += alpha * x[i];
}
template<class TYPE, class FUN>
void FunctionTable::axpy( const TYPE alpha, const Array<TYPE, FUN>& x, Array<TYPE, FUN>& y )
{
    if ( x.size() != y.size() )
        throw std::logic_error( "Array sizes do not match" );
    call_axpy( x.length(), alpha, x.data(), y.data() );
}


/********************************************************
 *  Multiply two arrays                                  *
 ********************************************************/
template<class TYPE>
inline void call_gemv( size_t M, size_t N, TYPE alpha, TYPE beta, const TYPE* A, const TYPE* x, TYPE* y );
template<>
inline void call_gemv<double>(
    size_t, size_t, double, double, const double*, const double*, double* )
{
    throw std::logic_error( "LapackWrappers not configured" );
}
template<>
inline void call_gemv<float>( size_t, size_t, float, float, const float*, const float*, float* )
{
    throw std::logic_error( "LapackWrappers not configured" );
}
template<class TYPE>
inline void call_gemv(
    size_t M, size_t N, TYPE alpha, TYPE beta, const TYPE* A, const TYPE* x, TYPE* y )
{
    for ( size_t i = 0; i < M; i++ )
        y[i] = beta * y[i];
    for ( size_t j = 0; j < N; j++ ) {
        for ( size_t i = 0; i < M; i++ )
            y[i] += alpha * A[i + j * M] * x[j];
    }
}
template<class TYPE>
inline void call_gemm(
    size_t M, size_t N, size_t K, TYPE alpha, TYPE beta, const TYPE* A, const TYPE* B, TYPE* C );
template<>
inline void call_gemm<double>( size_t, size_t, size_t, double, double, const double*, const double*, double* )
{
    throw std::logic_error( "LapackWrappers not configured" );
}
template<>
inline void call_gemm<float>( size_t, size_t, size_t, float, float, const float*, const float*, float* )
{
    throw std::logic_error( "LapackWrappers not configured" );
}
template<class TYPE>
inline void call_gemm(
    size_t M, size_t N, size_t K, TYPE alpha, TYPE beta, const TYPE* A, const TYPE* B, TYPE* C )
{
    for ( size_t i = 0; i < K * M; i++ )
        C[i] = beta * C[i];
    for ( size_t k = 0; k < K; k++ ) {
        for ( size_t j = 0; j < N; j++ ) {
            for ( size_t i = 0; i < M; i++ )
                C[i + k * M] += alpha * A[i + j * M] * B[j + k * N];
        }
    }
}
template<class TYPE, class FUN>
void FunctionTable::gemm( const TYPE alpha, const Array<TYPE, FUN>& a, const Array<TYPE, FUN>& b,
    const TYPE beta, Array<TYPE, FUN>& c )
{
    if ( a.ndim() == 2 && b.ndim() == 1 ) {
        if ( a.size( 1 ) != b.size( 0 ) )
            throw std::logic_error( "Inner dimensions must match" );
        call_gemv<TYPE>( a.size( 0 ), a.size( 1 ), alpha, beta, a.data(), b.data(), c.data() );
    } else if ( a.ndim() <= 2 && b.ndim() <= 2 ) {
        if ( a.size( 1 ) != b.size( 0 ) )
            throw std::logic_error( "Inner dimensions must match" );
        call_gemm<TYPE>(
            a.size( 0 ), a.size( 1 ), b.size( 1 ), alpha, beta, a.data(), b.data(), c.data() );
    } else {
        throw std::logic_error( "Not finished yet" );
    }
}
template<class TYPE, class FUN>
void FunctionTable::multiply(
    const Array<TYPE, FUN>& a, const Array<TYPE, FUN>& b, Array<TYPE, FUN>& c )
{
    if ( a.ndim() == 2 && b.ndim() == 1 ) {
        if ( a.size( 1 ) != b.size( 0 ) )
            throw std::logic_error( "Inner dimensions must match" );
        c.resize( a.size( 0 ) );
        call_gemv<TYPE>( a.size( 0 ), a.size( 1 ), 1, 0, a.data(), b.data(), c.data() );
    } else if ( a.ndim() <= 2 && b.ndim() <= 2 ) {
        if ( a.size( 1 ) != b.size( 0 ) )
            throw std::logic_error( "Inner dimensions must match" );
        c.resize( a.size( 0 ), b.size( 1 ) );
        call_gemm<TYPE>(
            a.size( 0 ), a.size( 1 ), b.size( 1 ), 1, 0, a.data(), b.data(), c.data() );
    } else {
        throw std::logic_error( "Not finished yet" );
    }
}


/********************************************************
 *  Check if two arrays are equal                        *
 ********************************************************/
template<class TYPE, class FUN>
inline typename std::enable_if<!std::is_floating_point<TYPE>::value, bool>::type
FunctionTableCompare( const Array<TYPE, FUN>& a, const Array<TYPE, FUN>& b, TYPE )
{
    bool pass = true;
    if ( a.size() != b.size() )
        throw std::logic_error( "Sizes of x and y do not match" );
    for ( size_t i = 0; i < a.length(); i++ )
        pass = pass && a( i ) == b( i );
    return pass;
}
template<class TYPE, class FUN>
inline typename std::enable_if<std::is_floating_point<TYPE>::value, bool>::type
FunctionTableCompare( const Array<TYPE, FUN>& a, const Array<TYPE, FUN>& b, TYPE tol )
{
    bool pass = true;
    if ( a.size() != b.size() )
        throw std::logic_error( "Sizes of x and y do not match" );
    for ( size_t i = 0; i < a.length(); i++ )
        pass = pass && ( std::abs( a( i ) - b( i ) ) < tol );
    return pass;
}
template<class TYPE, class FUN>
bool FunctionTable::equals( const Array<TYPE, FUN>& a, const Array<TYPE, FUN>& b, TYPE tol )
{
    return FunctionTableCompare( a, b, tol );
}


#endif
