/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

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
#include "common/UtilityMacros.h"

#include <algorithm>
#include <cstring>
#include <limits>
//#include <random>

/********************************************************
 *  Random number initialization                         *
 ********************************************************/
/*template<class TYPE> TYPE genRand();
template<class TYPE, class FUN>
inline void FunctionTable::rand( Array<TYPE, FUN> &x )
{
    for ( size_t i = 0; i < x.length(); i++ )
        x( i ) = genRand<TYPE>();
}
*/

/********************************************************
 *  Reduction                                            *
 ********************************************************/
template <class TYPE, class FUN, typename LAMBDA>
inline TYPE FunctionTable::reduce(LAMBDA &op, const Array<TYPE, FUN> &A,
                                  const TYPE &initialValue) {
    if (A.length() == 0)
        return TYPE();
    const TYPE *x = A.data();
    TYPE y = initialValue;
    for (size_t i = 0; i < A.length(); i++)
        y = op(x[i], y);
    return y;
}
template <class TYPE, class FUN, typename LAMBDA>
inline TYPE FunctionTable::reduce(LAMBDA &op, const Array<TYPE, FUN> &A,
                                  const Array<TYPE, FUN> &B,
                                  const TYPE &initialValue) {
    ARRAY_ASSERT(A.length() == B.length());
    if (A.length() == 0)
        return TYPE();
    const TYPE *x = A.data();
    const TYPE *y = B.data();
    TYPE z = initialValue;
    for (size_t i = 0; i < A.length(); i++)
        z = op(x[i], y[i], z);
    return z;
}

/********************************************************
 *  Unary transformation                                 *
 ********************************************************/
template <class TYPE, class FUN, typename LAMBDA>
inline void FunctionTable::transform(LAMBDA &fun, const Array<TYPE, FUN> &x,
                                     Array<TYPE, FUN> &y) {
    y.resize(x.size());
    const size_t N = x.length();
    for (size_t i = 0; i < N; i++)
        y(i) = fun(x(i));
}
template <class TYPE, class FUN, typename LAMBDA>
inline void FunctionTable::transform(LAMBDA &fun, const Array<TYPE, FUN> &x,
                                     const Array<TYPE, FUN> &y,
                                     Array<TYPE, FUN> &z) {
    if (x.size() != y.size())
        throw std::logic_error("Sizes of x and y do not match");
    z.resize(x.size());
    const size_t N = x.length();
    for (size_t i = 0; i < N; i++)
        z(i) = fun(x(i), y(i));
}

/********************************************************
 *  axpy                                                 *
 ********************************************************/
template <class TYPE>
void call_axpy(size_t N, const TYPE alpha, const TYPE *x, TYPE *y);
template <>
void call_axpy<float>(size_t N, const float alpha, const float *x, float *y);
template <>
void call_axpy<double>(size_t N, const double alpha, const double *x,
                       double *y);
template <class TYPE>
void call_axpy(size_t N, const TYPE alpha, const TYPE *x, TYPE *y) {
    for (size_t i = 0; i < N; i++)
        y[i] += alpha * x[i];
}
template <class TYPE, class FUN>
void FunctionTable::axpy(const TYPE alpha, const Array<TYPE, FUN> &x,
                         Array<TYPE, FUN> &y) {
    if (x.size() != y.size())
        throw std::logic_error("Array sizes do not match");
    call_axpy(x.length(), alpha, x.data(), y.data());
}

/********************************************************
 *  Multiply two arrays                                  *
 ********************************************************/
template <class TYPE>
void call_gemv(size_t M, size_t N, TYPE alpha, TYPE beta, const TYPE *A,
               const TYPE *x, TYPE *y);
template <>
void call_gemv<double>(size_t M, size_t N, double alpha, double beta,
                       const double *A, const double *x, double *y);
template <>
void call_gemv<float>(size_t M, size_t N, float alpha, float beta,
                      const float *A, const float *x, float *y);
template <class TYPE>
void call_gemv(size_t M, size_t N, TYPE alpha, TYPE beta, const TYPE *A,
               const TYPE *x, TYPE *y) {
    for (size_t i = 0; i < M; i++)
        y[i] = beta * y[i];
    for (size_t j = 0; j < N; j++) {
        for (size_t i = 0; i < M; i++)
            y[i] += alpha * A[i + j * M] * x[j];
    }
}
template <class TYPE>
void call_gemm(size_t M, size_t N, size_t K, TYPE alpha, TYPE beta,
               const TYPE *A, const TYPE *B, TYPE *C);
template <>
void call_gemm<double>(size_t M, size_t N, size_t K, double alpha, double beta,
                       const double *A, const double *B, double *C);
template <>
void call_gemm<float>(size_t M, size_t N, size_t K, float alpha, float beta,
                      const float *A, const float *B, float *C);
template <class TYPE>
void call_gemm(size_t M, size_t N, size_t K, TYPE alpha, TYPE beta,
               const TYPE *A, const TYPE *B, TYPE *C) {
    for (size_t i = 0; i < K * M; i++)
        C[i] = beta * C[i];
    for (size_t k = 0; k < K; k++) {
        for (size_t j = 0; j < N; j++) {
            for (size_t i = 0; i < M; i++)
                C[i + k * M] += alpha * A[i + j * M] * B[j + k * N];
        }
    }
}
template <class TYPE, class FUN>
void FunctionTable::gemm(const TYPE alpha, const Array<TYPE, FUN> &a,
                         const Array<TYPE, FUN> &b, const TYPE beta,
                         Array<TYPE, FUN> &c) {
    if (a.size(1) != b.size(0))
        throw std::logic_error("Inner dimensions must match");
    if (a.ndim() == 2 && b.ndim() == 1) {
        call_gemv<TYPE>(a.size(0), a.size(1), alpha, beta, a.data(), b.data(),
                        c.data());
    } else if (a.ndim() <= 2 && b.ndim() <= 2) {
        call_gemm<TYPE>(a.size(0), a.size(1), b.size(1), alpha, beta, a.data(),
                        b.data(), c.data());
    } else {
        throw std::logic_error("Not finished yet");
    }
}
template <class TYPE, class FUN>
void FunctionTable::multiply(const Array<TYPE, FUN> &a,
                             const Array<TYPE, FUN> &b, Array<TYPE, FUN> &c) {
    if (a.size(1) != b.size(0))
        throw std::logic_error("Inner dimensions must match");
    if (a.ndim() == 2 && b.ndim() == 1) {
        c.resize(a.size(0));
        call_gemv<TYPE>(a.size(0), a.size(1), 1, 0, a.data(), b.data(),
                        c.data());
    } else if (a.ndim() <= 2 && b.ndim() <= 2) {
        c.resize(a.size(0), b.size(1));
        call_gemm<TYPE>(a.size(0), a.size(1), b.size(1), 1, 0, a.data(),
                        b.data(), c.data());
    } else {
        throw std::logic_error("Not finished yet");
    }
}

/********************************************************
 *  Check if two arrays are equal                        *
 ********************************************************/
template <class TYPE, class FUN>
inline typename std::enable_if<std::is_integral<TYPE>::value, bool>::type
FunctionTableCompare(const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b,
                     TYPE) {
    bool pass = true;
    if (a.size() != b.size())
        throw std::logic_error("Sizes of x and y do not match");
    for (size_t i = 0; i < a.length(); i++)
        pass = pass && a(i) == b(i);
    return pass;
}
template <class TYPE, class FUN>
inline typename std::enable_if<std::is_floating_point<TYPE>::value, bool>::type
FunctionTableCompare(const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b,
                     TYPE tol) {
    bool pass = true;
    if (a.size() != b.size())
        throw std::logic_error("Sizes of x and y do not match");
    for (size_t i = 0; i < a.length(); i++)
        pass = pass && (std::abs(a(i) - b(i)) < tol);
    return pass;
}
template <class TYPE, class FUN>
bool FunctionTable::equals(const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b,
                           TYPE tol) {
    return FunctionTableCompare(a, b, tol);
}

/********************************************************
 *  Specialized Functions                                *
 ********************************************************/
template <class TYPE, class FUN, class ALLOC>
void FunctionTable::transformReLU(const Array<TYPE, FUN, ALLOC> &A,
                                  Array<TYPE, FUN, ALLOC> &B) {
    const auto &fun = [](const TYPE &a) {
        return std::max(a, static_cast<TYPE>(0));
    };
    transform(fun, A, B);
}

template <class TYPE, class FUN, class ALLOC>
void FunctionTable::transformAbs(const Array<TYPE, FUN, ALLOC> &A,
                                 Array<TYPE, FUN, ALLOC> &B) {
    B.resize(A.size());
    const auto &fun = [](const TYPE &a) { return std::abs(a); };
    transform(fun, A, B);
}
template <class TYPE, class FUN, class ALLOC>
void FunctionTable::transformTanh(const Array<TYPE, FUN, ALLOC> &A,
                                  Array<TYPE, FUN, ALLOC> &B) {
    B.resize(A.size());
    const auto &fun = [](const TYPE &a) { return tanh(a); };
    transform(fun, A, B);
}

template <class TYPE, class FUN, class ALLOC>
void FunctionTable::transformHardTanh(const Array<TYPE, FUN, ALLOC> &A,
                                      Array<TYPE, FUN, ALLOC> &B) {
    B.resize(A.size());
    const auto &fun = [](const TYPE &a) {
        return std::max(-static_cast<TYPE>(1.0),
                        std::min(static_cast<TYPE>(1.0), a));
    };
    transform(fun, A, B);
}

template <class TYPE, class FUN, class ALLOC>
void FunctionTable::transformSigmoid(const Array<TYPE, FUN, ALLOC> &A,
                                     Array<TYPE, FUN, ALLOC> &B) {
    B.resize(A.size());
    const auto &fun = [](const TYPE &a) { return 1.0 / (1.0 + exp(-a)); };
    transform(fun, A, B);
}

template <class TYPE, class FUN, class ALLOC>
void FunctionTable::transformSoftPlus(const Array<TYPE, FUN, ALLOC> &A,
                                      Array<TYPE, FUN, ALLOC> &B) {
    B.resize(A.size());
    const auto &fun = [](const TYPE &a) { return log1p(exp(a)); };
    transform(fun, A, B);
}

template <class TYPE, class FUN, class ALLOC>
TYPE FunctionTable::sum(const Array<TYPE, FUN, ALLOC> &A) {
    const auto &fun = [](const TYPE &a, const TYPE &b) { return a + b; };
    return reduce(fun, A, (TYPE)0);
}

template <class TYPE>
inline void FunctionTable::gemmWrapper(char, char, int, int, int, TYPE,
                                       const TYPE *, int, const TYPE *, int,
                                       TYPE, TYPE *, int) {
    ERROR("Not finished");
}

#endif
